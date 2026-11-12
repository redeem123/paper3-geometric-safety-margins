"""Reanalyze the saved-route data archive with NumPy and SciPy."""
from argparse import ArgumentParser
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import fisher_exact, mannwhitneyu

from clearance_profiles import from_route
from distance_bounds import Terrain, distance


ROOT = Path(__file__).resolve().parent
SCENARIOS = (101, 102, 103)
METHODS = ('AFSACPSO', 'PSO_Standard', 'SAEPSO', 'SAC_SAPSO', 'RLAM_PSO', 'LSHADE', 'GLEET')
FIGURE_EXAMPLE = (101, 'AFSACPSO', 1)


def read_csv(relative):
    with (ROOT / relative).open(newline='') as file:
        return list(csv.DictReader(file))


def run_key(row):
    return int(row['scenario']), row['method'], int(row['run'])


def holm(values):
    values = np.asarray(values)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1)))
    return adjusted


def verify_cost_statistics(groups):
    summaries = read_csv('summary.csv')
    assert len(summaries) == len(groups) == 21
    for row in summaries:
        group = groups[int(row['scenario']), row['method']]
        costs = [float(record['J']) for record in group]
        expected = [float(row[k]) for k in ('J_q1', 'J_median', 'J_q3')]
        assert len(costs) == int(row['runs']) == int(row['feasible']) == 30
        np.testing.assert_allclose(np.quantile(costs, [.25, .5, .75]), expected, atol=1e-12, rtol=0)

    tests = read_csv('comparisons.csv')
    assert len(tests) == 18
    raw_score, raw_feasibility = [], []
    for row in tests:
        scenario = int(row['scenario'])
        a, b = (groups[scenario, method] for method in ('AFSACPSO', row['comparator']))
        raw_score.append(mannwhitneyu(
            [float(r['score']) for r in a], [float(r['score']) for r in b],
            alternative='two-sided', method='asymptotic').pvalue)
        feasible = [sum(r['band'] == '3' for r in group) for group in (a, b)]
        raw_feasibility.append(fisher_exact([[n, 30 - n] for n in feasible]).pvalue)
    for name, values in [('score', raw_score), ('feasibility', raw_feasibility)]:
        np.testing.assert_allclose(values, [float(row[f'p_{name}']) for row in tests], rtol=1e-9, atol=0)
        np.testing.assert_allclose(holm(values), [float(row[f'holm_{name}']) for row in tests],
                                   rtol=1e-9, atol=0)


def verify_reversals(groups):
    rows = read_csv('profiles/reversal_contrasts.csv')
    by_pair = {(int(r['scenario']), r['methodA'], r['methodB']): r for r in rows}
    counts = []
    for scenario in SCENARIOS:
        for a, b in combinations(METHODS, 2):
            x, y = (groups[scenario, method] for method in (a, b))
            ax, by = [np.array([[float(r['fR']), float(r['fR']) / float(r['fL'])]
                                for r in group]) for group in (x, y)]
            differences = ax[:, None] - by[None, :]
            sign = np.where(np.abs(differences) <= 1e-7, 0, np.sign(differences))
            count = int(np.count_nonzero(sign[:, :, 0] * sign[:, :, 1] < 0))
            saved = by_pair[scenario, a, b]
            assert count == int(saved['reversal_pairs'])
            assert np.count_nonzero(np.any(sign == 0, axis=2)) == int(saved['tie_pairs'])
            counts.append(count / 900)
    assert len(counts) == len(rows) == 63
    return [round(x * 100, 2) for x in np.quantile(counts, [0, .25, .5, .75, 1])]


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--euclidean', choices=('saved', 'example', 'all'), default='saved',
                        help='verify saved point pairs or also recompute continuous bounds')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'SHA256.json').read_text())
    for name, digest in manifest.items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name

    rows = read_csv('verified_runs.csv')
    assert len(rows) == len({run_key(r) for r in rows}) == 630
    groups, vertices = defaultdict(list), defaultdict(list)
    for row in read_csv('data/returned_routes.csv'):
        vertices[run_key(row)].append((int(row['vertex']), [float(row[k]) for k in ('x', 'y', 'z')]))
    bounds = {run_key(row): row for row in read_csv('euclidean/route_distances.csv')}
    plan = json.loads((ROOT / 'euclidean/CONTINUOUS_DISTANCE_PLAN.json').read_text())
    envs = {s: loadmat(ROOT / f'data/maps/scenario{s}.mat', simplify_cells=True)['env'] for s in SCENARIOS}
    surfaces = {s: RegularGridInterpolator((env['terrainY'][:, 0], env['terrainX'][0]), env['terrainGrid'])
                for s, env in envs.items()}
    terrain_cache = {}
    sensitivity = {s: [0, 0, 0] for s in SCENARIOS}
    endpoint_counts = {s: 0 for s in SCENARIOS}
    example = None
    for row in rows:
        s, method, run = key = run_key(row)
        groups[s, method].append(row)
        assert int(row['seed']) == 1000 * run + 100 * s + METHODS.index(method) + 1
        ids, points = zip(*sorted(vertices[key]))
        assert ids == tuple(range(7))
        path = np.array(points)
        env = envs[s]
        np.testing.assert_array_equal(path[[0, -1]], [env['startPoint'], env['goalPoint']])
        profile = from_route(path, env)
        assert abs(profile.minimum - float(row['clearance'])) < 1e-8
        assert abs(profile.inverse_profile(64) - float(row['fR']) / float(row['fL'])) < 1e-7
        ends = path[[0, -1], 2] - surfaces[s](path[[0, -1]][:, [1, 0]])
        endpoint_counts[s] += int(abs(profile.minimum - ends.min()) <= 1e-8)

        bound = bounds[key]
        point_pair = json.loads((ROOT / f'euclidean/results/scenario{s}_{method}_run{run}.json').read_text())
        assert bound['rawSHA256'] == point_pair['rawSHA256'] == row['sha256']
        lower, upper = [float(bound[name]) for name in ('lower', 'upper')]
        assert (lower, upper) == (point_pair['lower'], point_pair['upper'])
        assert 0 <= upper - lower <= .001000000001
        q, u = np.array(point_pair['pathPoint']), np.array(point_pair['surfacePoint'])
        i, t = point_pair['segment'], point_pair['segmentFraction']
        assert 0 <= t <= 1
        np.testing.assert_allclose(q, path[i] + t * (path[i + 1] - path[i]), rtol=0, atol=1e-8)
        assert abs(surfaces[s](u[[1, 0]])[0] - u[2]) < 1e-7
        assert abs(np.linalg.norm(q - u) - upper) < 1e-8
        for index, threshold in enumerate((6, 7, 8)):
            assert upper < threshold or lower >= threshold
            sensitivity[s][index] += int(upper < threshold)
        if key == FIGURE_EXAMPLE:
            example = dict(vertical_at_route_point=float(q[2] - surfaces[s](q[[1, 0]])[0]),
                           lower=lower, upper=upper)
        if args.euclidean == 'all' or (args.euclidean == 'example' and key == FIGURE_EXAMPLE):
            if s not in terrain_cache:
                terrain_cache[s] = Terrain(env)
            recomputed = distance(path, terrain_cache[s], tolerance=plan['primaryAbsoluteGap'],
                                  max_pops=plan['maximumPopsPerPass'], initial_points=[u])
            if recomputed['lower'] < 8 <= recomputed['upper']:
                recomputed = distance(path, terrain_cache[s], tolerance=plan['nearThresholdGap'],
                                      max_pops=plan['maximumPopsPerPass'],
                                      initial_points=[recomputed['surfacePoint']])
            assert recomputed['converged']
            assert recomputed['lower'] <= upper + 1e-8 and recomputed['upper'] >= lower - 1e-8
            assert recomputed['upper'] < 8 if upper < 8 else recomputed['lower'] >= 8

    verify_cost_statistics(groups)
    print(json.dumps(dict(routes=len(rows), mean_median_groups=21, adjusted_score_tests=18,
                          thresholds=[6, 7, 8], counts_by_map=sensitivity,
                          endpoint_minima_by_map=endpoint_counts, reversal_quantiles_percent=verify_reversals(groups),
                          figure2b=example, euclidean_recomputation=args.euclidean), indent=2))


if __name__ == '__main__':
    main()
