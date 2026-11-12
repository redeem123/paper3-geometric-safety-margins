# Geometric safety margins in terrain-based UAV path planning

Research artifacts accompanying *Assessing geometric safety margins in
terrain-based UAV path planning*, by Trung-Kien Hoang, Thi Thanh-Hoa Lai and
Viet-Anh Ngo.

The study compares continuous vertical clearance with distance to a bilinear
terrain surface, then examines endpoint placement, exposure normalization and
the effect of imposing the spatial margin during search. It evaluates geometric
polylines on three transformed height fields; it does not certify physical flight.

## Data and reproducibility

Use the fixed release [v1.0.0](https://github.com/redeem123/paper3-geometric-safety-margins/releases/tag/v1.0.0), rather than a moving branch, to reproduce the
reported analyses. The release contains:

- [PAPER3_COMPLETE_ANALYSIS.zip](https://github.com/redeem123/paper3-geometric-safety-margins/releases/download/v1.0.0/PAPER3_COMPLETE_ANALYSIS.zip): the transformed terrain grids, all returned
  polylines and run-level outcomes from four distinct cohorts; independent
  evaluators and statistical analyses; numerical distance bounds and witnesses.
- [PAPER3_OPTIMIZER_REPLAY.zip](https://github.com/redeem123/paper3-geometric-safety-margins/releases/download/v1.0.0/PAPER3_OPTIMIZER_REPLAY.zip): the frozen search implementations, map
  transformations, deployment adapters, four actual pretrained model files,
  dependency information and root-entry-point replay commands.
- [SHA256SUMS.txt](https://github.com/redeem123/paper3-geometric-safety-margins/releases/download/v1.0.0/SHA256SUMS.txt): checksums of the release assets.

| Cohort | Search methods | Returned routes | Scientific purpose |
|---|---:|---:|---|
| Original fixed missions | 7 | 630 | Compare vertical and spatial assessment |
| Endpoint height intervention | 7 | 1,260 | Test endpoint limitation |
| Paired search constraints | 2 | 360 | Test enforcing the Euclidean margin |
| Waypoints, thresholds and additional missions | 2 | 540 | Examine sensitivity across 27 configurations |

These cohorts retain separate conditions, seed schedules and denominators.
They must not be pooled as a homogeneous sample of environments. Details and
verification commands are in each archive's README. Reanalyzing saved routes
requires Python, NumPy and SciPy; additional statistical scripts use pandas.
Replaying search requires MATLAB and, for GLEET, the recorded Python environment.
Fixed seeds do not imply cross-platform bitwise reproducibility.

The sensitivity cohort retains all nine randomly drawn missions. Five have an
endpoint already below the spatial threshold (cases 407, 408, 425–427); their
100 routes cannot meet that threshold. Separate summaries cover the other 440
routes and the 80 routes from the four random missions with admissible endpoints.
Geometric witness routes establish feasibility for all 22 endpoint-admissible
configurations. These witnesses are separate from optimizer outputs.

## Provenance

Terrain arrays derive from Darlan and colleagues' benchmark at commit
[b2525897bf8c381e13c6997e9fe08dfd877e1e6c](https://github.com/daisondarlan/UAV-Path-Planning-Benchmark/tree/b2525897bf8c381e13c6997e9fe08dfd877e1e6c).
The manuscript states the independent horizontal and relief transformations;
the evaluated problem is not the original multiobjective B-spline benchmark.

The fixed AFSACPSO controller and cellwise vertical-clearance calculation are
attributed to [Ngo et al.](https://doi.org/10.1016/j.swevo.2026.102599).
The common objective, new search cohorts and additional assessments in this
study are identified separately.

GLEET uses the upstream source at commit
[e7159273d4dc7ea78ef3ae3a701069cfd4d61191](https://github.com/MetaEvo/GLEET/tree/e7159273d4dc7ea78ef3ae3a701069cfd4d61191)
and its supplied `func1_10.pt` model. L-SHADE derives from Ryoji Tanabe's
[version 1.0.1 MATLAB implementation](https://ryojitanabe.github.io/publication).
Upstream authorship and notices remain with the archived code. This repository
does not impose a new blanket license on inherited components.

## Interpretation

The original 350/630 spatial-below-margin result is conditional on the original
three missions and the eight-model-unit comparison. A route below this margin
is not necessarily colliding. Numerical distance intervals and statistical
uncertainty intervals serve different purposes. Cost summaries describe the
implementations and fixed policies tested, rather than an intrinsic ranking of
controller architectures.
