# Paper 3 saved-route analysis

This archive supports reanalysis of 630 original routes and 2,160 additional routes in four separate cohorts.
It includes the three transformed bilinear maps, seven vertices per route,
run-level costs, summaries, fixed model identities, and continuous Euclidean
point pairs and distance bounds.

Run the following commands from the extracted directory with Python 3, NumPy
and SciPy installed.

    python reanalyze_saved_routes.py
    python reanalyze_saved_routes.py --euclidean example

The first command verifies the saved files and point pairs, recomputes vertical
minima, normalized reciprocal exposure, reversal counts, the 21 median/quartile
summaries and the two 18-test families. It also counts routes below each
Euclidean threshold using the recorded numerical bounds. The second command
additionally recomputes the Figure 2(b) continuous distance bounds. Pass
--euclidean all to recompute bounds for every saved route. The geometry code
uses floating-point allowances and produces numerical bounds.

- summary.csv contains the 21 group summaries and comparisons.csv contains
  the adjusted and unadjusted tests.
- data/returned_routes.csv uses model units and original vertex ordering.
  data/maps/ contains the six numerical environment fields needed for geometry.
- profiles/ contains route profiles, contrast counts, and archived pointwise
  bootstrap intervals. euclidean/ contains per-route bounds and point pairs.
- model_identities.json lists the four original weight files and their hashes.
  Replaying the original search also requires the frozen optimizer source,
  model binaries and MATLAB/PyTorch environments.
- SHA256.json verifies the archived files. The run CSV uses relative identities
  in the original-file column. Numeric records and raw-result hashes are retained.

Versioned release location: https://github.com/redeem123/paper3-geometric-safety-margins/releases/tag/v1.0.0; no DOI is assigned.

## Additional cohorts and revised analyses

    python reanalyze_additional_routes.py
    python reanalyze_additional_routes.py --euclidean examples
    python reanalyze_geometry_statistics.py
    python test_distance_bounds.py

The first command recomputes vertical clearance and all defined objective terms
for all 2,160 additional paths, verifies spatial witnesses and numerical intervals,
and verifies exact initial-population equality for the 180 constraint pairs.
It also recomputes sensitivity endpoint distances and 22 feasible geometric
existence witnesses. Five random missions have spatially inadmissible endpoints;
their 100 unavoidable below-margin paths are reported separately from the other
440 paths, including 80 from the four admissible random missions.
Use --euclidean all for a fresh continuous spatial assessment of all 900 routes
in the constraint and sensitivity cohorts. Endpoint routes were studied for
vertical binding, not spatial compliance. The original reanalyze_saved_routes.py
remains the entry point for the 630-route cohort.

Vertices are stored losslessly as decimal doubles. rawSHA256 identifies the
original native MAT record; those larger native optimizer records remain in the
research workspace. Each group retains its own endpoints and definitions. The
original_assessment folders preserve the original numerical receipts, including
original path identities where present; they are not executable local paths.

reanalyze_geometry_statistics.py additionally requires pandas and reproduces
trimmed minima, slope associations, all 63 cost contrasts, materiality-filtered
reversals, weight rescoring and the vertical-clearance profile bootstrap curves.
The weight-rescoring CSV retains the baseline under each component label;
deduplicating those four identical baseline settings yields nine distinct
weights per map and exposure definition, as reported in the manuscript.
reanalyze_surface_controls.py reproduces both triangulations, regridding and
spatial exposure examples; it can take substantial time. It uses the native
engine, included for Apple Silicon macOS. On another platform, build the C++17
source as a shared library named engine/libbilinear_distance.dylib using
`c++ -std=c++17 -O3 -shared -fPIC engine/bilinear_distance.cpp -o engine/libbilinear_distance.dylib`.
Rebuilding changes artifact
hashes: preserve the original package and record a new local manifest.

The separate PAPER3_OPTIMIZER_REPLAY.zip contains full optimization source,
map transformations, adapters and all four pretrained model binaries. Replaying
search requires MATLAB and the method-specific Python environment; the commands
above do not run an optimizer or train a policy.
