# Clearance-definition intervention

Authorization: the author accepted the proposed scientific reframing and controlled constraint comparison on12 November2026. This is a new cohort. Neither the630 original routes nor the running1260 endpoint runs are replaced or pooled with it.

## Scientific question

Determine how imposing an eight-unit Euclidean route-to-surface margin changes achieved separation, feasibility and geometric route costs under an otherwise identical UAV planning problem. This tests a modeling decision, not a new optimizer or controller. Lower observed objective in the Euclidean condition is possible for stochastic heuristics and will not be suppressed to impose an expected direction.

## Conditions and instances

Use the exact original adapted urban101, suburban102 and mountain103 grids, altitude ceilings and endpoints. All six endpoints independently have Euclidean lower bounds above8. A seven-vertex high-route witness satisfies both conditions on each map, establishing feasibility within the same15-dimensional path representation. Cases and witnesses are recorded before optimization in cases.json.

Vertical condition retains the original full-segment vertical margin8. Euclidean condition retains that vertical requirement and adds minimum distance8 to the original continuous bilinear graph. The four original objective components, equal weights, normalizers, waypoint count, bounds and complete search algorithms are identical across conditions. The objective still measures reciprocal vertical clearance, so changing the constraint does not silently change the cost definition. No architectural solid or physical vehicle is introduced.

PSO_Standard and L-SHADE represent two search families without learned controllers. Their implementations, parameters, initialization and boundary treatments are copied from the original frozen campaign. The sole score-wrapper intervention introduces spatial feasibility. Invalid routes and vertical-contact routes preserve their original handling. Other candidates use min(verticalMargin/8,EuclideanLower/8) for worst-margin violation. Original J is unchanged wherever defined. Numerical Euclidean lower bounds define conservative feasibility. Width is at most0.001; intervals straddling8 are refined to1e-7 and retained as conservative failures if still unresolved. Counts of unresolved calls and all convergence failures are retained. A convergence failure stops the task.

## Runs and budget

Smoke:800 search evaluations, run91,12 method/map/condition combinations. Gate:20300 evaluations, run92, the same12 combinations. Primary:100100 evaluations, runs1 through30,360 routes. Each algorithm receives the same budget in both conditions, including initialization. Additional spatial-query counts and wall time are reported separately. Stage cohorts are isolated and pilot runs never enter primary tests.

Seeds are3000000+(map-101)*100000+10*run+methodIndex, with methodIndex1 forPSO and2 forL-SHADE. Two conditions share each map/method/run seed under threefry-inversion-v1. Initial candidate populations are saved and must match exactly within pairs. Subsequent search trajectories may diverge. No cross-platform reproducibility claim follows from the seed alone.

## Outcomes and inference

Primary outcome is the paired difference in the probability that the returned route meets both modeled margins. Report the four paired success/failure counts and an exact two-sided McNemar test for each of six map/method blocks, with Holm correction over the six tests. Numerical threshold-ambiguous returned routes remain separately flagged and count conservatively as failures.

Secondary outcomes are paired differences in original J, length, squared turning, datum altitude, vertical exposure, Euclidean distance and wall time. Cost contrasts require a defined original J in both arms; report the number of retained pairs and excluded outcomes. A route-pair bootstrap with10000 resamples gives pointwise median-difference intervals, conditional on those pairs. No broad superiority or equal-training inference is drawn. Distance lower/upper bounds and their threshold classes are recomputed independently on every returned route. Failure rates remain visible even when cost comparisons condition on feasibility.

## Gates and resource policy

Component checks must confirm unchanged original J and baseline score, six admissible endpoints and three common feasible path witnesses. Geometry tests include independent analytical/minimization cases and agreement on the630 original saved routes. Smoke/gate require correct source identity, budgets, score/path agreement, unchanged terrain/endpoints/objective, exact paired initial populations and independent returned-route geometry. Gates do not select maps/methods or demand favorable performance. Failures stop expansion and retain evidence.

Use the isolated source root's run_comparison.m as the sole optimizer entry point. The live repository launcher is unchanged because the ongoing endpoint campaign pins its hash. The isolated root contains the same frozen algorithm/shared sources, an explicitly modified constraint wrapper and a dedicated dispatch branch. Record every source and binary hash before the first pilot.

At most two intervention workers run while the endpoint campaign remains active. Gate and primary expansion wait for successful independent verification of that campaign, then use at most12 single-thread MATLAB workers for the12 intervention tasks. Source, protocol and resume-state identity are checked on every launch. Historical results are never overwritten. A failed or uncertain task is inspected before any retry. Timing limits are reported as failures, not removed observations.
