# Endpoint-clearance follow-up: locked design

The author approved two additional endpoint pairs per map following the observation that 589/630 previous paths attain an endpoint-imposed clearance bound. These six follow-up instances were selected after the earlier results were known and before any follow-up optimization. The original 630 records and their original comparison family remain unchanged.

## Questions and cases

For each of the same three transformed terrain grids, the two new cases use start `(0, 200)` and goal `(200, 0)`. In the near condition each endpoint lies nine model units vertically above its own terrain; in the high condition it lies twenty units above. The two cases on a map have identical horizontal geometry, identical change in endpoint altitude, and identical direct endpoint distance. Their vertical anchors differ by eleven units. The original mission's horizontal diagonal differs, so these are a small targeted follow-up, not a random sample of terrain missions.

`cases.json` records exact endpoint values, terrain provenance and a separate seven-vertex feasibility witness for each case. The witnesses are only geometric feasibility checks; they are not used to initialize or guide any method. Each is certified by the independent cellwise bilinear-clearance routine before optimization.

Raising both endpoints changes their absolute altitude as well as their clearance. Because the objective includes mean altitude and the feasible paths can change, the near/high contrast represents the joint endpoint-height intervention on the fixed map and route direction. It does not isolate the effect of one objective term, establish flight safety, or by itself establish generalization across maps.

## Conditions and execution

Seven complete methods, their frozen checkpoints, native representations, repairs and operators, the four-term objective, clearance threshold eight, feasibility-first ordering, and search budget of 100,100 evaluations are the same as the original production study. Three fixed maps by two endpoint levels by seven methods by 30 independent runs gives 1,260 additional primary paths. Previously saved results are kept separate. The root `run_comparison.m` selects a new explicit suite that dispatches to the separately recorded runner and uses the original frozen V2 scientific source. The 714-file source manifest, root launcher, new driver, case definition, verifier and supervisor identities are checked.

Seeds are `1000000 + (caseId-201)*100000 + runId*10 + methodIndex`, where method indices are one-based and primary run IDs are 1–30. All 1,260 seeds are unique and disjoint from the 42 pilot seeds (run ID 91). Frozen MATLAB methods retain their original Threefry/Inversion profile, and GLEET retains its native NumPy/PyTorch sampled-action profile. Distinct seeds do not prove cross-platform or parallel reproducibility. GLEET receives the same integer seed in an integer JSON representation; this corrects a smoke-test serialization failure without changing the numeric seed or frozen algorithm.

The 42 pilot runs, one per case and method at 20,300 search evaluations, precede production. They activate the AFSACPSO temporal encoder and are excluded from primary analysis. Gate success means that all outputs, budgets, source/model identities and independently recomputed geometry/objective agree. A returned infeasible route is an experimental outcome, not a validation failure by itself. A separate preliminary 800-evaluation smoke test is also excluded; the first smoke was incomplete and remains archived with its six saved outputs.

The host scheduler starts up to 14 independent single-thread MATLAB workers with RAM and disk guards. It uses isolated per-case/per-method directories and atomic per-run saves. If a worker or verification fails, it stops the stage and preserves checkpoints. Production starts only after all 42 pilot paths pass independent verification. Runtime measurements under mixed machine load are descriptive only.

## Declared analysis

Within each case and method, report feasibility rate and Wilson interval; conditional on feasibility, report median and IQR of objective, proportion of paths attaining the endpoint-imposed clearance bound (tolerance `1e-8` model units), and median arc-length harmonic clearance. For each map and method, compare endpoint-binding counts between the near and high groups using a two-sided Fisher exact test, Holm-corrected over the 21 planned map–method comparisons. These tests are exploratory follow-up analyses. If either feasible sample is empty, mark the test untestable and retain a value of one for multiplicity accounting. Conditioning on feasibility must be discussed, especially when feasibility differs between groups.

Method outcomes are compared within the same case; raw objective differences across endpoint levels cannot be interpreted as an isolated effect of clearance because the altitude term also changes. The three original cases and six follow-up cases are not silently pooled into a single independent-population estimate.
