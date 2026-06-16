# Run from anywhere; paths are resolved relative to this file.
base_dir = normpath(joinpath(@__DIR__, "..", "results", "curvature_p1_3d"))

examples = [
    ("Elliptic_Torus_patch", (4e-2, 3e1)),
    ("Elliptic_Torus_face",  (4e-2, 3e1)),
    ("Decocube_patch",      (1e0, 2e2)),
    ("Decocube_face",       (1e0, 2e2)),
]

reference_rates = [1.0]

include(joinpath(@__DIR__, "convplotting.jl"))

for (example, y_limits) in examples
    process_example_data(base_dir, example, reference_rates, true, y_limits)
end
