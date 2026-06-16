# Run from anywhere; paths are resolved relative to this file.
base_dir = normpath(joinpath(@__DIR__, "..", "results", "curvature_3d_algoim"))

examples = [
    ("Elliptic_Torus_algoim_Q1_elem_Hexa",    (4e-2, 3e1)),
    ("Decocube_algoim_Q1_elem_Hexa",          (1e0, 2e2)),
    ("Elliptic_Torus_algoim_P2_elem_face_Tet", (1e-3, 1e0)),
    ("Decocube_algoim_P2_elem_face_Tet",       (1e0, 2e1)),
]

reference_rates = [1.0]

include(joinpath(@__DIR__, "convplotting.jl"))

for (example, y_limits) in examples
    process_example_data(base_dir, example, reference_rates, true, y_limits)
end
