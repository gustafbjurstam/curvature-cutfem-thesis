# Run from anywhere; paths are resolved relative to this file.
base_dir = normpath(joinpath(@__DIR__, "..", "results", "curvature_2d_algoim"))

examples = [
    ("Peanut_elem_Tris", (2e-6, 1e-1)),
    ("Star_elem_Tris",   (2e-5, 5e0)),
    ("Amoeba_elem_Tris", (1e-6, 3e-2)),
    ("Peanut_surf_Tris", (2e-6, 1e-1)),
    ("Star_surf_Tris",   (2e-5, 5e0)),
    ("Amoeba_surf_Tris", (1e-6, 3e-2)),
]

reference_rates = [1.0, 2.0]

include(joinpath(@__DIR__, "convplotting.jl"))

for (example, y_limits) in examples
    process_example_data(base_dir, example, reference_rates, true, y_limits)
end
