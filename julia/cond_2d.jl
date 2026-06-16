# Run from anywhere; paths are resolved relative to this file.
base_dir = normpath(joinpath(@__DIR__, "..", "results", "curvature_p1_2d"))

examples = [
    "Peanut_elem_Quads",
    "Star_elem_Quads",
    "Amoeba_elem_Quads",
    "Peanut_surf_Quads",
    "Star_surf_Quads",
    "Amoeba_surf_Quads",

    # Uncomment these if condition-number matrices were exported for the triangular runs.
    # "Peanut_elem_Tris",
    # "Star_elem_Tris",
    # "Amoeba_elem_Tris",
    # "Peanut_patch_Tris",
    # "Star_patch_Tris",
    # "Amoeba_patch_Tris",
]

include(joinpath(@__DIR__, "convplotting.jl"))
include(joinpath(@__DIR__, "condplotting.jl"))

for example in examples
    process_condition_number(base_dir, example, true)
end
