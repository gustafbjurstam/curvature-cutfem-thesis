# Run from anywhere; paths are resolved relative to this file.
base_dir = normpath(joinpath(@__DIR__, "..", "results", "curvature_2d_algoim"))

examples = [
    "Peanut_elem_Tris",
    "Star_elem_Tris",
    "Amoeba_elem_Tris",
    "Peanut_surf_Tris",
    "Star_surf_Tris",
    "Amoeba_surf_Tris",
]

include(joinpath(@__DIR__, "convplotting.jl"))
include(joinpath(@__DIR__, "condplotting.jl"))

for example in examples
    process_condition_number(base_dir, example, true)
end
