clear; clc; close all;

base_dir = '~/Work/repositories/SF250X/cutfem-workfiles/output_files/curvature_2d_algoim';
examples = [
            "Peanut_elem_Tris";
            "Star_elem_Tris";
            "Amoeba_elem_Tris";

            "Peanut_surf_Tris";
            "Star_surf_Tris";
            "Amoeba_surf_Tris";
            ];

for i = 1:length(examples)
    process_condition_number(base_dir, examples(i), false);
end