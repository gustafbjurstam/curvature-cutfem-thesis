clear; clc; close all;

base_dir = '~/Work/repositories/SF250X/cutfem-workfiles/output_files/curvature_p1_2d';
examples = [
            % "Peanut_elem_Quads";
            % "Star_elem_Quads";
            % "Amoeba_elem_Quads";
            % "Peanut_surf_Quads";
            % "Star_surf_Quads";
            % "Amoeba_surf_Quads";

            "Peanut_elem_Tris";
            "Star_elem_Tris";
            "Amoeba_elem_Tris";
            "Peanut_patch_Tris";
            "Star_patch_Tris";
            "Amoeba_patch_Tris";
            ];

for i = 1:length(examples)
    process_condition_number(base_dir, examples(i), true);
end