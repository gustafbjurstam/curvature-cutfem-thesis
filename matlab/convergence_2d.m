%% 
clc
clear
close all
% Convergence 2D
base_dir = '~/Work/repositories/SF250X/cutfem-workfiles/output_files/curvature_p1_2d';
examples = {
            % "Peanut_elem_Quads", [8e-3, 1.2e0];
            % "Star_elem_Quads", [5e-2, 1e1];
            % "Amoeba_elem_Quads", [1e-3, 5e-1];
            % "Peanut_surf_Quads", [8e-3, 1.2e0];
            % "Star_surf_Quads", [5e-2, 1e1];
            % "Amoeba_surf_Quads", [1e-3, 5e-1];

            "Peanut_elem_Tris", [5e-4, 1e0];
            "Star_elem_Tris", [1e-2, 1e1];
            "Amoeba_elem_Tris", [1e-4, 5e-1];
            "Peanut_patch_Tris", [5e-4, 1e0];
            "Star_patch_Tris", [1e-2, 1e1];
            "Amoeba_patch_Tris", [1e-4, 5e-1];
            "Peanut_face_only_Tris", [5e-4, 1e0];
            "Star_face_only_Tris", [1e-2, 1e1];
            "Amoeba_face_only_Tris", [1e-4, 5e-1];
            };


for i = 1:length(examples)
    process_example_data(base_dir, examples{i,1}, [1], true, examples{i,2});
end
