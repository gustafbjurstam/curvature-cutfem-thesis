%% 
clc
clear
close all
% Convergence 2D
base_dir = '~/Work/repositories/SF250X/cutfem-workfiles/output_files/curvature_2d_algoim/';
examples = {
            "Peanut_elem_Tris", [2e-6, 1e-1];
            "Star_elem_Tris", [2e-5, 5e0];
            "Amoeba_elem_Tris", [1e-6, 3e-2];
            "Peanut_surf_Tris", [2e-6, 1e-1];
            "Star_surf_Tris", [2e-5, 5e0];
            "Amoeba_surf_Tris", [1e-6, 3e-2];
            };


for i = 1:length(examples)
    process_example_data(base_dir, examples{i,1}, [1,2], true, examples{i,2});
end
