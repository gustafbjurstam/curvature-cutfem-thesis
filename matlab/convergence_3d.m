%% 
clc
clear
close all
% Convergence 3D
base_dir = '../results/curvature_p1_3d';
examples = {
            "Elliptic_Torus_patch", [4e-2, 3e1];
            "Elliptic_Torus_face", [4e-2, 3e1];
            % "Elliptic_Torus_patch_elem", [4e-2, 3e1];
            "Decocube_patch", [1e0, 2e2];
            "Decocube_face", [1e0, 2e2];
            % "Decocube_patch_elem", [1e0, 2e2];
            };


for i = 1:size(examples,1)
    process_example_data(base_dir, examples{i,1}, [1], true, examples{i,2});
end
