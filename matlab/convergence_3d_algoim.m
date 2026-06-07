clc
clear
close all
% Convergence 3D
base_dir = '../results/curvature_3d_algoim';
examples = {
            "Elliptic_Torus_algoim_q1_elem_Hexa", [4e-2, 3e1];
            "Decocube_algoim_q1_elem_Hexa", [1e0, 2e2];
            
            
            % "Elliptic_Torus_algoim_q1_patch_Hexa", [4e-2, 3e1];
            % "Decocube_algoim_q1_patch_Hexa", [1e0, 2e2];
            

            % "Elliptic_Torus_algoim_P2_patch_face_Tet", [4e-2, 3e1];
            % "Decocube_algoim_P2_patch_face_Tet", [1e0, 5e1];

            "Elliptic_Torus_algoim_P2_elem_face_Tet", [1e-3, 1e0];
            "Decocube_algoim_P2_elem_face_Tet", [1e0,2e1]
            };


for i = 1:size(examples,1)
    process_example_data(base_dir, examples{i,1}, [1,2], false, examples{i,2});
end
