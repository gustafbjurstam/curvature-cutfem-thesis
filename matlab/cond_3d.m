clear; clc; close all;

base_dir = '../results/curvature_p1_3d/';
examples = [
    "Elliptic_Torus_patch";
    "Elliptic_Torus_face";
    "Decocube_patch";
    "Decocube_face"
];

for i = 1:length(examples)
    process_condition_number(base_dir, examples(i), false);
end