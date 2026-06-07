clear; clc; close all;

base_dir = '~/Work/repositories/SF250X/cutfem-workfiles/output_files/curvature_p1_3d/';
examples = [
    "Elliptic_Torus_patch";
    "Elliptic_Torus_face";
    % "Decocube_patch";
    % "Decocube_face"
];

for i = 1:length(examples)
    process_condition_number(base_dir, examples(i), false);
end