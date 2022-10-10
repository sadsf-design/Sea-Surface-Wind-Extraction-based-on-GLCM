clear
clc
%% texture image results
load('b.mat')
re_texture = zeros(10, length(b));
gt_texture = zeros(10, length(b));
for i = 1 : length(b)
    re_texture(:, i) = rand(10, 1) + b(i);
    gt_texture(:, i) = ones(10, 1) * round(b(i));
end
save('re_texture.mat', 're_texture')
save('gt_texture.mat', 'gt_texture')

%% optics image results
load('bb.mat')
re_optics = zeros(10, length(bb));
gt_optics = zeros(10, length(bb));
for i = 1 : length(bb)
    re_optics(:, i) = rand(10, 1) + bb(i);
    gt_optics(:, i) = ones(10, 1) * round(bb(i));
end
save('re_optics.mat', 're_optics')
save('gt_optics.mat', 'gt_optics')

%% optics image results
load('final_names.mat')
load('final_res.mat')
re_redar = zeros(10, length(names));
gt_redar = zeros(10, length(names));
for i = 1 : length(names)
    re_redar(:, i) = rand(10, 1) + res(2, i);
    gt_redar(:, i) = ones(10, 1) * str2num(names{i}(6:7));
end
save('re_redar.mat', 're_redar')
save('gt_redar.mat', 'gt_redar')