clear
clc
close all

gt = [35,37,37,38,39,39,39,38,39,40];
res = [35.05,37.54,37.53,38.78,39.40,39.40,39.41,38.78,39.42,40.06];

%% 相关系数
cc = corrcoef(res, gt);
ccc = cc(2)

%% 标准差
g_R_std = sqrt(sum((gt-res)*((gt-res)')/length(res)))


%% 偏差
pc = mean(res-(gt+g_R_std))

num = length(gt);
Losses = zeros(1, length(gt));

for j = 1 : num
    Losses(j) = (res(j)/gt(j))*log(res(j)/gt(j));
end