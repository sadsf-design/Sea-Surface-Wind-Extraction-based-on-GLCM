clear
clc

repition = 10;
seqs = {['texture'], 15;
        ['optics'], 10;
        ['redar'], 40};
Losses = zeros(size(seqs, 1), 6);
%% results_calculation
for i = 1 : size(seqs, 1)
    load(['re_', seqs{i, 1}, '.mat'])
    load(['gt_', seqs{i, 1}, '.mat'])
    eval(['res_data = re_', seqs{i, 1}, ';'])
    eval(['gt_data = gt_', seqs{i, 1}, ';'])
    mr = sum(mean(res_data));
    vr = sum(std(res_data, 1, 1));
    res = res_data(:);
    gt = gt_data(:);
    num = length(gt);
%% weighted MSE Loss (W-MSE)
    lambda = rand;
    k = 1e-2;
    Loss = 0;
    for j = 1 : num
        weight = k/(exp(-k*(gt(j)-lambda)));
        Loss = Loss + 1/num * (weight*(res(j)-gt(j)).^2);
    end
    Losses(i, 1) = Loss;
%% Weighted Binary Cross Entropy (W-BCE)
    lambda = rand;
    Loss = 0;
    for j = 1 : num
        weight = k/(exp(-k*(gt(j)-lambda)));
        Loss = Loss - (weight*(gt(j)+1)*log(res(j)+1)+(1-weight)*((gt(j))*log(res(j))));
    end
    Losses(i, 2) = Loss;
%% Normalized Scanpath Saliency (NSS)
    fix = 3;
    Loss = 0;
    for j = 1 : num
        Loss = Loss + (1/num)*((res(j) - mr)/vr)*(weight*gt(j)).^fix;
    end
    Losses(i, 3) = Loss;
%% Pearson's Correlation Coefficient (CC)
    cc = cov(res, gt);
    Loss = cc(1,2)/(cc(1,1)*cc(2,2));
    Losses(i, 4) = Loss;
%% Kullback-Leibler Divergence (KLD)
    Loss = 0;
    for j = 1 : num
        Loss = Loss + (res(j))*log(res(j)/gt(j));
    end
    Losses(i, 5) = Loss;
%% Sqrt Exponential Absolute Eifference (SEAD)
    Loss = 0;
    for j = 1 : num
        Loss = Loss + (1/num)*(exp(abs(res(j)-gt(j)))-1);
    end
    Losses(i, 6) = sqrt(Loss);
end
