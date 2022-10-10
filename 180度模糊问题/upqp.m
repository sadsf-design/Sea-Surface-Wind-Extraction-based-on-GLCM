function P = upqp(P2, Q2, c)
load('bb.mat')
res1 = load('res1.mat');
res2 = load('res2.mat');

res11 = res1.results.reps;
res22 = res2.results.reps;
num = size(res11, 3);

[rc1, rc2] = size(P2);
[rc3, rc4] = size(Q2);

% 找列的最大值
[m, ~] = max(bast_p, [], 1);
[~, ci] = max(m);

% 找行的最大值
[m, ~] = max(bast_p, [], 2);
[~, ri] = max(m);

ru = ri - 1;
rd = size(bast_p, 1) - ri;

cu = ci - 1;
cd = size(bast_p, 2) - ci;

if numel(nnz(rand(rc1+rc2-rc3-rc4 + 2 , 1) > 0.5)) > 1
    r = floor(100+randn);
else
    r = floor(100+randn+randn);
end

best_p = zeros(200, 360);

ru_r = min(r - 1, ru);
rd_r = min(size(bast_p, 1) - ri, rd);

cu_r = min(c - 1, cu);
cd_r = min(size(bast_p, 2) - c, cd);

if r < 180
    best_p_temp = bast_p(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r);
    [rl, cl] = size(best_p_temp);
    best_p_temp2(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = best_p_temp.*mexResize(abs(res11(:, :, floor(num*c/360))), [rl, cl], 'auto');
    
    [m, ~] = max(best_p_temp2, [], 1);
    [~, ci] = max(m);
    
    [m, ~] = max(best_p_temp2, [], 2);
    [~, ri] = max(m);
    
    ru = ri - 1;
    rd = size(best_p_temp2, 1) - ri;

    cu = ci - 1;
    cd = size(best_p_temp2, 2) - ci;
    
    ru_r = min(r - 1, ru);
    rd_r = min(size(best_p_temp2, 1) - ri, rd);

    cu_r = min(c - 1, cu);
    cd_r = min(size(best_p_temp2, 2) - c, cd);
    
    best_p(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = best_p_temp2(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r);
else
    best_p_temp = flip(bast_p(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r), 2);
    [rl, cl] = size(best_p_temp);
    best_p_temp2(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = best_p_temp.*mexResize(abs(res22(:, :, floor(num*c/360))), [rl, cl], 'auto');
    
    [m, ~] = max(best_p_temp2, [], 1);
    [~, ci] = max(m);
    
    [m, ~] = max(best_p_temp2, [], 2);
    [~, ri] = max(m);
    
    ru = ri - 1;
    rd = size(best_p_temp2, 1) - ri;

    cu = ci - 1;
    cd = size(best_p_temp2, 2) - ci;
    
    ru_r = min(r - 1, ru);
    rd_r = min(size(best_p_temp2, 1) - ri, rd);

    cu_r = min(c - 1, cu);
    cd_r = min(size(best_p_temp2, 2) - c, cd);
    
    best_p(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = best_p_temp2(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r);
end
P = best_p;
end