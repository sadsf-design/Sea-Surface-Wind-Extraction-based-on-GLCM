clear
clc

load bb.mat

% 找列的最大值
[m, ~] = max(bast_p, [], 1);
[~, ci] = max(m);

% 找行的最大值
[m, ~] = max(bast_p, [], 2);
[~, ri] = max(m);

x = 1:size(bast_p,1);
y = 1:size(bast_p,2);
[X,Y] = meshgrid(x,y);
surf(X',Y',bast_p)
shading interp
axis off

ru = ri - 1;
rd = size(bast_p, 1) - ri;

cu = ci - 1;
cd = size(bast_p, 2) - ci;

best_p = zeros(200, 360);
r = 100;
c = 224;

ru_r = min(r - 1, ru);
rd_r = min(size(bast_p, 1) - ri, rd);

cu_r = min(c - 1, cu);
cd_r = min(size(bast_p, 2) - c, cd);

if r < 180
    best_p(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = bast_p(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r);
else
    best_p(r - ru_r:r + rd_r, c - cu_r:c + cd_r) = flip(bast_p(ri - ru_r:ri + rd_r, ci - cu_r:ci + cd_r), 2);
end
figure
surf(X',Y', best_p)
shading interp
axis off