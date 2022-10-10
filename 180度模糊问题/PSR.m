%
%   Reliable Patch Trackers: Robust Visual Tracking by Exploiting Reliable Patches
%
%   Yang Li, 2015
%   http://ihpdep.github.io
%
%   This is the research code of our RPT tracker. You can check the
%   details in our CVPR paper.

function [ psr ] = PSR( response, rate )
%PSR Summary of this function goes here
%   Detailed explanation goes here
% 取最大值
maxresponse = max(response(:));
sz = size(response);
% 计算系数k，4/图像面积
k = 4/(sz(1)*sz(2));%8 /2
% calculate the PSR
% 计算range，N*rate, N为图像的像素个数，其值为sz(1)*sz(2)
range = ceil(sqrt(numel(response))*rate);
% response=fftshift(response);
% 得到最大值的纵横坐标
[xx, yy] = find(response == maxresponse, 1);

idx = 1:sz(1);
idy = 1:sz(2);
% repmat(a, [r, c])函数的作用是将矩阵a，在行方向上复制r次，列方向上复制c次
idx = repmat(idx,[sz(2),1])' - xx;
idy = repmat(idy,[sz(1),1]) - yy;
% 经过运算，idx和idy都与response的大小相同
% idx表示的是各行各列的点与峰值点之间x轴方向（行）的距离
% idy表示的是各行各列的点与峰值点之间y轴方向（列）的距离

t = idx.^2 + idy.^2; % t表示各行各列的点与峰值点间的欧氏距离

delta = 1 - exp(-k*t); % 构造函数，该函数的作用是使得各函数值离峰值点越近，值越小，越远，越大

r = (maxresponse - response)./delta; % 继续构造，求出综合值的矩阵。

psr = min(r(:)); %取所有值中最小的值

% 
% idx = xx-range:xx+range;
% idy = yy-range:yy+range;
% idy(idy<1)=1;idx(idx<1)=1;
% idy(idy>size(response,2))=size(response,2);idx(idx>size(response,1))=size(response,1);
% response(idx,idy)=0;
% m = sum(response(:))/numel(response);
% d=sqrt(size(response(:),1)*var(response(:))/numel(response));
% psr =(maxresponse - m)/d ;

end



