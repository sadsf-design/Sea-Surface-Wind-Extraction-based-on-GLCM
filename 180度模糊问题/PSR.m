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
% ȡ���ֵ
maxresponse = max(response(:));
sz = size(response);
% ����ϵ��k��4/ͼ�����
k = 4/(sz(1)*sz(2));%8 /2
% calculate the PSR
% ����range��N*rate, NΪͼ������ظ�������ֵΪsz(1)*sz(2)
range = ceil(sqrt(numel(response))*rate);
% response=fftshift(response);
% �õ����ֵ���ݺ�����
[xx, yy] = find(response == maxresponse, 1);

idx = 1:sz(1);
idy = 1:sz(2);
% repmat(a, [r, c])�����������ǽ�����a�����з����ϸ���r�Σ��з����ϸ���c��
idx = repmat(idx,[sz(2),1])' - xx;
idy = repmat(idy,[sz(1),1]) - yy;
% �������㣬idx��idy����response�Ĵ�С��ͬ
% idx��ʾ���Ǹ��и��еĵ����ֵ��֮��x�᷽���У��ľ���
% idy��ʾ���Ǹ��и��еĵ����ֵ��֮��y�᷽���У��ľ���

t = idx.^2 + idy.^2; % t��ʾ���и��еĵ����ֵ����ŷ�Ͼ���

delta = 1 - exp(-k*t); % ���캯�����ú�����������ʹ�ø�����ֵ���ֵ��Խ����ֵԽС��ԽԶ��Խ��

r = (maxresponse - response)./delta; % �������죬����ۺ�ֵ�ľ���

psr = min(r(:)); %ȡ����ֵ����С��ֵ

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



