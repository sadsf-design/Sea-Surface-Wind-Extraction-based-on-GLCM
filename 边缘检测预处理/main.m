clear
close all
clc
% Parameters
Gamma = 0.98;
Phi = 200;
Epsilon = -0.1;
k = 1.6;
Sigma = 0.8;

fileID = fopen('D:/CoLSY/new_picture/list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);
% inputIm = imread('sampleImage.jpg');

root_name = 'D:/CoLSY/new_picture/';
save_name = 'D:/CoLSY/边缘检测预处理/边缘检测图像及结果/';

% 获取具体图片名称序列
names = raw_names{1};

thetas = zeros(1, length(names));

for ww = 1 : length(names)
    
inputIm = imread([root_name, names{ww}]);
inputIm = rgb2gray(inputIm);
inputIm = im2double(inputIm);

% Gauss Filters
gFilteredIm1 = imgaussfilt(inputIm, Sigma);
gFilteredIm2 = imgaussfilt(inputIm, Sigma * k);

differencedIm2 = gFilteredIm1 - (Gamma * gFilteredIm2);

x = size(differencedIm2,1);
y = size(differencedIm2,2);

% Extended difference of gaussians
for i=1:x
    for j=1:y
        if differencedIm2(i, j) < Epsilon
            differencedIm2(i, j) = 1;
        else
            differencedIm2(i, j) = 1 + tanh(Phi*(differencedIm2(i,j)));
        end
    end
end

% XDoG Filtered Image
% figure, imshow(mat2gray(differencedIm2));

XDOGFilteredImage = mat2gray(differencedIm2);

% take mean of XDoG Filtered image to use in thresholding operation
meanValue = mean2(XDOGFilteredImage);

x = size(XDOGFilteredImage,1);
y = size(XDOGFilteredImage,2);

% thresholding
for i=1:x
    for j=1:y
        if XDOGFilteredImage(i, j) <= meanValue
            XDOGFilteredImage(i, j) = 0.0;
        else 
            XDOGFilteredImage(i, j) = 1.0;
        end
    end
end
[~, amount, label] = connect4( XDOGFilteredImage );
% figure, imshow(label)
imwrite(mat2gray(label), [save_name, names{ww}]);

% [~, stack, label2] = connect4( label );



% index = find(amount(1,:)==0,1);
% [~, yi] = min(amount(1,:));
ddd = 20;
yi = size(label, 2) + (0:-1:-ddd);
id2 = [];
for j = 1 : ddd + 1
    if ~isempty(find(label(:,yi(j))==0, 1))
        id2 = [id2, find(label(:,yi(j))==0, 1)];
    end
end
% id2 = find(amount(2, :) == yi);

[~, yi2] = min(id2);

d2 = [id2(yi2);yi(yi2)];

% x1 = max(amount(1,:));
% in2 = [];
% th = 3;
% subs = max(amount(1,:)) + (0:-1:-(th-1));
% for i = 1 : th
%     in2 = [in2, find(amount(1,:)==subs(i))];
% end
% y1 = median(amount(2,in2(:))) - std(amount(2,in2(:)));

dddd = 3;
yi2 = size(label, 1) + (0:-1:-dddd);
id22 = [];

for i = 1 : dddd+1
    if ~isempty(find(label(yi2(i),:)==0, 1))
        id22 = [id22, find(label(yi2(i),:)==0, 1)];
    end
end

x1 = size(label, 1);
y1 = median(id2) - std(id2);


thetas(ww) = atand(abs(x1-d2(1))/abs(y1-d2(2)));

end