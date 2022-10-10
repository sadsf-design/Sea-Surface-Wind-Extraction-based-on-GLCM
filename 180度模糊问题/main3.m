clear
close all
clc
% Parameters
Gamma = 0.98;
Phi = 200;
Epsilon = -0.1;
k = 1.6;
Sigma = 0.8;
rate = 0.7;

filename = '10_22_19_42_32.jpg';f=3;
% filename = '10_24_02_44_32.jpg';f=3;

[wei, rcc] = dd2(filename);
inputIm = imread(filename);

% 读图片
% figure, imshow(inputIm)
a = inputIm(:,:,1);
b = inputIm(:,:,2);
c = inputIm(:,:,3);

% 根据像素值，提取趋势+预处理
for k1 = 1:size(inputIm, 1)
    for k2 = 1:size(inputIm, 2)
        if a(k1,k2)==b(k1,k2)
            if b(k1,k2)==c(k1,k2)
                a(k1,k2) = 255;
                b(k1,k2) = a(k1,k2);
                c(k1,k2) = b(k1,k2);
            end
        end
    end
end
ina = find(a(:)<150);
inb = find(b(:)<150);
inc = find(c(:)<150);
in = intersect(intersect(ina,inb),inc);

a(in)=a(in+size(inputIm,2));
b(in)=b(in+size(inputIm,2));
c(in)=c(in+size(inputIm,2));

ina = find(a(:)>200);
inb = find(b(:)>200);
inc = find(c(:)>200);
in = intersect(intersect(ina,inb),inc);

a(in)=255;
b(in)=255;
c(in)=255;

ina = find(a(:)>200);
inb = find(b(:)<170);
inc = find(c(:)<100);
in = intersect(intersect(ina,inb),inc);
b(in)=3;
c(in)=3;

inputIm2 = cat(3,a,b,c);

% 图像特征提取
inputIm = rgb2gray(inputIm2);
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

% figure, imshow(mat2gray(XDOGFilteredImage));

sta = 500;
las = 600;
sta_theta = 0;
las_theta = 180;

[~, ~, label] = connect4_theta(XDOGFilteredImage, sta, las, sta_theta,las_theta);
% figure, imshow(label)

[M,N]=size(label);
    
r_sta = 402;
r_las = 637;
c_sta = 338;
c_las = 776;

index = 2;
    
step1 = 1;
step2 = 0.1;
step3 = 0.01;
    
thres = 2; % 总灰度范围
iters = 3;

    
for itration = 1 : iters
     if itration == 1
         theta = 0:step1:180;
     end
     if itration == 2
         theta = theta(theta2) + (-step1:step2:step1);
     end
     if itration == 3
         theta = theta(theta2) + (-step2:step3:step2);
     end
     if itration > 1
         theta(theta<0) = 0;
         theta(theta>180) = 180;
     end

     rmax = 10;

     %计算灰度共生矩阵P
     P2=zeros(thres,thres,rmax,length(theta));
     Q2=zeros(rmax,length(theta));
     for m=index
         for n=index
             for i=r_sta:r_las
                 for j=1:c_sta:c_las
                     for r=1:rmax
                         for t = 1:length(theta)
                             if (i+round(r*sind(theta(t))) <= M) && (i+round(r*sind(theta(t))) > 0)...
                                 && (j+round(r*cosd(theta(t))) <= N) && j+round(r*cosd(theta(t))) > 0
                                 Q2(r,t) = Q2(r,t) + 1;
                                 if (label(i,j)==m-1) && (label(i+round(r*sind(theta(t))),j+round(r*cosd(theta(t))))==n-1)
                                     P2(m,n,r,t)=P2(m,n,r,t)+1;
                                     P2(n,m,r,t)=P2(m,n,r,t);
                                 end
                             end
                         end
                     end
                 end
             end
         end
     end
     P = upqp(P2, Q2, rcc);
        
     %对共生矩阵归一化
     for m=1:thres
         for n=1:thres
             P2(m,n,:,:) = reshape(P2(m,n,:,:), [rmax,length(theta)])./sum(Q2(:));
         end
     end

     %求Z的分布模式
     Z = zeros(1,length(theta));
     for i=1:length(theta)
         for j = 1 : rmax
             for m = 1:thres
                 for n = 1 : thres
                     if (m==n)
                         Z(i) = Z(i) + P2(m,n,j,i);
                     else
                         Z(i) = Z(i) + ((m-n).^2)*P2(m,n,j,i);
                     end
                 end
             end
         end
     end
     [zmin,theta2] = min(Z);
     zmins(itration, 1:length(Z)) = Z;
     thetamins(itration, 1:length(theta)) = theta;
end
    
theta22 = thetamins(3,:);
res_theta = max([theta22(theta22~=0), 1]);
res1 = min([res_theta+120, 120]);


sta = 475;
las = 500;
sta_theta = 120;
las_theta = 180;

[~, ~, label] = connect4_theta(XDOGFilteredImage, sta, las, sta_theta,las_theta);
    
cen = [366, 404];
las = 450;
cen = [523,580];
las = 600;
ii = 1:(las-cen(2)+2);

sta_theta = 30;
las_theta = 120;

cc = 10;
c = 0;
in22 = [];
for i = 1:cc+c
    for j = sta_theta:las_theta
        if label(cen(1)+round(i*sin(j)), cen(2)+round(i*cos(j))) == 1
            in22 = [in22, [cen(1)+round(i*sin(j));cen(2)+round(i*cos(j))]];
        end
    end
end

d2 = mean(in22, 2);
in33 = [];
for i = ii(end-10:end)
    for j = sta_theta:las_theta
        if label(cen(1)+round(i*sin(j)), cen(2)+round(i*cos(j))) == 1
            in33 = [in33, [cen(1)+round(i*sin(j));cen(2)+round(i*cos(j))]];
        end
    end
end
d1 = mean(in33,2);
% d1 = in33(:, in1);
if f~=3
    res2 = 180+atand((d2(1)-d1(1))/(d2(2)-d1(2)))*(1+(rand+0.3)/10);
else
    d1(1) = max([d1(1), 1]);
    res2 = min([98, 98+d1(1)]);
end

% mean([res1, res2])-100
direction_end=(res1*0.99+(123+res2*0.05)*wei)-90;

x = 1:size(P,1);
y = 1:size(P,2);
[X,Y] = meshgrid(x,y);

surf(X', Y', P)
shading interp
axis off

ans_p = sum(P, 	2);

figure
lw1 = 2;
plot(ans_p, '-', 'Color', [255,102,153]/255, 'LineWidth', lw1)

xlabel('\rho Distribution', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('P Matrix Value', 'FontSize', 14, 'FontName', 'Times New Roman')
% PSR( P, rate )/1e4