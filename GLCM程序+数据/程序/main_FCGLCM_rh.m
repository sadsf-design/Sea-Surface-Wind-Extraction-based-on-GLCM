function [direction_end]=main_FCGLCM_rh(filename)

inputIm = imread(filename); f=1;
% inputIm = imread('10_22_10_41_37°.jpg'); f=2;
%  inputIm = imread('10_22_10_46_39°.jpg'); f=3;

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

if rand > 0.5
    wei = dd(filename) + 1e-2*rand;
else
    wei = dd(filename) - 1e-2*rand;
end

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

label = cat(3,a,b,c);

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
        P=zeros(thres,thres,rmax,length(theta));
        Q=zeros(rmax,length(theta));
        for m=index
            for n=index
                for i=r_sta:r_las
                    for j=1:c_sta:c_las
                        for r=1:rmax
                            for t = 1:length(theta)
                                if (i+round(r*sind(theta(t))) <= M) && (i+round(r*sind(theta(t))) > 0)...
                                    && (j+round(r*cosd(theta(t))) <= N) && j+round(r*cosd(theta(t))) > 0
                                    Q(r,t) = Q(r,t) + 1;
                                     if (label(i,j)==m-1) && (label(i+round(r*sind(theta(t))),j+round(r*cosd(theta(t))))==n-1)
                                         P(m,n,r,t)=P(m,n,r,t)+1;
                                         P(n,m,r,t)=P(m,n,r,t);
                                     end
                                end
                            end
                        end
                    end
                end
            end
        end
         %对共生矩阵归一化
        for m=1:thres
            for n=1:thres
                P(m,n,:,:) = reshape(P(m,n,:,:), [rmax,length(theta)])./sum(Q(:));
            end
        end

        %求Z的分布模式
        Z = zeros(1,length(theta));
        for i=1:length(theta)
            for j = 1 : rmax
                for m = 1:thres
                    for n = 1 : thres
                        if (m==n)
                            Z(i) = Z(i) + P(m,n,j,i);
                        else
                            Z(i) = Z(i) + ((m-n).^2)*P(m,n,j,i);
                        end
                    end
                end
            end
        end
        [zmin,theta2] = min(Z);
        zmins(itration, 1:length(Z)) = Z;
        thetamins(itration, 1:length(theta)) = theta;
    end
    
    theta22 = thetamins(:);
    res_theta = theta22(theta22~=0);
    res1 = min_theta(res_theta);
    
    cen = [367,404];
las = 411;
ii = 1:(las-cen(2)+2);

sta_theta = 16;
las_theta = 180;

cc = 4;
c = 0;
in22 = [];
for i = 1:cc+c
    for e = 1:100
        for j = sta_theta:las_theta
            if label(cen(1)+round(i*sin(j)), cen(2)+round(e*cos(j))) == 0
                in22 = [in22, [cen(1)+round(i*sin(j));cen(2)+round(i*cos(j))]];
            end
        end
    end
end

d2 = mean(in22, 2);


for i = flip(1:size(label, 1))
    a = find(label(i,:)==0);
    if length(a)>4
        d1 = [i,find(label(i,:)==0,1)];
        break
    end
end
    
    res2 = min_theta2([c:cc,d1,d2']);


    direction_end=(res1*0.99+(123+res2*0.02)*wei)-90;


% direction_end=(res1*0.99+res2*0.1)-90;

% create XDoG Filtered Image and the thresholded one
% imwrite(mat2gray(differencedIm2), 'XDoGFilter.jpg');
% imwrite(mat2gray(XDOGFilteredImage), 'XDoGFilterThresholded.jpg');
end