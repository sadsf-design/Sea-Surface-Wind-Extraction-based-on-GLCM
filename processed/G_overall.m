% 基于GLCM纹理特征提取，d=1,θ=0°、45°、90°、135°共四个矩阵
% 所用图像灰度级均为256
clear
clc
% 图像路径
root_name = 'F:/new_picture/';
fid = fopen('list.txt','r');
bb = textscan(fid,'%s');
names = bb{1};
% 读取图像
for s = 1 : length(names)
    try load([root_name, names{s}(1:end-4), '.mat'])
        if ~isempty(load([root_name, names{s}(1:end-4), '.mat']))
            continue;
        end
    catch
        disp(['Start from image ',names{s}])
    end
    file_name = [root_name, names{s}];
    Image = imread(file_name);
    [M,N]=size(Image(:,:,1));
    Gray=rgb2gray(Image);
    thres = 16;
    rmax = floor(sqrt(M.^2+N.^2)/2);
    theta = [0,45,90,135];
    %将原始图像灰度级压缩，将Gray量化成16级
        for i=1:M
            for j=1:N
                for n=1:max(Gray(:)+1)/thres
                    if (n-1)*thres<=Gray(i,j) && Gray(i,j)<=(n-1)*thres+thres-1
                        Gray(i,j)=n-1;
                    end
                end
            end
        end
        % figure
        % imshow(Gray2)
        %计算四个共生矩阵P，取距离d=1，θ分别为0、45、90、135
        P=zeros(thres,thres,rmax,length(theta));
        Q=zeros(rmax,length(theta));
        for m=1:thres
            for n=1:thres
                for i=1:M
                    for j=1:N
                        for r=1:rmax
                            for t = 1:length(theta)
                                if (i+round(r*sind(theta(t))) <= M) && (i+round(r*sind(theta(t))) > 0)...
                                        && (j+round(r*cosd(theta(t))) <= N) && j+round(r*cosd(theta(t))) > 0
                                    Q(r,t) = Q(r,t) + 1;
                                    if (Gray(i,j)==m-1) && (Gray(i+round(r*sind(theta(t))),j+round(r*cosd(theta(t))))==n-1)
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
                    Z(i) = Z(i) + ((m-n).^2)*P(m,n,j,i);
                end
            end
        end
    end
    [zmin,theta2] = min(Z);
    fprintf('The minimum value of z'' in %s is %f, and the corresponding theta is %d degree.', names{s}, zmin, theta(theta2))
    save([root_name, names{s}(1:end-4), '.mat'], 'Z')
end