clear
clc
rect = imread('矩形.jpg');
% imshow(rect(75:790, 165:970, :))
% imwrite(rect(75:790, 165:970, :), '矩形_processed.jpg')
rect2 = rect(75:790, 165:970, :);
Image = rect2;

[M,N]=size(Image(:,:,1));
%将各颜色分量转化为灰度
Gray=rgb2gray(Image);
thres = max(Gray(:))-min(Gray(:))+1;
flag = 0;
if flag == 1
    load('P.mat')
    load('Q.mat')
else
    P = [];
    Q = [];
end
rmax = floor(sqrt(M.^2+N.^2)/2);
theta = [0,45,90,135];
%将原始图像灰度级压缩，将Gray量化成thres级
% imshow(Gray)
% Gray2=Gray;
if isempty(P) || isempty(Q)
%     for i=1:M
%         for j=1:N
%             for n=1:max(Gray(:)+1)/thres
%                 if (n-1)*thres<=Gray(i,j) && Gray(i,j)<=(n-1)*thres+thres-1
%                     Gray(i,j)=n-1;
%                 end
%             end
%         end
%     end
    % figure
    % imshow(Gray2)
    %计算四个共生矩阵P，取距离d=1，θ分别为0、45、90、135
    P=zeros(16,16,rmax,length(theta));
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
fprintf('The minimum value of z'' is %f, and the corresponding theta is %d degree.', zmin, theta(theta2))
xs = (M/2)*ones(1,M/2+1);
ys = M/2-floor(M/4):M/2+floor(M/4);
% imshow(Gray, [])
% hold on
% gray2 = bitmapplot(ys,xs,Gray);
imshow(Image)
hold on 
plot(xs,ys,'*b','LineWidth',4)
% plot(ys,xs,'*b','LineWidth',4)
hold off