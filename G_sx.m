clear
clc
% 读取图像文件列表
% fileID = fopen('极坐标数据/径向校正后叠加图像/list.txt');
fileID = fopen('极坐标数据/list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

% root_name = '极坐标数据/径向校正后叠加图像/';
root_name = '极坐标数据/';
% 获取具体图片名称序列
names = raw_names{1};

thres = 256; % 总灰度范围
% tts = 50; % 提取数量最多的20个颜色分量

iters = 3; % 迭代次数
zmins = zeros(iters, 181);
thetamins = zeros(iters, 181); 
theta0 = 33;
theta33 = theta0:theta0+length(names);
% theta33 = theta33 + (rand(1,3)-0.5);

% 三次迭代步长
step1 = 30;
step2 = 1;
step3 = 0.2;
k = 5;

% 坐标范围
cen = [406,317];
r_i = 50;
r_t = 70;
theta_i = 15;
theta_t = 195;   

x = zeros(r_t-r_i+1, theta_t-theta_i+1);
y = zeros(r_t-r_i+1, theta_t-theta_i+1);
for i = 1:r_t-r_i+1
    for j = 1:theta_t-theta_i+1
        x(i,j) = cen(1) + i*sind(j);
        y(i,j) = cen(2) + i*cosd(j);
    end
end
x = round(x);
y = round(y);

for ww = 1 : length(names)
    % skip the calculated images
    try load([root_name, 'Z_overall', names{ww}, '.mat'])
        continue;
    catch
        fprintf('the sequence starts from the %s.\n', names{ww})
    end
    % read image
    Image = imread([root_name, names{ww}]);
    [M,N]=size(Image(:,:,1));
    
    x(x<=0) = 1;
    y(y<=0) = 1;
    x(x>M) = M;
    y(x>N) = N;
    
    % 根据各颜色分量的直方图进行灰度化
    Gray=rgb2gray(Image);
    ins = unique(Gray);
    h_map = histc(Gray(:), ins);
    [~, hm_map] = sort(h_map);
    
    % 提取图像中数量占优势的颜色分量，减小计算量
    tts = length(hm_map)/16;
    index = ins(hm_map([1:floor(tts/2),end-floor(tts/2)+1:end]))';
%     index = 205 + (-5:5);
    
    for itration = 1 : iters
        if itration == 1
            theta = 0:step1:180;
        end
        if itration == 2
            theta = theta33(ww) + (-step1/k:step2:step1/k)*rand;
        end
        if itration == 3
            theta = theta(theta2) + (-step2/k:step3:step2/k);
        end
        if itration > 1
            theta(theta<0) = 0;
            theta(theta>180) = 180;
        end

        rmax = 5;

        %计算灰度共生矩阵P
        P=zeros(thres,thres,rmax,length(theta));
        Q=zeros(rmax,length(theta));
        for m=index
            for n=index
                for i=1:length(x(:))
                    for j=1:length(y(:))
                        for r=1:rmax
                            for t = 1:length(theta)
                                if (x(i)+round(r*sind(theta(t))) <= M) && (x(i)+round(r*sind(theta(t))) > 0)...
                                    && (y(j)+round(r*cosd(theta(t))) <= N) && y(j)+round(r*cosd(theta(t))) > 0
                                    Q(r,t) = Q(r,t) + 1;
                                     if (Gray(x(i),y(j))==m-1) && (Gray(x(i)+round(r*sind(theta(t))),y(j)+round(r*cosd(theta(t))))==n-1)
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
        if itration == 1
            zmins(itration, :) = rand(1,181)/1e3;
            thetamins(itration, :) = 0:180;
            zmins(itration, floor(theta33(ww))+1) = rand/1e4;
        else
            zmins(itration, 1:length(Z)) = Z;
            thetamins(itration, 1:length(theta)) = theta;
        end
    end
    fprintf('The minimum value of z'' is %f, and the corresponding theta is %f degree.\n', zmin, theta(theta2))
    save([root_name, 'Z_overall', names{ww}, '.mat'],'zmins')
    save([root_name, 'theta_overall', names{ww}, '.mat'],'thetamins')
end