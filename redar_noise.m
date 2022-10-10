clear
clc
% 读取图像文件列表
% fileID = fopen('new_picture/list.txt');
fileID = fopen('纹理图像库/list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

root_name = '纹理图像库/';
% 获取具体图片名称序列
names = raw_names{1};
root_name2 = 'noise_images/';
branch_name1 = 'additive_noise/';
branch_name2 = 'productive_noise/';
% 瑞利噪声参数
a_ray = -0.2;
b_ray = 0.03;
% 椒盐噪声参数
a_pepper = 0.05;
b_pepper = 0.05;

for i = 1 : length(names)
    % 读取图像
    Image = double(imread([root_name, names{i}]));
    % 载入大小信息
    [M, N, ~] = size(Image);
    % 生成瑞利噪声
    n_rayleigh = a_ray + (-b_ray .* log(1 - rand(M,N))).^0.5;
    % 两种噪声模型
    additive_ray_Image = Image + n_rayleigh;
    imwrite(additive_ray_Image, [root_name2, branch_name1, names{i}(1:end-4), 'rayleigh.jpg'])
    productive_ray_Image = Image.*(1 + n_rayleigh);
    imwrite(productive_ray_Image, [root_name2, branch_name2, names{i}(1:end-4), 'rayleigh.jpg'])
    % 生成椒盐噪声
    x = rand(M,N);
    g_sp = Image;
    g_sp(x<=a_pepper) = 0;
    g_sp(x > a_pepper & x<(a_pepper+b_pepper)) = 1;
    % 两种噪声模型
    imwrite(g_sp, [root_name2, branch_name1, names{i}(1:end-4), 'pepper.jpg'])
    imwrite(g_sp, [root_name2, branch_name2, names{i}(1:end-4), 'pepper.jpg'])
    % 生成白噪声
    n_white = wn(M,N);
    wh_additive_Image = Image + n_white;
    imwrite(wh_additive_Image, [root_name2, branch_name1, names{i}(1:end-4), 'wn.jpg'])
    wh_productive_Image = Image .* (1 + n_white);
    imwrite(wh_productive_Image, [root_name2, branch_name2, names{i}(1:end-4), 'wn.jpg'])
end