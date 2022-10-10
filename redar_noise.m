clear
clc
% ��ȡͼ���ļ��б�
% fileID = fopen('new_picture/list.txt');
fileID = fopen('����ͼ���/list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

root_name = '����ͼ���/';
% ��ȡ����ͼƬ��������
names = raw_names{1};
root_name2 = 'noise_images/';
branch_name1 = 'additive_noise/';
branch_name2 = 'productive_noise/';
% ������������
a_ray = -0.2;
b_ray = 0.03;
% ������������
a_pepper = 0.05;
b_pepper = 0.05;

for i = 1 : length(names)
    % ��ȡͼ��
    Image = double(imread([root_name, names{i}]));
    % �����С��Ϣ
    [M, N, ~] = size(Image);
    % ������������
    n_rayleigh = a_ray + (-b_ray .* log(1 - rand(M,N))).^0.5;
    % ��������ģ��
    additive_ray_Image = Image + n_rayleigh;
    imwrite(additive_ray_Image, [root_name2, branch_name1, names{i}(1:end-4), 'rayleigh.jpg'])
    productive_ray_Image = Image.*(1 + n_rayleigh);
    imwrite(productive_ray_Image, [root_name2, branch_name2, names{i}(1:end-4), 'rayleigh.jpg'])
    % ���ɽ�������
    x = rand(M,N);
    g_sp = Image;
    g_sp(x<=a_pepper) = 0;
    g_sp(x > a_pepper & x<(a_pepper+b_pepper)) = 1;
    % ��������ģ��
    imwrite(g_sp, [root_name2, branch_name1, names{i}(1:end-4), 'pepper.jpg'])
    imwrite(g_sp, [root_name2, branch_name2, names{i}(1:end-4), 'pepper.jpg'])
    % ���ɰ�����
    n_white = wn(M,N);
    wh_additive_Image = Image + n_white;
    imwrite(wh_additive_Image, [root_name2, branch_name1, names{i}(1:end-4), 'wn.jpg'])
    wh_productive_Image = Image .* (1 + n_white);
    imwrite(wh_productive_Image, [root_name2, branch_name2, names{i}(1:end-4), 'wn.jpg'])
end