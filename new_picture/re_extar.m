clear
clc
% ��ȡͼ���ļ��б�
fileID = fopen('list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

% ��ȡ����ͼƬ��������
names = raw_names{1};

% �������
theta_name = 'theta_overall';
z_name = 'Z_overall';
tail = '.mat';

res = zeros(2, length(names));

for i = 1 : length(names)
    load([theta_name, names{i}, tail])
    load([z_name, names{i}, tail])
    z_temp = zmins(3, :);
    z_temp(z_temp == 0) = [];
    [val, index] = min(z_temp);
    res(1, i) = val;
    res(2, i) = thetamins(3, index);
end

save('marked_pictures/final_res.mat', 'res')