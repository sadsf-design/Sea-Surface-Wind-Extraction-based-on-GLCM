clear
clc
name1 = 'theta_overall';
name2 = 'Z_overall';

% ��ȡͼ���ļ��б�
fileID = fopen('list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

% ��ȡ����ͼƬ��������
names = raw_names{1};
rs = 55:575;
cs = 145:665;
root_names = 'processed_imgs/';

% �������
res = zeros(2, length(names));

for i = 1 : length(names)
    load([name1, names{i}, '.mat'])
    load([name2, names{i}, '.mat'])
    
    theta_temp = thetamins(3, :);
    theta_temp(theta_temp == 0) = [];
    
    [m, n] = min(theta_temp);
    res(1, i) = m;
    res(2, i) = zmins(3, n);
end
save([root_names, names{i}, '.mat'], 'res')