clear
clc
% ��ȡͼ���ļ��б�
fileID = fopen('list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

% ��ȡ����ͼƬ��������
names = raw_names{1};
rs = 55:575;
cs = 145:665;
root_names = 'processed_imgs/';

for i = 1 : length(names)
    img_temp = imread(names{i});
    processed_img = img_temp(rs, cs, :);
    imwrite(processed_img, [root_names, names{i}])
end