clear
clc
% 读取图像文件列表
fileID = fopen('list.txt');
raw_names = textscan(fileID, '%s');
fclose(fileID);

% 获取具体图片名称序列
names = raw_names{1};
rs = 55:575;
cs = 145:665;
root_names = 'processed_imgs/';

for i = 1 : length(names)
    img_temp = imread(names{i});
    processed_img = img_temp(rs, cs, :);
    imwrite(processed_img, [root_names, names{i}])
end