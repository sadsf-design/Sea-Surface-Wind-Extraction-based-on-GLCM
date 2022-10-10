clear
clc
fid = fopen('list.txt','r');
bb = textscan(fid,'%s');
names = bb{1};
root_name = 'F:/方框图叠加平均/';
new_root_name = 'F:/new_picture/';
for i = 1 : length(names)
    file_name = [root_name, names{i}];
    raw_data = imread(file_name);
    processed_data = raw_data(75:790, 165:970, :);
    new_name = [new_root_name, names{i}];
    imwrite(processed_data,new_name)
end