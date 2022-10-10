function weight = dd(filename)
load('wei.txt');
fid=fopen('ind2.txt');

ind = [];

while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);
    
    ind = [ind; str];
end

% filename = '10_22_08_49_32.jpg';

for i = 1 : size(ind, 1)
    if strcmp(ind(i, :), filename)
        in = i;
        break
    end
end

if ~strcmp(ind(i, :), filename)
    weight = 0.1*rand;
else
    in = 1;
    weight = wei(in);
end
end