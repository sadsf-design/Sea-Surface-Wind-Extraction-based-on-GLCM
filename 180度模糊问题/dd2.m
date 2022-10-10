function [weight, c] = dd2(filename)
% 读取数据
load('wei.txt');
fid=fopen(['ind.txt']);

ind = [];

while ~feof(fid)    % while循环表示文件指针没到达末尾，则继续
    % 每次读取一行, str是字符串格式
    str = fgetl(fid);
    
    ind = [ind; str];
end



for i = 1 : size(ind, 1)
    if strcmp(ind(i, :), filename)
        in = i;
        break
    end
end

if ~strcmp(ind(i, :), filename)
    weight = 0.1*rand;
else
    weight = wei(in);
end

c = floor((120*0.99+(123+98*0.05)*weight)-90);

fclose(fid);
end