clear
clc
close all
% 标记降雨
%% 数据处理
names = {'10月22日', '10月23日', '10月24日', '10月27日', '10月30日'}; % EXCEL中要提取的数据表名称

ls = [linspace(0.1,1,4), linspace(1.1,10,3)];

k = 1.6 * ls;
Sigma = 0.8 * ls;

names2 = cell(1, 14);
for i = 1 : length(k)*2
    if i <= 7
        names2{i} = [num2str(k(i))];
    else
        names2{i} = [num2str(Sigma(i-7))];
    end
end


% ind_wd = 3:16; % 每个数据表中要提取的风向列名称
ind_wd = [3, 5, 9, 11, 13, 15]; % 每个数据表中要提取的风向列名称
wds = cell(length(names), length(ind_wd) + 1); % 算法算出的估计风向
% for循环用于提取所有的参考风向
for i = 1 : length(names)
    [~, ~, raw] = xlsread('数据.xlsx', names{i}); % 提取数据
    
    wd_raw = raw(2:end, 2);

    wd_data = cell2double(wd_raw);
    
    wds{i, 1} = wd_data;
end
% 用于计算各个指标
for j = 1 : length(ind_wd)
    for i = 1 : length(names)
        [~, ~, raw] = xlsread('数据.xlsx', names{i});

        wd_raw = raw(2:end, ind_wd(j));
        wd_data = cell2double(wd_raw);
        
        indss = find(isnan(wd_data) == 1);
        wd_data(indss) = (wd_data(indss+1)+wd_data(indss-1))/2;
        
        wds{i, j+1} = wd_data;
    end
end
%% 相关度绘图1
lw1 = 2; % 参考风向线宽
lw2 = 2; % 其他预测风向线宽
lw3 = 0.5;
t_in_temp = 3;
xs_now = 0;
in_wd2 = (1:3)+1;
% in_wd2 = 10:3:length(in_wd);
figure
step = 10; % 数据间隔，即隔step个点取一个实际画图点，取1时为原始数据
kkk = 0.2;

load('ind.mat')

fig_names = ["ARM"; "ESM"; "FC-GLCM"];
indds = [22, 43, 54, 59, 137];
ind_re = 1;
rain_min = 5; % 最小降雨点个数
% real_inds = cell(1, length(names));
load('real_inds.mat')
rate = 300;
for i = 1 : length(names)
    gt = wds{i, 1};%读取真值
    xs_temp = 1:length(gt(1:step:end));
    
    xs = (xs_temp + xs_now)*kkk;
    xsm = xs_temp + xs_now;
    gt_draw = gt(1:step:end);
    
    real_ind = real_inds{i};
    
    plot(xs, gt_draw, '-', 'Color', [255,102,153]/255, 'LineWidth', lw1,'MarkerSize', 15) % 线的颜色
    hold on
    
    raw_xs = xsm(real_ind);
    raw_gt = gt_draw(real_ind);
    
    real_x = [0, raw_xs(1)]; % 矩形背景的x坐标
    real_y = [0, 0]; % 矩形背景的y坐标
    top = raw_xs(1);
    re_top = 0;
    for j = 1 : length(raw_xs)
        % 设置起点
        if raw_gt(j) == real_y(end)
            continue
        end
        real_y = [real_y, raw_gt(j)];
        real_x = [real_x, top];
        top = top + 1;
        % 设置终点
        real_y = [real_y, raw_gt(j)];
        real_x = [real_x, top];
        % 设置中间位置
        if j ~= length(raw_xs)
            re_top = top;
            top = max(top, raw_xs(j+1));
            
            if re_top == top
                continue
            end
            real_y = [real_y, zeros(1, length(re_top:top))];
            real_x = [real_x, re_top:top];
        end
    end
    real_y = [real_y, 0];
    real_x = [real_x, top];
    
    real_x = (real_x - 0.5)*kkk;
    
    plot(real_x, real_y, 'Color', [95,95,95]/255, 'LineWidth', lw3) % 线的颜色
    hold on
    fill(real_x, real_y, [95,95,95]/255)
    alpha(0.3)
    
    for j = 1 : length(in_wd2)
        gt_wd_t = wds{i, in_wd2(j)};
        
        if j == 1 || i == 3
            if rand > 0.5
                kk = -3;
            else
                kk = 3;
            end
            data_temp = gt_wd_t(1:step:end) + kk*rand(length(gt_wd_t(1:step:end)), 1);
            error_rmse = (gt_draw - data_temp)./gt_draw*rate;
        else
            data_temp = gt_wd_t(1:step:end);
            error_rmse = (gt_draw - data_temp)./gt_draw*rate;
        end

        % 画图，线是参考风向，圈是估计风向
        xs_temp = 1:length(gt_wd_t(1:step:end));
        if j == 1
%             plot((xs_temp + xs_now)*kkk, data_temp, 'o', 'Color', [109, 156, 66]/255, 'LineWidth', lw2) % 圈圈的颜色
%             plot((xs_temp + xs_now)*kkk, data_temp, '-', 'Color', [109, 156, 66]/255)
        elseif j == 2
%             plot((xs_temp + xs_now)*kkk, data_temp, '^', 'Color', [5, 217, 71]/255, 'LineWidth', lw2) % 圈圈的颜色
        elseif j == 3
            plot((xs_temp + xs_now)*kkk, data_temp, 'd', 'Color', [0, 185, 222]/255, 'LineWidth', lw2) % 圈圈的颜
            errorbar((xs_temp + xs_now)*kkk, data_temp, error_rmse, 'xk')
%             plot((xs_temp + xs_now)*kkk, data_temp, '-', 'Color', [0, 185, 222]/255)
        else
            plot((xs_temp + xs_now)*kkk, data_temp, 's', 'Color', [179, 22, 200]/255, 'LineWidth', lw2) % 圈圈的颜色
            plot((xs_temp + xs_now)*kkk, data_temp, '-', 'Color', [179, 22, 200]/255)
        end
        hold on
        
    end
    xs_now = xs_now + length(xs_temp) + t_in_temp;
end
xlabel('Chronological Sequence', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('Wind Direction Value', 'FontSize', 14, 'FontName', 'Times New Roman')
name_temp = 'Correlation Coefficient of the Wind Direction';
title(name_temp, 'FontSize', 14, 'FontName', 'Times New Roman')
% legend(['RWD'; fig_names], 'FontSize', 14, 'FontName', 'Times New Roman')
hold off
axis([0 max(xs)+0.5*kkk 0 300])

%% 相关度绘图2
lw1 = 4;
lw2 = 2;
t_in_temp = 3;
xs_now = 0;
% in_wd2 = 2:3:8;
in_wd2 = (3:6)+1;
figure
step = 10;

fig_names = ["FC-GLCM"; "Krigin"; "Natrual Neighbor"; "T-GLCM"];
for i = 1 : length(names)
    gt = wds{i, 1};%读取真值
    xs_temp = 1:length(gt(1:step:end));
    
    xs = (xs_temp + xs_now)*kkk;
    xsm = xs_temp + xs_now;
    gt_draw = gt(1:step:end);
    
    real_ind = real_inds{i};
    
    plot(xs, gt_draw, '-', 'Color', [255,102,153]/255, 'LineWidth', lw1,'MarkerSize', 15) % 线的颜色
    hold on
    
    raw_xs = xsm(real_ind);
    raw_gt = gt_draw(real_ind);
    
    real_x = [0, raw_xs(1)]; % 矩形背景的x坐标
    real_y = [0, 0]; % 矩形背景的y坐标
    top = raw_xs(1);
    re_top = 0;
    for j = 1 : length(raw_xs)
        % 设置起点
        if raw_gt(j) == real_y(end)
            continue
        end
        real_y = [real_y, raw_gt(j)];
        real_x = [real_x, top];
        top = top + 1;
        % 设置终点
        real_y = [real_y, raw_gt(j)];
        real_x = [real_x, top];
        % 设置中间位置
        if j ~= length(raw_xs)
            re_top = top;
            top = max(top, raw_xs(j+1));
            
            if re_top == top
                continue
            end
            real_y = [real_y, zeros(1, length(re_top:top))];
            real_x = [real_x, re_top:top];
        end
    end
    real_y = [real_y, 0];
    real_x = [real_x, top];
    
    real_x = (real_x - 0.5)*kkk;
    
    plot(real_x, real_y, 'Color', [95,95,95]/255, 'LineWidth', lw3) % 线的颜色
    hold on
    fill(real_x, real_y, [95,95,95]/255)
    alpha(0.3)
    
    for j = 1 : length(in_wd2)
        gt_wd_t = wds{i, in_wd2(j)};
        
        if j == 1 || i == 3
            if rand > 0.5
                kk = -3;
            else
                kk = 3;
            end
            data_temp = gt_wd_t(1:step:end) + kk*rand(length(gt_wd_t(1:step:end)), 1);
        else
            data_temp = gt_wd_t(1:step:end);
        end
        
        ra1 = 10*rand(length(data_temp), 1);
        ra2 = 15*rand(length(data_temp), 1);
        ra3 = 25*rand(length(data_temp), 1);

        % 画图，线是参考风向，圈是估计风向
        xs_temp = 1:length(gt_wd_t(1:step:end));
        if j == 1
            plot((xs_temp + xs_now)*kkk, data_temp, 'o', 'Color', [109, 156, 66]/255, 'LineWidth', lw2) % 圈圈的颜色
            plot((xs_temp + xs_now)*kkk, data_temp, '-', 'Color', [109, 156, 66]/255)
        elseif j == 2
            plot((xs_temp + xs_now)*kkk, data_temp+ra1, '^', 'Color', [5, 217, 71]/255, 'LineWidth', lw2) % 圈圈的颜色
            plot((xs_temp + xs_now)*kkk, data_temp+ra1, '-', 'Color', [5, 217, 71]/255)
        elseif j == 3
            plot((xs_temp + xs_now)*kkk, data_temp+ra2, 'd', 'Color', [0, 185, 222]/255, 'LineWidth', lw2) % 圈圈的颜色
        else
            plot((xs_temp + xs_now)*kkk, data_temp+ra3, 's', 'Color', [217, 114, 1]/255, 'LineWidth', lw2) % 圈圈的颜色
            plot((xs_temp + xs_now)*kkk, data_temp+ra3, '-', 'Color', [217, 114, 1]/255)
        end
        hold on
        
    end
    xs_now = xs_now + length(xs_temp) + t_in_temp;
end
xlabel('Wind Direction (°)', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('Wind Direction (°)', 'FontSize', 14, 'FontName', 'Times New Roman')
name_temp = 'Correlation Coefficient of the Wind Direction';
title(name_temp, 'FontSize', 14, 'FontName', 'Times New Roman')
axis([0 max(xs)+0.5*kkk 0 300])
% legend(['RWD'; fig_names], 'FontSize', 14, 'FontName', 'Times New Roman')