clear
clc
close all

% 读取数据
raw_p_xjd = load('P_Xjd.mat');
raw_p_djd = load('P_djd.mat');
load('P_ans_good2.mat')

% 真正读取数据
p_xjd = raw_p_xjd.P;
p_djd = raw_p_djd.P;

lw1 = 2;
% 画phi的图像
figure
plot(sum(p_xjd, 1), '-', 'Color', [255, 102,153]/255, 'LineWidth', lw1)
hold on
plot(ans_p/4, '--', 'Color', [74, 182, 230]/255, 'LineWidth', lw1)

xlabel('Typical Huge Wind Direction Distribution', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('P Matrix Value', 'FontSize', 14, 'FontName', 'Times New Roman')
hold off
draw_h = ans_p/4;
axis([0 360 0 max(draw_h(:))+2e3])

% 画rho的图像
figure
plot(sum(p_xjd, 2), '-', 'Color', [255, 102,153]/255, 'LineWidth', lw1)
hold on
plot(sum(p_xjd, 1), '--', 'Color', [74, 182, 230]/255, 'LineWidth', lw1)

xlabel('Typical Trival Wind Direction Distribution', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('P Matrix Value', 'FontSize', 14, 'FontName', 'Times New Roman')
hold off
draw_h = sum(p_xjd, 2);
axis([0 360 0 max(draw_h(:))+2e3])