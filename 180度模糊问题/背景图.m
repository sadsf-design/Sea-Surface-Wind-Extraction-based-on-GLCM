clear
clc 
close all

x = 0:0.01:3;
t = 4500 - 500*x.^2;
xs = 1:length(x);

lw1 = 2;

y = t+50*randn(1, length(t));

plot(xs, y, '-', 'Color', [255,102,153]/255, 'LineWidth', lw1)
fill([0, 0, xs, xs(end)], [0, y(1), y, 0], [132, 132, 132]/255)

xlabel('Pixel Distribution', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('Pixel Value of the Image', 'FontSize', 14, 'FontName', 'Times New Roman')

axis([0, 300, 0, 4500])

figure
x1 = 4200;
x2 = 3200;
x3 = 2800;
ts = [];
for i = 1 : 3
    k = (x1 - x2)/50;
    
    x = 1:0.1:6;
    t = x1 - k*x.^2;
    t = t + 50*randn(1, length(t));
    
    ts = [ts, t];
    
    k = (x2 - x3)/4;
    
    x = 1:0.1:2;
    t = x2 - k*x.^2;
    t = t + 50*randn(1, length(t));
    ts = [ts, t];
    
    if i < 4
        x1 = x3;
        x2 = x1 - 1000;
        x3 = x2 - 400;
    end
end

plot(1:length(ts), ts, '-', 'Color', [255,102,153]/255, 'LineWidth', lw1)
fill([0, 0, 1:length(ts), length(ts)], [0, ts(1), ts, 0], [132, 132, 132]/255)
axis([0, length(ts), 0, 4200])

xlabel('Pixel Distribution', 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel('Pixel Value of the Image', 'FontSize', 14, 'FontName', 'Times New Roman')