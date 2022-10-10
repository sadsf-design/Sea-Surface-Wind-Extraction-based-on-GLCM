clear
clc
close all

a = 0;    %高斯噪声均值
b = 0.08; %高斯噪声标准差

num = 200;

t = 1:200;
t2 = 10 : 209;
t3 = -99:100;

gaussian = a + b .* randn(1, num);
gaussian = sort(gaussian)+sort(gaussian).*0.2.*randn(1,200);




a_6 = sort(randraw('tukeylambda', -0.06, [1, 200]));
% figure, plot(((a_6+8)/25.4545)*33, 'bo', 'MarkerFaceColor', 'b')
figure, plot(t, (gaussian+0.3)*35, 'bo', 'MarkerFaceColor', 'b')
hold on
% p1 = 0.0334;
% p2 = -3.514;
% plot(t, ((p1*t+p2+8)/25.4545)*33, 'r-')
a1 = 0.017;
b1 = 0.022;
c1 = -0.026;
d1 = -0.022;
y = a1*exp(b1*t3) + c1*exp(d1*t3);
plot(t3, y, 'r-')
hold off

a_26 = sort(randraw('tukeylambda', -0.26, [1, 200]));
figure, plot(a_26, 'bo', 'MarkerFaceColor', 'b')
p1 = 0.04336;
p2 = -4.65;
hold on
plot(t, p1*t+p2, 'r-')
hold off