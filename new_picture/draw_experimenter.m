% draw_experimenter
clear
clc
a = imread('1140£¨37+180£©.jpg');

imshow(a)
hold on
xs = ones(1, 80);
ys = 1 : 80;

ss = [xs; ys];
theta = 60;
rot = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];

sx = rot*ss;

plot(sx(1,:)+size(a,2)/2, sx(2,:)+size(a,1)/2,'*k','LineWidth',4)
hold off