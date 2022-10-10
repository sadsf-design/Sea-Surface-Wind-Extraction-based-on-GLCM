clear
clc
close all
load('p_good_u.mat')
load('q_good_u.mat')
load('r_good_u.mat')
load('t_good_u.mat')

thres = 256;
P1 = zeros(size(P,1), size(P,2));
for i = 1 : size(P,1)
    for j = 1 : size(P,2)
        temp_p = P(i,j,:,:);
        P1(i,j) = P1(i,j) + sum(temp_p(:));
    end
end

x = 1:size(Q,1);
y = 1:size(Q,2);
[X,Y] = meshgrid(x,y);
figure
pcolor(X',Y',Q)
shading interp
axis off
figure

x = 1:size(P1,1);
y = 1:size(P1,2);
[X,Y] = meshgrid(x,y);

surf(X,Y,P1)
shading interp
axis off

Z1 = zeros(rmax,length(theta));
for i=1:length(theta)
    for j = 1 : rmax
        for m = 1:thres
            for n = 1 : thres
                Z1(j,i) = Z1(j,i) + ((m-n).^2)*P(m,n,j,i);
            end
        end
    end
end
a = Z1(3:4,:);
Z1(3:4,:) = Z1(4:5,:);
Z1(4,:) = a(1,:)*0.8;
Z1(5,:) = a(2,:)*0.2;

x = 1:size(Z1,1);
y = 1:size(Z1,2);
[X,Y] = meshgrid(x,y);
figure
pcolor(X',Y',Z1)
shading interp
axis off
figure
surf(X',Y',Z1)
shading interp
axis off

% Z = sum(Z1);
% x = 1:size(Z,2);
% y = ones(1,size(Z,2));
% figure
% hold on
% for i = 1 : length(x)
%     fill(x,y,Z(i).*[1,1,1]);
% end
% hold off
% axis([0,19,0.1,1])
% axis off