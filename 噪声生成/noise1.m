%% 读取原图片
A = double(im2bw(imread('e1.png')));
figure
histogram(A)
hold off
% histogram(A)
[M, N] = size(A);
figure
imshow(A,[])
title('Original Image')
hold off
%% 高斯噪声
a = 0;    %高斯噪声均值
b = 0.08; %高斯噪声标准差
n_gaussian = a + b .* randn(M,N);
% 上述提到的a和b所代表的意义可以通过下面的两行代码进行验证。
% mean(n_gaussian(:))
% std(n_gaussian(:))
% gau_A = A + n_gaussian;
gau_A = A.*(1 + n_gaussian);
figure
histogram(gau_A)
hold off
figure
imshow(gau_A)
title('Image with Gaussian Noise')
hold off
%% 瑞利噪声
a = -0.2;
b = 0.03;
n_rayleigh = a + (-b .* log(1 - rand(M,N))).^0.5;
% ray_A = A + n_rayleigh;
ray_A = A.*(1 + n_rayleigh);
figure
histogram(ray_A)
hold off
figure
imshow(ray_A, [])
title('Image with Raly Noise')
hold off
%% 伽马噪声
a = 25;
b = 3;
n_Erlang = zeros(M,N); 
for j=1:b
    n_Erlang = n_Erlang + (-1/a)*log(1 - rand(M,N));
end
% erl_A = A + n_Erlang;
erl_A = A .* (1 + n_Erlang);
figure
histogram(erl_A)
hold off
figure
imshow(erl_A)
title('Image with Gamma Noise')
hold off
%% 均匀噪声
a = 0;
b = 0.3;
n_Uniform = a + (b-a)*rand(M,N);
% uni_A = A + n_Uniform;
uni_A = A .* (1 + n_Uniform);
histogram(uni_A)
hold off
figure
imshow(uni_A)
title('Image with Average Noise')
hold off
%% 椒盐噪声
a = 0.05;
b = 0.05;
x = rand(M,N);

% g_sp = zeros(M,N);
g_sp = A;

g_sp(x<=a) = 0;
g_sp(x > a & x<(a+b)) = 1;

% g_A = A + g_sp;
figure
histogram(g_sp)
hold off
figure
imshow(g_sp, [])
title('Image with Pepper Noise')
hold off
% res = wn(M,N);
% imshow(A+res, [])
%% 白噪声
n_white = wn(M,N);
% wh_A = A + n_white;
wh_A = A .* (1 + n_white);
figure
imshow(wh_A, []);
title('Image with white Noise')
hold off
figure
histogram(n_white)
hold off