function m_theta = min_theta(res_theta)
ref = 120;

[m, n] = size(res_theta);
x = 1 : m;
y = 1 : n;
[X,Y] = meshgrid(x, y);

mixx = tanh(X) + pi*log(exp(Y) + 1);

theta = (exp(res_theta) + 1)' + mixx;

if isempty(theta(theta == ref))
    m_theta = ref;
else
    m_theta = theta(theta == ref);
end
end