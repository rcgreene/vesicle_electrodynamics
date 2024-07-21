function x = lobatto_nodes(n)

x = sin(linspace(-pi/2, pi/2, n));
x(1) = -1;
x(end) = 1;
x_resid = zeros(1, n - 2);
x_deriv = zeros(1, n - 2);
tmp = zeros(1, n - 2);

for k = 1:40
    x_resid(1:end) = 1;
    tmp(1:end) = x(2:end - 1);
    for i = 1:n-1
        x_deriv(1:end) = ((2*i + 1)*x(2:end - 1).*tmp - i*x_resid)/(i + 1);
        x_resid(1:end) = tmp;
        tmp(1:end) = x_deriv;
    end
    x_resid(1:end) = x_resid - ((2*n  + 1)*x(2:end - 1).*tmp - n*x_resid)/(n + 1);
    x_resid(1:end) = -x_resid/(2*n + 1);
    if k < 25
        x(2:end - 1) = x(2:end - 1) - .1*x_resid./x_deriv;
    else
        x(2:end - 1) = x(2:end - 1) - x_resid./x_deriv;
    end
end

end