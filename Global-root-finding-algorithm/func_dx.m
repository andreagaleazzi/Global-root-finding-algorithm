function y = func_dx(x)
% y = 2*cos(x).^2.*sin(x)-sin(x).^3;
% y = 2 .* x.^2;
% y = sec(x).^2;
% y = 2.*x.^3.*cos(x).*sin(x) + 3.*x.^2.*sin(x).^2;
y = cos(x) - x.*sin(x) + 2*x.^3.*cos(x).*sin(x) + 3*x.^2.*sin(x).^2;
end