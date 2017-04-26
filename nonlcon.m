function [c, ceq] = nonlcon(x)
a = 0.2;
b = 20;
lambda_t = 2*pi/3;
c = ones(6,1).*a*exp(-b*(x(1)-lambda_t)^2)-x(5);
ceq = [];
