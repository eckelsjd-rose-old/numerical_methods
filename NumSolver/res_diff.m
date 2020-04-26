function ret = res_diff(x)

% computes the Jacobian of the given functions:
% R1(x1,x2) = x1^2 + x2^2 - 4;
% R2(x1,x2) = x1x2 + x2 - 1;
ret = [2*x(1) 2*x(2); x(2) x(1)+1];