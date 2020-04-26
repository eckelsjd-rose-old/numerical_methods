function ret = approxJ(x,coeff)

% used for solving a general cubic function
% x is passed in as a column vector of initial guesses [a,b,c]'
a = x(1);
b = x(2);
c = x(3);
% where a,b,c are from the factored form of the std form cubic function:
%
% std form: s^3 + As^2 + Bs + C = 0
% factored: (s-a)(s-b)(s-c)     = 0

% computes the approximate Jacobian of the residual functions:
% R1(a,b,c) = -(a+b+c) - A
% R2(a,b,c) = ab + bc + ac - B
% R3(a,b,c) = -abc - c

n = size(x,1);
ret = zeros(n);
R = @residual2;
for i = 1:n
    x_new = x;
    x_new(i) = x(i)*1.01;
    ret(:,i) = (R(x_new,coeff) - R(x,coeff)) / (0.01*x(i));
end