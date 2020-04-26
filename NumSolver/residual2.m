function res = residual2(x,coeff)

% used for solving a general cubic function
% x is passed in as a column vector of initial guesses [a,b,c]'
% coeff is passed in as a column vector of the 3 coefficients of a cubic
% funtion in standard form [A,B,C]'

% a = x(1);
% b = x(2);
% c = x(3);
% where a,b,c are from the factored form of the std form cubic function:
%
% std form: s^3 + As^2 + Bs + C = 0
% factored: (s-a)(s-b)(s-c)     = 0

res = zeros(3,1);
A = coeff(1);
B = coeff(2);
C = coeff(3);
res(1,1) = -(x(1) + x(2) + x(3)) - A;
res(2,1) = x(1)*x(2) + x(2)*x(3) + x(1)*x(3) - B;
res(3,1) = -x(1)*x(2)*x(3) - C;