% computes the error using sum(yfit - y)^4
% goal is to minimize the error when calculated in this way
% Parameters:
% x: [m,c,k] values passed in from fminunc
% t: time data
% f: force data
% y: displacement data
function [err] = error_obj(x, t, f, y)

    m = x(1);     % mass [kg]
    c = x(2);     % damping ratio 
    k = x(3);     % spring constant [N/m]
    wn = sqrt(k/m);
    zeta = c/(2*sqrt(m*k));
    yfit = (f./k).*(1-(1/(sqrt(1-zeta^2))).*exp(-zeta*wn.*t).*sin(wn*sqrt(1-zeta^2).*t+acos(zeta)));
    
    e = yfit - y; % error vector
    err = sum(e.^2);
end

