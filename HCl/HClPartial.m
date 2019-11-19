function [resX, resY, resZ] = HClPartial(eq, y, gamma, b)
% This function calculates partial equilibrium for one gamma
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   y is point to calculate function
%   gamma is 1-by-n vector gamma
%   b is vector of ballances
%
%Outputs:
%   resX, resY, resZ is tree coordinates of partial equilibrium point


    % y contains H2, Cl2 and HCl only
    c = [y(1), b(1) - 2 * y(1) - y(3), y(2), b(2) - 2 * y(2) - y(3), y(3)];
    %    H2    H                       Cl2    Cl                     HCl
    
    % Correction for small negatives
    c(c < 0) = 0;
    
    % Search partial equilibrium
    % Define the minimal and maximal values of t in c+t*gamma from the
    % condition of positivity of c
    ind = gamma < 0;
    tmp = - c(ind) ./ gamma(ind);
    ma = min(tmp);
    ind = gamma > 0;
    tmp = c(ind) ./ gamma(ind);
    mi = -min(tmp);

    options = optimset('TolX', 1.e-7);

    chi = fminbnd(@(x) H(eq, c + x * gamma), mi, ma, options);

    resX = y(1) + chi * gamma(1);
    resY = y(2) + chi * gamma(3);
    resZ = y(3) + chi * gamma(5);
end
