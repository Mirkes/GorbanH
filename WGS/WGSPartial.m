function res = WGSPartial(eq, y, gamma, b)
% This function calculates partial equilibrium for one gamma
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   y is point to calculate function
%   gamma is 1-by-n vector gamma
%   b is vector of ballances
%
%Outputs:
%   res is two coordinates of partial equilibrium point

    % y contains H2O and CO only
    c = [y(1), b(1) - y(1), y(2), b(2) - y(2), 0,  0];
    %    H2O    H2          CO    CO2          red Ox
    c(6) = b(3) - c(1) - c(3) - 2 * c(4); % red
    c(5) = b(4) - c(6);                   % Ox

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

    res = y + chi * gamma([1, 3]);
end
