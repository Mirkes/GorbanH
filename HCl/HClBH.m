function res = HClBH(eq, y, b)
% This function calculates Boltzmann’s H function
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   y is point to calculate function
%   b is vector of ballances
%
%Outputs:
%   res is value of function


    % y contains H2, Cl2 and HCl only
    c = [y(1), b(1) - 2 * y(1) - y(3), y(2), b(2) - 2 * y(2) - y(3), y(3)];
    %    H2    H                       Cl2    Cl                     HCl

    ind = c > 0;
    res = sum(c(ind) .* (log(c(ind) ./ eq(ind)) -1));
end
