function res = WGSBH(eq, y, b)
% This function calculates Boltzmann’s H function
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   y is point to calculate function
%   b is vector of ballances
%
%Output:
%   res is value of function

    % y contains H2O and CO only
    c = [y(1), b(1) - y(1), y(2), b(2) - y(2), 0,  0];
    %    H2O    H2          CO    CO2          red Ox
    c(6) = b(3) - c(1) - c(3) - 2 * c(4); % red
    c(5) = b(4) - c(6);                   % Ox
    
    ind = c > 0;
    res = sum(c(ind) .* (log(c(ind) ./ eq(ind)) -1));
end