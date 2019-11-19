function res = H(eq, c)
% This function calculates Boltzmann’s H function
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   c is point to calculate function
%
%Output:
%   res is value of function

    ind = c > 0;
    res = sum(c(ind) .* (log(c(ind) ./ eq(ind)) -1));
end