function res = WGSGH(eq, y, Gamma, b)
% This function calculates Gorban’s H function
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   y is point to calculate function
%   Gamma is m-by-n matrix with one vector gamma in each row
%   b is vector of ballances
%
%Outputs:
%   res is value of function


    % y contains H2O and CO only
    c = [y(1), b(1) - y(1), y(2), b(2) - y(2), 0,  0];
    %    H2O    H2          CO    CO2          red Ox
    c(6) = b(3) - c(1) - c(3) - 2 * c(4); % red
    c(5) = b(4) - c(6);                   % Ox

    % For each gamma from Gamma search partial equilibrium
    gh = zeros(size(Gamma, 1), 1);
    for k = 1:size(Gamma, 1)
        gamma = Gamma(k, :);
        % Define the minimal and maximal values of t in c+t*gamma from the
        % condition of positivity of c
        % Minimal value
        ind = gamma ~= 0;
        tmp = -c(ind) ./ gamma(ind);
        ind = tmp < 0;
        ma = min(tmp(~ind))*0.99999;
        mi = max(tmp(ind))*0.99999;
        [~, gh(k)] = fminbnd(@(x) H(eq, c + x * gamma), mi, ma);
    end
    res = max(gh);
end
