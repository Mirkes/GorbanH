function res = GH(eq, c, Gamma)
% This function calculates Gorban’s H function
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   c is point to calculate function
%   Gamma is m-by-n matrix with one vector gamma in each row
%
%Outputs:
%   res is value of function

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
