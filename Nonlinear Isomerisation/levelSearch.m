function dots = levelSearch(eq, h, Gamma)
%levelSearch searchs all points with the same level of Gorban's function
%
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   h is required level
%   Gamma is m-by-n matrix with one vector gamma in each row
%
%Output:
%   dots is K-by-n matrix of points with the GH(dots) == h

    nDot = 200;
    tmp = linspace(0.01, 1, nDot)';
    zers = zeros(nDot, 1);
    dots = bsxfun(@minus, [[tmp, 1 - tmp, zers]; [1 - tmp, zers, tmp]; [zers, tmp, 1 - tmp]], eq);
    for k = 1:nDot * 3
        chi = fminbnd(@(x) (h - GH(eq, eq + x * dots(k, :), Gamma))^2, 0, 0.9999);
        dots(k, :) = eq + chi * dots(k, :);
    end
end
