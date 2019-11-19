function dots = levelSearchH(eq, h)
%levelSearchH searchs all points with the same level of Boltzmann's H function
%
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   h is required level
%
%Output:
%   dots is K-by-n matrix of points with the GH(dots) == h

    nDot = 200;
    tmp = linspace(0.01, 1, nDot)';
    zers = zeros(nDot, 1);
    dots = bsxfun(@minus, [[tmp, 1 - tmp, zers]; [1 - tmp, zers, tmp]; [zers, tmp, 1 - tmp]], eq);
    for k = 1:nDot * 3
        chi = fminbnd(@(x) (h - H(eq, eq + x * dots(k, :)))^2, 0, 0.9999);
        dots(k, :) = eq + chi * dots(k, :);
    end
end
