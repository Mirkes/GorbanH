function dots = HClLevelSearch(eq, h, Gamma, dots, b)
%levelSearch searchs all points with the same level of Gorban's function
%
%Inputs:
%   eq is 1-by-n vector of equilibrium point
%   h is required level
%   Gamma is m-by-n matrix with one vector gamma in each row
%   dots is set of points (directions) to find points with level h
%   b is vector of ballances
%
%Output:
%   dots is K-by-n matrix of points with the GH(dots) == h

    nDot = size(dots, 1);
    dots = dots - eq([1, 3, 5]);
    
    for k = 1:nDot
        chi = fminbnd(@(x) (h - HClGH(eq, eq([1, 3, 5]) + x * dots(k, :), Gamma, b))^2, 0.0000001, 0.999999);
        dots(k, :) = eq([1, 3, 5]) + chi * dots(k, :);
    end
end
