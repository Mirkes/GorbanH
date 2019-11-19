function [res, eq] = integrWGSDetail(param, coeq, tt, b)
% Calculate result of integration of ODE system during 0.59 sec. 
%
% Inputs:
%   param is vector of four elements:
%       param(1) is b(4): total number of catalyst
%       param(2) is equilibrium value of free catalyst
%       param(3) is reaction rate constants for the first inverse reaction
%       param(4) is reaction rate constants for the second inverse reaction
%   coeq is fraction of CO which is NOT coverted into CO2 in eqilibrium
%   b is vector of four elements with values of balances b_H, b_C, b_O, and b_r
%
% Outputs:
%   res is 200-by-7 matrix with time in the first column and concentration
%       of six substances in other 6 columns.
%   eq is equilibrium point.
%

    % Check the physical conditions
    if any(param <= 0)
        res = Inf;
        return
    end
    
    % Calculate equilibrium
    b(4) = param(1);
    eq = zeros(6, 1);
    eq(5) = param(2);
    eq(3) = coeq * b(2);
    eq(4) = b(2) - eq(3);
    eq(6) = b(4) - eq(5);
    eq(2) = eq(4) + eq(6);
    eq(1) = b(1) - eq(2);
    
    % Check the physical conditions
    if any(eq <= 0)
        res = Inf;
        return
    end
    % Calcu;ate reaction rates constants
    kM = param(3:4);
    kP = [kM(1) * eq(2) * eq(6) / (eq(1) * eq(5)),...
        kM(2) * eq(4) * eq(5) / (eq(3) * eq(6))];
    
    % Initial point
    c0 = [0.9999 * b(1), 0.9999 * b(2)];
    
    % Parameters
    opts = odeset('Reltol',1e-13,'AbsTol',1e-14); %,'Stats','on');
    
    % Integrate
    [t, c] = ode113(@(ttt, y) WGSODE(ttt, y, kP, kM, b), linspace(0, tt, 590), c0, opts);
    
    % Restore all substances and form output
    c6 = (b(3) - 2 * b(2)) - c(:, 1) + c(:, 2);
    res = [t, c(:, 1), b(1) - c(:, 1), c(:, 2), b(2) - c(:, 2), b(4) - c6, c6];
end