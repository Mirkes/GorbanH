% HCl system study. Data preparation
% 
% Since modified system is not rigid (see kp below) integration requires
% less than a second. 

% Define parameters
% Equilibrium
eq = [0.2, 0.2, 0.25, 0.1, 0.4];

% Initial state for three concentrations (H2, Cl2, HCl)
c0 = [0.4998, 0.4999, 0.0001];

% Balances
b = [1, 1];

modified = 3;
% Reaction rate constants for direct reactions
if modified == 1
    kp = [5, 10, 2, 1];
    suf = 'Mod';
    tt = 15;
    ts = 1;
    tss = 0.02;
elseif modified == 2
    kp = [1, 2, 5, 10];
    suf = 'Mod2';
    tt = 15;
    ts = 1;
    tss = 0.02;
elseif modified == 3
    kp = [0, 1, 1, 1];
    suf = 'Mod3';
    tt = 40;
    ts = 4;
    tss = 0.5;
end

% Reverse reaction rate calculation
km = [kp(1) * eq(1) / (eq(2) ^ 2),...
      kp(2) * eq(3) / (eq(4) ^ 2),...
      kp(3) * eq(2) * eq(3) / (eq(4) * eq(5)),...
      kp(4) * eq(1) * eq(4) / (eq(2) * eq(5))];

% Integration accuracy (non rigid system can be integrated with high
% accuracy
opts = odeset('Reltol',1e-13,'AbsTol',1e-14);

% Time to stop

% Inital fragment with small time steps
tic;
[t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b),...
    linspace(0, ts, 10001), c0, opts); %#ok<ASGLU>

fprintf('Initial fragment is completed');
toc
save([suf, 'TimStart.mat'], 'c', 't');

% Shorter inital fragment with small time steps
tic;
[t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b),...
    linspace(0, tss, 1001), c0, opts); %#ok<ASGLU>

fprintf('Initial fragment is completed');
toc
save([suf, 'TimStartStart.mat'], 'c', 't');

% Integrate
tic;
[t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b), linspace(0, tt, 1001), c0, opts);

fprintf('Full trajectory with final time %g is completed', tt);
toc

% Save current result
save([suf, 'Tim.mat'], 'c', 't');
