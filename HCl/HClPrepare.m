% HCl system study. Data preparation
% 
% Since system is very rigid (see kp below) integration can require a lot
% of time.
% Calculated data are presented in repository and required following time:
%
%                                               Calculation time
% Stage	Time from	Time to     Time step       seconds hours
% 1     0           1.e-10      1.e-13          16109    4.47
% 2     1.e-10      1.5e-10     5.e-12           8981    2.49
% 3     1.5e-10     2.25e-10	7.5e-12         11905    3.31
% 4     2.25e-10	3.375e-10	1.125e-11       17652    4.90
% 5     3.375e-10	5.0625e-10	1.6875e-11      25882    7.19
% 6     5.0625e-10	7.5938e-10	2.53125e-11     42293   11.75
% 7     7.5938e-10	1.1391e-9	3.796875e-11	55770   15.49

% Define parameters
% Maximal number of calculation stages (see table above)
maxStages = 7;

% Equilibrium
eq = [0.198, 0.004, 0.1995, 0.001, 0.6];

% Initial state for three concentrations (H2, Cl2, HCl)
c0 = [0.4998, 0.4999, 0.0001];

% Balances
b = [1, 1];

% Reaction rate constants for direct reactions
kp = [1.e16, 1.e16, 1.7e11, 1.59e8];

% Reverse reaction rate calculation
km = [kp(1) * eq(1) / (eq(2) ^ 2),...
      kp(2) * eq(3) / (eq(4) ^ 2),...
      kp(3) * eq(2) * eq(3) / (eq(4) * eq(5)),...
      kp(4) * eq(1) * eq(4) / (eq(2) * eq(5))];

% Form very detailed initial fragment for figures with switches. Reuired
% time is less than a second.
if  ~exist('timStart.mat','file')
    % Integration accuracy (short interval can be calculated with high
    % accuracy
    opts = odeset('Reltol',1e-13,'AbsTol',1e-14);

    % Inital fragment with small time steps
    tic;
    [t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b),...
        linspace(0, 1.e-19, 1001), c0, opts);
    
    fprintf('Initial fragment is completed');
    toc
    save('timStart.mat', 'c', 't');
end

% Integration accuracy (we reduce accuracy to decrease time)
opts = odeset('Reltol',1e-8,'AbsTol',1e-8);

% Time to stop
tt =1.e-10;


% Search of current state
n = 1;
while exist(['tim',num2str(n),'.mat'],'file')
    fprintf('Step %d is completed \n', n);
    n = n + 1;
end

if n == 1
    % Integrate
    tic;
    [t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b), linspace(0, tt, 1001), c0, opts);
    
    fprintf('Initial fragment with final time %g is completed', tt);
    toc
    
    % Save current result
    save(['tim',num2str(n),'.mat'], 'c', 't');
    n = n + 1;
else
    load(['tim',num2str(n - 1),'.mat']);
end

if n >= maxStages
    return
end

% Main cycle
while true
    tt = t(end);
    c0 = c(end, :);

    % Integrate
    tic;
    [t, c] = ode113(@(ttt, y) HClODE(ttt, y, kp, km, b), linspace(tt, tt * 1.5, 11), c0, opts);
    
    fprintf('Step %d with start time %g is completed', n, tt);
    toc
    
    % Save current result
    save(['tim',num2str(n),'.mat'], 'c', 't');
    n = n + 1;
end
