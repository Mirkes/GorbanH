% Model reaction figures
% To form figures for randomly generated equilibria it is necessary to
% set required number of equilibrias into variable nRep in line 29 and
% remove file EqSel.mat from the work directory.

% What to draw
drawDirections = true;         % Draw triangle and projections
drawGorbanLevels = true;       % Draw triangle with Gorban's H level sets
drawBoltzmannLevels = true;    % Draw triangle with Boltzmann's H level sets
drawTrajectory = true;         % Draw triangle with trajecory and graphs of 
                               % Gorban's and Bolzmann's H 


% Define name of folder to save figures
dirName = 'FiguresDot/';

% Figure format
drive = '-dpng'; % For png images
% drive = '-depsc'; % For esp images
% drive = '-dpdf'; % For pdf format

% Define dot to illustrate the calculation of Gorban's H function
c0 = [0.6, 0.4, 0.06];
c0 = c0 / sum(c0);

% Form Gamma with stoichiometric vectors only
Gamma = [-1,  1,  0;
          0, -1,  1; 
          1,  0, -1];
nReac = size(Gamma, 1);

if  ~exist('EqSel.mat','file')
    % Number of points 0 for my selection, positive integer for random
    % generation
    nRep = 0;
    if nRep == 0
        % My selected equilibria
        eqs = [1/3, 1/3, 1/3; 0.13, 0.29, 0.58; 0.36, 0.07, 0.57];
    else
        % Randomly generated equilibriums
        eqs = rand(nRep, 3);
        % Renormalise to the unit sum
        for k = 1:nRep
            eqs(k, :) = eqs(k, :) / sum(eqs(k, :));
        end
    end
    save('EqSel', 'eqs');
else
    % Load previously prepared equilibria
    load('EqSel');
end
nRep = size(eqs, 1);

if drawDirections
    % Draw Triangle with dot and directions
    for k = 1:nRep
        DrawTriangle(eqs(k, :), c0);
        saveFigures(sprintf('%sMod2Direct%03d', dirName, k), drive);
        close();
    end
end

if drawGorbanLevels
    % Draw triangles with Gorban's H levels
    nLev = 10;
    res = zeros(600, 3, nLev);
    for kk = 1:nRep
        % Identify the furthest (with greatest GH value) vertix
        % Test values in the almost vertices
        gh = [GH(eqs(kk, :), [0.99990, 0.00005, 0.00005], Gamma),...
              GH(eqs(kk, :), [0.00005, 0.99990, 0.00005], Gamma),...
              GH(eqs(kk, :), [0.00005, 0.00005, 0.99990], Gamma)];
        [~, ind] = max(gh);
        ind = [1, 2, 3] == ind;
        base = linspace(eqs(kk, ind), 0.9999, nLev + 1);
        base = base(2:end);
        dot = zeros(nLev, 3);
        dot(:, ind) = base;
        dot(:, ~ind) = repmat(eqs(kk, ~ind) / sum(eqs(kk, ~ind)), nLev, 1) .* repmat((1 - base)',1, 2);
        for k = 1:nLev
            res(:, :, k) = levelSearch(eqs(kk, :), GH(eqs(kk, :), dot(k, :), Gamma), Gamma);
        end
        DrawTriangleLevels(eqs(kk, :), res);
        saveFigures(sprintf('%sMod2Levels%03d', dirName, kk), drive);
        close();
    end
end

if drawBoltzmannLevels
    % Draw triangles with Boltzmann's H levels
    nLev = 10;
    res = zeros(600, 3, nLev);
    for kk = 1:nRep
        % Identify the furthest (with greatest GH value) vertix
        % Test values in the almost vertices
        gh = [H(eqs(kk, :), [0.99990, 0.00005, 0.00005]),...
              H(eqs(kk, :), [0.00005, 0.99990, 0.00005]),...
              H(eqs(kk, :), [0.00005, 0.00005, 0.99990])];
        [~, ind] = max(gh);
        ind = [1, 2, 3] == ind;
        base = linspace(eqs(kk, ind), 0.9999, nLev + 1);
        base = base(2:end);
        dot = zeros(nLev, 3);
        dot(:, ind) = base;
        dot(:, ~ind) = repmat(eqs(kk, ~ind) / sum(eqs(kk, ~ind)), nLev, 1) .* repmat((1 - base)',1, 2);
        for k = 1:nLev
            res(:, :, k) = levelSearchH(eqs(kk, :), H(eqs(kk, :), dot(k, :)));
        end
        DrawTriangleLevels(eqs(kk, :), res);
        saveFigures(sprintf('%sMod2LevelsBH%03d', dirName, kk), drive);
        close();
    end
end

% Specify font size for figures
fontSize = 20;

% Define origin of trajectory
c0 = [0.9999, 0.00005, 0.00005];
c0 = c0 / sum(c0);

% Define reaction rate constants of direct reactions: 2 sets of constants
% for each equilibrium.
tmp = 1 / 3 - 0.001;
kPs = [0.1, 0.2, 0.3;  tmp,  tmp, tmp;
       0.5, 0.6, 0.1;  0.5,  0.6, 0.1;
       0.2, 0.5, 0.1;  0.05, 0.1, 0.853];
% Specify values of reaction rate constant for the thirs reverse reaction:
% 2 values of constant for each equilibrium. 
kMs = [0.6, 0.001, 1.1, 10, 0.5, 2];
       
% Colours
cols = ['r', 'm', 'b'];

% Stop time for each case
stopTime = [10, 12, 4, 1, 8, 2];
stopTimeD = [10, 10, 5, 6, 12, 4];

% Convergence level to select most interesting initial part
cLevel = -0.99;

% Number of reaction reate set
rr = 0;
% Equilibrium number
for kE = 1:nRep
    eq = eqs(kE, :);
    % Reaction rate set number
    for kS = 1:2
        % Actual reaction rate constant set number
        rr = rr + 1;

        % Complex balance
        % Get reaction rate constants for direct reactions
        kP = kPs(rr, :);
        % Get reaction rate constant for the third revers reaction
        kM(3) = kMs(rr);
        % calculate two other reaction rate constants of reverse reactions
        kM(1) = (kM(3) * eq(1) + kP(1) * eq(1) - kP(3) * eq(3)) / eq(2);
        kM(2) = (kM(3) * eq(1) + kP(2) * eq(2) - kP(3) * eq(3)) / eq(3);
        
        % Parameters
        opts = odeset('Reltol',1e-13,'AbsTol',1e-14); %,'Stats','on');
        
        % Integrate
        [t, c] = ode113(@(tt, y) modelODE(tt, y, kP, kM), linspace(0, stopTime(rr), 100), c0([1, 2]), opts);
        
        % Calculate complete set of coordinates
        c = [c(:, 1), c(:, 2), 1 - c(:, 1) - c(:, 2)];
        nDots = size(c, 1);
        
        % Draw trajectory
        DrawTriangleLevels(eq, c, true);
        saveFigures(sprintf('%sMod2Traject%d_%d', dirName, kE, kS), drive);
        close();
        
        % Calculate H*gamma.
        % Memory allocation
        HGam = zeros(nDots, nReac);
        HB = zeros(nDots, 1);
        for k = 1:nDots
            HGam(k, :) = HG(eq, c(k, :), Gamma);
            % Calculate Boltzmann's H
            HB(k) = H(eq, c(k, :));
        end
        HGo = max(HGam, [], 2);
        
        % Form graphs
        figure;
        for k = 1:3
            plot(t, HGam(:, k), cols(k), 'Linewidth', 1.5);
            hold on;
        end
        set(gca, "FontSize", fontSize);
        plot(t, HGo, ':k', 'Linewidth', 4);
        xlim([0, t(end)]);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(A_3\leftrightarrow A_1)', '\it H_{\Gamma}');
        saveFigures(sprintf('%sMod2TrajectH%d_%d', dirName, kE, kS), drive);
        close();
        
        % Conv in this case is half of total number of points
        conv = 50;

        % Form graphs
        figure;
        for k = 1:3
            plot(t(1:conv), HGam(1:conv, k), cols(k), 'Linewidth', 1.5);
            hold on;
        end
        set(gca, "FontSize", fontSize);
        plot(t(1:conv), HGo(1:conv), ':k', 'Linewidth', 4);
        xlim([0, t(conv + 1)]);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(A_3\leftrightarrow A_1)', '\it H_{\Gamma}');
        saveFigures(sprintf('%sMod2TrajectHT%d_%d', dirName, kE, kS), drive);
        close();

        figure;
        plot(t(1:conv), HGo(1:conv), 'r', 'Linewidth', 1.5);
        hold on;
        plot(t(1:conv), HB(1:conv), 'b', 'Linewidth', 1.5);
        xlim([0, t(conv + 1)]);
        set(gca, "FontSize", fontSize);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\Gamma}', '\it H');
        saveFigures(sprintf('%sMod2TrajectHTB%d_%d', dirName, kE, kS), drive);
        close();

        % Detailed balance
        % Get reaction rate constants for direct reactions
        kP = kPs(rr, :);
        % Get reaction rate constant for the third revers reaction
        % Calculate all reaction rate constants of reverse reactions
        kM(1) = kP(1) * eq(1) / eq(2);
        kM(2) = kP(2) * eq(2) / eq(3);
        kM(3) = kP(3) * eq(3) / eq(1);
        
        % Parameters
        opts = odeset('Reltol',1e-13,'AbsTol',1e-14); %,'Stats','on');
        
        % Integrate
        [t, c] = ode113(@(tt, y) modelODE(tt, y, kP, kM), linspace(0, stopTimeD(rr), 100), c0([1, 2]), opts);
        
        % Calculate complete set of coordinates
        c = [c(:, 1), c(:, 2), 1 - c(:, 1) - c(:, 2)];
        nDots = size(c, 1);
        
        % Draw trajectory
        DrawTriangleLevels(eq, c, true);
        saveFigures(sprintf('%sMod2TrajectDet%d_%d', dirName, kE, kS), drive);
        close();
        
        % Calculate H*gamma.
        % Memory allocation
        HGam = zeros(nDots, nReac);
        HB = zeros(nDots, 1);
        for k = 1:nDots
            HGam(k, :) = HG(eq, c(k, :), Gamma);
            % Calculate Boltzmann's H
            HB(k) = H(eq, c(k, :));
        end
        HGo = max(HGam, [], 2);
        
        % Form graphs
        figure;
        for k = 1:3
            plot(t, HGam(:, k), cols(k), 'Linewidth', 1.5);
            hold on;
        end
        set(gca, "FontSize", fontSize);
        plot(t, HGo, ':k', 'Linewidth', 4);
        xlim([0, t(end)]);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(A_3\leftrightarrow A_1)', '\it H_{\Gamma}');
        saveFigures(sprintf('%sMod2TrajectDetH%d_%d', dirName, kE, kS), drive);
        close();
        
        % Conv in this case is half of total number of points
        conv = 50;

        % Form graphs
        figure;
        for k = 1:3
            plot(t(1:conv), HGam(1:conv, k), cols(k), 'Linewidth', 1.5);
            hold on;
        end
        set(gca, "FontSize", fontSize);
        plot(t(1:conv), HGo(1:conv), ':k', 'Linewidth', 4);
        xlim([0, t(conv + 1)]);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(A_3\leftrightarrow A_1)', '\it H_{\Gamma}');
        saveFigures(sprintf('%sMod2TrajectDetHT%d_%d', dirName, kE, kS), drive);
        close();

        figure;
        plot(t(1:conv), HGo(1:conv), 'r', 'Linewidth', 1.5);
        hold on;
        plot(t(1:conv), HB(1:conv), 'b', 'Linewidth', 1.5);
        xlim([0, t(conv + 1)]);
        set(gca, "FontSize", fontSize);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\Gamma}', '\it H');
        saveFigures(sprintf('%sMod2TrajectDetHTB%d_%d', dirName, kE, kS), drive);
        close();
    end
end
