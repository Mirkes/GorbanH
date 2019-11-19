% Model reaction figures
% To form figures for randomly generated equilibria it is necessary to
% set required number of equilibrias into variable nRep in line 20

% Define name of folder to save figures
dirName = 'FiguresDot/';

% Figure format
drive = '-dpng'; % For png images
% drive = '-depsc'; % For esp images
% drive = '-dpdf'; % For pdf format

% Specify font size for figures
fontSize = 20;

% Colours
cols = ['r', 'm', 'b'];

% Define dot to illustrate the calculation of Gorban's H function
c0 = [0.6, 0.4, 0.06];
c0 = c0 / sum(c0);

% Form Gamma with stoichiometric vectors only
Gamma = [-1,  1, 0;
          0, -1, 1; 
         -2,  1, 1];
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

% Draw Triangle with dot and directions
for k = 1:nRep
    DrawTriangle(eqs(k, :), c0);
    saveFigures(sprintf('%sdirect%03d', dirName, k), drive);
    close();
end

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
    saveFigures(sprintf('%slevels%03d', dirName, kk), drive);
    close();
end

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
    saveFigures(sprintf('%slevelsBH%03d', dirName, kk), drive);
    close();
end

% Define origin of trajectory
c0 = [0.9999, 0.00005, 0.00005];
c0 = c0 / sum(c0);

% Define sum of reaction rates to present four versions of models with
% different reaction rates for each equilibrium
kSumS = [1, 1, 1; 10, 1, 1; 10, 5, 1; 10, 1, 5];

% Stop time
stopTime = 10;

% Form time span for integration
ts = linspace(0, stopTime, stopTime * 100 + 1);

% Convergence level to select most interesting initial part
cLevel = -0.99;

% Equilibrium number
for kE = 1:nRep
    eq = eqs(kE, :);

    % Sum number
    for kS = 1:size(kSumS, 1)
        kSum = kSumS(kS, :);
        
        % Calculate fraction of reaction rate constants
        kFrac = zeros(1, nReac);
        for r = 1:nReac
            indP = Gamma(r, :) > 0;
            indM = Gamma(r, :) < 0;
            kFrac(r) = prod(eq(indP) .^ Gamma(r, indP))/prod(eq(indM) .^ (-Gamma(r, indM)));
        end
        
        % Calculate all reaction rate constants
        kM = kSum ./ (kFrac + 1);
        kP = kFrac .* kM;
        
        % Parameters
        opts = odeset('Reltol',1e-13,'AbsTol',1e-14); %,'Stats','on');
        
        % Integrate
        [t, c] = ode113(@(tt, y) modelODE(tt, y, kP, kM), ts, c0([1, 3]), opts);
        
        % Calculate complete set of coordinates
        c = [c(:, 1), 1 - c(:, 1) - c(:, 2), c(:, 2)];
        nDots = size(c, 1);
        
        % Draw trajectory
        DrawTriangleLevels(eq, c);
        saveFigures(sprintf('%straject%d_%d', dirName, kE, kS), drive);
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
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(2A_1\leftrightarrow A_{2}+A_3 )', '\it H_{\Gamma}');
        saveFigures(sprintf('%strajectH%d_%d', dirName, kE, kS), drive);
        close();
        
        % Serch point of convergence
        conv = find(HB < cLevel, 1);
        if isempty(conv)
            conv = length(HB);
        end
        conv = find(t >= (ceil(t(conv) / 0.5) * 0.5), 1);
        
        % Form graphs
        figure;
        for k = 1:3
            plot(t(1:conv), HGam(1:conv, k), cols(k), 'Linewidth', 1.5);
            hold on;
        end
        set(gca, "FontSize", fontSize);
        plot(t(1:conv), HGo(1:conv), ':k', 'Linewidth', 4);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\gamma}(A_1\leftrightarrow A_2 )', '\it H_{\gamma}(A_2\leftrightarrow A_3 )',...
            '\it H_{\gamma}(2A_1\leftrightarrow A_{2}+A_3 )', '\it H_{\Gamma}');
        saveFigures(sprintf('%strajectHT%d_%d', dirName, kE, kS), drive);
        close();

        figure;
        plot(t(1:conv), HGo(1:conv), 'r', 'Linewidth', 1.5);
        hold on;
        plot(t(1:conv), HB(1:conv), 'b', 'Linewidth', 1.5);
        set(gca, "FontSize", fontSize);
        xlabel('Time');
        ylabel('Function value');
        legend('\it H_{\Gamma}', '\it H');
        saveFigures(sprintf('%strajectHTB%d_%d', dirName, kE, kS), drive);
        close();
    end
end
