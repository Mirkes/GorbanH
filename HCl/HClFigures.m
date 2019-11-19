% HCl system figures
% Please run HClModPrepare for mdified system and HClPrepare for original
% system to calculate data to draw trajectories and other trajectory based
% graphs.
% 
% All figures have prefixes according to value of suf variable below.

% Usage of modified system:
modified = 1;

% Suffix for file names and names of substances
subst = {'A_1', 'A_2', 'A_3', 'A_4', 'A_5'};
if modified == 0
    suf = '';
    eq = [0.198, 0.004, 0.1995, 0.001, 0.6];
    subst = {'H_2', 'H', 'Cl_2', 'Cl', 'HCl'};
elseif modified == 1
    suf = 'Mod';
    eq = [0.2, 0.2, 0.25, 0.1, 0.4];
elseif modified == 2
    suf = 'Mod2';
    eq = [0.2, 0.2, 0.25, 0.1, 0.4];
elseif modified == 3
    suf = 'Mod3';
    eq = [0.2, 0.2, 0.25, 0.1, 0.4];
end

% What to draw
drawConcentrations = true;     % Concentrations versus time
drawHGB = true;                % H_Gamma and H versus time
drawPolyhedron = true;         % Draw polyhedron
drawTrajectory = true;         % Draw trajectory in polyhedron
drawHGLevels = true;           % Draw level set for Gorban's H function
drawHBLevels = true;           % Draw level set for Boltzman's H function
drawPartialEquilibria = true;  % Draw surfaces of partial equilibria
drawLevelSetSection = true;    % Draw section of partial equilibria

% parameters for level sets
drawEdges = 'none';

% Figure format
% drive = '-dpng'; % For png images
drive = '-depsc'; % For esp images
% drive = '-dpdf'; % For pdf format

% Define constants

% Level to draw level set
lvls = [-0.5, -0.7, -0.9, -1.1];

% Balances
b = [1, 1];

% Define name of folder to save figures
dirName = 'FiguresDot/';

% Specify font size for figures
fontSize = 20;

% Line width
lineWidth = 1.5;

% Form Gamma with stoichiometric vectors only
Gamma = [-1,  2,  0,  0, 0;
          0,  0, -1,  2, 0;
          0, -1, -1,  1, 1;
         -1,  1,  0, -1, 1];
nReact = size(Gamma, 1);

% Colours
cols = [0.635, 0.078, 0.184;...
        0.000, 0.000, 1.000;...
        0.460, 0.674, 0.188;...
        1.000, 0.000, 1.000];

% Vertices
vert = [  0,   0, 1;
          0,   0, 0;
          0, 0.5, 0;
        0.5,   0, 0;
        0.5, 0.5, 0];

legSubst = strcat({'\it '}, subst);
legReac = {[subst{1},'\leftrightarrow 2',subst{2}],...
           [subst{3},'\leftrightarrow 2',subst{4}],...
           [subst{2},'+',subst{3},'\leftrightarrow ',subst{5},'+',subst{4}],...
           [subst{4},'+',subst{1},'\leftrightarrow ',subst{5},'+',subst{2}]};
legHGam = [strcat({'\it H_{\gamma}('}, legReac, {')'}), '\it H_\Gamma'];
    
if modified > 0
    % Load data for HCl reaction from file ModTimX.mat (see file HClPrepare.m)
    load([suf, 'Tim.mat']);
    tt = t;
    % Calculate all concentrations because of files contains values for H2,
    % Cl2, and HCl only.
    cc = [c(:, 1), b(1) - 2 * c(:, 1) - c(:, 3), c(:, 2), b(2) - 2 * c(:, 2) - c(:, 3), c(:, 3)];
    nDot = size(cc, 1);

    % Load data for short initial fragment of graphs
    load([suf, 'TimStartStart.mat']);
    % Calculate all concentrations because of files contains values for H2,
    % Cl2, and HCl only.
    cs = [c(:, 1), b(1) - 2 * c(:, 1) - c(:, 3), c(:, 2), b(2) - 2 * c(:, 2) - c(:, 3), c(:, 3)];
    %    H2       H                              Cl2     Cl                            HCl
    ts = t;
    nDotss = size(cs, 1);

    % Load data for initial fragment of graphs
    load([suf, 'TimStart.mat']);
    % Calculate all concentrations because of files contains values for H2,
    % Cl2, and HCl only.
    c = [c(:, 1), b(1) - 2 * c(:, 1) - c(:, 3), c(:, 2), b(2) - 2 * c(:, 2) - c(:, 3), c(:, 3)];
    %    H2       H                              Cl2     Cl                            HCl
    nDots = size(c, 1);
elseif modified == 0
    % Load data for HCl reaction (files timX.mat, where X is number of stage
    % (see file HClPrepare.m)
    tim = [];
    y = [];
    for k = 1:6
        load(['tim',num2str(k),'.mat'])
        tim = [tim; t]; %#ok<AGROW>
        y = [y; c]; %#ok<AGROW>
    end

    % Calculate all concentrations because of files contains values for H2,
    % Cl2, and HCl only.
    cc = [y(:, 1), b(1) - 2 * y(:, 1) - y(:, 3), y(:, 2), b(2) - 2 * y(:, 2) - y(:, 3), y(:, 3)];
    %    H2       H                              Cl2     Cl                            HCl
    tt = tim;
    clear y tim
    nDot = size(cc, 1);

    % Load data for initial fragment of graphs
    load('timStart.mat');
    % Calculate all concentrations because of files contains values for H2,
    % Cl2, and HCl only.
    c = [c(:, 1), b(1) - 2 * c(:, 1) - c(:, 3), c(:, 2), b(2) - 2 * c(:, 2) - c(:, 3), c(:, 3)];
    %    H2       H                              Cl2     Cl                            HCl
    nDots = size(c, 1);
    cs = c;
    ts = t;
    nDotss = nDots;
    
end

if drawConcentrations
    % Form graph c(t) for full data
    figure;
    plot(tt, cc(:, 1), 'b-', 'Linewidth', lineWidth);
    hold on;
    plot(tt, cc(:, 2), 'b--', 'Linewidth', lineWidth);
    plot(tt, cc(:, 3), 'r--', 'Linewidth', lineWidth);
    plot(tt, cc(:, 4), 'r-', 'Linewidth', lineWidth);
    plot(tt, cc(:, 5), 'm-', 'Linewidth', lineWidth);
    set(gca, "FontSize", fontSize);
    legend(legSubst, 'Location', 'eastoutside');
    xlabel('Time');
    ylabel('Concentration');
    set(gcf,'Position',[100 100 600 420]);
    saveFigures(sprintf('%sHCl%sConc', dirName, suf), drive);
    close();
    % Form graph c(t) for start fragment
    figure;
    plot(t, c(:, 1), 'b-', 'Linewidth', lineWidth);
    hold on;
    plot(t, c(:, 2), 'b--', 'Linewidth', lineWidth);
    plot(t, c(:, 3), 'r--', 'Linewidth', lineWidth);
    plot(t, c(:, 4), 'r-', 'Linewidth', lineWidth);
    plot(t, c(:, 5), 'm-', 'Linewidth', lineWidth);
    
    set(gca, "FontSize", fontSize);
    legend(legSubst, 'Location', 'eastoutside');
    xlabel('Time');
    ylabel('Concentration');
    set(gcf,'Position',[100 100 600 420]);
    saveFigures(sprintf('%sHCl%sConcStart', dirName, suf), drive);
    close();
end

if drawHGB
    % The first graph: H_gamma for long time
    % Calculate partial equilibrium H_gamma functions, Gorban's and Boltzmann's
    % H functions.
    options = optimset('TolX', 1.e-4);
    HGam = zeros(nDot, nReact);
    HB = zeros(nDot, 1);
    for k = 1:nDot
        HGam(k, :) = HG(eq, cc(k, :), Gamma, options);
        % Calculate Boltzmann's H
        HB(k) = H(eq, cc(k, :));
    end
    HGo = max(HGam, [], 2);
    figure;
    for k = 1:nReact
        plot(tt, HGam(:, k), 'Color', cols(k, :), 'Linewidth', lineWidth);
        hold on;
    end
    set(gca, "FontSize", fontSize);
    plot(tt, HGo, ':k', 'Linewidth', 4);
    xlabel('Time');
    ylabel('Function value');
    legend(legHGam, 'Location', 'eastoutside');
    set(gcf,'Position',[100 100 855 420]);
    saveFigures(sprintf('%sHCl%sH', dirName, suf), drive);
    close();

    % The second graph: GH and BH for long time
    figure;
    plot(tt, HGo, 'r', 'Linewidth', lineWidth);
    hold on;
    plot(tt, HB, 'b', 'Linewidth', lineWidth);
    set(gca, "FontSize", fontSize);
    xlabel('Time');
    ylabel('Function value');
    legend('\it H_\Gamma', '\it H');
    saveFigures(sprintf('%sHCl%sHGHB', dirName, suf), drive);
    close();

    % The third graph. Very short time
    % Calculate partial equilibrium H_gamma functions, Gorban's and Boltzmann's
    % H functions.
    options = optimset('TolX', 1.e-5);
    HGam = zeros(nDotss, nReact);
    HB = zeros(nDotss, 1);
    for k = 1:nDotss
        HGam(k, :) = HG(eq, cs(k, :), Gamma, options);
        % Calculate Boltzmann's H
        HB(k) = H(eq, cs(k, :));
    end
    HGo = max(HGam, [], 2);
    figure;
    for k = 1:nReact
        plot(ts, HGam(:, k), 'Color', cols(k, :), 'Linewidth', lineWidth);
        hold on;
    end
    set(gca, "FontSize", fontSize);
    plot(ts, HGo, ':k', 'Linewidth', 4);
    xlabel('Time');
    ylabel('Function value');
    legend(legHGam, 'Location', 'eastoutside');
    set(gcf,'Position',[100 100 855 420]);
    saveFigures(sprintf('%sHCl%sHStart', dirName, suf), drive);
    close();

    % The second graph: GH and BH for very short time
    figure;
    plot(ts, HGo, 'r', 'Linewidth', lineWidth);
    hold on;
    plot(ts, HB, 'b', 'Linewidth', lineWidth);
    set(gca, "FontSize", fontSize);
    xlabel('Time');
    ylabel('Function value');
    legend('\it H_\Gamma', '\it H');
    saveFigures(sprintf('%sHCl%sHGHBStart', dirName, suf), drive);
    close();
end

% Number of splitting per coordinate
nSplits = 100;

if drawPolyhedron
    figure;
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    saveFigures(sprintf('%sHCl%sPolyhedron', dirName, suf), drive);
    close();
end

if drawTrajectory
    figure;
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    plot3(cc(:, 1), cc(:, 3), cc(:, 5), 'r', 'Linewidth', lineWidth); 
    saveFigures(sprintf('%sHCl%sTrajectory', dirName, suf), drive);
    close();
end

% Create set of points and triangulation if necessary
if drawHGLevels || drawHBLevels
    % create sets of dots
    vX = linspace(0, 0.5, nSplits);
    vY = linspace(0, 0.5, nSplits);
    vZ = linspace(0,   1, nSplits);
    % base
    dotsX = repmat(vX, 1, nSplits)';
    dotsY = repmat(vY, nSplits, 1);
    dotsY = dotsY(:);
    dotsZ = zeros(nSplits * nSplits, 1);
    % matrix of dot numbers
    tmp = reshape(1:nSplits * nSplits, nSplits, nSplits);
    tri = [reshape(tmp(1:end - 1, 1:end - 1), [], 1),...
        reshape(tmp(1:end - 1, 2:end), [], 1),...
        reshape(tmp(2:end, 1:end - 1), [], 1);...
        reshape(tmp(1:end - 1, 2:end), [], 1),...
        reshape(tmp(2:end, 1:end - 1), [], 1),...
        reshape(tmp(2:end, 2:end), [], 1)];
    
    % Face c1,c5
    dotsX1 = repmat(vX, 1, nSplits)';
    dotsY1 = zeros(nSplits * nSplits, 1);
    dotsZ1 = repmat(vZ, nSplits, 1);
    dotsZ1 = dotsZ1(:);
    
    % Find indeces of inappropriate points
    ind = (2 * dotsX1 + dotsZ1) > 1;
    
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    
    tri1 = delaunay(dotsX1, dotsZ1);
    
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
    
    % Face c3,c5
    dotsX1 = zeros(nSplits * nSplits, 1);
    dotsY1 = repmat(vY,     1, nSplits)';
    dotsZ1 = repmat(vZ, nSplits,     1);
    dotsZ1 = dotsZ1(:);
    
    % Find indeces of inappropriate points
    ind = (2 * dotsY1 + dotsZ1) > 1;
    
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    
    tri1 = delaunay(dotsY1, dotsZ1);
    
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
    
    % Face 1 of c1, c3, c5
    dotsX1 = repmat(vX, 1, nSplits)';
    dotsY1 = repmat(vY, nSplits, 1);
    dotsY1 = dotsY1(:);
    dotsZ1 = repmat(fliplr(vZ), 1, nSplits)';
    
    % Find indeces of inappropriate points
    ind = (2 * dotsY1 + dotsZ1) > 1;
    
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    
    tri1 = delaunay(dotsX1, dotsY1);
    
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
    
    % Face 2 of c1, c3, c5
    dotsX1 = repmat(vX, 1, nSplits)';
    dotsY1 = repmat(vY, nSplits, 1);
    dotsY1 = dotsY1(:);
    dotsZ1 = repmat(fliplr(vZ), nSplits, 1);
    dotsZ1 = dotsZ1(:);
    
    % Find indeces of inappropriate points
    ind = (2 * dotsX1 + dotsZ1) > 1;
    
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    
    tri1 = delaunay(dotsX1, dotsY1);
    
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
end

if drawHGLevels
    for lvl = lvls
        % Search levelset
        res = HClLevelSearch(eq, lvl, Gamma, [dotsX, dotsY, dotsZ], b);
        resLev = zeros(size(res, 1), 1);
        for kk = 1:size(res, 1)
            resLev(kk) = HClGH(eq, res(kk, :), Gamma, b);
        end
        % Index of non fitted points (Borders)
        ind = abs(resLev - lvl) > 0.01 * abs(lvl);
        % Create colours and change coloutr for non-fitted points
        col = 0.4 * ones(size(res, 1), 1);
        col(ind) = 0;
        % Draw surface
        figure;
        tsr = trisurf(tri, res(:, 1), res(:, 2), res(:, 3), col, 'EdgeColor', drawEdges);
        hold on
        drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
        title(sprintf('\\it H_\\Gamma=%4.1f',lvl));
        lt = camlight('heardligth','local');
        tsr.DiffuseStrength = 0.9;
        saveFigures(strrep(sprintf('%sHCl%sGLevel%4.1fPos0', dirName, suf, lvl), '.', '_'), drive);
        view(-27, 25);
        camlight(lt, 'heardligth','local');
        saveFigures(strrep(sprintf('%sHCl%sGLevel%4.1fPos1', dirName, suf, lvl), '.', '_'), drive);
        close();
    end
end

if drawHBLevels
    for lvl = lvls
        % Search levelset
        res = HClLevelSearchH(eq, lvl, [dotsX, dotsY, dotsZ], b);
        resLev = zeros(size(res, 1), 1);
        for kk = 1:size(res, 1)
            resLev(kk) = HClBH(eq, res(kk, :), b);
        end
        % Index of non fitted points (Borders)
        ind = abs(resLev - lvl) > 0.01 * abs(lvl);
        % Create colours and change colour for non-fitted points
        col = 0.4 * ones(size(res, 1), 1);
        col(ind) = 0;
        % Draw surface
        figure;
        tsr = trisurf(tri, res(:, 1), res(:, 2), res(:, 3), col, 'EdgeColor', drawEdges);
        hold on
        drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
        title(sprintf('\\it H=%4.1f',lvl));
        lt = camlight('heardligth','local');
        tsr.DiffuseStrength = 0.9;
        saveFigures(strrep(sprintf('%sHCl%sBLevel%4.1fPos0', dirName, suf, lvl), '.', '_'), drive);
        view(-27, 25);
        camlight(lt, 'heardligth','local');
        saveFigures(strrep(sprintf('%sHCl%sBLevel%4.1fPos1', dirName, suf, lvl), '.', '_'), drive);
        close();
    end
end

if drawPartialEquilibria
    nSplits = 20;
    % create sets of dots
    vX = linspace(0, 0.5, nSplits);
    vY = linspace(0, 0.5, nSplits);
    vZ = linspace(0,   1, nSplits);

    % The first reaction
    gamma = Gamma(1, :);
    % Form surface prototype
    % Face c3,c5
    dotsX = zeros(nSplits * nSplits, 1);
    dotsY = repmat(vY, 1, nSplits)';
    dotsZ = repmat(vZ, nSplits, 1);
    dotsZ = dotsZ(:);
    % Find indeces of inappropriate points
    ind = (2 * dotsY + dotsZ) > 1;
    % Remove inappropriate points
    dotsX(ind) = [];
    dotsY(ind) = [];
    dotsZ(ind) = [];
    % Create triangulation
    tri = delaunay(dotsY, dotsZ);
    % Search minimum for all points.
    for k = 1:size(dotsX, 1)
        [dotsX(k), dotsY(k), dotsZ(k)] = HClPartial(eq, [dotsX(k), dotsY(k), dotsZ(k)], gamma, b);
    end
    figure;
    trisurf(tri, dotsX, dotsY, dotsZ, 'FaceColor', 'c');
    hold on
    if drawTrajectory
        plot3(cc(:, 1), cc(:, 3), cc(:, 5), 'r', 'Linewidth', lineWidth); 
    end
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    title(['\it ', legReac{1}]);
    saveFigures(sprintf('%sHCl%sPartial_1', dirName, suf), drive);
    close();
    
    % The second reaction
    gamma = Gamma(2, :);
    % Face c1,c5
    % Form surface prototype
    dotsX = repmat(vX, 1, nSplits)';
    dotsY = zeros(nSplits * nSplits, 1);
    dotsZ = repmat(vZ, nSplits, 1);
    dotsZ = dotsZ(:);
    % Find indeces of inappropriate points
    ind = (2 * dotsX + dotsZ) > 1;
    % Remove inappropriate points
    dotsX(ind) = [];
    dotsY(ind) = [];
    dotsZ(ind) = [];
    % Create triangulation
    tri = delaunay(dotsX, dotsZ);
    % Search minimum for all points.
    for k = 1:size(dotsX, 1)
        [dotsX(k), dotsY(k), dotsZ(k)] = HClPartial(eq, [dotsX(k), dotsY(k), dotsZ(k)], gamma, b);
    end
    figure;
    trisurf(tri, dotsX, dotsY, dotsZ, 'FaceColor', 'c');
    hold on
    if drawTrajectory
        plot3(cc(:, 1), cc(:, 3), cc(:, 5), 'r', 'Linewidth', lineWidth); 
    end
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    title(['\it ', legReac{2}]);
    view(227, 27);
    saveFigures(sprintf('%sHCl%sPartial_2', dirName, suf), drive);
    close();
    
    % The third reaction
    gamma = Gamma(3, :);
    % Face c1,c5
    % Form surface prototype
    dotsX = repmat(vX, 1, nSplits)';
    dotsY = zeros(nSplits * nSplits, 1);
    dotsZ = repmat(vZ, nSplits, 1);
    dotsZ = dotsZ(:);
    % Find indeces of inappropriate points
    ind = (2 * dotsX + dotsZ) > 1;
    % Remove inappropriate points
    dotsX(ind) = [];
    dotsY(ind) = [];
    dotsZ(ind) = [];
    % Create triangulation
    tri = delaunay(dotsX, dotsZ);
    % Face 1 of c1, c3, c5
    dotsX1 = repmat(vX, 1, nSplits)';
    dotsY1 = repmat(vY, nSplits, 1);
    dotsY1 = dotsY1(:);
    dotsZ1 = repmat(fliplr(vZ), 1, nSplits)';
    % Find indeces of inappropriate points
    ind = (2 * dotsY1 + dotsZ1) > 1;
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    tri1 = delaunay(dotsX1, dotsY1);
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
    % Search minimum for all points.
    for k = 1:size(dotsX, 1)
        [dotsX(k), dotsY(k), dotsZ(k)] = HClPartial(eq, [dotsX(k), dotsY(k), dotsZ(k)], gamma, b);
    end
    figure;
    trisurf(tri, dotsX, dotsY, dotsZ, 'FaceColor', 'c');
    hold on
    if drawTrajectory
        plot3(cc(:, 1), cc(:, 3), cc(:, 5), 'r', 'Linewidth', lineWidth); 
    end
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    title(['\it ', legReac{3}]);
    view(69.7, 48.1334);
    saveFigures(sprintf('%sHCl%sPartial_3', dirName, suf), drive);
    close();

    % The fourth reaction
    gamma = Gamma(4, :);
    % Face c1,c5
    % Form surface prototype
    dotsX = zeros(nSplits * nSplits, 1);
    dotsY = repmat(vY, 1, nSplits)';
    dotsZ = repmat(vZ, nSplits, 1);
    dotsZ = dotsZ(:);
    % Find indeces of inappropriate points
    ind = (2 * dotsY + dotsZ) > 1;
    % Remove inappropriate points
    dotsX(ind) = [];
    dotsY(ind) = [];
    dotsZ(ind) = [];
    % Create triangulation
    tri = delaunay(dotsY, dotsZ);
    % Face 2 of c1, c3, c5
    dotsX1 = repmat(vX, 1, nSplits)';
    dotsY1 = repmat(vY, nSplits, 1);
    dotsY1 = dotsY1(:);
    dotsZ1 = repmat(fliplr(vZ), nSplits, 1);
    dotsZ1 = dotsZ1(:);
    % Find indeces of inappropriate points
    ind = (2 * dotsX1 + dotsZ1) > 1;
    % Remove inappropriate points
    dotsX1(ind) = [];
    dotsY1(ind) = [];
    dotsZ1(ind) = [];
    % Create triangulation
    tri1 = delaunay(dotsX1, dotsY1);
    % Append new points into data
    tri = [tri; tri1 + size(dotsX, 1)];
    dotsX = [dotsX; dotsX1];
    dotsY = [dotsY; dotsY1];
    dotsZ = [dotsZ; dotsZ1];
    % Search minimum for all points.
    for k = 1:size(dotsX, 1)
        [dotsX(k), dotsY(k), dotsZ(k)] = HClPartial(eq, [dotsX(k), dotsY(k), dotsZ(k)], gamma, b);
    end
    figure;
    trisurf(tri, dotsX, dotsY, dotsZ, 'FaceColor', 'c');
    hold on
    if drawTrajectory
        plot3(cc(:, 1), cc(:, 3), cc(:, 5), 'r', 'Linewidth', lineWidth); 
    end
    drawPolyhedronF(vert, lineWidth, fontSize, eq, legSubst);
    title(['\it ', legReac{4}]);
    view(191.8334, 40.1334);
    saveFigures(sprintf('%sHCl%sPartial_4', dirName, suf), drive);
    close();
end

if drawLevelSetSection
    % Calcualte 9 figures and draw them on 9 subplots
    figure;
    % level
    lvl = -0.9;
    for k = 1:9
        HCl = 0.1 * k;
        subplot(3, 3, k);
        siz = (1 - HCl) / 2;
        rectangle('Position', [0, 0, siz, siz]);
        axis([0, 0.5, 0, 0.5]);
        axis square;
        hold on;
        
        % 1. Search equilibrium in this plane
        ceq = zeros(1, 5);
        ceq(5) = HCl;
        r = roots([4 * eq(1), -eq(2) .^ 2 - 4 * eq(1) * (b(1) - ceq(5)), eq(1) * ((b(1) - ceq(5)) .^ 2)]);
        r(r > siz) = [];
        r(r < 0) = [];
        ceq(1) = r(1);
        ceq(2) = b(1) - ceq(5) - 2 * ceq(1);
        r = roots([4 * eq(3), -eq(4) .^ 2 - 4 * eq(3) * (b(2) - ceq(5)), eq(3) * ((b(2) - ceq(5)) .^ 2)]);
        r(r > siz) = [];
        r(r < 0) = [];
        ceq(3) = r(1);
        ceq(4) = b(2) - ceq(5) - 2 * ceq(3);
        % 2. Calculate H(ceq). If H(ceq)>lvl then stop
        if H(eq, ceq) < lvl
            % 2. Form border rectangle
            n = 100;
            zers = zeros(n, 1);
            sizes = siz * ones(n, 1);
            last = HCl * ones(n, 1);
            vec = linspace(0, siz, n)';
            revVec = flipud(vec);
            % Gorban's H
            borders = [zers, vec, last; vec, sizes, last; sizes, revVec, last; revVec, zers, last];
            nDot = size(borders, 1);
            borders = borders - ceq([1, 3, 5]);
            for kk = 1:nDot
                chi = fminbnd(@(x) (lvl - HClGH(eq, ceq([1, 3, 5]) + x * borders(kk, :), Gamma, b))^2, 0.0000001, 0.999999);
                borders(kk, :) = ceq([1, 3, 5]) + chi * borders(kk, :);
            end
            plot(borders(:, 1), borders(:, 2), 'b-', 'Linewidth', lineWidth);
            
            % Boltzmann's H
            borders = [zers, vec, last; vec, sizes, last; sizes, revVec, last; revVec, zers, last];
            nDot = size(borders, 1);
            borders = borders - ceq([1, 3, 5]);
            for kk = 1:nDot
                chi = fminbnd(@(x) (lvl - HClBH(eq, ceq([1, 3, 5]) + x * borders(kk, :), b))^2, 0.0000001, 0.999999);
                borders(kk, :) = ceq([1, 3, 5]) + chi * borders(kk, :);
            end
            plot(borders(:, 1), borders(:, 2), 'm-', 'Linewidth', lineWidth);
        end
        %Draw equilibrium point
        plot3(ceq(1), ceq(3), ceq(5), 'r*', 'Linewidth', lineWidth);
        rectangle('Position', [0, 0, siz, siz]);
        title(sprintf('\\it HCl=%4.2f', HCl));
    end
    set(gcf, 'Position', [2197, 792, 922, 910]);
    saveFigures(sprintf('%sHCl%sLevelsSections', dirName, suf), drive);
    close();
end

function drawPolyhedronF(vert, lineWidth, fontSize, eq, subst)
% drawPolyhedronF depict reaction polyhedron.
% Inputs:
%   vert is matrix with coordinates of polyhedron vertices.
%   lineWidth is required width of lines.
%   fontSize is required size of font.
%   eq is equilibrium point.
    plot3(vert([2, 3, 5, 4, 2], 1), vert([2, 3, 5, 4, 2], 2),...
        vert([2, 3, 5, 4, 2], 3), 'k-', 'Linewidth', lineWidth);
    hold on
    set(gca, "FontSize", fontSize);
    for k = 2:5
        plot3(vert([1, k], 1), vert([1, k], 2), vert([1, k], 3),...
            'k-', 'Linewidth', lineWidth);
    end
    % Depict equilibrium
    plot3(eq(1), eq(3), eq(5), 'r*', 'Linewidth', lineWidth);
    xlabel(subst{1});
    ylabel(subst{3});
    zlabel(subst{5});
    axis([0, 0.5, 0, 0.5, 0, 1]);
    axis equal
    view(45.1667, 31.0667);
end