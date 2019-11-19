% General constants
% Variable modified must be
% 0 for original data
% 1 for the first modified system
modified = 1;

% Call WGSPrepare to get app required information. Variable modified MUST
% be set before
WGSPrepare;

% Specify what to draw
drawConcentrations = true;     % Concentrations versus time
drawHGB = true;                % H_Gamma and H versus time
drawPartial = true;            % Draw lines of partial equilibria
drawTrajectory = true;         % Trajectory in the reaction polygon
drawLevelSets = true;          % Level sets for Gorban's and Boltzmann's H functions
drawLevelSetsAndPartial = true;% Add partial equilibria to Level sets for Gorban's and Boltzmann's H functions

% Define name of folder to save figures
dirName = 'FiguresDot/';

% Figure format
% drive = '-dpng'; % For png images
drive = '-depsc'; % For esp images
% drive = '-dpdf'; % For pdf format

% Decoration
fontSize = 20;

% Colours
cols = ['r', 'b'];

% Form legends strings
labReact1 = ['\it ',subst{1},'+',subst{5},'\leftrightarrow ',subst{2},'+',subst{6}];
labReact2 = ['\it ',subst{3},'+',subst{6},'\leftrightarrow ',subst{4},'+',subst{5}];
legReact1 = ['\it H_{\gamma}(',subst{1},'+',subst{5},'\leftrightarrow ',subst{2},'+',subst{6},')'];
legReact2 = ['\it H_{\gamma}(',subst{3},'+',subst{6},'\leftrightarrow ',subst{4},'+',subst{5},')'];
legSubst = strcat({'\it '}, subst);

% Equilibrium
eq = zeros(1, 6);
eq(5) = param(2);
eq(3) = coeq * b(2);
eq(4) = b(2) - eq(3);
eq(6) = b(4) - eq(5);
eq(2) = eq(4) + eq(6);
eq(1) = b(1) - eq(2);

% Initial point
c0 = [0.9999 * b(1), 0.9999 * b(2)];

% Form Gamma with stoichiometric vectors only
Gamma = [-1, 1,  0, 0, -1,  1; 
          0, 0, -1, 1,  1, -1];

% Line width
lineWidth = 1.5;

if drawConcentrations
    % Draw graphs of concentrations over time
    figure;
    plot(cFull(:, 1), cFull(:, 2), 'b-', 'Linewidth', lineWidth); 
    set(gca, "FontSize", fontSize);
    hold on;
    plot(cFull(:, 1), cFull(:, 3), 'b--', 'Linewidth', lineWidth); 
    plot(cFull(:, 1), cFull(:, 4), 'r-', 'Linewidth', lineWidth); 
    plot(cFull(:, 1), cFull(:, 5), 'r--', 'Linewidth', lineWidth); 
    plot(cFull(:, 1), cFull(:, 6), 'm-', 'Linewidth', lineWidth); 
    plot(cFull(:, 1), cFull(:, 7), 'm--', 'Linewidth', lineWidth); 
    legend(legSubst,  'Location', 'northoutside', 'Orientation','horizontal');
    xlabel('Time');
    ylabel('Concentration');
    set(gcf,'Position',[100 100 800 500]);
    saveFigures(sprintf('%sWGS%sConcentrations', dirName, suf), drive);
    close();

    % Draw graphs of concentrations over logarithm of time
    figure;
    semilogx(cShort(:, 1), cShort(:, 2), 'b-', 'Linewidth', lineWidth); 
    set(gca, "FontSize", fontSize);
    hold on;
    semilogx(cShort(:, 1), cShort(:, 3), 'b-', 'Linewidth', lineWidth); 
    semilogx(cShort(:, 1), cShort(:, 4), 'r-', 'Linewidth', lineWidth); 
    semilogx(cShort(:, 1), cShort(:, 5), 'r--', 'Linewidth', lineWidth); 
    semilogx(cShort(:, 1), cShort(:, 6), 'm-', 'Linewidth', lineWidth); 
    semilogx(cShort(:, 1), cShort(:, 7), 'm--', 'Linewidth', lineWidth); 
    legend(legSubst,  'Location', 'northoutside', 'Orientation','horizontal');
    xlabel('Time (log)');
    ylabel('Concentration');
    set(gcf,'Position',[100 100 800 500]);
    saveFigures(sprintf('%sWGS%sConcentrationsLog', dirName, suf), drive);
    close();
end

% Form graphs of H_\gamma and H_Gamma
if drawHGB
    % Calculate HG and HB for cReact
    nDots = size(cReact, 1);
    nGam = size(Gamma, 1);
    % Memory allocation
    HGam = zeros(nDots, nGam);
    HB = zeros(nDots, 1);
    for k = 1:nDots
        HGam(k, :) = HG(eq, cReact(k, 2:end), Gamma);
        % Calculate Boltzmann's H
        HB(k) = H(eq, cReact(k, 2:end));
    end
    HGo = max(HGam, [], 2);

    % Complete graph
    figure;
    for k = 1:nGam
        plot(cReact(:, 1), HGam(:, k), cols(k), 'Linewidth', 1.5);
        hold on;
    end
    set(gca, "FontSize", fontSize);
    plot(cReact(:, 1), HGo, ':k', 'Linewidth', 4);
    xlabel('Time');
    ylabel('Function value');
    xlim([0, cReact(end, 1)]);
    legend(legReact1, legReact2, '\it H_\Gamma');
    saveFigures(sprintf('%sWGS%sH', dirName, suf), drive);
    close();

    %Form graph of GH and BH
    figure;
    plot(cReact(:, 1), HGo, 'r', 'Linewidth', 1.5);
    hold on;
    plot(cReact(:, 1), HB, 'b', 'Linewidth', 1.5);
    set(gca, "FontSize", fontSize);
    xlabel('Time');
    ylabel('Function value');
    xlim([0, cReact(end, 1)]);
    legend('\it H_\Gamma', '\it H');
    saveFigures(sprintf('%sWGS%sHGHB', dirName, suf), drive);
    close();
    
    % Calculate HG and HB for cShort
    nDots = size(cShort, 1);
    nGam = size(Gamma, 1);
    % Memory allocation
    HGam = zeros(nDots, nGam);
    HB = zeros(nDots, 1);
    for k = 1:nDots
        HGam(k, :) = HG(eq, cShort(k, 2:end), Gamma);
        % Calculate Boltzmann's H
        HB(k) = H(eq, cShort(k, 2:end));
    end
    HGo = max(HGam, [], 2);

    % Form graph for the time beginning
    figure;
    for k = 1:nGam
        plot(cShort(:, 1), HGam(:, k), cols(k), 'Linewidth', 1.5);
        hold on;
    end
    set(gca, "FontSize", fontSize);
    plot(cShort(:, 1), HGo, ':k', 'Linewidth', 4);
    xlabel('Time');
    ylabel('Function value');
    xlim([0, cShort(end, 1)]);
    legend(legReact1, legReact2, '\it H_\Gamma');
    saveFigures(sprintf('%sWGS%sHStart', dirName, suf), drive);
    close();

    %Form graph of GH and BH
    figure;
    plot(cShort(:, 1), HGo, 'r', 'Linewidth', 1.5);
    hold on;
    plot(cShort(:, 1), HB, 'b', 'Linewidth', 1.5);
    set(gca, "FontSize", fontSize);
    xlabel('Time');
    ylabel('Function value');
    xlim([0, cShort(end, 1)]);
    legend('\it H_\Gamma', '\it H');
    saveFigures(sprintf('%sWGS%sHGHBStart', dirName, suf), drive);
    close();
end

% Form list of border points to calculate level sets and partial equilibria
% lines.
% Form set of points on border of reaction polygon
borders = [linspace(0, 1, 100)', linspace(0, 1, 100)';... % Main diagonal
           linspace(1, 1-b(4), 50)', ones(50, 1);... % Top side
           linspace(1-b(4), 0, 100)', linspace(1-b(4), 0, 100)' + b(4);... % The second diagonal
           zeros(50, 1),  linspace(b(4), 0, 50)'];

if drawPartial || drawLevelSetsAndPartial       
    % memory allocation
    part1 = [linspace(0, 1, 100)', linspace(0, 1, 100)'];
    % For each gamma separately. It is simple to do by hands
    % The first gamma
    for k = 1:size(part1, 1)
        part1(k, :) = WGSPartial(eq, part1(k, :), Gamma(1, :), b);
    end
    % memory allocation
    part2 = [linspace(0, 1, 100)', linspace(0, 1, 100)'];
    for k = 1:size(part1, 1)
        part2(k, :) = WGSPartial(eq, part2(k, :), Gamma(2, :), b);
    end
end

if drawPartial
    drawPolygon(true, eq, b, lineWidth, fontSize, subst);
    % define horizontal size of arrow
    has = 0.1;
    plot(part1(:, 1), part1(:, 2), 'm-', 'Linewidth', lineWidth);
    % Draw label for partial equilibrium
    ps = size(part1, 1) - 20;
    quiver(part1(ps, 1) - has, part1(ps, 2), has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
    text(part1(ps, 1) - has, part1(ps, 2) +0.01, labReact1, "FontSize", fontSize, 'HorizontalAlignment', 'right');
    
    % The second gamma
    plot(part2(:, 1), part2(:, 2), 'b-', 'Linewidth', lineWidth);
    % Draw label for partial equilibrium
    ps = size(part2, 1) - 13;
    quiver(part2(ps, 1) - has, part2(ps, 2), has, 0, 'b', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
    text(part2(ps, 1) - has, part2(ps, 2) +0.01, labReact2, "FontSize", fontSize, 'HorizontalAlignment', 'right');

    % DRaw trajectory
    if drawTrajectory
        % Draw trajectory
        plot(cFull(:, 2), cFull(:, 4), 'r-', 'Linewidth', lineWidth);
    end
    
    saveFigures(sprintf('%sWGS%sPartial', dirName, suf), drive);
    close();
end

if drawTrajectory
    drawPolygon(true, eq, b, lineWidth, fontSize, subst);
    % Draw trajectory
    plot(cFull(:, 2), cFull(:, 4), 'r-', 'Linewidth', lineWidth);
    saveFigures(sprintf('%sWGS%sTrajectory', dirName, suf), drive);
    close();
end

if drawLevelSets
    % Number of H levels
    nLev = 20;

    % The furthest point is point (1, 1) (c(1), c(3))
    % Calculate vector from this point to equilibrium
    dir = [1, 1] - eq([1, 3]);

    % Subtract equilibrium
    borders = borders - eq([1, 3]);

    % The first figure: Gorban's H
    figure;
    % Calculate each level set and depict it
    for k = 1:nLev
        lev = WGSGH(eq, eq([1, 3]) + dir * k / (nLev + 1), Gamma, b);
        res = WGSlevelSearch(eq, lev, Gamma, borders, b);
        % Depict found set
        plot(res(:, 1), res(:, 2), 'b-', 'Linewidth', lineWidth);
        hold on;
    end

    if drawLevelSetsAndPartial
        % Draw partial equilibria if requested
        % The first reaction
        has = 0.1;
        plot(part1(:, 1), part1(:, 2), 'm-', 'Linewidth', lineWidth);
        % Draw label for partial equilibrium
        ps = size(part1, 1) - 20;
        quiver(part1(ps, 1) - has, part1(ps, 2), has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(part1(ps, 1) - has, part1(ps, 2) +0.01, labReact1, "FontSize", fontSize, 'HorizontalAlignment', 'right');
        % The second gamma
        plot(part2(:, 1), part2(:, 2), 'm-', 'Linewidth', lineWidth);
        % Draw label for partial equilibrium
        ps = size(part2, 1) - 13;
        quiver(part2(ps, 1) - has, part2(ps, 2), has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(part2(ps, 1) - has, part2(ps, 2) +0.01, labReact2, "FontSize", fontSize, 'HorizontalAlignment', 'right');
    end
    
    % Add polygon and decorations
    drawPolygon(false, eq, b, lineWidth, fontSize, subst);
    saveFigures(sprintf('%sWGS%sGLevels', dirName, suf), drive);
    close();

    % The second figure: Boltzmann's H
    figure;
    % Calculate each level set and depict it
    for k = 1:nLev
        lev = WGSBH(eq, eq([1, 3]) + dir * k / (nLev + 1), b);
        res = WGSlevelSearchH(eq, lev, borders, b);
        % Depict found set
        plot(res(:, 1), res(:, 2), 'b-', 'Linewidth', lineWidth);
        hold on;
    end
    
    if drawLevelSetsAndPartial
        % Draw partial equilibria if requested
        % The first reaction
        has = 0.1;
        plot(part1(:, 1), part1(:, 2), 'm-', 'Linewidth', lineWidth);
        % Draw label for partial equilibrium
        ps = size(part1, 1) - 20;
        quiver(part1(ps, 1) - has, part1(ps, 2), has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(part1(ps, 1) - has, part1(ps, 2) +0.01, labReact1, "FontSize", fontSize, 'HorizontalAlignment', 'right');
        % The second gamma
        plot(part2(:, 1), part2(:, 2), 'm-', 'Linewidth', lineWidth);
        % Draw label for partial equilibrium
        ps = size(part2, 1) - 13;
        quiver(part2(ps, 1) - has, part2(ps, 2), has, 0, 'm', 'MaxHeadSize', 0.9, 'Linewidth', 1.5);
        text(part2(ps, 1) - has, part2(ps, 2) +0.01, labReact2, "FontSize", fontSize, 'HorizontalAlignment', 'right');
    end

    % Add polygon and decorations
    drawPolygon(false, eq, b, lineWidth, fontSize, subst);
    saveFigures(sprintf('%sWGS%sBLevels', dirName, suf), drive);
    close();
end


% Function to draw reaction polygon
function drawPolygon(newFig, eq, b, lineWidth, fontSize, subs)
% drawPolygon depicts reaction polygon and equilibrium
% Inputs:
%   newFig is boolean to coltrol figure creation
%   eq is equilibrium
%   b is array of balances
%   lineWidth is width of lines
%   fontSize is size of font

    % Create figure if requested
    if newFig
        figure;
    end
    % Draw polygon
    plot([0, 1], [0, 1], 'k-', 'Linewidth', lineWidth);
    set(gca, "FontSize", fontSize);
    hold on;
    plot([0, 1 - b(4)], [b(4), 1], 'k-', 'Linewidth', lineWidth);
    plot([0, 0], [0, b(4)], 'k-', 'Linewidth', lineWidth);
    plot([1-b(4), 1], [1, 1], 'k-', 'Linewidth', lineWidth);

    % Depict equilibrium
    plot(eq(1), eq(3), '*r', 'Linewidth', 1.5);

    axis square;
    xlabel(['\it ', subs{1}]);
    ylabel(['\it ', subs{3}]);
end
