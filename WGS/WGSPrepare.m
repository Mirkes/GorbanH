% Prepare all data for WGS system
% Value of modified MUST be set before start of this script. Acceptable
% values of this variable are described in script WGSFigures.

if modified == 0
    % Original system
    % Suffix for file names
    suf = '';
    % Substances
    subst = {'H_2O', 'H_2', 'CO', 'CO_2', 'red', 'Ox'};
    % Define known parameters
    % Fraction of non converted CO in equilibrium is 1 - degree of conversion
    % (87%)
    coeq = 1 - 0.87;
    % Time of going through reactor
    tt = 0.59;
    % Fraction of non converted CO at time tt is 1 - degree of conversion (70%)
    cott = 1 - 0.7;
    % Balances:
    %   1 for H2 (arbitrary choice)
    %   b(H2) for C (proportion of H2O and CO is 1:1)
    %   b(H2) + b(CO) (H2O and CO contains one atom of oxygen each)
    %   b(4) is unknown. Expected that this value is essentially less than b(H2).
    b = [1, 1, 2, 0.1];

    % Set of parameters
    %       param(1) is b(4): total number of catalyst
    %       param(2) is equilibrium value of free catalyst
    %       param(3) is reaction rate constants for the first inverse reaction
    %       param(4) is reaction rate constants for the second inverse reaction

    % Parameters can be saved in the WGSParamFound.mat file or can be found by
    % Nelde-Mead search. This search is not very quick and can require several
    % hours.

    if  ~exist('WGSParamFound.mat','file')
        % There is no previously found parameters. Implement paramaters search
        % Initial set
        param = [0.1, 0.001, 100, 100];
        % Manually check of reasonability/plausability of results for specified parameters
        % res = integrWGS(param, coeq, cott, tt, b);

        % Search optimal values of parameters
        % To observe process of search
        options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);
        [x, fval, exitflag, output] = fminsearch(@(x) integrWGS(x, coeq, cott, tt, b),...
                                                 param, options);
        % Round param for reasonable precision
        param = round(x, 5, 'significant');
        % To manual check the rounded parameters quality
        % res = integrWGS(x, coeq, cott, tt, b);
        % res1 = integrWGS(param, coeq, cott, tt, b);

        % Save calculated parameters
        save('WGSParamFound.mat', 'param');
    else
        load('WGSParamFound.mat');
    end
    
    % Apply found value of b(4)
    b(4) = param(1);
    
    % If necessary then calculate tragectories with different length of time
    % interval for different purposes.
    % This procedure reqires several munutes only.
    % Time is the first column of cFull

    % To draw complete trajectory
    if  ~exist('WGSResTrajectory.mat','file')
        % Calculate trajectory with 590 equidistant steps in time for all variables
        [cFull, eq] = integrWGSDetail(param, coeq, tt * 3, b);
        save('WGSResTrajectory.mat', 'cFull');
    else
        load('WGSResTrajectory.mat');
    end

    % To draw trajectory in the reactor
    if  ~exist('WGSRes.mat','file')
        % Calculate trajectory with 590 equidistant steps in time for all variables
        [cReact, eq] = integrWGSDetail(param, coeq, tt, b);
        save('WGSRes.mat', 'cReact');
    else
        load('WGSRes.mat');
    end

    % To draw initial process in details
    if  ~exist('WGSResShortTime.mat','file')
        % Calculate trajectory with 590 equidistant steps in time for all variables
        [cShort, eq] = integrWGSDetail(param, coeq, tt * 1.e-5, b);
        save('WGSResShortTime.mat', 'cShort');
    else
        load('WGSResShortTime.mat');
    end
elseif modified == 1
    % Modified system
    % Suffix for file names
    suf = 'Mod';
    % Substances
    subst = {'A_1', 'A_2', 'A_3', 'A_4', 'A_5', 'A_6'};
    % All parameters are known
    % Set of parameters
    %       param(1) is b(4): total number of catalyst
    %       param(2) is equilibrium value of free catalyst
    %       param(3) is reaction rate constants for the first inverse reaction
    %       param(4) is reaction rate constants for the second inverse reaction
    param = [0.5, 0.25, 1, 1];

    % Balances:
    %   1 for H2 (arbitrary choice)
    %   b(H2) for C (proportion of H2O and CO is 1:1)
    %   b(H2) + b(CO) (H2O and CO contains one atom of oxygen each)
    %   Found in calculations.
    b = [1, 1, 2, param(1)];

    % Fraction of non converted CO in equilibrium is 1 - degree of conversion
    % (87%)
    coeq = 0.5;

    tt = 3.2;

    % If necessary then calculate tragectories with different length of time
    % interval for different purposes.
    % This procedure reqires several munutes only.
    % Time is the first column of cFull

    % To draw complete trajectory
    % Calculate trajectory with 590 equidistant steps in time for all variables
    [cFull, eq] = integrWGSDetail(param, coeq, tt * 3, b);
    save(sprintf('WGS%sResTrajectory.mat', suf), 'cFull');

    % To draw trajectory in the reactor
    % Calculate trajectory with 590 equidistant steps in time for all variables
    [cReact, eq] = integrWGSDetail(param, coeq, tt, b);
    save(sprintf('WGS%sRes.mat', suf), 'cReact');

    % To draw initial process in details
    % Calculate trajectory with 590 equidistant steps in time for all variables
    [cShort, eq] = integrWGSDetail(param, coeq, tt * 0.5, b);
    save(sprintf('WGS%sResShortTime.mat', suf), 'cShort');
end