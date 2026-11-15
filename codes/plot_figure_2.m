%% four_panel_p62_simulation.m
%
% Four simulations:
%
% A. Normal cells + ERV gene:
%       p1 = 10000, p2 = 0.04, p3 = 0.05,
%       p56 = 0.1, p62 = 1
%
% B. Normal cells + conventional gene:
%       p1 = 10000, p2 = 0.04, p3 = 0.05,
%       p56 = 0.1, p62 = 0
%
% C. Cancer cells + conventional gene:
%       p1 = 0, p2 = 0.5, p3 = 0.05,
%       p56 = 10, p62 = 0
%
% D. Cancer cells + ERV gene:
%       p1 = 0, p2 = 0.5, p3 = 0.05,
%       p56 = 10, p62 = 1
%
% Required files:
%   main_test_herv.m
%   tumor_immune_interaction_herv.m
%
% The ODE function must have the form:
%   dydt = tumor_immune_interaction_herv(t,y,par)

clear;
clc;
close all;

%% Check required files

if exist('main_test_herv.m','file') ~= 2
    error(['main_test_herv.m was not found. ', ...
        'Place this script in the same folder as the model files.']);
end

if exist('tumor_immune_interaction_herv.m','file') ~= 2
    error(['tumor_immune_interaction_herv.m was not found. ', ...
        'Place this script in the same folder as the model files.']);
end

%% Initialize the original model

% main_test_herv.m should create:
%   par : parameter structure
%   y0  : 36-dimensional initial state vector
run('main_test_herv.m');

if ~exist('par','var')
    error('main_test_herv.m did not create the parameter structure par.');
end

if ~exist('y0','var')
    error('main_test_herv.m did not create the initial state y0.');
end

if numel(y0) ~= 36
    error('The model must contain 36 state variables.');
end

y0 = y0(:);

% Preserve all original parameter values, including p58.
parOriginal = par;
y0Original = y0;

%% General settings

% Normal-cell reference population in the absence of immune killing:
%
% T_ref = p1 / (p3 - p2)
%
% Here:
% T_ref = 10000/(0.05-0.04) = 1,000,000
normalReference = 10000/(0.05 - 0.04);

% Low expression level.
normalP56 = 0.1;

% Numerical tumor-clearance threshold.
clearanceThreshold = 1e-4;

% Long simulation for normal-cell dynamics.
tspanNormal = [0 10000];

% Cancer-cell simulation interval.
tspanCancer = [0 100];

%% Solver settings

% The model contains variables with very different scales.
% Component-wise tolerances improve numerical robustness.
absoluteTolerance = 1e-10 * max(ones(36,1),abs(y0Original));

solverOptions = odeset( ...
    'RelTol',1e-7, ...
    'AbsTol',absoluteTolerance, ...
    'MaxStep',1);

%% Define four simulation conditions

condition(1).label = 'A  Normal cells + ERV gene';
condition(1).p1 = 10000;
condition(1).p2 = 0.04;
condition(1).p3 = 0.05;
condition(1).p56 = normalP56;
condition(1).p62 = 1;
condition(1).tspan = tspanNormal;
condition(1).initialPopulation = normalReference;
condition(1).isNormal = true;

condition(2).label = 'B  Normal cells + conventional gene';
condition(2).p1 = 10000;
condition(2).p2 = 0.04;
condition(2).p3 = 0.05;
condition(2).p56 = normalP56;
condition(2).p62 = 0;
condition(2).tspan = tspanNormal;
condition(2).initialPopulation = normalReference;
condition(2).isNormal = true;

condition(3).label = 'C  Cancer cells + conventional gene';
condition(3).p1 = 0;
condition(3).p2 = 0.5;
condition(3).p3 = 0.05;
condition(3).p56 = 10;
condition(3).p62 = 0;
condition(3).tspan = tspanCancer;
condition(3).initialPopulation = 100;
condition(3).isNormal = false;

condition(4).label = 'D  Cancer cells + ERV gene';
condition(4).p1 = 0;
condition(4).p2 = 0.5;
condition(4).p3 = 0.05;
condition(4).p56 = 10;
condition(4).p62 = 1;
condition(4).tspan = tspanCancer;
condition(4).initialPopulation = 100;
condition(4).isNormal = false;

%% Run simulations

Results = struct([]);

for i = 1:4

    % Start from the original parameter set.
    parCase = parOriginal;

    % Important:
    % p58 is deliberately preserved from the original model.
    parCase.p58 = parOriginal.p58;

    % p62 is the only gene-type switch.
    parCase.p62 = condition(i).p62;

    % Set cellular and transcription parameters.
    parCase.p1 = condition(i).p1;
    parCase.p2 = condition(i).p2;
    parCase.p3 = condition(i).p3;
    parCase.p56 = condition(i).p56;

    % Copy the original initial condition.
    y0Case = y0Original;

    % Set the initial cell population.
    y0Case(1) = condition(i).initialPopulation;

    % Solve the model.
    [t,Y] = ode15s( ...
        @(t,y) tumor_immune_interaction_herv(t,y,parCase), ...
        condition(i).tspan, ...
        y0Case, ...
        solverOptions);

    % Numerical validation.
    if any(~isfinite(Y(:)))
        error('NaN or Inf detected in condition %s.', ...
            condition(i).label);
    end

    if any(Y(:,1) < -1e-8)
        warning('Negative cell values detected in condition %s.', ...
            condition(i).label);
    end

    % Store the result.
Results(i).label = condition(i).label;
Results(i).p1 = condition(i).p1;
Results(i).p2 = condition(i).p2;
Results(i).p3 = condition(i).p3;
Results(i).p56 = condition(i).p56;
Results(i).p58 = parCase.p58;
Results(i).p62 = condition(i).p62;

% 保存该条件是否为正常细胞模拟
Results(i).isNormal = condition(i).isNormal;

Results(i).t = t;
Results(i).Y = Y;
Results(i).par = parCase;

    % Cell population.
    Results(i).T = Y(:,1);

    % RNA per cell.
    validPopulation = Results(i).T > 10*clearanceThreshold;
    Results(i).RNAperCell = nan(size(Results(i).T));

    Results(i).RNAperCell(validPopulation) = ...
        max(Y(validPopulation,31),0) ./ ...
        Results(i).T(validPopulation);

    % Simulation summary.
    Results(i).initialPopulation = Results(i).T(1);
    Results(i).finalPopulation = Results(i).T(end);
    Results(i).peakPopulation = max(Results(i).T);

    Results(i).intrinsicGrowth = ...
        condition(i).p2 - condition(i).p3;
end

%% Print numerical summary

fprintf('\n');
fprintf('====================================================\n');
fprintf('Four-condition p62 simulation summary\n');
fprintf('====================================================\n');

for i = 1:4

    fprintf('\n%s\n',Results(i).label);
    fprintf('  p1  = %.6g\n',Results(i).p1);
    fprintf('  p2  = %.6g\n',Results(i).p2);
    fprintf('  p3  = %.6g\n',Results(i).p3);
    fprintf('  p56 = %.6g\n',Results(i).p56);
    fprintf('  p58 = %.6g\n',Results(i).p58);
    fprintf('  p62 = %.6g\n',Results(i).p62);
    fprintf('  Initial cell population = %.6g\n', ...
        Results(i).initialPopulation);
    fprintf('  Final cell population   = %.6g\n', ...
        Results(i).finalPopulation);
    fprintf('  Peak cell population    = %.6g\n', ...
        Results(i).peakPopulation);
end

fprintf('\n');
fprintf('Normal-cell reference population = %.6g\n', ...
    normalReference);
fprintf('Cancer intrinsic growth rate     = %.6g\n', ...
    0.5 - 0.05);
fprintf('====================================================\n');

%% Create four-panel figure

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 19 15]);

tileLayout = tiledlayout(fig,2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Colors
ervColor = [0.82 0.25 0.12];
conventionalColor = [0.10 0.32 0.65];

for i = 1:4

    ax = nexttile(tileLayout,i);

    T = Results(i).T;
    t = Results(i).t;

    if Results(i).p62 == 1
        curveColor = ervColor;
    else
        curveColor = conventionalColor;
    end

    if Results(i).isNormal
        plottedPopulation = T;
    else
        % Avoid log(0) in case the model reaches exactly zero.
        plottedPopulation = max(T,clearanceThreshold/10);
    end

    if Results(i).isNormal

        plot(ax,t,plottedPopulation, ...
            'Color',curveColor, ...
            'LineWidth',2);

        yline(ax,1,':', ...
            'Color',[0.35 0.35 0.35], ...
            'LineWidth',1);

        ylim(ax,'auto');

    else

        semilogy(ax,t,plottedPopulation, ...
            'Color',curveColor, ...
            'LineWidth',2);

        yline(ax,clearanceThreshold,':', ...
            'Color',[0.35 0.35 0.35], ...
            'LineWidth',1);

    end

    xlabel(ax,'Time (model units)');
    ylabel(ax,'Cell population');

    title(ax,Results(i).label, ...
        'FontName','Arial', ...
        'FontWeight','bold', ...
        'FontSize',11);

    set(ax, ...
        'FontName','Arial', ...
        'FontSize',10, ...
        'LineWidth',0.8, ...
        'Box','off', ...
        'TickDir','out');

    grid(ax,'on');
    ax.GridAlpha = 0.15;

    if Results(i).isNormal
        xlim(ax,[0 10000]);
        ylabel(ax,'Normal cells / reference population');
    else
        xlim(ax,[0 100]);
        ylabel(ax,'Malignant cells, log scale');
    end
end

title(tileLayout, ...
    'Cell-population dynamics under conventional and ERV-associated expression', ...
    'FontName','Arial', ...
    'FontWeight','bold', ...
    'FontSize',12);

%% Export figure and data

outputFolder = 'four_panel_p62_results';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

exportgraphics(fig, ...
    fullfile(outputFolder,'Figure_p62_four_panel.pdf'), ...
    'ContentType','vector');

exportgraphics(fig, ...
    fullfile(outputFolder,'Figure_p62_four_panel.png'), ...
    'Resolution',600);

save(fullfile(outputFolder,'Figure_p62_simulation_data.mat'), ...
    'Results', ...
    'normalReference', ...
    'clearanceThreshold');

fprintf('\nFigure and simulation data saved in:\n');
fprintf('  %s\n',fullfile(pwd,outputFolder));