%% Equilibria, stability, and threshold-induced growth
% Model:
%   dy/dt = k1*y*z/(z+k2) - c1*y + c3
%   dz/dt = k3*y - c2*z
%
% State vector: u = [y; z]
%
% Outputs:
%   equilibrium_threshold_demo.pdf  (vector graphics)
%   equilibrium_threshold_demo.png  (600 dpi)
%
% MATLAB R2020a or newer recommended.
% Local functions are provided at the end of this script.

clear;
close all;
clc;

%% 1. Parameters

p.k1 = 2.0;
p.c1 = 1.0;
p.k2 = 1.0;
p.k3 = 1.0;
p.c2 = 1.0;

assert(p.k1 > p.c1, 'This example requires k1 > c1.');

a = p.c2 / p.k3;

C3SN = a * p.k2 * ...
    (sqrt(p.k1) - sqrt(p.k1 - p.c1))^2;

zSN = p.k2 * (sqrt(p.k1 / (p.k1 - p.c1)) - 1);
ESN = [a*zSN; zSN];

c3Low  = 0.10;
c3High = 0.25;

assert(c3Low > 0 && c3Low < C3SN);
assert(c3High > C3SN);

%% 2. Equilibria below the threshold

% A*z^2 + B*z + C = 0
A = a * (p.k1 - p.c1);
B = c3Low - a * p.c1 * p.k2;
C = c3Low * p.k2;

Delta = B^2 - 4*A*C;

zStable = (-B - sqrt(Delta)) / (2*A);
zSaddle = (-B + sqrt(Delta)) / (2*A);

EStable = [a*zStable; zStable];
ESaddle = [a*zSaddle; zSaddle];

JStable = modelJacobian(EStable, p);
JSaddle = modelJacobian(ESaddle, p);

lambdaStable = sort(eig(JStable));
lambdaSaddle = sort(eig(JSaddle));

assert(all(lambdaStable < 0), ...
    'The lower equilibrium should be a stable node.');

assert(lambdaSaddle(1) < 0 && lambdaSaddle(2) > 0, ...
    'The upper equilibrium should be a saddle.');

fprintf('\n--- Threshold and equilibria ---\n');
fprintf('C3_SN = %.10f\n', C3SN);
fprintf('E_SN  = (%.10f, %.10f)\n', ESN(1), ESN(2));

fprintf('\nc3_low = %.4f,  c3_low/C3_SN = %.4f\n', ...
    c3Low, c3Low/C3SN);
fprintf('Stable node: (%.10f, %.10f)\n', ...
    EStable(1), EStable(2));
fprintf('Eigenvalues: %.8f, %.8f\n', ...
    lambdaStable(1), lambdaStable(2));

fprintf('\nSaddle: (%.10f, %.10f)\n', ...
    ESaddle(1), ESaddle(2));
fprintf('Eigenvalues: %.8f, %.8f\n', ...
    lambdaSaddle(1), lambdaSaddle(2));

fprintf('\nc3_high = %.4f,  c3_high/C3_SN = %.4f\n', ...
    c3High, c3High/C3SN);
fprintf('No positive equilibrium exists for c3_high.\n\n');

%% 3. Numerical settings

odeOpts = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'MaxStep', 0.10);

% Phase portraits are restricted to a finite display window.
% Stopping at its boundary only clips the displayed trajectory.
boxMax = [1.35; 1.35];

phaseOpts = odeset(odeOpts, ...
    'Events', @(t,u) stopAtBox(t,u,boxMax));

Tphase = 80;
Tlow   = 50;
Thigh  = 45;

% These three states lie in the attraction basin of the stable node.
ICstable = [ ...
    0.04, 0.04;
    0.35, 0.35;
    0.65, 0.65];

% This state lies above the saddle and produces unbounded growth.
ICescape = [1.00; 1.00];

% Above the threshold, include the old stable equilibrium as an IC.
IChigh = [ ...
    EStable.';
    0.04, 0.04;
    0.45, 0.45;
    0.95, 0.95];

%% 4. Plot style

% Colorblind-friendly colors
blue       = [0.000, 0.447, 0.698];
orange     = [0.835, 0.369, 0.000];
green      = [0.000, 0.620, 0.451];
purple     = [0.800, 0.475, 0.655];
darkGray   = [0.25, 0.25, 0.25];
lightGray  = [0.80, 0.82, 0.84];

trajectoryColors = [blue; green; purple];

fig = figure( ...
    'Color', 'w', ...
    'Units', 'centimeters', ...
    'Position', [2, 2, 19, 15.5], ...
    'Renderer', 'painters');

tl = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% 5. Panel (a): phase portrait below threshold

ax1 = nexttile(tl, 1);
hold(ax1, 'on');

drawDirectionField(ax1, p, c3Low, boxMax, lightGray);

% Attracting trajectories
for j = 1:size(ICstable,1)
    u0 = ICstable(j,:).';

    [~,u] = ode45( ...
        @(t,u) modelRHS(t,u,p,c3Low), ...
        [0,Tphase], u0, phaseOpts);

    plot(ax1, u(:,1), u(:,2), ...
        'Color', trajectoryColors(j,:), ...
        'LineWidth', 1.65, ...
        'HandleVisibility', 'off');

    plot(ax1, u0(1), u0(2), 'o', ...
        'Color', trajectoryColors(j,:), ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', 4, ...
        'HandleVisibility', 'off');
end

% Escaping trajectory below the threshold
[~,uEscape] = ode45( ...
    @(t,u) modelRHS(t,u,p,c3Low), ...
    [0,Tphase], ICescape, phaseOpts);

hEscape = plot(ax1, uEscape(:,1), uEscape(:,2), ...
    'Color', orange, 'LineWidth', 1.8);

plot(ax1, ICescape(1), ICescape(2), 'o', ...
    'Color', orange, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', 4, ...
    'HandleVisibility', 'off');

% Approximate the saddle's stable manifold:
% initialize near the saddle along its stable eigenvector and
% integrate backward in time.
[V,D] = eig(JSaddle);
[~,stableIndex] = min(diag(D));

vStable = real(V(:,stableIndex));
vStable = vStable / norm(vStable);

epsilon = 1e-7;

for branch = [-1,1]
    u0 = ESaddle + branch*epsilon*vStable;

    [~,uManifold] = ode45( ...
        @(t,u) modelRHS(t,u,p,c3Low), ...
        [0,-30], u0, phaseOpts);

    plot(ax1, uManifold(:,1), uManifold(:,2), '--', ...
        'Color', darkGray, ...
        'LineWidth', 1.35, ...
        'HandleVisibility', 'off');
end

hManifold = plot(ax1, nan, nan, '--', ...
    'Color', darkGray, 'LineWidth', 1.35);

hStable = plot(ax1, EStable(1), EStable(2), 'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', blue, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.8);

hSaddle = plot(ax1, ESaddle(1), ESaddle(2), 'd', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', orange, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.8);

xlabel(ax1, '$y$', 'Interpreter', 'latex');
ylabel(ax1, '$z$', 'Interpreter', 'latex');

title(ax1, sprintf( ...
    '(a) Below threshold: $c_3=%.2f$', c3Low), ...
    'Interpreter', 'latex');

xlim(ax1, [0,boxMax(1)]);
ylim(ax1, [0,boxMax(2)]);
pbaspect(ax1, [1,1,1]);

legend(ax1, [hStable,hSaddle,hManifold,hEscape], ...
    {'Stable node', 'Saddle', ...
     'Stable manifold', 'Escaping trajectory'}, ...
    'Location', 'northwest', ...
    'FontSize', 8, ...
    'Box', 'off');

styleAxes(ax1);

%% 6. Panel (b): convergence below threshold

ax2 = nexttile(tl, 2);
hold(ax2, 'on');

tLow = linspace(0,Tlow,1200);

for j = 1:size(ICstable,1)
    [t,u] = ode45( ...
        @(t,u) modelRHS(t,u,p,c3Low), ...
        tLow, ICstable(j,:).', odeOpts);

    % Same color = same initial condition
    % Solid = y; dashed = z
    plot(ax2, t, u(:,1), '-', ...
        'Color', trajectoryColors(j,:), ...
        'LineWidth', 1.6, ...
        'HandleVisibility', 'off');

    plot(ax2, t, u(:,2), '--', ...
        'Color', trajectoryColors(j,:), ...
        'LineWidth', 1.6, ...
        'HandleVisibility', 'off');
end

hEq = yline(ax2, EStable(1), ':', ...
    'Color', darkGray, 'LineWidth', 1.3);

hY = plot(ax2, nan, nan, '-', ...
    'Color', darkGray, 'LineWidth', 1.6);

hZ = plot(ax2, nan, nan, '--', ...
    'Color', darkGray, 'LineWidth', 1.6);

xlabel(ax2, '$t$', 'Interpreter', 'latex');
ylabel(ax2, '$y(t),\,z(t)$', 'Interpreter', 'latex');

title(ax2, '(b) Convergence to the stable node', ...
    'Interpreter', 'latex');

xlim(ax2, [0,Tlow]);
ylim(ax2, [0,0.70]);

legend(ax2, [hY,hZ,hEq], ...
    {'$y(t)$', '$z(t)$', '$y^*=z^*$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'FontSize', 9, ...
    'Box', 'off');

styleAxes(ax2);

%% 7. Panel (c): phase portrait above threshold

ax3 = nexttile(tl, 3);
hold(ax3, 'on');

drawDirectionField(ax3, p, c3High, boxMax, lightGray);

highColors = [orange; blue; green; purple];

for j = 1:size(IChigh,1)
    u0 = IChigh(j,:).';

    [~,u] = ode45( ...
        @(t,u) modelRHS(t,u,p,c3High), ...
        [0,Tphase], u0, phaseOpts);

    plot(ax3, u(:,1), u(:,2), ...
        'Color', highColors(j,:), ...
        'LineWidth', 1.7, ...
        'HandleVisibility', 'off');

    plot(ax3, u0(1), u0(2), 'o', ...
        'Color', highColors(j,:), ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', 4, ...
        'HandleVisibility', 'off');
end

% This is an initial condition, not an equilibrium at c3High.
hOldEq = plot(ax3, EStable(1), EStable(2), 's', ...
    'MarkerSize', 6.5, ...
    'MarkerFaceColor', orange, ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.8);

xlabel(ax3, '$y$', 'Interpreter', 'latex');
ylabel(ax3, '$z$', 'Interpreter', 'latex');

title(ax3, sprintf( ...
    '(c) Above threshold: $c_3=%.2f$', c3High), ...
    'Interpreter', 'latex');

xlim(ax3, [0,boxMax(1)]);
ylim(ax3, [0,boxMax(2)]);
pbaspect(ax3, [1,1,1]);

legend(ax3, hOldEq, ...
    {'Initial state: former stable node'}, ...
    'Location', 'northwest', ...
    'FontSize', 8, ...
    'Box', 'off');

styleAxes(ax3);

%% 8. Panel (d): growth after increasing c3 across the threshold

% Initial condition = stable equilibrium for c3Low.
% At t = 0, c3 is changed to c3High and then held constant.
tHigh = linspace(0,Thigh,1600);

[tGrowth,uGrowth] = ode45( ...
    @(t,u) modelRHS(t,u,p,c3High), ...
    tHigh, EStable, odeOpts);

ax4 = nexttile(tl, 4);
hold(ax4, 'on');
set(ax4, 'YScale', 'log');

hGrowthY = plot(ax4, tGrowth, uGrowth(:,1), '-', ...
    'Color', blue, 'LineWidth', 1.9);

hGrowthZ = plot(ax4, tGrowth, uGrowth(:,2), '--', ...
    'Color', orange, 'LineWidth', 1.9);

xlabel(ax4, '$t$', 'Interpreter', 'latex');
ylabel(ax4, '$y(t),\,z(t)$ (log scale)', ...
    'Interpreter', 'latex');

title(ax4, '(d) Growth after crossing the threshold', ...
    'Interpreter', 'latex');

xlim(ax4, [0,Thigh]);

legend(ax4, [hGrowthY,hGrowthZ], ...
    {'$y(t)$', '$z(t)$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northwest', ...
    'FontSize', 9, ...
    'Box', 'off');

styleAxes(ax4);

%% 9. Export

drawnow;

exportgraphics(fig, 'equilibrium_threshold_demo.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

exportgraphics(fig, 'equilibrium_threshold_demo.png', ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

fprintf('\nFigures saved:\n');
fprintf('  equilibrium_threshold_demo.pdf\n');
fprintf('  equilibrium_threshold_demo.png\n');

%% Local functions
% Keep all local functions at the end of this script.

function du = modelRHS(~, u, p, c3)
% Right-hand side of the ODE system.
% u(1) = y; u(2) = z.

    y = u(1);
    z = u(2);

    dy = p.k1*y*z/(z + p.k2) - p.c1*y + c3;
    dz = p.k3*y - p.c2*z;

    du = [dy; dz];
end

function J = modelJacobian(u, p)
% Jacobian with respect to the state variables [y; z].
% The constant source c3 does not contribute to the Jacobian.

    y = u(1);
    z = u(2);

    J = [ ...
        p.k1*z/(z + p.k2) - p.c1, ...
        p.k1*p.k2*y/(z + p.k2)^2; ...
        p.k3, -p.c2 ...
    ];
end

function [value, isterminal, direction] = stopAtBox(~, u, boxMax)
% Stop integration when the trajectory reaches the boundary
% of the displayed phase-plane window:
%   0 <= y <= boxMax(1)
%   0 <= z <= boxMax(2)
%
% direction = 0 allows detection in either integration direction.
% This is useful for backward integration of the stable manifold.

    value = [ ...
        u(1); ...
        u(2); ...
        boxMax(1) - u(1); ...
        boxMax(2) - u(2) ...
    ];

    isterminal = ones(4,1);
    direction  = zeros(4,1);
end

function drawDirectionField(ax, p, c3, boxMax, arrowColor)
% Draw a normalized direction field.
% Arrow lengths indicate direction only, not the speed of motion.

    nGrid = 19;

    % Place arrows slightly inside the plotting window.
    yGrid = linspace(0.025*boxMax(1), 0.975*boxMax(1), nGrid);
    zGrid = linspace(0.025*boxMax(2), 0.975*boxMax(2), nGrid);

    [Y,Z] = meshgrid(yGrid, zGrid);

    dY = p.k1.*Y.*Z./(Z + p.k2) - p.c1.*Y + c3;
    dZ = p.k3.*Y - p.c2.*Z;

    speed = hypot(dY,dZ);

    % Avoid division by zero if a grid point is an equilibrium.
    speed(speed == 0) = 1;

    arrowLength = 0.65 * min( ...
        yGrid(2) - yGrid(1), ...
        zGrid(2) - zGrid(1));

    U = arrowLength .* dY ./ speed;
    V = arrowLength .* dZ ./ speed;

    quiver(ax, Y, Z, U, V, 0, ...
        'Color', arrowColor, ...
        'LineWidth', 0.65, ...
        'MaxHeadSize', 0.8, ...
        'HandleVisibility', 'off');
end

function styleAxes(ax)
% Consistent formatting for all panels.

    set(ax, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'LineWidth', 0.9, ...
        'TickDir', 'out', ...
        'TickLength', [0.015, 0.015], ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'XColor', [0.15, 0.15, 0.15], ...
        'YColor', [0.15, 0.15, 0.15]);

    grid(ax, 'on');

    set(ax, ...
        'GridAlpha', 0.12, ...
        'GridLineStyle', '-', ...
        'MinorGridAlpha', 0.06);

    ax.Title.FontSize = 11;
    ax.Title.FontWeight = 'normal';
    ax.XLabel.FontSize = 11;
    ax.YLabel.FontSize = 11;
end