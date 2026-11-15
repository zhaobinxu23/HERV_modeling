%% Model3_HERV_ageing_figure2.m
%
% Model 3:
% Age-dependent HERV-antibody dynamics
%
% This implementation follows the equations and parameter values given in
% the supplied manuscript.
%
% State variables:
%   V  : HERV burden
%   A  : HERV-specific antibody
%   C  : HERV-antibody immune complex
%   C1 : self-antigen-antibody immune complex
%   S  : self-antigen
%
% Model equations:
%
% dV/dt  = k1*(1 + t/alpha)*V - k2*V*A + k9
%
% dA/dt  = -k2*V*A ...
%          + k4*(1 - t/beta)*C ...
%          + k6*C1 ...
%          - k5*S*A ...
%          - k7*A
%
% dC/dt  = k2*V*A - k3*C
%
% dC1/dt = k5*S*A - k3*C1
%
% dS/dt  = -k5*S*A + k8
%
% The code is intended to provide a publication-quality visualization of
% the model behaviour rather than to reproduce the original graphical
% layout pixel by pixel.

clear;
close all;
clc;

%% ================================================================
% 1. Model parameters from the manuscript
% =================================================================

par = struct();

% HERV-related parameters
par.k1    = 0.1;
par.k2    = 1e-13;
par.k3    = 0.5;
par.k4    = 1.0;

% Self-antigen and antibody parameters
par.k5    = 1e-17;
par.k6    = 0.5;
par.k7    = 1e-3;
par.k8    = 1e10;

% Age-dependent parameters
par.alpha = 200;
par.beta  = 1000;

% Constant HERV burst/input term
par.k9 = 1e4;

%% ================================================================
% 2. Initial conditions from the manuscript
% =================================================================

% State order:
% y = [V, A, C, C1, S]

V0  = 1e10;
A0  = 1e12;
C0  = 2e9;
C10 = 2e10;
S0  = 1e15;

y0 = [V0; A0; C0; C10; S0];

%% ================================================================
% 3. Time interval
% =================================================================

t_start = 0;
t_end   = 500;

% Dense output grid for smooth publication-quality curves
t_eval = linspace(t_start, t_end, 5001);

%% ================================================================
% 4. Numerical integration
% =================================================================

% The variables span several orders of magnitude. ode15s is therefore
% used as a robust solver for the potentially stiff system.
%
% NonNegative prevents small negative concentrations due only to numerical
% round-off and does not change the model equations.

solver_options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', [1e-1, 1e-1, 1e-1, 1e-1, 1e2], ...
    'NonNegative', 1:5, ...
    'MaxStep', 0.5);

[t, y] = ode15s( ...
    @(t, y) model3_rhs(t, y, par), ...
    t_eval, ...
    y0, ...
    solver_options);

%% ================================================================
% 5. Extract variables
% =================================================================

V  = y(:,1);
A  = y(:,2);
C  = y(:,3);
C1 = y(:,4);
S  = y(:,5);

%% ================================================================
% 6. Time-dependent coefficients
% =================================================================

% Effective HERV replication coefficient:
% k1_eff(t) = k1*(1+t/alpha)
k1_eff = par.k1 .* (1 + t ./ par.alpha);

% Effective antibody regeneration coefficient:
% k4_eff(t) = k4*(1-t/beta)
k4_eff = par.k4 .* (1 - t ./ par.beta);

%% ================================================================
% 7. Characteristic points for figure annotations
% =================================================================

% Maximum HERV-specific antibody level
[A_max, i_Amax] = max(A);
t_Amax = t(i_Amax);

% Minimum HERV burden
[V_min, i_Vmin] = min(V);
t_Vmin = t(i_Vmin);

% Maximum HERV-antibody complex concentration
[C_max, i_Cmax] = max(C);
t_Cmax = t(i_Cmax);

% Values at the beginning and end of the simulation
V_initial = V(1);
V_final   = V(end);

A_initial = A(1);
A_final   = A(end);

C_initial = C(1);
C_final   = C(end);

%% ================================================================
% 8. Publication-quality figure
% =================================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [80, 80, 1450, 900]);

tl = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% ----------------------------------------------------------------
% Panel A: self-antigen
% -----------------------------------------------------------------

axA = nexttile(tl, 1);

plot(axA, t, S ./ 1e15, ...
    'Color', [0.00, 0.42, 0.70], ...
    'LineWidth', 2.4);

xlabel(axA, 'Age-associated time, t');
ylabel(axA, 'Self-antigen, S (\times10^{15})');

title(axA, '\bfA\rm  Self-antigen dynamics');

xlim(axA, [0, 500]);
grid(axA, 'on');
box(axA, 'off');
format_axis(axA);

text(axA, 20, 0.90 * max(S ./ 1e15), ...
    'Continuous replenishment of self-antigen', ...
    'FontSize', 9, ...
    'Color', [0.20, 0.20, 0.20], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

%% ----------------------------------------------------------------
% Panel B: HERV burden, full trajectory plus late-stage inset
% -----------------------------------------------------------------

axB = nexttile(tl, 2);

plot(axB, t, V ./ 1e10, ...
    'Color', [0.80, 0.12, 0.08], ...
    'LineWidth', 2.5);

xlabel(axB, 'Age-associated time, t');
ylabel(axB, 'HERV burden, V (\times10^{10})');

title(axB, '\bfB\rm  HERV burden');

xlim(axB, [0, 500]);
grid(axB, 'on');
box(axB, 'off');
format_axis(axB);

% Main-panel annotation
text(axB, 18, 0.88 * max(V ./ 1e10), ...
    'Age-dependent HERV amplification', ...
    'FontSize', 9, ...
    'Color', [0.55, 0.08, 0.06], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

% Inset: late-stage HERV expansion
posB = axB.Position;

axB_inset = axes( ...
    'Parent', fig, ...
    'Position', [ ...
        posB(1) + 0.53 * posB(3), ...
        posB(2) + 0.48 * posB(4), ...
        0.39 * posB(3), ...
        0.38 * posB(4)]);

late = t >= 350;

plot(axB_inset, t(late), V(late) ./ 1e10, ...
    'Color', [0.80, 0.12, 0.08], ...
    'LineWidth', 1.8);

xlabel(axB_inset, 't', 'FontSize', 8);
ylabel(axB_inset, 'V', 'FontSize', 8);

title(axB_inset, 'Late phase', ...
    'FontSize', 8.5, ...
    'FontWeight', 'bold');

xlim(axB_inset, [350, 500]);
grid(axB_inset, 'on');
box(axB_inset, 'on');
format_inset_axis(axB_inset);

%% ----------------------------------------------------------------
% Panel C: HERV-specific antibody
% -----------------------------------------------------------------

axC = nexttile(tl, 3);

plot(axC, t, A ./ 1e12, ...
    'Color', [0.00, 0.30, 0.72], ...
    'LineWidth', 2.5);

hold(axC, 'on');

plot(axC, t_Amax, A_max ./ 1e12, 'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', [0.00, 0.30, 0.72], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.0);

xlabel(axC, 'Age-associated time, t');
ylabel(axC, 'HERV-specific antibody, A (\times10^{12})');

title(axC, '\bfC\rm  HERV-specific antibody dynamics');

xlim(axC, [0, 500]);
grid(axC, 'on');
box(axC, 'off');
format_axis(axC);

text(axC, t_Amax + 14, A_max ./ 1e12, ...
    sprintf('Maximum at t = %.1f', t_Amax), ...
    'FontSize', 9, ...
    'Color', [0.00, 0.20, 0.52], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

text(axC, 20, 0.12 * max(A ./ 1e12), ...
    'Initial compensation followed by late decline', ...
    'FontSize', 9, ...
    'Color', [0.00, 0.20, 0.52], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

%% ----------------------------------------------------------------
% Panel D: HERV-antibody complex, full trajectory plus late inset
% -----------------------------------------------------------------

axD = nexttile(tl, 4);

plot(axD, t, C ./ 1e9, ...
    'Color', [0.47, 0.12, 0.65], ...
    'LineWidth', 2.5);

xlabel(axD, 'Age-associated time, t');
ylabel(axD, 'HERV-antibody complex, C (\times10^{9})');

title(axD, '\bfD\rm  HERV-antibody immune-complex dynamics');

xlim(axD, [0, 500]);
grid(axD, 'on');
box(axD, 'off');
format_axis(axD);

text(axD, 18, 0.88 * max(C ./ 1e9), ...
    'Antigen-antibody engagement', ...
    'FontSize', 9, ...
    'Color', [0.33, 0.07, 0.48], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

% Inset: late-stage immune-complex dynamics
posD = axD.Position;

axD_inset = axes( ...
    'Parent', fig, ...
    'Position', [ ...
        posD(1) + 0.53 * posD(3), ...
        posD(2) + 0.48 * posD(4), ...
        0.39 * posD(3), ...
        0.38 * posD(4)]);

plot(axD_inset, t(late), C(late) ./ 1e9, ...
    'Color', [0.47, 0.12, 0.65], ...
    'LineWidth', 1.8);

xlabel(axD_inset, 't', 'FontSize', 8);
ylabel(axD_inset, 'C', 'FontSize', 8);

title(axD_inset, 'Late phase', ...
    'FontSize', 8.5, ...
    'FontWeight', 'bold');

xlim(axD_inset, [350, 500]);
grid(axD_inset, 'on');
box(axD_inset, 'on');
format_inset_axis(axD_inset);

%% ----------------------------------------------------------------
% Overall title and explanatory annotation
% -----------------------------------------------------------------

title(tl, ...
    ['Figure 2. Model 3 predicts age-associated loss of ' ...
     'HERV-specific immune surveillance'], ...
    'FontName', 'Arial', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

annotation(fig, 'textbox', ...
    [0.16, 0.005, 0.68, 0.032], ...
    'String', ...
    ['Ageing increases HERV replication through k_1(1+t/\alpha) ' ...
     'while reducing immune-complex-mediated antibody regeneration ' ...
     'through k_4(1-t/\beta).'], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', 'Arial', ...
    'FontSize', 9.5, ...
    'Color', [0.20, 0.20, 0.20], ...
    'EdgeColor', 'none');

%% ================================================================
% 9. Export files
% =================================================================

exportgraphics(fig, ...
    'Model3_Figure2_HERV_ageing.png', ...
    'Resolution', 600);

exportgraphics(fig, ...
    'Model3_Figure2_HERV_ageing.pdf', ...
    'ContentType', 'vector');

savefig(fig, ...
    'Model3_Figure2_HERV_ageing.fig');

save('Model3_Figure2_HERV_ageing_data.mat', ...
    't', 'V', 'A', 'C', 'C1', 'S', ...
    'k1_eff', 'k4_eff', ...
    'par', 'y0');

%% ================================================================
% 10. Numerical summary
% =================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('Model 3: HERV ageing dynamics\n');
fprintf('====================================================\n');

fprintf('Initial HERV burden:          %.6e\n', V_initial);
fprintf('Final HERV burden:            %.6e\n', V_final);

fprintf('Initial antibody level:       %.6e\n', A_initial);
fprintf('Maximum antibody level:       %.6e\n', A_max);
fprintf('Time of antibody maximum:     %.4f\n', t_Amax);
fprintf('Final antibody level:         %.6e\n', A_final);

fprintf('Initial immune complex:       %.6e\n', C_initial);
fprintf('Maximum immune complex:       %.6e\n', C_max);
fprintf('Time of complex maximum:      %.4f\n', t_Cmax);
fprintf('Final immune complex:         %.6e\n', C_final);

fprintf('Effective k1 at t = 0:        %.6f\n', k1_eff(1));
fprintf('Effective k1 at t = 500:      %.6f\n', k1_eff(end));

fprintf('Effective k4 at t = 0:        %.6f\n', k4_eff(1));
fprintf('Effective k4 at t = 500:      %.6f\n', k4_eff(end));

fprintf('====================================================\n');

%% ================================================================
% Local function: Model 3 ODE system
% =================================================================

function dydt = model3_rhs(t, y, par)

    % State variables
    V  = y(1);
    A  = y(2);
    C  = y(3);
    C1 = y(4);
    S  = y(5);

    % Age-dependent coefficients specified in the manuscript
    k1_t = par.k1 * (1 + t / par.alpha);
    k4_t = par.k4 * (1 - t / par.beta);

    % Model 3 equations
    dV = ...
        k1_t * V ...
        - par.k2 * V * A ...
        + par.k9;

    dA = ...
        - par.k2 * V * A ...
        + k4_t * C ...
        + par.k6 * C1 ...
        - par.k5 * S * A ...
        - par.k7 * A;

    dC = ...
        par.k2 * V * A ...
        - par.k3 * C;

    dC1 = ...
        par.k5 * S * A ...
        - par.k3 * C1;

    dS = ...
        - par.k5 * S * A ...
        + par.k8;

    dydt = [dV; dA; dC; dC1; dS];

end

%% ================================================================
% Local function: main-axis formatting
% =================================================================

function format_axis(ax)

    ax.FontName = 'Arial';
    ax.FontSize = 10.5;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.TickLength = [0.015, 0.015];
    ax.Layer = 'top';
    ax.GridAlpha = 0.14;
    ax.MinorGridAlpha = 0.06;

end

%% ================================================================
% Local function: inset-axis formatting
% =================================================================

function format_inset_axis(ax)

    ax.FontName = 'Arial';
    ax.FontSize = 8;
    ax.LineWidth = 0.7;
    ax.TickDir = 'out';
    ax.TickLength = [0.025, 0.025];
    ax.Layer = 'top';
    ax.GridAlpha = 0.16;

end