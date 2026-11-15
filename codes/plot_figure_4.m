%% ========================================================================
%  Model3_Figure3_vaccine.m
%
%  Model 3: HERV-antibody dynamics with HERV-antigen vaccination
%
%  This script simulates the vaccine-extended Model 3 and generates a
%  publication-quality Figure 3.
%
%  State variables:
%       V   : endogenous retrovirus/HERV burden
%       V1  : vaccine antigen
%       A   : antibody concentration
%       C   : HERV-antibody immune complex
%       C1  : self-antigen-antibody immune complex
%       C2  : vaccine-antigen-antibody immune complex
%       S   : self-antigen concentration
%
%  State-vector order:
%       y = [V, V1, A, C, C1, C2, S]
%
%  The parameter values are taken from the supplied attachment:
%       V1(0) = 1e13
%       k10   = 0.05
%
%  The plotting layout is redesigned for clearer SCI-style presentation.
%  It does not attempt pixel-by-pixel reproduction of the original figure.
%
%  MATLAB version recommended: R2020a or later.
% ========================================================================

clear;
close all;
clc;

%% ========================================================================
% 1. Model parameters
% ========================================================================

par = struct();

% HERV and antibody-interaction parameters
par.k1     = 0.1;       % basal HERV amplification coefficient
par.k2     = 1e-13;     % antigen-antibody binding coefficient
par.k3     = 0.5;       % immune-complex clearance coefficient
par.k4     = 1.0;       % immune-complex-mediated antibody regeneration
par.k5     = 1e-17;     % self-antigen-antibody binding coefficient
par.k6     = 0.5;       % C1-mediated antibody regeneration
par.k7     = 1e-3;      % antibody degradation coefficient
par.k8     = 1e10;      % self-antigen replenishment rate

% Age-dependent parameters
par.alpha  = 200;       % HERV amplification timescale
par.beta   = 1000;      % antibody-regeneration timescale

% Continuous HERV burst/input term
par.k9     = 1e4;

% Vaccine-specific parameter from the attachment
par.k10    = 0.05;      % vaccine-antigen clearance coefficient

%% ========================================================================
% 2. Initial conditions
% ========================================================================

% Initial conditions from the attachment
V0  = 1e10;             % initial HERV burden
V10 = 1e13;             % initial vaccine-antigen dose
A0  = 1e12;             % initial antibody concentration
C0  = 2e9;              % initial HERV-antibody complex
C10 = 2e10;             % initial self-antigen-antibody complex
C20 = 0;                % initial vaccine-antibody complex
S0  = 1e15;             % initial self-antigen concentration

y0 = [V0; V10; A0; C0; C10; C20; S0];

%% ========================================================================
% 3. Simulation time
% ========================================================================

t_start = 0;
t_end   = 500;

% Dense output grid for smooth curves and high-resolution export
t_eval = linspace(t_start, t_end, 5001).';

%% ========================================================================
% 4. Numerical solver settings
% ========================================================================

% Because the model contains variables spanning many orders of magnitude,
% ode15s is used as a robust solver for this potentially stiff system.
%
% The NonNegative option prevents numerical round-off from producing small
% negative concentrations. It does not alter the intended reaction model.

solver_options = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e2], ...
    'NonNegative', 1:7, ...
    'MaxStep', 0.5);

%% ========================================================================
% 5. Solve the vaccine model
% ========================================================================

[t, y] = ode15s( ...
    @(t, y) model3_vaccine_rhs(t, y, par), ...
    t_eval, ...
    y0, ...
    solver_options);

%% ========================================================================
% 6. Extract state variables
% ========================================================================

V  = y(:,1);
V1 = y(:,2);
A  = y(:,3);
C  = y(:,4);
C1 = y(:,5);
C2 = y(:,6);
S  = y(:,7);

%% ========================================================================
% 7. Effective age-dependent coefficients
% ========================================================================

k1_eff = par.k1 .* (1 + t ./ par.alpha);
k4_eff = par.k4 .* (1 - t ./ par.beta);

%% ========================================================================
% 8. Derived quantities for analysis and annotations
% ========================================================================

% Total antigen burden relevant to antibody binding
antigen_total = V + V1;

% Total immune-complex burden involving HERV or vaccine antigen
effector_complex = C + C2;

% Antibody-regeneration contribution through HERV and vaccine complexes
antibody_regeneration = k4_eff .* (C + C2);

% Characteristic points
[A_max, i_Amax] = max(A);
t_Amax = t(i_Amax);

[V_max, i_Vmax] = max(V);
t_Vmax = t(i_Vmax);

[V1_max, i_V1max] = max(V1);
t_V1max = t(i_V1max);

[C2_max, i_C2max] = max(C2);
t_C2max = t(i_C2max);

% Values at selected time points
idx0   = 1;
idx100 = find(t >= 100, 1, 'first');
idx250 = find(t >= 250, 1, 'first');
idx500 = length(t);

%% ========================================================================
% 9. Print numerical summary
% ========================================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('Model 3: HERV-antibody dynamics with vaccine antigen\n');
fprintf('================================================================\n');
fprintf('Simulation interval:             t = %.1f to %.1f\n', t_start, t_end);
fprintf('\n');
fprintf('Initial HERV burden V(0):        %.6e\n', V(1));
fprintf('Final HERV burden V(500):        %.6e\n', V(end));
fprintf('\n');
fprintf('Initial vaccine antigen V1(0):   %.6e\n', V1(1));
fprintf('Final vaccine antigen V1(500):  %.6e\n', V1(end));
fprintf('Maximum vaccine antigen:         %.6e\n', V1_max);
fprintf('Time of vaccine-antigen maximum: %.4f\n', t_V1max);
fprintf('\n');
fprintf('Initial antibody A(0):           %.6e\n', A(1));
fprintf('Maximum antibody:                %.6e\n', A_max);
fprintf('Time of antibody maximum:        %.4f\n', t_Amax);
fprintf('Final antibody A(500):           %.6e\n', A(end));
fprintf('\n');
fprintf('Maximum HERV-antibody complex:   %.6e\n', max(C));
fprintf('Maximum vaccine complex C2:      %.6e\n', C2_max);
fprintf('Time of maximum vaccine complex: %.4f\n', t_C2max);
fprintf('\n');
fprintf('k1_eff(0):                       %.6f\n', k1_eff(1));
fprintf('k1_eff(500):                     %.6f\n', k1_eff(end));
fprintf('k4_eff(0):                       %.6f\n', k4_eff(1));
fprintf('k4_eff(500):                     %.6f\n', k4_eff(end));
fprintf('================================================================\n');

%% ========================================================================
% 10. Create publication-quality Figure 3
% ========================================================================

fig = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [80, 60, 1500, 930]);

layout = tiledlayout(fig, 2, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% ------------------------------------------------------------------------
% Panel A: HERV burden
% -------------------------------------------------------------------------

axA = nexttile(layout, 1);

plot(axA, t, V ./ 1e10, ...
    'Color', [0.80, 0.12, 0.08], ...
    'LineWidth', 2.3);

xlabel(axA, 'Age-associated time, t');
ylabel(axA, 'HERV burden, V (\times10^{10})');
title(axA, '\bfA\rm  HERV burden');

xlim(axA, [0, 500]);
grid(axA, 'on');
box(axA, 'off');
format_main_axis(axA);

text(axA, 18, 0.88 * max(V ./ 1e10), ...
    'Age-dependent HERV amplification', ...
    'FontSize', 9, ...
    'Color', [0.55, 0.08, 0.06], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

%% ------------------------------------------------------------------------
% Panel B: vaccine antigen
% -------------------------------------------------------------------------

axB = nexttile(layout, 2);

plot(axB, t, V1 ./ 1e13, ...
    'Color', [0.10, 0.50, 0.28], ...
    'LineWidth', 2.3);

xlabel(axB, 'Age-associated time, t');
ylabel(axB, 'Vaccine antigen, V_1 (\times10^{13})');
title(axB, '\bfB\rm  Vaccine-antigen decay');

xlim(axB, [0, 500]);
grid(axB, 'on');
box(axB, 'off');
format_main_axis(axB);

text(axB, 18, 0.82 * max(V1 ./ 1e13), ...
    'V_1(0) = 10^{13}', ...
    'FontSize', 9, ...
    'Color', [0.06, 0.33, 0.18], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

%% ------------------------------------------------------------------------
% Panel C: antibody concentration
% -------------------------------------------------------------------------

axC = nexttile(layout, 3);

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
ylabel(axC, 'Antibody, A (\times10^{12})');
title(axC, '\bfC\rm  Antibody dynamics');

xlim(axC, [0, 500]);
grid(axC, 'on');
box(axC, 'off');
format_main_axis(axC);

text(axC, t_Amax + 12, A_max ./ 1e12, ...
    sprintf('Peak: t = %.1f', t_Amax), ...
    'FontSize', 9, ...
    'Color', [0.00, 0.20, 0.52], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

%% ------------------------------------------------------------------------
% Panel D: HERV-antibody and vaccine-antibody complexes
% -------------------------------------------------------------------------

axD = nexttile(layout, 4);

plot(axD, t, C ./ 1e9, ...
    'Color', [0.49, 0.16, 0.63], ...
    'LineWidth', 2.2, ...
    'DisplayName', 'HERV-antibody complex, C');
hold(axD, 'on');

plot(axD, t, C2 ./ 1e9, ...
    'Color', [0.90, 0.47, 0.05], ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Vaccine-antibody complex, C_2');

xlabel(axD, 'Age-associated time, t');
ylabel(axD, 'Immune complexes (\times10^{9})');
title(axD, '\bfD\rm  Antigen-antibody complexes');

xlim(axD, [0, 500]);
grid(axD, 'on');
box(axD, 'off');
format_main_axis(axD);
legend(axD, 'Location', 'best', 'Box', 'off');

%% ------------------------------------------------------------------------
% Panel E: total antigen and total effector complex
% -------------------------------------------------------------------------

axE = nexttile(layout, 5);

yyaxis(axE, 'left');
plot(axE, t, antigen_total ./ 1e13, ...
    'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 2.0, ...
    'DisplayName', 'V + V_1');
ylabel(axE, 'Total antigen, V+V_1 (\times10^{13})');

axE.YColor = [0.20, 0.20, 0.20];

yyaxis(axE, 'right');
plot(axE, t, effector_complex ./ 1e9, ...
    'Color', [0.10, 0.58, 0.55], ...
    'LineWidth', 2.0, ...
    'DisplayName', 'C+C_2');
ylabel(axE, 'Effector complexes, C+C_2 (\times10^{9})');

axE.YColor = [0.10, 0.45, 0.43];
xlabel(axE, 'Age-associated time, t');
title(axE, '\bfE\rm  Antigenic load and immune-complex formation');

xlim(axE, [0, 500]);
grid(axE, 'on');
box(axE, 'off');
format_main_axis(axE);
legend(axE, 'Location', 'best', 'Box', 'off');

%% ------------------------------------------------------------------------
% Panel F: self-antigen and age-dependent coefficients
% -------------------------------------------------------------------------

axF = nexttile(layout, 6);

yyaxis(axF, 'left');
plot(axF, t, S ./ 1e15, ...
    'Color', [0.00, 0.42, 0.70], ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Self-antigen S');
ylabel(axF, 'Self-antigen, S (\times10^{15})');
axF.YColor = [0.00, 0.30, 0.58];

yyaxis(axF, 'right');
plot(axF, t, k4_eff, '--', ...
    'Color', [0.82, 0.20, 0.18], ...
    'LineWidth', 1.9, ...
    'DisplayName', 'k_{4,eff}(t)');
ylabel(axF, 'Effective regeneration coefficient');
axF.YColor = [0.65, 0.10, 0.08];

xlabel(axF, 'Age-associated time, t');
title(axF, '\bfF\rm  Self-antigen and immune-regeneration capacity');

xlim(axF, [0, 500]);
grid(axF, 'on');
box(axF, 'off');
format_main_axis(axF);
legend(axF, 'Location', 'best', 'Box', 'off');

%% ------------------------------------------------------------------------
% Overall title and explanatory annotation
% -------------------------------------------------------------------------

title(layout, ...
    ['Figure 3. Vaccine-associated enhancement of HERV-specific ' ...
     'immune-complex formation in Model 3'], ...
    'FontName', 'Arial', ...
    'FontSize', 15.5, ...
    'FontWeight', 'bold');

annotation(fig, 'textbox', ...
    [0.12, 0.005, 0.76, 0.032], ...
    'String', ...
    ['The vaccine antigen V_1 is initialized at 10^{13}. ' ...
     'It is cleared with k_{10}=0.05 and binds antibody to form C_2, ' ...
     'thereby contributing to immune-complex-mediated antibody regeneration.'], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', 'Arial', ...
    'FontSize', 9.2, ...
    'Color', [0.20, 0.20, 0.20], ...
    'EdgeColor', 'none');

%% ========================================================================
% 11. Export figure and numerical data
% ========================================================================

exportgraphics(fig, ...
    'Figure3_Model3_vaccine.png', ...
    'Resolution', 600);

exportgraphics(fig, ...
    'Figure3_Model3_vaccine.pdf', ...
    'ContentType', 'vector');

savefig(fig, 'Figure3_Model3_vaccine.fig');

save('Figure3_Model3_vaccine_data.mat', ...
    't', 'y', ...
    'V', 'V1', 'A', 'C', 'C1', 'C2', 'S', ...
    'antigen_total', 'effector_complex', ...
    'antibody_regeneration', ...
    'k1_eff', 'k4_eff', ...
    'par', 'y0');

fprintf('\nFiles exported successfully:\n');
fprintf('  Figure3_Model3_vaccine.png\n');
fprintf('  Figure3_Model3_vaccine.pdf\n');
fprintf('  Figure3_Model3_vaccine.fig\n');
fprintf('  Figure3_Model3_vaccine_data.mat\n');

%% ========================================================================
% Local function 1: Model 3 vaccine-extended ODE system
% ========================================================================

function dydt = model3_vaccine_rhs(t, y, par)

    % ---------------------------------------------------------------------
    % State variables
    % ---------------------------------------------------------------------
    V  = y(1);     % HERV burden
    V1 = y(2);     % vaccine antigen
    A  = y(3);     % antibody
    C  = y(4);     % HERV-antibody complex
    C1 = y(5);     % self-antigen-antibody complex
    C2 = y(6);     % vaccine-antibody complex
    S  = y(7);     % self-antigen

    % ---------------------------------------------------------------------
    % Age-dependent coefficients
    % ---------------------------------------------------------------------
    % Age-associated HERV amplification
    k1_t = par.k1 * (1 + t / par.alpha);

    % Age-associated immune-regeneration capacity
    k4_t = par.k4 * (1 - t / par.beta);

    % ---------------------------------------------------------------------
    % Model equations
    % ---------------------------------------------------------------------

    % HERV burden
    dV = ...
        k1_t * V ...
        - par.k2 * V * A ...
        + par.k9;

    % Vaccine antigen
    % The vaccine pool is reduced by intrinsic clearance and antibody
    % binding. The initial vaccine dose is represented by V1(0)=1e13.
    dV1 = ...
        - par.k10 * V1 ...
        - par.k2 * V1 * A;

    % Antibody
    % Antibody is consumed by both HERV and vaccine-antigen binding,
    % regenerated through C and C2, regenerated through C1, consumed by
    % self-antigen binding, and removed by degradation.
    dA = ...
        - par.k2 * V1 * A ...
        - par.k2 * V * A ...
        + k4_t * (C + C2) ...
        + par.k6 * C1 ...
        - par.k5 * S * A ...
        - par.k7 * A;

    % HERV-antibody complex
    dC = ...
        par.k2 * V * A ...
        - par.k3 * C;

    % Self-antigen-antibody complex
    dC1 = ...
        par.k5 * S * A ...
        - par.k3 * C1;

    % Vaccine-antibody complex
    dC2 = ...
        par.k2 * V1 * A ...
        - par.k3 * C2;

    % Self-antigen
    dS = ...
        - par.k5 * S * A ...
        + par.k8;

    dydt = [dV; dV1; dA; dC; dC1; dC2; dS];

end

%% ========================================================================
% Local function 2: Main-axis formatting
% ========================================================================

function format_main_axis(ax)

    ax.FontName = 'Arial';
    ax.FontSize = 9.5;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.TickLength = [0.015, 0.015];
    ax.Layer = 'top';
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.06;

end