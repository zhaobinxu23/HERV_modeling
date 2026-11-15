function sensitivity_analysis_herv()
% 参数敏感性分析主脚本
% 输出：sensitivity_results/ 文件夹中的 .mat、.csv、.png

close all; clc;
output_dir = 'sensitivity_results';
if ~exist(output_dir,'dir'), mkdir(output_dir); end

%% 1. 基线设置（癌症条件）
[par0, y0, tspan, opts] = setup_baseline();

%% 2. 全局 OAT 扫描
global_params = {'p2','p3','p4','p5','p7','p8','p9','p10','p11','p12', ...
    'p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23', ...
    'p24','p25','p26','p27','p28','p29','p30','p31','p32','p33','p34', ...
    'p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45', ...
    'p46','p47','p48','p49','p50','p51','p52','p53','p56','p57','p58', ...
    'p59','p61','p62','p63'};

fprintf('==== Global OAT sensitivity ====\n');
global_results = run_global_oat(par0, y0, tspan, opts, global_params);

%% 3. 重点参数扫描
fprintf('\n==== Focused parameter scans ====\n');
focused_results = run_focused_scans(par0, y0, tspan, opts);

%% 4. 画图与保存
plot_global_sensitivity(global_results, output_dir);
plot_focused_scans(focused_results, output_dir);
plot_phase_diagrams(par0, y0, tspan, opts, output_dir);
save_results(global_results, focused_results, output_dir);

fprintf('\nAll done. Results saved in: %s\n', fullfile(pwd, output_dir));
end

%% ================= 基线 =================
function [par, y0, tspan, opts] = setup_baseline()
par = struct();
par.p1=0; par.p2=0.5; par.p3=0.05; par.p4=10; par.p5=1e20;
par.p6=1e-9; par.p7=1e8; par.p8=1e-12; par.p9=0.1; par.p10=0.9;
par.p11=1e9; par.p12=9*1.22e10; par.p13=0.1; par.p14=1e-9;
par.p15=1.0; par.p16=0.01; par.p17=1e-10; par.p18=1;
par.p19=1e8; par.p20=0.01; par.p21=1e1; par.p22=1e1;
par.p23=0.01; par.p24=0.01; par.p25=1e7; par.p26=1e8;
par.p27=0.01; par.p28=0.5; par.p29=1e1; par.p30=1;
par.p31=1; par.p32=1; par.p33=0.1; par.p34=1e4;
par.p35=1e9; par.p36=9*1.22e10; par.p37=1e-11; par.p38=1.1;
par.p39=1e-10; par.p40=2; par.p41=0.01; par.p42=1e12;
par.p43=0.01; par.p44=0.01; par.p45=0.9; par.p46=0.01;
par.p47=1; par.p48=0.01; par.p49=1; par.p50=1;
par.p51=1e5; par.p52=50; par.p53=0.01; par.p54=0;
par.p55=0; par.p56=10; par.p57=1e2; par.p58=1e2;
par.p59=1e2; par.p60=0; par.p61=1e8; par.p62=1; par.p63=1e4;
par.Tu_cutoff = 1e-4;

y0 = zeros(36,1);
y0(1)=100; y0(5)=1e10; y0(7)=1e6; y0(8)=9999*1e6;
y0(10)=1e10; y0(11)=1e6; y0(12)=9999*1e6;
y0(19)=1e10; y0(23)=1e10; y0(24)=9999*1e10;
y0(25)=1e6; y0(28)=1e10;

tspan = [0 100];
opts = odeset('RelTol',1e-6, 'AbsTol', 1e-8*max(ones(36,1),abs(y0)), ...
              'NonNegative', 1:36, 'MaxStep', 5);
end

%% ================= 全局 OAT =================
function results = run_global_oat(par0, y0, tspan, opts, params)
multipliers = [0.1, 0.3, 1, 3, 10];
n_params = numel(params); n_mult = numel(multipliers);
peakTu = nan(n_params, n_mult);
logPeak = nan(n_params, n_mult);

for i = 1:n_params
    pname = params{i};
    base_val = par0.(pname);
    fprintf('  [%d/%d] %s (base=%.3g)\n', i, n_params, pname, base_val);
    for j = 1:n_mult
        par = par0;
        par.(pname) = base_val * multipliers(j);
        if par.(pname) <= 0, par.(pname) = 1e-30; end
        pt = simulate_peak_Tu(par, y0, tspan, opts);
        peakTu(i,j) = pt;
        if isnan(pt) || pt <= 0
            logPeak(i,j) = NaN;
        else
            logPeak(i,j) = log10(max(pt, 1e-6));
        end
    end
end

S = nan(n_params,1);
for i = 1:n_params
    if ~isnan(logPeak(i,1)) && ~isnan(logPeak(i,5))
        S(i) = (logPeak(i,5) - logPeak(i,1)) / 2;
    end
end
results.params = params;
results.multipliers = multipliers;
results.peakTu = peakTu;
results.logPeak = logPeak;
results.S = S;
end

%% ================= 重点扫描 =================
function results = run_focused_scans(par0, y0, tspan, opts)
results = struct();

results.p56 = scan_param(par0,'p56',[0.1,0.3,1,3,10,30,100,300,1000],y0,tspan,opts);
results.p57 = scan_param(par0,'p57',[1,10,30,100,300,1000,3000],y0,tspan,opts);
results.p63 = scan_param(par0,'p63',[1e2,1e3,1e4,1e5,1e6,1e7],y0,tspan,opts);
results.p62 = scan_param(par0,'p62',[0,1],y0,tspan,opts);
results.p51 = scan_param(par0,'p51',[1e3,1e4,1e5,1e6,1e7],y0,tspan,opts);
results.p29 = scan_param(par0,'p29',[0.1,0.3,1,3,10,30,100],y0,tspan,opts);

% p0：统一缩放 p21,p22,p29,p51,p52
p0_vals = [0.1,0.3,1,3,10,30];
p0_peak = nan(size(p0_vals));
for i = 1:numel(p0_vals)
    par = apply_immune_capacity(par0, p0_vals(i));
    p0_peak(i) = simulate_peak_Tu(par, y0, tspan, opts);
end
results.p0.vals = p0_vals;
results.p0.peakTu = p0_peak;
results.p0.pname = 'p0';

results.p8 = scan_param(par0,'p8',[1e-14,1e-13,1e-12,1e-11,1e-10,1e-9],y0,tspan,opts);
results.p4 = scan_param(par0,'p4',[0.01,0.1,1,3,10,30,100,300,1000],y0,tspan,opts);
results.p5 = scan_param(par0,'p5',[1e14,1e16,1e18,1e20,1e22,1e24,1e26],y0,tspan,opts);
results.p7 = scan_param(par0,'p7',[1e2,1e4,1e6,1e8,1e10,1e12,1e14],y0,tspan,opts);
results.p1 = scan_param(par0,'p1',[0,1,10,100,1000,10000],y0,tspan,opts);
end

function out = scan_param(par0, pname, vals, y0, tspan, opts)
peak = nan(size(vals));
for i = 1:numel(vals)
    par = par0;
    par.(pname) = vals(i);
    peak(i) = simulate_peak_Tu(par, y0, tspan, opts);
    fprintf('    %s = %.3g, peak Tu = %.3g\n', pname, vals(i), peak(i));
end
out.vals = vals;
out.peakTu = peak;
out.pname = pname;
end

function par = apply_immune_capacity(par0, cap)
par = par0;
par.p21 = par0.p21 * cap;
par.p22 = par0.p22 * cap;
par.p29 = par0.p29 * cap;
par.p51 = par0.p51 * cap;
par.p52 = par0.p52 * cap;
end

%% ================= 单次模拟 =================
function peakTu = simulate_peak_Tu(par, y0, tspan, opts)
try
    [~, Y] = ode15s(@(t,y) tumor_immune_interaction_herv_sens(t,y,par), ...
                    tspan, y0, opts);
    if isempty(Y) || any(~isfinite(Y(:)))
        peakTu = NaN; return;
    end
    peakTu = max(Y(:,1));
    if peakTu < par.Tu_cutoff, peakTu = par.Tu_cutoff; end
catch ME
    warning('Simulation failed: %s', ME.message);
    peakTu = NaN;
end
end

%% ================= 画图 =================
function plot_global_sensitivity(res, output_dir)
    S = res.S(:);                 % 列向量
    params = res.params(:);       % 列 cell
    [~, idx] = sort(abs(S), 'descend', 'MissingPlacement','last');

    S_sorted      = S(idx);
    names_sorted  = params(idx);

    fig = figure('Color','w','Position',[100 100 900 1000]);
    barh(S_sorted);
    set(gca,'YTick',1:numel(names_sorted), ...
        'YTickLabel',names_sorted, ...
        'TickLabelInterpreter','none','FontSize',9);
    xlabel('S = d log_{10}(peak Tu) / d log_{10}(\theta)');
    title('Global OAT sensitivity of peak tumor burden');
    grid on;
    xline(0,'k--');
    saveas(fig, fullfile(output_dir,'Fig_S1_global_sensitivity.png'));
end

function plot_focused_scans(res, output_dir)
fig = figure('Color','w','Position',[100 100 1400 1000]);

ax = subplot(3,3,1); plot_logx(ax,res.p56,'p56','HERV transcription');
ax = subplot(3,3,2); plot_logx(ax,res.p57,'p57','Translation coefficient');
ax = subplot(3,3,3); plot_logx(ax,res.p63,'p63','Feedback half-saturation');
ax = subplot(3,3,4); plot_bar_cat(ax,res.p62,'p62','Reverse-transcription switch');
ax = subplot(3,3,5); plot_logx(ax,res.p51,'p51','B-cell expansion');
ax = subplot(3,3,6); plot_logx(ax,res.p29,'p29','CD8+ T-cell recovery');
ax = subplot(3,3,7); plot_logx(ax,res.p0,'p0','Immune competence');
ax = subplot(3,3,8); plot_logx(ax,res.p8,'p8','Antibody/BCR affinity');
ax = subplot(3,3,9); plot_logx(ax,res.p4,'p4','Max killing rate');

saveas(fig, fullfile(output_dir,'Fig_S2_focused_scans_part1.png'));

fig2 = figure('Color','w','Position',[100 100 900 400]);
ax = subplot(1,2,1); plot_logx(ax,res.p5,'p5','Antibody killing threshold');
ax = subplot(1,2,2); plot_logx(ax,res.p7,'p7','TCR killing threshold');
saveas(fig2, fullfile(output_dir,'Fig_S3_focused_scans_part2.png'));
end

function plot_logx(ax, res, pname, label)
semilogx(ax, res.vals, res.peakTu, '-o', 'LineWidth', 2, 'MarkerSize', 5);
xlabel(ax, sprintf('%s (%s)', pname, label));
ylabel(ax, 'Peak tumor burden');
grid on; set(ax,'FontSize',9);
end

function plot_bar_cat(ax, res, pname, label)
bar(ax, 1:numel(res.vals), res.peakTu);
set(ax,'XTick',1:numel(res.vals), ...
    'XTickLabel', arrayfun(@(x) sprintf('%s=%d',pname,x),res.vals,'UniformOutput',false));
xlabel(ax,label); ylabel(ax,'Peak tumor burden'); grid on;
end

function plot_phase_diagrams(par0, y0, tspan, opts, output_dir)
p56_vals = [0.1,0.3,1,3,10,30,100];
p4_vals  = [0.1,1,10,100];
p5_vals  = [1e16,1e18,1e20,1e22,1e24];
p7_vals  = [1e4,1e6,1e8,1e10,1e12];

% p56 vs p4
peak = nan(numel(p4_vals),numel(p56_vals));
for i=1:numel(p4_vals)
    for j=1:numel(p56_vals)
        par=par0; par.p56=p56_vals(j); par.p4=p4_vals(i);
        peak(i,j)=simulate_peak_Tu(par,y0,tspan,opts);
    end
end
fig=figure('Color','w'); imagesc(log10(peak)); colorbar;
set(gca,'XTick',1:numel(p56_vals),'XTickLabel',p56_vals, ...
    'YTick',1:numel(p4_vals),'YTickLabel',p4_vals);
xlabel('p56'); ylabel('p4'); title('log10(peak Tu): p56 vs p4');
saveas(fig,fullfile(output_dir,'Fig_S4_phase_p56_p4.png'));

% p56 vs p5
peak = nan(numel(p5_vals),numel(p56_vals));
for i=1:numel(p5_vals)
    for j=1:numel(p56_vals)
        par=par0; par.p56=p56_vals(j); par.p5=p5_vals(i);
        peak(i,j)=simulate_peak_Tu(par,y0,tspan,opts);
    end
end
fig=figure('Color','w'); imagesc(log10(peak)); colorbar;
set(gca,'XTick',1:numel(p56_vals),'XTickLabel',p56_vals, ...
    'YTick',1:numel(p5_vals),'YTickLabel',p5_vals);
xlabel('p56'); ylabel('p5'); title('log10(peak Tu): p56 vs p5');
saveas(fig,fullfile(output_dir,'Fig_S5_phase_p56_p5.png'));

% p56 vs p7
peak = nan(numel(p7_vals),numel(p56_vals));
for i=1:numel(p7_vals)
    for j=1:numel(p56_vals)
        par=par0; par.p56=p56_vals(j); par.p7=p7_vals(i);
        peak(i,j)=simulate_peak_Tu(par,y0,tspan,opts);
    end
end
fig=figure('Color','w'); imagesc(log10(peak)); colorbar;
set(gca,'XTick',1:numel(p56_vals),'XTickLabel',p56_vals, ...
    'YTick',1:numel(p7_vals),'YTickLabel',p7_vals);
xlabel('p56'); ylabel('p7'); title('log10(peak Tu): p56 vs p7');
saveas(fig,fullfile(output_dir,'Fig_S6_phase_p56_p7.png'));
end

%% ================= 保存 =================
function save_results(global_results, focused_results, output_dir)
    save(fullfile(output_dir,'sensitivity_results.mat'), ...
         'global_results','focused_results');

    % 强制统一为列向量，避免 table 行数不一致
    params_col = global_results.params(:);
    S_col      = global_results.S(:);

    % 防御性检查：长度必须一致
    assert(numel(params_col) == numel(S_col), ...
        'params 和 S 长度不一致：%d vs %d', numel(params_col), numel(S_col));

    T = table(params_col, S_col, ...
        'VariableNames', {'Parameter','SensitivityIndex'});

    writetable(T, fullfile(output_dir,'global_sensitivity.csv'));
end