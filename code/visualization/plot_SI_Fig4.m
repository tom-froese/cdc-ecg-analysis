function plot_SI_Fig4()
% PLOT_SI_FIG4 - Supplementary Figure 4: Kaplan-Meier survival curves by
%   CDC-deviation tertile, overall and sex-stratified
%
% Three-panel figure:
%   a: Overall KM curves by CDC-deviation tertile (with 95% CI)
%   b: Female KM curves by CDC-deviation tertile
%   c: Male KM curves by CDC-deviation tertile
%
% Colour palette matches Figure 2 (plot_Fig2.m):
%   Green  [0.30 0.65 0.40]  = Near 1/e  (T1, nearest third)
%   Amber  [0.85 0.65 0.13]  = Moderate  (T2, middle third)
%   Red    [0.80 0.30 0.25]  = Far from 1/e (T3, farthest third)
%
% Group sample sizes are shown in the legend. Survival probability on the
% y-axis allows assessment of numbers at risk at each follow-up time.
%
% Data source: survival_curve_results.mat (via analyze_survival_curves.m)
%
% Nature Aging formatting: double-column width (183 mm), bold lowercase
%   panel labels outside axes, 7-pt tick labels, log-log 95% CI bands.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();

    %% ================================================================
    %  LOAD PRECOMPUTED RESULTS
    %  ================================================================

    S = load(fullfile(paths.results, 'survival_curve_results.mat'));

    km_overall   = S.km_overall;
    km_by_sex    = S.km_by_sex;
    risk_times   = S.risk_times;
    at_risk_all  = S.at_risk_overall;
    at_risk_sex  = S.at_risk_sex;
    tert_labels  = S.tert_labels;
    sex_labels   = S.sex_labels;
    tau          = S.rmst_tau;

    %% ================================================================
    %  COLOUR SCHEME (matches Figure 2)
    %  ================================================================

    col_near = [0.30 0.65 0.40];   % Green  — nearest to 1/e
    col_mod  = [0.85 0.65 0.13];   % Amber  — moderate deviation
    col_far  = [0.80 0.30 0.25];   % Red    — farthest from 1/e
    colors = [col_near; col_mod; col_far];

    % Lighter shades for CI bands
    ci_alpha = 0.12;

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================

    fig_w_cm = 18.3;   % double-column width
    fig_h_cm = 22.0;   % three stacked panels with spacing for labels

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 1 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters
    leg_fs   = 7;    % legend (matches tick/p-value font)

    %% ================================================================
    %  PANEL (a): Overall KM curves
    %  ================================================================

    ax1 = subplot(3, 1, 1);
    hold on; box on;

    h_lines = plot_km_panel(km_overall, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    t1 = title(sprintf('Overall (N = %s; %s deaths)', ...
        format_comma(S.n_total), format_comma(S.n_deceased)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');
    t1.Units = 'normalized'; t1.Position(2) = t1.Position(2) + 0.04;

    ylim([0.88 1.005]);
    xlim([0 tau]);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.88:0.02:1.00);

    % Legend (positioned in layout section)
    leg_strs = cell(3, 1);
    for ti = 1:3
        leg_strs{ti} = sprintf('%s (n=%s)', ...
            tert_labels{ti}, format_comma(km_overall(ti).n));
    end
    leg1 = legend(h_lines, leg_strs, ...
        'FontSize', leg_fs, 'Box', 'off');
    leg1.ItemTokenSize = [12 8];

    % Panel label
    text(-0.08, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): Female KM curves
    %  ================================================================

    ax2 = subplot(3, 1, 2);
    hold on; box on;

    h_lines_f = plot_km_panel(km_by_sex(1).km, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    t2 = title(sprintf('%s (n = %s; %s deaths)', ...
        sex_labels{1}, format_comma(km_by_sex(1).n), format_comma(km_by_sex(1).d)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');
    t2.Units = 'normalized'; t2.Position(2) = t2.Position(2) + 0.04;

    ylim([0.88 1.005]);
    xlim([0 tau]);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.88:0.02:1.00);

    % Legend (positioned in layout section)
    leg_strs_f = cell(3, 1);
    for ti = 1:3
        leg_strs_f{ti} = sprintf('%s (n=%s)', ...
            tert_labels{ti}, format_comma(km_by_sex(1).km(ti).n));
    end
    leg2 = legend(h_lines_f, leg_strs_f, ...
        'FontSize', leg_fs, 'Box', 'off');
    leg2.ItemTokenSize = [12 8];

    % Panel label
    text(-0.08, 1.06, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (c): Male KM curves
    %  ================================================================

    ax3 = subplot(3, 1, 3);
    hold on; box on;

    h_lines_m = plot_km_panel(km_by_sex(2).km, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    t3 = title(sprintf('%s (n = %s; %s deaths)', ...
        sex_labels{2}, format_comma(km_by_sex(2).n), format_comma(km_by_sex(2).d)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');
    t3.Units = 'normalized'; t3.Position(2) = t3.Position(2) + 0.04;

    ylim([0.88 1.005]);
    xlim([0 tau]);
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.88:0.02:1.00);

    % Legend (positioned in layout section)
    leg_strs_m = cell(3, 1);
    for ti = 1:3
        leg_strs_m{ti} = sprintf('%s (n=%s)', ...
            tert_labels{ti}, format_comma(km_by_sex(2).km(ti).n));
    end
    leg3 = legend(h_lines_m, leg_strs_m, ...
        'FontSize', leg_fs, 'Box', 'off');
    leg3.ItemTokenSize = [12 8];

    % Panel label
    text(-0.08, 1.06, '\bfc', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    %  Three stacked panels, evenly spaced.

    % Panel positions: [left, bottom, width, height]
    panel_w = 0.50;
    panel_h = 0.22;
    left    = 0.10;
    gap     = 0.10;
    panel_bottoms = [0.72, 0.40, 0.08];

    set(ax1, 'Position', [left, panel_bottoms(1), panel_w, panel_h]);
    set(ax2, 'Position', [left, panel_bottoms(2), panel_w, panel_h]);
    set(ax3, 'Position', [left, panel_bottoms(3), panel_w, panel_h]);

    % Legends — right side, one per panel, left-aligned, vertically centered
    leg_left = 0.62;
    legs = [leg1, leg2, leg3];
    for k = 1:3
        set(legs(k), 'Units', 'normalized');
        auto_pos = get(legs(k), 'Position');
        leg_y = panel_bottoms(k) + (panel_h - auto_pos(4)) / 2;
        set(legs(k), 'Position', ...
            [leg_left, leg_y, auto_pos(3), auto_pos(4)]);
    end

    % Log-rank p-values — right side, below each legend
    p_vals = [S.logrank_omnibus.p, km_by_sex(1).logrank_p, km_by_sex(2).logrank_p];
    for k = 1:3
        if p_vals(k) < 0.001
            pk_str = 'Log-rank \itp\rm < 0.001';
        else
            pk_str = sprintf('Log-rank \\itp\\rm = %.3f', p_vals(k));
        end
        pk_y = panel_bottoms(k) + 0.01;
        annotation('textbox', [leg_left, pk_y, 0.36, 0.03], ...
            'String', pk_str, 'FontSize', ax_fs, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'FitBoxToText', 'on');
    end

    % Number-at-risk tables removed — full n per group shown in legend,
    % and y-axis (survival probability) allows assessment at each time point.

    out_pdf = fullfile(paths.figures, 'SI_Fig4_survival_curves.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig4_survival_curves.png');
    out_fig = fullfile(paths.figures, 'SI_Fig4_survival_curves.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 4 saved.\n');

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================

    fprintf('\n--- Summary for SI Fig 4 legend ---\n');
    fprintf('N = %s patients (%s deaths, %.1f%%)\n', ...
        format_comma(S.n_total), format_comma(S.n_deceased), ...
        100 * S.n_deceased / S.n_total);
    fprintf('Median follow-up: %.2f years\n', S.median_followup);
    fprintf('Tertile boundaries: |CDC - 1/e| = %.4f / %.4f\n', ...
        S.tertile_bounds(1), S.tertile_bounds(2));

    fprintf('\nOverall survival by tertile:\n');
    for ti = 1:3
        fprintf('  %-25s  N=%s, Deaths=%d (%.2f%%), RMST(%dyr)=%.3f\n', ...
            tert_labels{ti}, format_comma(km_overall(ti).n), ...
            km_overall(ti).d, 100 * km_overall(ti).d / km_overall(ti).n, ...
            tau, km_overall(ti).rmst);
    end
    fprintf('RMST difference (T1-T3): %.3f years\n', S.rmst_diff_t1_t3);

    fprintf('\nLog-rank tests:\n');
    fprintf('  Omnibus: chi2=%.2f, p=%.2e\n', ...
        S.logrank_omnibus.chi2, S.logrank_omnibus.p);
    fprintf('  T1 vs T3: chi2=%.2f, p=%.2e\n', ...
        S.logrank_t1_vs_t3.chi2, S.logrank_t1_vs_t3.p);

    fprintf('\nSex-stratified log-rank (omnibus):\n');
    for si = 1:2
        fprintf('  %s: chi2=%.2f, p=%.2e (N=%s)\n', ...
            sex_labels{si}, km_by_sex(si).logrank_chi2, ...
            km_by_sex(si).logrank_p, format_comma(km_by_sex(si).n));
    end
end


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function h_lines = plot_km_panel(km_struct, colors, ci_alpha, n_groups)
% PLOT_KM_PANEL - Plot KM step functions with CI bands for all groups
%
%   Renders in reverse order (T3 first) so that T1 (nearest) is on top.

    h_lines = gobjects(n_groups, 1);

    for ti = n_groups:-1:1
        t = km_struct(ti).t;
        s = km_struct(ti).s;
        lo = km_struct(ti).lo;
        hi = km_struct(ti).hi;
        col = colors(ti, :);

        % Step-function versions for CI bands
        [t_step, lo_step] = stairs(t, lo);
        [~, hi_step] = stairs(t, hi);

        % CI band
        fill([t_step; flipud(t_step)], [lo_step; flipud(hi_step)], ...
            col, 'FaceAlpha', ci_alpha, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % KM step function
        h_lines(ti) = stairs(t, s, 'Color', col, 'LineWidth', 1.5);
    end
end



function save_large_figure(fig, out_pdf, out_png, out_fig, w_cm, h_cm)
% SAVE_LARGE_FIGURE - Save figure without crashing on large scatter data
%
% For figures with many graphic objects, MATLAB's painters renderer and
% savefig serialize every element, consuming gigabytes of memory.
%
% Strategy:
%   PNG: exportgraphics (raster, always safe)
%   PDF: exportgraphics with ContentType 'image' (raster-in-PDF wrapper,
%        avoids the painters memory explosion while producing a PDF that
%        embeds at 300 dpi)
%   FIG: skipped (the .fig format stores all graphic objects and can
%        itself become multi-GB)

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [w_cm h_cm]);
    set(fig, 'PaperPosition', [0 0 w_cm h_cm]);

    % PNG — raster, always safe
    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf('  Saved: %s (raster, 300 dpi)\n', out_png);

    % PDF — raster-in-PDF (avoids painters memory explosion)
    exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 300);
    fprintf('  Saved: %s (raster-in-PDF, 300 dpi)\n', out_pdf);

    % FIG — skip for large figures
    fprintf('  Skipped: %s (too many graphic objects for .fig format)\n', out_fig);

    close(fig);
end


function s = format_comma(n)
% FORMAT_COMMA - Format integer with thousands separators
    s = num2str(n);
    if n >= 1000
        idx = length(s) - 2;
        while idx > 1
            s = [s(1:idx-1) ',' s(idx:end)];
            idx = idx - 3;
        end
    end
end
