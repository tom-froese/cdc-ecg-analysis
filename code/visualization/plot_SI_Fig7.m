function plot_SI_Fig7()
% PLOT_SI_FIG7 - Supplementary Figure 7: Kaplan-Meier survival curves by
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
% Each panel includes a number-at-risk table below the axes, consistent
% with standard survival analysis reporting (CONSORT guidelines).
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
    fig_h_cm = 20.0;   % tall to accommodate three panels + at-risk tables

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 1 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters
    leg_fs   = 6;    % legend
    risk_fs  = 5.5;  % at-risk table

    %% ================================================================
    %  PANEL (a): Overall KM curves
    %  ================================================================

    ax1 = subplot(3, 1, 1);
    hold on; box on;

    h_lines = plot_km_panel(km_overall, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('Overall (N = %s; %s deaths)', ...
        format_comma(S.n_total), format_comma(S.n_deceased)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim([0.90 1.005]);
    xlim([0 tau]);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.90:0.02:1.00);

    % Legend
    leg_strs = cell(3, 1);
    for ti = 1:3
        leg_strs{ti} = sprintf('%s (n=%s)', ...
            tert_labels{ti}, format_comma(km_overall(ti).n));
    end
    leg = legend(h_lines, leg_strs, ...
        'Location', 'southwest', 'FontSize', leg_fs, 'Box', 'off');
    leg.ItemTokenSize = [12 8];

    % Log-rank p-value
    add_logrank_annotation(S.logrank_omnibus.p, ax_fs);

    % Panel label
    text(-0.08, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    % Number-at-risk table
    add_risk_table(ax1, risk_times, at_risk_all, tert_labels, ...
                   colors, risk_fs, tau);

    hold off;

    %% ================================================================
    %  PANEL (b): Female KM curves
    %  ================================================================

    ax2 = subplot(3, 1, 2);
    hold on; box on;

    h_lines_f = plot_km_panel(km_by_sex(1).km, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('%s (n = %s; %d deaths)', ...
        sex_labels{1}, format_comma(km_by_sex(1).n), km_by_sex(1).d), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim([0.90 1.005]);
    xlim([0 tau]);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.90:0.02:1.00);

    % Log-rank p-value
    add_logrank_annotation(km_by_sex(1).logrank_p, ax_fs);

    % Panel label
    text(-0.08, 1.06, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    % Number-at-risk table
    add_risk_table(ax2, risk_times, at_risk_sex{1}, tert_labels, ...
                   colors, risk_fs, tau);

    hold off;

    %% ================================================================
    %  PANEL (c): Male KM curves
    %  ================================================================

    ax3 = subplot(3, 1, 3);
    hold on; box on;

    h_lines_m = plot_km_panel(km_by_sex(2).km, colors, ci_alpha, 3);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('%s (n = %s; %d deaths)', ...
        sex_labels{2}, format_comma(km_by_sex(2).n), km_by_sex(2).d), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim([0.90 1.005]);
    xlim([0 tau]);
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], ...
        'YTick', 0.90:0.02:1.00);

    % Log-rank p-value
    add_logrank_annotation(km_by_sex(2).logrank_p, ax_fs);

    % Panel label
    text(-0.08, 1.06, '\bfc', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    % Number-at-risk table
    add_risk_table(ax3, risk_times, at_risk_sex{2}, tert_labels, ...
                   colors, risk_fs, tau);

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    %  Three stacked panels with room below each for the at-risk table.

    % Panel positions: [left, bottom, width, height]
    % Leave extra space at the bottom of each panel for the risk table
    panel_w = 0.82;
    panel_h = 0.22;
    left    = 0.12;
    gap     = 0.105;  % gap includes risk table space

    set(ax1, 'Position', [left, 0.72, panel_w, panel_h]);
    set(ax2, 'Position', [left, 0.39, panel_w, panel_h]);
    set(ax3, 'Position', [left, 0.06, panel_w, panel_h]);

    out_pdf = fullfile(paths.figures, 'SI_Fig7_survival_curves.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig7_survival_curves.png');
    out_fig = fullfile(paths.figures, 'SI_Fig7_survival_curves.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 7 saved.\n');

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================

    fprintf('\n--- Summary for SI Fig 7 legend ---\n');
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


function add_logrank_annotation(p_val, fs)
% ADD_LOGRANK_ANNOTATION - Place log-rank p-value in lower-left corner

    if isnan(p_val), return; end

    if p_val < 0.001
        p_str = 'Log-rank \itp\rm < 0.001';
    else
        p_str = sprintf('Log-rank \\itp\\rm = %.3f', p_val);
    end

    text(0.97, 0.06, p_str, 'Units', 'normalized', ...
        'FontSize', fs, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);
end


function add_risk_table(ax, risk_times, at_risk, tert_labels, colors, fs, tau)
% ADD_RISK_TABLE - Number-at-risk table below KM axes
%
%   Displays only times within [0, tau] and thins labels for readability.

    % Select time points to display (thin if > 8)
    valid = risk_times <= tau;
    disp_times = risk_times(valid);
    disp_risk  = at_risk(:, valid);

    if length(disp_times) > 8
        step = ceil(length(disp_times) / 8);
        idx = [1:step:length(disp_times), length(disp_times)];
        idx = unique(idx);
        disp_times = disp_times(idx);
        disp_risk  = disp_risk(:, idx);
    end

    % Get axes position in normalised figure coordinates
    ax_pos = get(ax, 'Position');  % [left, bottom, width, height]

    n_groups = size(at_risk, 1);
    row_h = 0.012;    % height of each row in normalised units
    table_top = ax_pos(2) - 0.015;  % start just below axes

    % Short labels for the risk table
    short_labels = {'T1 (Near)', 'T2 (Mod.)', 'T3 (Far)'};
    if n_groups > length(short_labels)
        short_labels = tert_labels;
    end

    for ti = 1:n_groups
        y_row = table_top - (ti - 1) * row_h;

        % Row label
        annotation('textbox', ...
            [ax_pos(1) - 0.09, y_row - row_h/2, 0.08, row_h], ...
            'String', short_labels{ti}, ...
            'FontSize', fs, 'FontWeight', 'bold', ...
            'Color', colors(ti, :) * 0.8, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FitBoxToText', 'off');

        % Numbers
        for ri = 1:length(disp_times)
            x_frac = ax_pos(1) + ax_pos(3) * disp_times(ri) / tau;
            annotation('textbox', ...
                [x_frac - 0.02, y_row - row_h/2, 0.04, row_h], ...
                'String', format_comma(disp_risk(ti, ri)), ...
                'FontSize', fs, ...
                'Color', [0.3 0.3 0.3], ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FitBoxToText', 'off');
        end
    end

    % "No. at risk" header
    annotation('textbox', ...
        [ax_pos(1) - 0.09, table_top + row_h * 0.3, 0.08, row_h], ...
        'String', 'No. at risk', ...
        'FontSize', fs, 'FontAngle', 'italic', ...
        'Color', [0.4 0.4 0.4], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off');
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
