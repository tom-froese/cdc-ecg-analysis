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
% Uses manual axes positioning and data-coordinate annotations to
% avoid exportgraphics displacement bugs.
%
% Data source: survival_curve_results.mat (via analyze_survival_curves.m)
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

    col_near = [0.30 0.65 0.40];
    col_mod  = [0.85 0.65 0.13];
    col_far  = [0.80 0.30 0.25];
    colors = [col_near; col_mod; col_far];
    ci_alpha = 0.12;

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================

    fig_w_cm = 18.3;
    fig_h_cm = 28.0;   % taller to give at-risk tables clearance

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 0.5 fig_w_cm fig_h_cm]);

    % Font sizes
    ax_fs    = 9;
    lab_fs   = 10;
    title_fs = 11;
    panel_fs = 13;
    leg_fs   = 8;
    risk_fs  = 7;
    lr_fs    = 9;

    km_lw    = 2.0;

    % Panel layout: [left, bottom, width, height]
    % Compact panel height; large gaps for x-axis labels + at-risk table
    panel_w = 0.82;
    panel_h = 0.17;
    left    = 0.12;

    pos_a = [left, 0.78, panel_w, panel_h];
    pos_b = [left, 0.46, panel_w, panel_h];
    pos_c = [left, 0.14, panel_w, panel_h];

    % Consistent y-axis across all panels
    y_lims = [0.88 1.005];
    y_ticks = 0.88:0.02:1.00;

    %% ================================================================
    %  PANEL (a): Overall KM curves
    %  ================================================================

    ax1 = axes('Position', pos_a);
    hold on; box on;

    h_lines = plot_km_panel(km_overall, colors, ci_alpha, 3, km_lw);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('Overall (N = %s; %s deaths)', ...
        format_comma(S.n_total), format_comma(S.n_deceased)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim(y_lims); xlim([0 tau]);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], 'YTick', y_ticks);

    % Legend
    leg_strs = cell(3, 1);
    for ti = 1:3
        leg_strs{ti} = sprintf('%s (n=%s)', ...
            tert_labels{ti}, format_comma(km_overall(ti).n));
    end
    leg1 = legend(h_lines, leg_strs, ...
        'Location', 'southwest', 'FontSize', leg_fs, 'Box', 'on');
    set(leg1, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
    leg1.ItemTokenSize = [15 10];

    add_logrank_data(ax1, S.logrank_omnibus.p, tau, y_lims, lr_fs);

    axes(ax1);
    xl = xlim; yl = ylim;
    text(xl(1) - 0.08*diff(xl), yl(2) + 0.03*diff(yl), '\bfa', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    add_risk_table(ax1, risk_times, at_risk_all, tert_labels, ...
                   colors, risk_fs, tau);
    hold off;

    %% ================================================================
    %  PANEL (b): Female KM curves
    %  ================================================================

    ax2 = axes('Position', pos_b);
    hold on; box on;

    plot_km_panel(km_by_sex(1).km, colors, ci_alpha, 3, km_lw);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('%s (n = %s; %s deaths)', ...
        sex_labels{1}, format_comma(km_by_sex(1).n), ...
        format_comma(km_by_sex(1).d)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim(y_lims); xlim([0 tau]);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], 'YTick', y_ticks);

    add_logrank_data(ax2, km_by_sex(1).logrank_p, tau, y_lims, lr_fs);

    axes(ax2);
    xl = xlim; yl = ylim;
    text(xl(1) - 0.08*diff(xl), yl(2) + 0.03*diff(yl), '\bfb', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    add_risk_table(ax2, risk_times, at_risk_sex{1}, tert_labels, ...
                   colors, risk_fs, tau);
    hold off;

    %% ================================================================
    %  PANEL (c): Male KM curves
    %  ================================================================

    ax3 = axes('Position', pos_c);
    hold on; box on;

    plot_km_panel(km_by_sex(2).km, colors, ci_alpha, 3, km_lw);

    xlabel('Follow-up (years)', 'FontSize', lab_fs);
    ylabel('Survival probability', 'FontSize', lab_fs);
    title(sprintf('%s (n = %s; %s deaths)', ...
        sex_labels{2}, format_comma(km_by_sex(2).n), ...
        format_comma(km_by_sex(2).d)), ...
        'FontSize', title_fs, 'FontWeight', 'bold');

    ylim(y_lims); xlim([0 tau]);
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015], 'YTick', y_ticks);

    add_logrank_data(ax3, km_by_sex(2).logrank_p, tau, y_lims, lr_fs);

    axes(ax3);
    xl = xlim; yl = ylim;
    text(xl(1) - 0.08*diff(xl), yl(2) + 0.03*diff(yl), '\bfc', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    add_risk_table(ax3, risk_times, at_risk_sex{2}, tert_labels, ...
                   colors, risk_fs, tau);
    hold off;

    %% ================================================================
    %  SAVE
    %  ================================================================

    out_pdf = fullfile(paths.figures, 'SI_Fig7_survival_curves.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig7_survival_curves.png');
    out_fig = fullfile(paths.figures, 'SI_Fig7_survival_curves.fig');

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf('  Saved: %s (raster, 300 dpi)\n', out_png);

    exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 300);
    fprintf('  Saved: %s (raster-in-PDF, 300 dpi)\n', out_pdf);

    fprintf('  Skipped: %s (too many graphic objects)\n', out_fig);
    fprintf('\nSupplementary Figure 7 saved.\n');

    %% ================================================================
    %  CONSOLE SUMMARY
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

function h_lines = plot_km_panel(km_struct, colors, ci_alpha, n_groups, lw)
    h_lines = gobjects(n_groups, 1);
    for ti = n_groups:-1:1
        t = km_struct(ti).t;
        s = km_struct(ti).s;
        lo = km_struct(ti).lo;
        hi = km_struct(ti).hi;
        col = colors(ti, :);

        [t_step, lo_step] = stairs(t, lo);
        [~, hi_step] = stairs(t, hi);

        fill([t_step; flipud(t_step)], [lo_step; flipud(hi_step)], ...
            col, 'FaceAlpha', ci_alpha, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h_lines(ti) = stairs(t, s, 'Color', col, 'LineWidth', lw);
    end
end


function add_logrank_data(ax, p_val, tau, y_lims, fs)
    if isnan(p_val), return; end
    axes(ax);
    if p_val < 0.001
        p_str = 'Log-rank \itp\rm < 0.001';
    else
        p_str = sprintf('Log-rank \\itp\\rm = %.3f', p_val);
    end
    x_pos = tau * 0.97;
    y_pos = y_lims(1) + 0.06 * diff(y_lims);
    text(x_pos, y_pos, p_str, ...
        'FontSize', fs, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.5 0.5 0.5], ...
        'Margin', 2);
end


function add_risk_table(ax, risk_times, at_risk, tert_labels, colors, fs, tau)
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

    ax_pos = get(ax, 'Position');
    n_groups = size(at_risk, 1);
    row_h = 0.012;
    % Push well below axes to clear x-axis tick labels and xlabel
    table_top = ax_pos(2) - 0.050;

    short_labels = {'T1 (Near)', 'T2 (Mod.)', 'T3 (Far)'};
    if n_groups > length(short_labels)
        short_labels = tert_labels;
    end

    for ti = 1:n_groups
        y_row = table_top - (ti - 1) * row_h;

        annotation('textbox', ...
            [ax_pos(1) - 0.10, y_row - row_h/2, 0.09, row_h], ...
            'String', short_labels{ti}, ...
            'FontSize', fs, 'FontWeight', 'bold', ...
            'Color', colors(ti, :) * 0.8, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FitBoxToText', 'off');

        for ri = 1:length(disp_times)
            x_frac = ax_pos(1) + ax_pos(3) * disp_times(ri) / tau;
            annotation('textbox', ...
                [x_frac - 0.025, y_row - row_h/2, 0.05, row_h], ...
                'String', format_comma(disp_risk(ti, ri)), ...
                'FontSize', fs, ...
                'Color', [0.3 0.3 0.3], ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FitBoxToText', 'off');
        end
    end

    annotation('textbox', ...
        [ax_pos(1) - 0.10, table_top + row_h * 0.4, 0.09, row_h], ...
        'String', 'No. at risk', ...
        'FontSize', fs, 'FontAngle', 'italic', ...
        'Color', [0.4 0.4 0.4], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off');
end


function s = format_comma(n)
    s = num2str(n);
    if n >= 1000
        idx = length(s) - 2;
        while idx > 1
            s = [s(1:idx-1) ',' s(idx:end)];
            idx = idx - 3;
        end
    end
end