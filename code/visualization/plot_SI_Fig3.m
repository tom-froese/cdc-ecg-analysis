function plot_SI_Fig3()
% PLOT_SI_FIG3 - Supplementary Figure 3: Bland-Altman agreement between
%   automatic and manual annotation pipelines
%
% Four-panel figure:
%   a  R-peak timing error (ms) by beat index
%   b  T-end timing error (ms) vs manual CDC
%   c  CDC agreement (beat-level Bland-Altman)
%   d  CDC agreement (subject-level median Bland-Altman)
%
% Data source: validation_results.mat (via analyze_gold_standard_validation.m)
%
% Nature Aging formatting: double-column (183 mm), bold lowercase
%   panel labels outside axes, 7-pt tick labels, no sgtitle.
%
% Figure export uses print() with -painters renderer and savefig(),
% which are stable across MATLAB versions and avoid the rendering-engine
% crashes that exportgraphics can trigger with large scatter plots.
%
% Tom Froese, OIST Embodied Cognitive Science Unit, March 2026

    paths = config();

    %% ================================================================
    %  LOAD DATA
    %  ================================================================

    mat_file = fullfile(paths.results, 'validation_results.mat');
    if ~exist(mat_file, 'file')
        error(['validation_results.mat not found.\n' ...
               'Run analyze_gold_standard_validation() first.']);
    end

    S = load(mat_file, 'results');
    a = S.results.ludb;
    b = S.results.qtdb;

    fprintf('\nSupplementary Figure 3 — data loaded:\n');
    fprintf('  LUDB: %d matched R-peaks, %d matched T-ends, %d subjects\n', ...
        a.n_matched_r, a.n_matched_t, a.subj_n);
    fprintf('  QTDB: %d matched R-peaks, %d matched T-ends, %d subjects\n', ...
        b.n_matched_r, b.n_matched_t, b.subj_n);

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================

    % Nature double-column: 183 mm wide.  Aspect ~ 16:9 for 4 panels.
    fig_w_mm = 183;
    fig_h_mm = 110;
    fig_w_in = fig_w_mm / 25.4;
    fig_h_in = fig_h_mm / 25.4;

    fig = figure('Units', 'inches', 'Position', [1 1 fig_w_in fig_h_in], ...
                 'Color', 'w', 'Visible', 'on');

    % Font sizes for final print dimensions (Nature minimum: 5 pt)
    ax_fs  = 9;    % axis tick labels
    lab_fs = 12;    % axis labels
    ttl_fs = 12;   % panel titles
    ann_fs = 10;    % annotation text (bias, MAE, r)
    leg_fs = 10;    % legend

    % Colours: LUDB blue, QTDB orange (consistent with SI Fig 4)
    c_ludb = [0.20 0.47 0.76];
    c_qtdb = [0.85 0.40 0.25];
    ms = 4;   % marker size (points)

    %% ================================================================
    %  PANEL a: R-peak timing error
    %  ================================================================

    ax1 = subplot(2,2,1); hold on;
    n1 = length(a.r_err_ms);
    n2 = length(b.r_err_ms);
    if n1 > 0
        scatter((1:n1)', a.r_err_ms, ms, c_ludb, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'LUDB');
    end
    if n2 > 0
        scatter(n1+(1:n2)', b.r_err_ms, ms, c_qtdb, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'QTDB');
    end
    d_r = [a.r_err_ms; b.r_err_ms];
    draw_ba_lines(d_r);
    legend('LUDB', 'QTDB', 'Location', 'southwest', 'FontSize', leg_fs, 'Box', 'off');
    ylabel('R-peak error (ms)', 'FontSize', lab_fs);
    xlabel('Beat index', 'FontSize', lab_fs);
    title('a  R-peak (beat-level)', 'FontWeight', 'bold', 'FontSize', ttl_fs);
    set(gca, 'FontSize', ax_fs, 'Box', 'on');
    annotate_panel(ax1, d_r, [], ann_fs, 0.05);

    %% ================================================================
    %  PANEL b: T-end timing error vs manual CDC
    %  ================================================================

    ax2 = subplot(2,2,2); hold on;
    if ~isempty(a.t_err_ms) && length(a.cdc_man) == length(a.t_err_ms)
        scatter(a.cdc_man, a.t_err_ms, ms, c_ludb, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'LUDB');
    end
    if ~isempty(b.t_err_ms) && length(b.cdc_man) == length(b.t_err_ms)
        scatter(b.cdc_man, b.t_err_ms, ms, c_qtdb, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'QTDB');
    end
    d_t = [a.t_err_ms; b.t_err_ms];
    draw_ba_lines(d_t);
    xlabel('Manual CDC', 'FontSize', lab_fs);
    ylabel('T-end error (ms)', 'FontSize', lab_fs);
    title('b  T-end (beat-level)', 'FontWeight', 'bold', 'FontSize', ttl_fs);
    set(gca, 'FontSize', ax_fs, 'Box', 'on');
    annotate_panel(ax2, d_t, [], ann_fs, 0.05);

    %% ================================================================
    %  PANEL c: CDC Bland-Altman (beat-level)
    %  ================================================================

    ax3 = subplot(2,2,3); hold on;
    cm = [a.cdc_man; b.cdc_man];
    ca = [a.cdc_aut; b.cdc_aut];
    if length(cm) == length(ca) && ~isempty(cm)
        mx = (cm + ca) / 2;
        df = ca - cm;
        n1c = length(a.cdc_man);
        if n1c > 0
            scatter(mx(1:n1c), df(1:n1c), ms, c_ludb, 'filled', ...
                'MarkerFaceAlpha', 0.4, 'DisplayName', 'LUDB');
        end
        if ~isempty(b.cdc_man)
            scatter(mx(n1c+1:end), df(n1c+1:end), ms, c_qtdb, 'filled', ...
                'MarkerFaceAlpha', 0.4, 'DisplayName', 'QTDB');
        end
        draw_ba_lines(df);
        r_val = corr(cm, ca);
    else
        df = []; r_val = NaN;
    end
    xlabel('Mean CDC', 'FontSize', lab_fs);
    ylabel('\DeltaCDC', 'FontSize', lab_fs);
    title('c  CDC (beat-level)', 'FontWeight', 'bold', 'FontSize', ttl_fs);
    legend('Location', 'southwest', 'FontSize', leg_fs, 'Box', 'off');
    set(gca, 'FontSize', ax_fs, 'Box', 'on');
    annotate_panel(ax3, df, r_val, ann_fs, 0.025);

    %% ================================================================
    %  PANEL d: CDC Bland-Altman (subject-level median)
    %  ================================================================

    ax4 = subplot(2,2,4); hold on;
    sm = [a.subj_cdc_manual; b.subj_cdc_manual];
    sa = [a.subj_cdc_auto;   b.subj_cdc_auto];
    if length(sm) == length(sa) && ~isempty(sm)
        mx2 = (sm + sa) / 2;
        df2 = sa - sm;
        n1s = length(a.subj_cdc_manual);
        if n1s > 0
            scatter(mx2(1:n1s), df2(1:n1s), 20, c_ludb, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', 'LUDB');
        end
        if ~isempty(b.subj_cdc_manual)
            scatter(mx2(n1s+1:end), df2(n1s+1:end), 20, c_qtdb, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'DisplayName', 'QTDB');
        end
        draw_ba_lines(df2);
        r_val2 = corr(sm, sa);
    else
        df2 = []; r_val2 = NaN;
    end
    xlabel('Mean median CDC', 'FontSize', lab_fs);
    ylabel('\DeltaCDC', 'FontSize', lab_fs);
    title('d  CDC (subject-level median)', 'FontWeight', 'bold', 'FontSize', ttl_fs);
    set(gca, 'FontSize', ax_fs, 'Box', 'on');
    annotate_panel(ax4, df2, r_val2, ann_fs, 0.01);

    %% ================================================================
    %  EXPORT
    %  ================================================================

    out_pdf = fullfile(paths.figures, 'SI_Fig3_validation.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig3_validation.png');
    out_fig = fullfile(paths.figures, 'SI_Fig3_validation.fig');

    fig_w_cm = fig_w_mm / 10;
    fig_h_cm = fig_h_mm / 10;
    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 3 saved.\n');
end


%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================

function draw_ba_lines(d)
% DRAW_BA_LINES - Mean bias line + 95% limits of agreement
    if isempty(d), return; end
    yline(mean(d), 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(mean(d) + 1.96*std(d), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(mean(d) - 1.96*std(d), 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    yline(0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'HandleVisibility', 'off');
end

function annotate_panel(ax, d, r_val, fs, y_frac)
% ANNOTATE_PANEL - In-plot text: bias, MAE, and correlation
%   y_frac: fraction of y-range from top (0.05 = 5% from top)
    if isempty(d), return; end
    if nargin < 5, y_frac = 0.05; end
    axes(ax);
    xl = xlim; yl = ylim;
    x0 = xl(1) + 0.03 * diff(xl);
    y_top = yl(2) - y_frac * diff(yl);
    dy = 0.10 * diff(yl);
    
    b = mean(d); m = mean(abs(d));
    if abs(b) < 1
        text(x0, y_top,      sprintf('Bias: %+.4f', b), 'FontSize', fs, 'VerticalAlignment', 'top');
        text(x0, y_top - dy,  sprintf('MAE: %.4f', m),   'FontSize', fs, 'VerticalAlignment', 'top');
    else
        text(x0, y_top,      sprintf('Bias: %+.1f ms', b), 'FontSize', fs, 'VerticalAlignment', 'top');
        text(x0, y_top - dy,  sprintf('MAE: %.1f ms', m),   'FontSize', fs, 'VerticalAlignment', 'top');
    end
    if ~isempty(r_val) && ~isnan(r_val)
        text(x0, y_top - 2*dy, sprintf('r = %.3f', r_val), 'FontSize', fs, 'VerticalAlignment', 'top');
    end
end

function save_large_figure(fig, out_pdf, out_png, out_fig, w_cm, h_cm)
% SAVE_LARGE_FIGURE - Save figure with many graphic objects without crashing
%
% For figures with large scatter plots (>10k points), MATLAB's painters
% renderer and savefig serialize every point as a vector element, consuming
% gigabytes of memory and often hanging the process.
%
% Strategy:
%   PNG: exportgraphics (raster, always safe)
%   PDF: exportgraphics with ContentType 'image' (raster-in-PDF wrapper,
%        avoids the painters memory explosion while producing a PDF file
%        that embeds at 300 dpi — sufficient for main figures)
%   FIG: skipped for large figures (the .fig format stores all graphic
%        objects and can itself become multi-GB)
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
end
