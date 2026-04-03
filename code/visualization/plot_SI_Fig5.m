function plot_SI_Fig5()
% PLOT_SI_FIG5 - Supplementary Figure 5: CDC distributions in large-scale
%   databases with algorithmic or hybrid annotation
%
% Four-panel figure:
%   a: Fantasia        — Healthy Controls (n=40, verified volunteers)
%   b: Autonomic Aging — Healthy Controls (n~1,100, verified volunteers)
%   c: PTB             — Healthy Control vs Pathological (manual T-end)
%   d: PTB-XL          — Clinically Normal vs Pathological (ECGDeli)
%
% X-axis: ΔCDC from theoretical optimum (1/e ≈ 0.3679).
% Subject-level median CDC (computed in analyze_large_scale.m with
% uniform quality filters applied).
%
% Colour palette matches main figures (plot_Fig1.m):
%   Blue    [0.20 0.55 0.85]  = Healthy Control (verified volunteers)
%   Green   [0.25 0.70 0.35]  = Clinically Normal (hospital patients, normal ECG)
%   Red     [0.85 0.25 0.20]  = Pathological
%
% Group classification follows the hierarchical model:
%   Fantasia, Autonomic Aging, PTB 'healthy' → HC (Blue)
%   PTB-XL 'healthy'                         → CN (Green)
%   PTB, PTB-XL other                        → Pathological (Red)
%
% Data source: large_scale_results.mat (via analyze_large_scale.m)
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD PRECOMPUTED RESULTS
    %  ================================================================
    S = load(fullfile(paths.results, 'large_scale_results.mat'), 'results');
    results = S.results;

    %% ================================================================
    %  COLOUR PALETTE (matches main figures)
    %  ================================================================
    col_hc   = [0.20 0.55 0.85];   % blue  — healthy controls
    col_cn   = [0.25 0.70 0.35];   % green — clinically normal
    col_path = [0.85 0.25 0.20];   % red   — pathological

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================
    %  120 mm width, four-panel stacked layout.

    fig_w_cm = 12.0;
    fig_h_cm = 22.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 1 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters
    leg_fs   = 6;    % legend

    % Histogram and KDE parameters
    edges = (0.20:0.01:0.65) - inv_e;   % ΔCDC bin edges
    x_limits = [0.20 - inv_e, 0.65 - inv_e];
    kde_pts = 500;

    % Panel labels
    panel_letters = {'a', 'b', 'c', 'd'};

    %% ================================================================
    %  PANEL (a): Fantasia — Healthy Controls
    %  ================================================================
    ax1 = subplot(4, 1, 1);
    hold on; box on;

    fant = results.fantasia;
    d_fant = fant.all_ratios - inv_e;
    mode_fant_d = fant.mode_all - inv_e;

    histogram(d_fant, edges, 'FaceColor', col_hc, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    [f, x] = ksdensity(d_fant, 'NumPoints', kde_pts);
    plot(x, f, 'Color', col_hc * 0.7, 'LineWidth', 1.8);

    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);
    if ~isnan(mode_fant_d)
        plot(mode_fant_d * [1 1], [0 yl(2)], '--', ...
             'Color', col_hc * 0.7, 'LineWidth', 1.0);
    end

    ylabel('Density', 'FontSize', lab_fs);
    title('Fantasia: healthy volunteers (ages 21–85)', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    ph = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg1 = legend(ph, {sprintf('Healthy controls (n=%d, \\Delta=%+.3f)', ...
           fant.n_all, mode_fant_d)}, ...
           'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
    set(leg1, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(x_limits); grid on;
    xlabel('');  
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabel', []);

    text(0.03, 0.88, {'Database R-peaks,', 'tangent T-end'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);
    text(-0.12, 1.06, ['\bf' panel_letters{1}], 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): Autonomic Aging — Healthy Controls
    %  ================================================================
    ax2 = subplot(4, 1, 2);
    hold on; box on;

    aa = results.autonomic_aging;
    d_aa = aa.all_ratios - inv_e;
    mode_aa_d = aa.mode_all - inv_e;

    histogram(d_aa, edges, 'FaceColor', col_hc, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    [f, x] = ksdensity(d_aa, 'NumPoints', kde_pts);
    plot(x, f, 'Color', col_hc * 0.7, 'LineWidth', 1.8);

    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);
    if ~isnan(mode_aa_d)
        plot(mode_aa_d * [1 1], [0 yl(2)], '--', ...
             'Color', col_hc * 0.7, 'LineWidth', 1.0);
    end

    ylabel('Density', 'FontSize', lab_fs);
    title('Autonomic Aging: healthy volunteers (ages 18–92)', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    pa = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg2 = legend(pa, {sprintf('Healthy controls (n=%s, \\Delta=%+.3f)', ...
           format_comma(aa.n_all), mode_aa_d)}, ...
           'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
    set(leg2, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(x_limits); grid on;
    xlabel('');  
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabel', []);

    text(0.03, 0.88, {'Fully automatic', '(Pan-Tompkins + tangent)'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);
    text(-0.12, 1.06, ['\bf' panel_letters{2}], 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (c): PTB — Healthy Control vs Pathological
    %  ================================================================
    ax3 = subplot(4, 1, 3);
    hold on; box on;

    ptb = results.ptb;

    if ptb.n_hc >= 3
        d_hc = ptb.hc_ratios - inv_e;
        mode_hc_d = ptb.mode_hc - inv_e;
        histogram(d_hc, edges, 'FaceColor', col_hc, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
        [f, x] = ksdensity(d_hc, 'NumPoints', kde_pts);
        plot(x, f, 'Color', col_hc * 0.7, 'LineWidth', 1.8);
    end

    d_path = ptb.path_ratios - inv_e;
    mode_path_d = ptb.mode_path - inv_e;
    histogram(d_path, edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');
    [f, x] = ksdensity(d_path, 'NumPoints', kde_pts);
    plot(x, f, 'Color', col_path * 0.7, 'LineWidth', 1.8);

    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    if ptb.n_hc >= 3 && ~isnan(mode_hc_d)
        plot(mode_hc_d * [1 1], [0 yl(2)], '--', ...
             'Color', col_hc * 0.7, 'LineWidth', 1.0);
    end
    plot(mode_path_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_path * 0.7, 'LineWidth', 1.0);

    ylabel('Density', 'FontSize', lab_fs);
    title('PTB: ages 17–87', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    % Legend
    leg_items = gobjects(2, 1);
    leg_strs = cell(2, 1);
    if ptb.n_hc >= 3 && ~isnan(ptb.mode_hc)
        leg_items(1) = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Healthy control (n=%d, \\Delta=%+.3f)', ...
                              ptb.n_hc, mode_hc_d);
    else
        leg_items(1) = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Healthy control (n=%d)', ptb.n_hc);
    end
    leg_items(2) = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg_strs{2} = sprintf('Pathological (n=%d, \\Delta=%+.3f)', ...
                           ptb.n_path, mode_path_d);
    leg3 = legend(leg_items, leg_strs, 'Location', 'northeast', ...
           'FontSize', leg_fs, 'Box', 'on');
    set(leg3, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(x_limits); grid on;
    xlabel('');  
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabel', []);

    add_p_annotation(ptb.p_value, ax_fs);
    text(0.03, 0.88, {'Manual T-end (5 referees),', 'Pan-Tompkins R-peak'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);
    text(-0.12, 1.06, ['\bf' panel_letters{3}], 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (d): PTB-XL — Clinically Normal vs Pathological
    %  ================================================================
    ax4 = subplot(4, 1, 4);
    hold on; box on;

    ptbxl = results.ptbxl;

    d_cn = ptbxl.cn_ratios - inv_e;
    mode_cn_d = ptbxl.mode_cn - inv_e;
    d_path = ptbxl.path_ratios - inv_e;
    mode_path_d = ptbxl.mode_path - inv_e;

    histogram(d_cn, edges, 'FaceColor', col_cn, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(d_path, edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    [f, x] = ksdensity(d_cn, 'NumPoints', kde_pts);
    plot(x, f, 'Color', col_cn * 0.7, 'LineWidth', 1.8);
    [f, x] = ksdensity(d_path, 'NumPoints', kde_pts);
    plot(x, f, 'Color', col_path * 0.7, 'LineWidth', 1.8);

    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);
    plot(mode_cn_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_cn * 0.7, 'LineWidth', 1.0);
    plot(mode_path_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_path * 0.7, 'LineWidth', 1.0);

    % X-axis label strictly on this bottom panel
    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', ...
           'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    title('PTB-XL: ages 2–90', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    % Legend
    ph = patch(NaN, NaN, col_cn, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg4 = legend([ph, pp], { ...
        sprintf('Clinically normal (n=%s, \\Delta=%+.3f)', ...
                format_comma(ptbxl.n_cn), mode_cn_d), ...
        sprintf('Pathological (n=%s, \\Delta=%+.3f)', ...
                format_comma(ptbxl.n_path), mode_path_d)}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
    set(leg4, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(x_limits); grid on;
    set(ax4, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    add_p_annotation(ptbxl.p_value, ax_fs);
    text(0.03, 0.88, {'ECGDeli automatic', '(R-peak + T-end)'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);
    text(-0.12, 1.06, ['\bf' panel_letters{4}], 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    % Four evenly spaced panels
    panel_h = 0.19;
    gap = 0.045;
    left = 0.14;
    width = 0.80;
    bottom_start = 0.06;

    % We explicitly array the axes in order (a, b, c, d) 
    % to guarantee 'ax1' gets the top slot and 'ax4' gets the bottom slot.
    all_axes = [ax1, ax2, ax3, ax4];
    
    for k = 1:4
        set(all_axes(k), 'Position', [left, bottom_start + (4-k)*(panel_h + gap), width, panel_h]);
    end

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'SI_Fig5_large_scale.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig5_large_scale.png');
    out_fig = fullfile(paths.figures, 'SI_Fig5_large_scale.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nSupplementary Figure 5 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================
    fprintf('\n--- Summary for SI Fig 5 legend ---\n');
    fprintf('Fantasia:        HC n=%d (mode=%.3f, dCDC=%+.4f)\n', ...
            fant.n_all, fant.mode_all, fant.mode_all - inv_e);
    fprintf('Autonomic Aging: HC n=%s (mode=%.3f, dCDC=%+.4f)\n', ...
            format_comma(aa.n_all), aa.mode_all, aa.mode_all - inv_e);
    fprintf('PTB:             HC n=%d, Path n=%d, p=%.2e\n', ...
            ptb.n_hc, ptb.n_path, ptb.p_value);
    fprintf('PTB-XL:          CN n=%s, Path n=%s, p=%.2e\n', ...
            format_comma(ptbxl.n_cn), format_comma(ptbxl.n_path), ptbxl.p_value);
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function add_p_annotation(p_val, fs)
    if isnan(p_val), return; end
    if p_val < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3f', p_val);
    end
    text(0.97, 0.08, p_str, 'Units', 'normalized', ...
         'FontSize', fs, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
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

