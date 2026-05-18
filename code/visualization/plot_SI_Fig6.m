function plot_SI_Fig6()
% PLOT_SI_FIG6 - Supplementary Figure 6: CDC distributions in gold-standard
%   databases with fully manual annotation
%
% Two-panel figure:
%   a: LUDB — Healthy Controls vs Patients (pathological ECG) (full manual annotation)
%   b: QTDB — Patients (non-pathological ECG) vs Patients (pathological ECG) vs Patients (sudden death) (full manual)
%
% X-axis: ΔCDC from theoretical optimum (1/e ≈ 0.3679).
% Subject-level median CDC (computed in analyze_gold_standard.m).
%
% Colour palette matches main figures (plot_Fig1.m):
%   Blue    [0.20 0.55 0.85]  = Healthy Control (verified healthy volunteers)
%   Green   [0.25 0.70 0.35]  = Patients (non-pathological ECG; hospital patients with normal ECG)
%   Red     [0.85 0.25 0.20]  = Patients (pathological ECG)
%   Crimson [0.55 0.00 0.15]  = Patients (sudden death) (QTDB only)
%
% LUDB "Healthy" are genuine healthy volunteers (Kalyakulina et al. 2020,
% IEEE Access: "healthy volunteers and patients of the Nizhny Novgorod
% City Hospital No 5"). → Blue (Healthy Control).
%
% QTDB "healthy" are MIT-BIH Normal Sinus Rhythm subjects: hospital
% referrals to the Arrhythmia Laboratory at Beth Israel Hospital who
% were found to have no significant arrhythmias. → Green (Patients with non-pathological ECG).
%
% Data source: gold_standard_results.mat (via analyze_gold_standard.m)
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD PRECOMPUTED RESULTS
    %  ================================================================
    S = load(fullfile(paths.results, 'gold_standard_results.mat'), 'results');
    results = S.results;

    %% ================================================================
    %  COLOUR PALETTE (matches main figures)
    %  ================================================================
    col_hc   = [0.20 0.55 0.85];   % blue    — healthy controls (LUDB)
    col_npe   = [0.25 0.70 0.35];   % green   — Patients (non-pathological ECG) (QTDB healthy)
    col_pe = [0.85 0.25 0.20];   % red     — pathological
    col_sd   = [0.55 0.00 0.15];   % crimson — patients (sudden death) (QTDB only)

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================
    %  120 mm width, stacked two-panel layout.

    fig_w_cm = 18.3;
    fig_h_cm = 14.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 3 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters (a, b)
    leg_fs   = 7;    % legend (matches tick/p-value font)

    % Histogram and KDE parameters
    edges = (0.20:0.01:0.60) - inv_e;   % ΔCDC bin edges
    x_limits = [0.20 - inv_e, 0.60 - inv_e];
    kde_pts = 500;

    %% ================================================================
    %  PANEL (a): LUDB — Healthy Controls vs Patients (pathological ECG)
    %  ================================================================
    ax1 = subplot(2, 1, 1);
    hold on; box on;

    ludb = results.ludb;

    d_hc   = ludb.healthy_ratios - inv_e;
    d_pe = ludb.pe_ratios - inv_e;

    mode_hc_d   = ludb.mode_healthy - inv_e;
    mode_pe_d = ludb.mode_pe - inv_e;

    % Histograms
    histogram(d_hc, edges, 'FaceColor', col_hc, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(d_pe, edges, 'FaceColor', col_pe, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    % KDE overlays
    [f_hc, x_hc] = ksdensity(d_hc, 'NumPoints', kde_pts);
    plot(x_hc, f_hc, 'Color', col_hc * 0.7, 'LineWidth', 1.8);
    [f_pe, x_pe] = ksdensity(d_pe, 'NumPoints', kde_pts);
    plot(x_pe, f_pe, 'Color', col_pe * 0.7, 'LineWidth', 1.8);

    % 1/e reference line (ΔCDC = 0)
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    plot(mode_hc_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_hc * 0.7, 'LineWidth', 1.0);
    plot(mode_pe_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_pe * 0.7, 'LineWidth', 1.0);

    ylabel('Density', 'FontSize', lab_fs);
    t1 = title('LUDB: healthy volunteers and cardiac patients', 'FontSize', title_fs, 'FontWeight', 'bold');
    t1.Units = 'normalized'; t1.Position(2) = t1.Position(2) + 0.04;

    % Legend
    ph = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_pe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg1 = legend([ph, pp], { ...
        sprintf('Healthy controls (n=%d, \\Delta=%+.3f)', ludb.n_healthy, mode_hc_d), ...
        sprintf('Patients (pathological ECG) (n=%d, \\Delta=%+.3f)', ludb.n_pe, mode_pe_d)}, ...
        'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabel', []);

    % Annotation method note (bullet list)
    text(0.03, 0.92, {'\bullet Manual R-peaks', '\bullet Manual T-end'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45], ...
         'VerticalAlignment', 'top');

    % Panel label (Nature: bold lowercase, outside axes)
    text(-0.12, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): QTDB — Patients (non-pathological ECG) vs Patients (pathological ECG) vs Patients (sudden death)
    %  ================================================================
    ax2 = subplot(2, 1, 2);
    hold on; box on;

    qtdb = results.qtdb;

    d_npe   = qtdb.npe_ratios - inv_e;
    d_pe = qtdb.pe_ratios - inv_e;
    d_sd   = qtdb.fatal_ratios - inv_e;

    mode_npe_d   = qtdb.mode_npe - inv_e;
    mode_pe_d = qtdb.mode_pe - inv_e;
    mode_sd_d   = qtdb.mode_fatal - inv_e;

    % Histograms — layered back to front (widest distribution first)
    histogram(d_sd, edges, 'FaceColor', col_sd, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.40, 'Normalization', 'pdf');
    histogram(d_pe, edges, 'FaceColor', col_pe, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.40, 'Normalization', 'pdf');
    if length(d_npe) >= 3
        histogram(d_npe, edges, 'FaceColor', col_npe, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    end

    % KDE overlays
    if length(d_npe) >= 3
        [f_npe, x_npe] = ksdensity(d_npe, 'NumPoints', kde_pts);
        plot(x_npe, f_npe, 'Color', col_npe * 0.7, 'LineWidth', 1.8);
    end
    [f_pe, x_pe] = ksdensity(d_pe, 'NumPoints', kde_pts);
    plot(x_pe, f_pe, 'Color', col_pe * 0.7, 'LineWidth', 1.8);
    [f_sd, x_sd] = ksdensity(d_sd, 'NumPoints', kde_pts);
    plot(x_sd, f_sd, 'Color', col_sd * 0.8, 'LineWidth', 1.8);

    % 1/e reference line
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    if length(d_npe) >= 3 && ~isnan(mode_npe_d)
        plot(mode_npe_d * [1 1], [0 yl(2)], '--', ...
             'Color', col_npe * 0.7, 'LineWidth', 1.0);
    end
    plot(mode_pe_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_pe * 0.7, 'LineWidth', 1.0);
    plot(mode_sd_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_sd * 0.8, 'LineWidth', 1.0);

    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', ...
           'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    t2 = title('QTDB: cardiac and sudden death patients', 'FontSize', title_fs, 'FontWeight', 'bold');
    t2.Units = 'normalized'; t2.Position(2) = t2.Position(2) + 0.04;

    % Legend — three groups
    leg_items = gobjects(3, 1);
    leg_strs  = cell(3, 1);

    if length(d_npe) >= 3 && ~isnan(mode_npe_d)
        leg_items(1) = patch(NaN, NaN, col_npe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Patients (non-pathological ECG) (n=%d, \\Delta=%+.3f)', ...
                              qtdb.n_npe, mode_npe_d);
    else
        leg_items(1) = patch(NaN, NaN, col_npe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Patients (non-pathological ECG) (n=%d)', qtdb.n_npe);
    end
    leg_items(2) = patch(NaN, NaN, col_pe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg_strs{2} = sprintf('Patients (pathological ECG) (n=%d, \\Delta=%+.3f)', ...
                           qtdb.n_pe, mode_pe_d);
    leg_items(3) = patch(NaN, NaN, col_sd, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg_strs{3} = sprintf('Patients (sudden death) (n=%d, \\Delta=%+.3f)', ...
                           qtdb.n_fatal, mode_sd_d);

    leg2 = legend(leg_items, leg_strs, ...
           'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    % Annotation method note (bullet list)
    text(0.03, 0.92, {'\bullet Manual R-peaks', '\bullet Manual T-end'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45], ...
         'VerticalAlignment', 'top');

    % Panel label
    text(-0.12, 1.06, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    set(ax1, 'Position', [0.10  0.56  0.48  0.37]);
    set(ax2, 'Position', [0.10  0.10  0.48  0.37]);
    % Legends — left-edge aligned, auto-width
    leg_left = 0.60;
    legs = [leg1, leg2];
    panel_bottoms = [0.56, 0.10];
    panel_h = 0.37;
    for k = 1:2
        set(legs(k), 'Units', 'normalized');
        auto_pos = get(legs(k), 'Position');
        leg_y = panel_bottoms(k) + (panel_h - auto_pos(4)) / 2;
        set(legs(k), 'Position', ...
            [leg_left, leg_y, auto_pos(3), auto_pos(4)]);
    end

    % p-value annotations below each legend
    if ludb.p_value < 0.001
        p1_str = 'Wilcoxon rank-sum p < 0.001';
    else
        p1_str = sprintf('Wilcoxon rank-sum p = %.3f', ludb.p_value);
    end
    annotation('textbox', [leg_left 0.64 0.38 0.05], 'String', p1_str, ...
        'FontSize', ax_fs, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'FitBoxToText', 'on');

    if qtdb.p_kruskal < 0.001
        p2_str = 'Kruskal-Wallis p < 0.001';
    else
        p2_str = sprintf('Kruskal-Wallis p = %.3f', qtdb.p_kruskal);
    end
    annotation('textbox', [leg_left 0.16 0.38 0.05], 'String', p2_str, ...
        'FontSize', ax_fs, 'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'FitBoxToText', 'on');

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'SI_Fig6_gold_standard.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig6_gold_standard.png');
    out_fig = fullfile(paths.figures, 'SI_Fig6_gold_standard.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nSupplementary Figure 6 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================
    fprintf('\n--- Summary for SI Fig 6 legend ---\n');
    fprintf('LUDB: Healthy controls n=%d (mode=%.3f, dCDC=%+.4f), ', ...
            ludb.n_healthy, ludb.mode_healthy, ludb.mode_healthy - inv_e);
    fprintf('Patients (pathological ECG) n=%d (mode=%.3f, dCDC=%+.4f)\n', ...
            ludb.n_pe, ludb.mode_pe, ludb.mode_pe - inv_e);
    fprintf('  Wilcoxon rank-sum p = %.2e\n', ludb.p_value);

    fprintf('QTDB: Patients (non-pathological ECG) n=%d (mode=%.3f), ', ...
            qtdb.n_npe, qtdb.mode_npe);
    fprintf('Patients (pathological ECG) n=%d (mode=%.3f), ', ...
            qtdb.n_pe, qtdb.mode_pe);
    fprintf('Patients (sudden death) n=%d (mode=%.3f)\n', ...
            qtdb.n_fatal, qtdb.mode_fatal);
    fprintf('  Kruskal-Wallis p = %.2e\n', qtdb.p_kruskal);
end


