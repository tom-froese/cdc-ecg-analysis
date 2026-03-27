function plot_SI_Fig4()
% PLOT_SI_FIG4 - Supplementary Figure 4: CDC distributions in gold-standard
%   databases with fully manual annotation
%
% Two-panel figure:
%   a: LUDB — Healthy Controls vs Pathological (full manual annotation)
%   b: QTDB — Clinically Normal vs Pathological vs Sudden Death (full manual)
%
% X-axis: ΔCDC from theoretical optimum (1/e ≈ 0.3679).
% Subject-level median CDC (computed in analyze_gold_standard.m).
%
% Colour palette matches main figures (plot_Fig1.m):
%   Blue    [0.20 0.55 0.85]  = Healthy Control (verified healthy volunteers)
%   Green   [0.25 0.70 0.35]  = Clinically Normal (hospital patients, normal ECG)
%   Red     [0.85 0.25 0.20]  = Pathological
%   Crimson [0.55 0.00 0.15]  = Sudden Death (QTDB only)
%
% LUDB "Healthy" are genuine healthy volunteers (Kalyakulina et al. 2020,
% IEEE Access: "healthy volunteers and patients of the Nizhny Novgorod
% City Hospital No 5"). → Blue (Healthy Control).
%
% QTDB "healthy" are MIT-BIH Normal Sinus Rhythm subjects: hospital
% referrals to the Arrhythmia Laboratory at Beth Israel Hospital who
% were found to have no significant arrhythmias. → Green (Clinically Normal).
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
    col_cn   = [0.25 0.70 0.35];   % green   — clinically normal (QTDB healthy)
    col_path = [0.85 0.25 0.20];   % red     — pathological
    col_sd   = [0.55 0.00 0.15];   % crimson — sudden death (QTDB only)

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================
    %  120 mm width, stacked two-panel layout.

    fig_w_cm = 12.0;
    fig_h_cm = 14.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 3 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters (a, b)
    leg_fs   = 6.5;  % legend

    % Histogram and KDE parameters
    edges = (0.20:0.01:0.60) - inv_e;   % ΔCDC bin edges
    x_limits = [0.20 - inv_e, 0.60 - inv_e];
    kde_pts = 500;

    %% ================================================================
    %  PANEL (a): LUDB — Healthy Controls vs Pathological
    %  ================================================================
    ax1 = subplot(2, 1, 1);
    hold on; box on;

    ludb = results.ludb;

    d_hc   = ludb.healthy_ratios - inv_e;
    d_path = ludb.patient_ratios - inv_e;

    mode_hc_d   = ludb.mode_healthy - inv_e;
    mode_path_d = ludb.mode_patient - inv_e;

    % Histograms
    histogram(d_hc, edges, 'FaceColor', col_hc, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(d_path, edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    % KDE overlays
    [f_hc, x_hc] = ksdensity(d_hc, 'NumPoints', kde_pts);
    plot(x_hc, f_hc, 'Color', col_hc * 0.7, 'LineWidth', 1.8);
    [f_path, x_path] = ksdensity(d_path, 'NumPoints', kde_pts);
    plot(x_path, f_path, 'Color', col_path * 0.7, 'LineWidth', 1.8);

    % 1/e reference line (ΔCDC = 0)
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    plot(mode_hc_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_hc * 0.7, 'LineWidth', 1.0);
    plot(mode_path_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_path * 0.7, 'LineWidth', 1.0);

    ylabel('Density', 'FontSize', lab_fs);
    title('LUDB: full manual annotation', 'FontSize', title_fs, 'FontWeight', 'bold');

    % Legend
    ph = patch(NaN, NaN, col_hc, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    legend([ph, pp], { ...
        sprintf('Healthy controls (n=%d, \\Delta=%+.3f)', ludb.n_healthy, mode_hc_d), ...
        sprintf('Pathological (n=%d, \\Delta=%+.3f)', ludb.n_patient, mode_path_d)}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02], 'XTickLabel', []);

    % p-value annotation
    add_p_annotation(ludb.p_value, ax_fs);

    % Annotation method note
    text(0.03, 0.88, 'Manual R-peaks, manual T-end', ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);

    % Panel label (Nature: bold lowercase, outside axes)
    text(-0.12, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): QTDB — Clinically Normal vs Pathological vs Sudden Death
    %  ================================================================
    ax2 = subplot(2, 1, 2);
    hold on; box on;

    qtdb = results.qtdb;

    d_cn   = qtdb.normal_ratios - inv_e;
    d_path = qtdb.pathological_ratios - inv_e;
    d_sd   = qtdb.fatal_ratios - inv_e;

    mode_cn_d   = qtdb.mode_normal - inv_e;
    mode_path_d = qtdb.mode_pathological - inv_e;
    mode_sd_d   = qtdb.mode_fatal - inv_e;

    % Histograms — layered back to front (widest distribution first)
    histogram(d_sd, edges, 'FaceColor', col_sd, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.40, 'Normalization', 'pdf');
    histogram(d_path, edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.40, 'Normalization', 'pdf');
    if length(d_cn) >= 3
        histogram(d_cn, edges, 'FaceColor', col_cn, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    end

    % KDE overlays
    if length(d_cn) >= 3
        [f_cn, x_cn] = ksdensity(d_cn, 'NumPoints', kde_pts);
        plot(x_cn, f_cn, 'Color', col_cn * 0.7, 'LineWidth', 1.8);
    end
    [f_path, x_path] = ksdensity(d_path, 'NumPoints', kde_pts);
    plot(x_path, f_path, 'Color', col_path * 0.7, 'LineWidth', 1.8);
    [f_sd, x_sd] = ksdensity(d_sd, 'NumPoints', kde_pts);
    plot(x_sd, f_sd, 'Color', col_sd * 0.8, 'LineWidth', 1.8);

    % 1/e reference line
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    if length(d_cn) >= 3 && ~isnan(mode_cn_d)
        plot(mode_cn_d * [1 1], [0 yl(2)], '--', ...
             'Color', col_cn * 0.7, 'LineWidth', 1.0);
    end
    plot(mode_path_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_path * 0.7, 'LineWidth', 1.0);
    plot(mode_sd_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_sd * 0.8, 'LineWidth', 1.0);

    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', ...
           'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    title('QTDB: full manual annotation', 'FontSize', title_fs, 'FontWeight', 'bold');

    % Legend — three groups
    leg_items = gobjects(3, 1);
    leg_strs  = cell(3, 1);

    if length(d_cn) >= 3 && ~isnan(mode_cn_d)
        leg_items(1) = patch(NaN, NaN, col_cn, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Clinically normal (n=%d, \\Delta=%+.3f)', ...
                              qtdb.n_normal, mode_cn_d);
    else
        leg_items(1) = patch(NaN, NaN, col_cn, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        leg_strs{1} = sprintf('Clinically normal (n=%d)', qtdb.n_normal);
    end
    leg_items(2) = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg_strs{2} = sprintf('Pathological (n=%d, \\Delta=%+.3f)', ...
                           qtdb.n_pathological, mode_path_d);
    leg_items(3) = patch(NaN, NaN, col_sd, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg_strs{3} = sprintf('Sudden death (n=%d, \\Delta=%+.3f)', ...
                           qtdb.n_fatal, mode_sd_d);

    legend(leg_items, leg_strs, 'Location', 'northeast', ...
           'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    % Kruskal-Wallis p-value
    if qtdb.p_kruskal < 0.001
        p_str = 'Kruskal-Wallis p < 0.001';
    else
        p_str = sprintf('Kruskal-Wallis p = %.3f', qtdb.p_kruskal);
    end
    text(0.97, 0.08, p_str, 'Units', 'normalized', ...
         'FontSize', ax_fs, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    % Annotation method note
    text(0.03, 0.88, 'Manual R-peaks, manual T-end', ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);

    % Panel label
    text(-0.12, 1.06, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    set(ax1, 'Position', [0.14  0.56  0.80  0.37]);
    set(ax2, 'Position', [0.14  0.10  0.80  0.37]);

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'SI_Fig4_gold_standard.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig4_gold_standard.png');
    out_fig = fullfile(paths.figures, 'SI_Fig4_gold_standard.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nSupplementary Figure 4 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================
    fprintf('\n--- Summary for SI Fig 4 legend ---\n');
    fprintf('LUDB: Healthy controls n=%d (mode=%.3f, dCDC=%+.4f), ', ...
            ludb.n_healthy, ludb.mode_healthy, ludb.mode_healthy - inv_e);
    fprintf('Pathological n=%d (mode=%.3f, dCDC=%+.4f)\n', ...
            ludb.n_patient, ludb.mode_patient, ludb.mode_patient - inv_e);
    fprintf('  Wilcoxon rank-sum p = %.2e\n', ludb.p_value);

    fprintf('QTDB: Clinically normal n=%d (mode=%.3f), ', ...
            qtdb.n_normal, qtdb.mode_normal);
    fprintf('Pathological n=%d (mode=%.3f), ', ...
            qtdb.n_pathological, qtdb.mode_pathological);
    fprintf('Sudden death n=%d (mode=%.3f)\n', ...
            qtdb.n_fatal, qtdb.mode_fatal);
    fprintf('  Kruskal-Wallis p = %.2e\n', qtdb.p_kruskal);
end


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function add_p_annotation(p_val, fs)
% ADD_P_ANNOTATION - Place a formatted p-value in the lower right
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
