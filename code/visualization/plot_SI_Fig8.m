function plot_SI_Fig8()
% PLOT_SI_FIG8 - Supplementary Figure 8: CODE-15% CDC distributions and
%   age trajectories
%
% Three-panel figure:
%   a: ΔCDC distribution — Patients (non-pathological ECG) vs Patients (pathological ECG)
%   b: ΔCDC vs Age by clinical group (scatter + OLS trend)
%   c: RR interval vs Age by clinical group (scatter + OLS trend)
%
% Colour palette matches main figures (plot_Fig1.m):
%   Green   [0.25 0.70 0.35]  = Patients (non-pathological ECG)
%   Red     [0.85 0.25 0.20]  = Patients (pathological ECG)
%
% Data source: code15_results.mat (via analyze_code15.m)
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD PRECOMPUTED RESULTS
    %  ================================================================
    S = load(fullfile(paths.results, 'code15_results.mat'), ...
             'all_data', 'mode_npe', 'mode_pe', 'p_ranksum', ...
             'n_npe', 'n_pe', 'inv_e');
    all_data = S.all_data;
    all_data.dCDC = all_data.CDC - inv_e;
    all_data.RR   = 60 ./ all_data.HR;   % bpm → seconds

    %% ================================================================
    %  COLOUR PALETTE (matches main figures)
    %  ================================================================
    col_npe   = [0.25 0.70 0.35];   % green — Patients (non-pathological ECG)
    col_pe = [0.85 0.25 0.20];   % red   — pathological

    %% Group masks
    is_npe   = all_data.Group == 'NonPathECG';
    is_pe = all_data.Group == 'PathECG';

    n_npe   = S.n_npe;
    n_pe = S.n_pe;
    n_total = n_npe + n_pe;

    mode_npe_d   = S.mode_npe - inv_e;
    mode_pe_d = S.mode_pe - inv_e;

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================
    fig_w_cm = 18.3;   % double-column width
    fig_h_cm = 16.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 fig_w_cm fig_h_cm]);

    % Font sizes (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters
    leg_fs   = 7;    % legend (matches tick/p-value font)

    % Scatter aesthetics
    dot_size  = 2;
    alpha_dot = 0.04;

    % Shared age limits
    age_lim = [15 97];

    %% ================================================================
    %  PANEL (a): ΔCDC distribution
    %  ================================================================
    ax1 = subplot(2, 2, [1 2]);
    hold on; box on;

    edges = (0.17:0.01:0.65) - inv_e;
    x_limits = [0.17 - inv_e, 0.65 - inv_e];
    kde_pts = 500;

    % Histograms
    histogram(all_data.dCDC(is_npe), edges, 'FaceColor', col_npe, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(all_data.dCDC(is_pe), edges, 'FaceColor', col_pe, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    % KDE overlays
    [f_npe, x_npe] = ksdensity(all_data.dCDC(is_npe), 'NumPoints', kde_pts);
    plot(x_npe, f_npe, 'Color', col_npe * 0.7, 'LineWidth', 1.8);
    [f_pe, x_pe] = ksdensity(all_data.dCDC(is_pe), 'NumPoints', kde_pts);
    plot(x_pe, f_pe, 'Color', col_pe * 0.7, 'LineWidth', 1.8);

    % 1/e reference line
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    plot(mode_npe_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_npe * 0.7, 'LineWidth', 1.0);
    plot(mode_pe_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_pe * 0.7, 'LineWidth', 1.0);

    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', 'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    t1 = title(sprintf('CODE-15%%: clinical ECG patients (N = %s)', format_comma(n_total)), ...
          'FontSize', title_fs, 'FontWeight', 'bold');
    t1.Units = 'normalized'; t1.Position(2) = t1.Position(2) + 0.04;

    % Legend (positioned in layout section)
    ph = patch(NaN, NaN, col_npe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_pe, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg1 = legend([ph, pp], { ...
        sprintf('Patients (non-pathological ECG)\n(n=%s, \\Delta=%+.3f)', format_comma(n_npe), mode_npe_d), ...
        sprintf('Patients (pathological ECG)\n(n=%s, \\Delta=%+.3f)', format_comma(n_pe), mode_pe_d)}, ...
        'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    % Annotation method (bullet list)
    text(0.03, 0.92, {'\bullet Pan-Tompkins R-peaks', '\bullet Tangent T-end'}, ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45], ...
         'VerticalAlignment', 'top');

    text(-0.12, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): ΔCDC vs Age
    %  ================================================================
    ax2 = subplot(2, 2, 3);
    hold on; box on;

    groups = {'NonPathECG', 'PathECG'};
    colors = [col_npe; col_pe];
    masks  = {is_npe, is_pe};

    h_lines = gobjects(2, 1);

    % Plot back to front: pathological first, then NPE on top
    for g = [2, 1]
        idx = masks{g};

        % OLS trend line with CI
        mdl = fitlm(all_data.Age(idx), all_data.dCDC(idx));
        xfit = linspace(age_lim(1), age_lim(2), 300);
        [yfit, ci] = predict(mdl, xfit', 'Alpha', 0.05);

        fill([xfit fliplr(xfit)], [ci(:,1)' fliplr(ci(:,2)')], colors(g,:), ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h_lines(g) = plot(xfit, yfit, 'Color', colors(g,:) * 0.75, 'LineWidth', 1.8);

        scatter(all_data.Age(idx), all_data.dCDC(idx), dot_size, colors(g,:), ...
                'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    end

    % Reference line
    yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    text(age_lim(2) - 1, -0.008, 'Optimal (1/\ite\rm)', ...
        'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 2);

    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);

    legend(h_lines, { ...
        sprintf('Patients (non-pathological ECG) (n=%s)', format_comma(n_npe)), ...
        sprintf('Patients (pathological ECG) (n=%s)', format_comma(n_pe))}, ...
        'Location', 'northwest', 'FontSize', leg_fs, 'Box', 'off');

    xlim(age_lim);
    ylim([-0.15 0.28]);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    text(-0.16, 1.06, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (c): RR interval vs Age
    %  ================================================================
    ax3 = subplot(2, 2, 4);
    hold on; box on;

    h_lines2 = gobjects(2, 1);

    for g = [2, 1]
        idx = masks{g};

        mdl = fitlm(all_data.Age(idx), all_data.RR(idx));
        xfit = linspace(age_lim(1), age_lim(2), 300);
        [yfit, ci] = predict(mdl, xfit', 'Alpha', 0.05);

        fill([xfit fliplr(xfit)], [ci(:,1)' fliplr(ci(:,2)')], colors(g,:), ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h_lines2(g) = plot(xfit, yfit, 'Color', colors(g,:) * 0.75, 'LineWidth', 1.8);

        scatter(all_data.Age(idx), all_data.RR(idx), dot_size, colors(g,:), ...
                'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    end

    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('Median RR interval (s)', 'FontSize', lab_fs);

    legend(h_lines2, { ...
        sprintf('Patients (non-pathological ECG) (n=%s)', format_comma(n_npe)), ...
        sprintf('Patients (pathological ECG) (n=%s)', format_comma(n_pe))}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'off');

    xlim(age_lim);
    ylim([0.3 1.8]);   % cap at 1.8 s (>99.5th percentile; excludes 20 extreme outliers)
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    text(-0.16, 1.06, '\bfc', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    % Top panel narrowed for legend; bottom two panels side by side
    panel_a_pos = [0.10  0.58  0.50  0.36];
    set(ax1, 'Position', panel_a_pos);
    set(ax2, 'Position', [0.10  0.08  0.38  0.40]);
    set(ax3, 'Position', [0.58  0.08  0.38  0.40]);

    % Legend — right side of panel (a), vertically centered
    leg_left_a = 0.62;
    set(leg1, 'Units', 'normalized');
    auto_pos = get(leg1, 'Position');
    leg_y = panel_a_pos(2) + (panel_a_pos(4) - auto_pos(4)) / 2;
    set(leg1, 'Position', [leg_left_a, leg_y, auto_pos(3), auto_pos(4)]);

    % p-value below legend
    if S.p_ranksum < 0.001
        p1_str = 'Wilcoxon rank-sum \itp\rm < 0.001';
    else
        p1_str = sprintf('Wilcoxon rank-sum \\itp\\rm = %.3f', S.p_ranksum);
    end
    annotation('textbox', [leg_left_a, leg_y - 0.03, 0.36, 0.03], ...
        'String', p1_str, 'FontSize', ax_fs, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on');

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'SI_Fig8_code15.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig8_code15.png');
    out_fig = fullfile(paths.figures, 'SI_Fig8_code15.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 8 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================
    fprintf('\n--- Summary for SI Fig 8 legend ---\n');
    fprintf('N = %s patients (%s NPE, %s PE)\n', ...
            format_comma(n_total), format_comma(n_npe), format_comma(n_pe));
    fprintf('NPE mode: %.3f (dCDC=%+.4f)\n', S.mode_npe, mode_npe_d);
    fprintf('PE mode: %.3f (dCDC=%+.4f)\n', S.mode_pe, mode_pe_d);
    fprintf('Wilcoxon p = %.2e\n', S.p_ranksum);

    % OLS slopes for legend
    fprintf('\nOLS age slopes (for legend):\n');
    for g = 1:2
        idx = masks{g};
        mdl = fitlm(all_data.Age(idx), all_data.dCDC(idx));
        coeffs = mdl.Coefficients;
        slope = coeffs.Estimate(2);
        p_slope = coeffs.pValue(2);
        fprintf('  %-20s  dCDC: %+.5f/yr (p=%.2e)\n', groups{g}, slope, p_slope);
    end
    for g = 1:2
        idx = masks{g};
        mdl = fitlm(all_data.Age(idx), all_data.RR(idx));
        coeffs = mdl.Coefficients;
        slope = coeffs.Estimate(2);
        p_slope = coeffs.pValue(2);
        fprintf('  %-20s  RR:   %+.5f s/yr (p=%.2e)\n', groups{g}, slope, p_slope);
    end
end


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================


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
%        that embeds at 300 dpi — sufficient for SI figures)
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
