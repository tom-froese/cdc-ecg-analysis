function plot_SI_Fig6()
% PLOT_SI_FIG6 - Supplementary Figure 6: CODE-15% CDC distributions and
%   age trajectories
%
% Three-panel figure:
%   a: ΔCDC distribution — Clinically Normal vs Pathological
%   b: ΔCDC vs Age by clinical group (scatter + OLS trend)
%   c: RR interval vs Age by clinical group (scatter + OLS trend)
%
% Colour palette matches main figures (plot_Fig1.m):
%   Green   [0.25 0.70 0.35]  = Clinically Normal
%   Red     [0.85 0.25 0.20]  = Pathological
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
             'all_data', 'mode_cn', 'mode_path', 'p_ranksum', ...
             'n_cn', 'n_path', 'inv_e');
    all_data = S.all_data;
    all_data.dCDC = all_data.CDC - inv_e;
    all_data.RR   = 60 ./ all_data.HR;   % bpm → seconds

    %% ================================================================
    %  COLOUR PALETTE (matches main figures)
    %  ================================================================
    col_cn   = [0.25 0.70 0.35];   % green — clinically normal
    col_path = [0.85 0.25 0.20];   % red   — pathological

    %% Group masks
    is_cn   = all_data.Group == 'ClinicallyNormal';
    is_path = all_data.Group == 'Pathological';

    n_cn   = S.n_cn;
    n_path = S.n_path;
    n_total = n_cn + n_path;

    mode_cn_d   = S.mode_cn - inv_e;
    mode_path_d = S.mode_path - inv_e;

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
    leg_fs   = 6;    % legend

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

    edges = (0.20:0.01:0.65) - inv_e;
    x_limits = [0.20 - inv_e, 0.65 - inv_e];
    kde_pts = 500;

    % Histograms
    histogram(all_data.dCDC(is_cn), edges, 'FaceColor', col_cn, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(all_data.dCDC(is_path), edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    % KDE overlays
    [f_cn, x_cn] = ksdensity(all_data.dCDC(is_cn), 'NumPoints', kde_pts);
    plot(x_cn, f_cn, 'Color', col_cn * 0.7, 'LineWidth', 1.8);
    [f_path, x_path] = ksdensity(all_data.dCDC(is_path), 'NumPoints', kde_pts);
    plot(x_path, f_path, 'Color', col_path * 0.7, 'LineWidth', 1.8);

    % 1/e reference line
    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', 1.8);

    % Mode dashed lines
    plot(mode_cn_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_cn * 0.7, 'LineWidth', 1.0);
    plot(mode_path_d * [1 1], [0 yl(2)], '--', ...
         'Color', col_path * 0.7, 'LineWidth', 1.0);

    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', 'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    title('CODE-15%: Clinically Normal vs Pathological', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    % Legend
    ph = patch(NaN, NaN, col_cn, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    legend([ph, pp], { ...
        sprintf('Clinically normal (n=%s, \\Delta=%+.3f)', format_comma(n_cn), mode_cn_d), ...
        sprintf('Pathological (n=%s, \\Delta=%+.3f)', format_comma(n_path), mode_path_d)}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'off');

    xlim(x_limits); grid on;
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    % p-value
    add_p_annotation(S.p_ranksum, ax_fs);

    % N and method annotations
    text(0.03, 0.88, sprintf('N = %s patients', format_comma(n_total)), ...
         'Units', 'normalized', 'FontSize', ax_fs, 'Color', [0.3 0.3 0.3]);
    text(0.03, 0.78, 'Fully automatic (Pan-Tompkins + tangent)', ...
         'Units', 'normalized', 'FontSize', ax_fs - 1, ...
         'FontAngle', 'italic', 'Color', [0.45 0.45 0.45]);

    text(-0.08, 1.06, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): ΔCDC vs Age
    %  ================================================================
    ax2 = subplot(2, 2, 3);
    hold on; box on;

    groups = {'ClinicallyNormal', 'Pathological'};
    colors = [col_cn; col_path];
    masks  = {is_cn, is_path};

    h_lines = gobjects(2, 1);

    % Plot back to front: pathological first, then CN on top
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
    text(age_lim(2) - 1, 0.004, 'Optimal (1/\ite\rm)', ...
        'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);

    legend(h_lines, { ...
        sprintf('Clinically Normal (n=%s)', format_comma(n_cn)), ...
        sprintf('Pathological (n=%s)', format_comma(n_path))}, ...
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
        sprintf('Clinically Normal (n=%s)', format_comma(n_cn)), ...
        sprintf('Pathological (n=%s)', format_comma(n_path))}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'off');

    xlim(age_lim);
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    text(-0.16, 1.06, '\bfc', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================
    % Top panel spans full width; bottom two panels side by side
    set(ax1, 'Position', [0.10  0.58  0.85  0.36]);
    set(ax2, 'Position', [0.10  0.08  0.38  0.40]);
    set(ax3, 'Position', [0.58  0.08  0.38  0.40]);

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'SI_Fig6_code15.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig6_code15.png');
    out_fig = fullfile(paths.figures, 'SI_Fig6_code15.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 6 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================
    fprintf('\n--- Summary for SI Fig 6 legend ---\n');
    fprintf('N = %s patients (%s CN, %s Path)\n', ...
            format_comma(n_total), format_comma(n_cn), format_comma(n_path));
    fprintf('CN mode: %.3f (dCDC=%+.4f)\n', S.mode_cn, mode_cn_d);
    fprintf('Path mode: %.3f (dCDC=%+.4f)\n', S.mode_path, mode_path_d);
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
