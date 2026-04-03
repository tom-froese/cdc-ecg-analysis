function plot_SI_Fig6()
% PLOT_SI_FIG6 - Supplementary Figure 6: CODE-15% CDC distributions and
%   age trajectories
%
% Three-panel figure:
%   a: ΔCDC distribution — Clinically Normal vs Pathological (full width)
%   b: ΔCDC vs Age by clinical group (half width, bottom left)
%   c: RR interval vs Age by clinical group (half width, bottom right)
%
% Colour palette matches main figures:
%   Green   [0.25 0.70 0.35]  = Clinically Normal
%   Red     [0.85 0.25 0.20]  = Pathological
%
% Uses manual subplot positioning (not tiledlayout) to avoid
% exportgraphics coordinate-displacement bugs with text annotations.
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
    %  COLOUR PALETTE
    %  ================================================================
    col_cn   = [0.25 0.70 0.35];
    col_path = [0.85 0.25 0.20];

    %% Group masks
    is_cn   = all_data.Group == 'ClinicallyNormal';
    is_path = all_data.Group == 'Pathological';

    n_cn   = S.n_cn;
    n_path = S.n_path;
    n_total = n_cn + n_path;

    mode_cn_d   = S.mode_cn - inv_e;
    mode_path_d = S.mode_path - inv_e;

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================

    fig_w_cm = 18.3;
    fig_h_cm = 18.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 fig_w_cm fig_h_cm]);

    % Font sizes
    ax_fs    = 10;
    lab_fs   = 11;
    title_fs = 12;
    panel_fs = 14;
    leg_fs   = 9;
    ann_fs   = 10;

    % Scatter aesthetics
    dot_size  = 2;
    alpha_dot = 0.04;
    line_w    = 1.5;

    age_lim = [15 97];

    % Histogram parameters
    edges = (0.20:0.01:0.65) - inv_e;
    x_limits = [0.20 - inv_e, 0.65 - inv_e];
    kde_pts = 500;

    %% ================================================================
    %  PANEL (a): ΔCDC distribution — full width, top half
    %  ================================================================
    %  Position: [left bottom width height] in normalised figure coords
    ax1 = axes('Position', [0.09 0.57 0.86 0.37]);
    hold on; box on;

    histogram(all_data.dCDC(is_cn), edges, 'FaceColor', col_cn, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.55, 'Normalization', 'pdf');
    histogram(all_data.dCDC(is_path), edges, 'FaceColor', col_path, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.45, 'Normalization', 'pdf');

    [f_cn, x_cn] = ksdensity(all_data.dCDC(is_cn), 'NumPoints', kde_pts);
    plot(x_cn, f_cn, 'Color', col_cn * 0.7, 'LineWidth', line_w);
    [f_path, x_path] = ksdensity(all_data.dCDC(is_path), 'NumPoints', kde_pts);
    plot(x_path, f_path, 'Color', col_path * 0.7, 'LineWidth', line_w);

    yl = ylim;
    plot([0 0], [0 yl(2)], 'k-', 'LineWidth', line_w);
    plot(mode_cn_d * [1 1], [0 yl(2)], '--', 'Color', col_cn * 0.7, 'LineWidth', 1.2);
    plot(mode_path_d * [1 1], [0 yl(2)], '--', 'Color', col_path * 0.7, 'LineWidth', 1.2);

    xlabel('\DeltaCDC from optimal (1/\ite\rm \approx 0.368)', 'FontSize', lab_fs);
    ylabel('Density', 'FontSize', lab_fs);
    title('CODE-15%: Clinically Normal vs Pathological (ages 17–100)', ...
          'FontSize', title_fs, 'FontWeight', 'bold');

    ph = patch(NaN, NaN, col_cn, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    pp = patch(NaN, NaN, col_path, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    leg1 = legend([ph, pp], { ...
        sprintf('Clinically normal (n=%s, \\Delta=%+.3f)', format_comma(n_cn), mode_cn_d), ...
        sprintf('Pathological (n=%s, \\Delta=%+.3f)', format_comma(n_path), mode_path_d)}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
    set(leg1, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(x_limits); grid on;
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, 'TickDir', 'out');

    % Annotations using data coordinates (avoids exportgraphics displacement)
    axes(ax1);
    xl = xlim; yl = ylim;
    x0 = xl(1) + 0.02 * diff(xl);

    % p-value (bottom right)
    if S.p_ranksum < 0.001, p_str = '\itp\rm < 0.001';
    else, p_str = sprintf('\\itp\\rm = %.3f', S.p_ranksum); end
    text(xl(2) - 0.02*diff(xl), yl(1) + 0.05*diff(yl), p_str, ...
         'FontSize', ann_fs, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    % N and method (top left)
    text(x0, yl(2) - 0.05*diff(yl), ...
         sprintf('N = %s patients', format_comma(n_total)), ...
         'FontSize', ann_fs, 'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'top');
    text(x0, yl(2) - 0.15*diff(yl), ...
         {'Fully automatic', '(Pan-Tompkins + tangent)'}, ...
         'FontSize', ann_fs - 1, 'FontAngle', 'italic', ...
         'Color', [0.45 0.45 0.45], 'VerticalAlignment', 'top');

    % Panel label
    text(xl(1) - 0.08*diff(xl), yl(2) + 0.04*diff(yl), '\bfa', ...
         'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    hold off;

    %% ================================================================
    %  PANEL (b): ΔCDC vs Age — bottom left
    %  ================================================================
    ax2 = axes('Position', [0.09 0.08 0.40 0.38]);
    hold on; box on;

    colors = [col_cn; col_path];
    masks  = {is_cn, is_path};
    h_lines = gobjects(2, 1);

    for g = [2, 1]
        idx = masks{g};
        mdl = fitlm(all_data.Age(idx), all_data.dCDC(idx));
        xfit = linspace(age_lim(1), age_lim(2), 300);
        [yfit, ci] = predict(mdl, xfit', 'Alpha', 0.05);

        fill([xfit fliplr(xfit)], [ci(:,1)' fliplr(ci(:,2)')], colors(g,:), ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h_lines(g) = plot(xfit, yfit, 'Color', colors(g,:) * 0.75, 'LineWidth', line_w);

        scatter(all_data.Age(idx), all_data.dCDC(idx), dot_size, colors(g,:), ...
                'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    end

    yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);

    leg2 = legend(h_lines, { ...
        sprintf('Clinically Normal (n=%s)', format_comma(n_cn)), ...
        sprintf('Pathological (n=%s)', format_comma(n_path))}, ...
        'Location', 'northwest', 'FontSize', leg_fs, 'Box', 'on');
    set(leg2, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(age_lim);
    ylim([-0.15 0.28]);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, 'TickDir', 'out');

    % Optimal label (data coordinates)
    axes(ax2);
    text(age_lim(2) - 1, 0.006, 'Optimal (1/\ite\rm)', ...
        'FontSize', ann_fs - 1, 'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    % Panel label
    xl2 = xlim; yl2 = ylim;
    text(xl2(1) - 0.15*diff(xl2), yl2(2) + 0.04*diff(yl2), '\bfb', ...
         'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    hold off;

    %% ================================================================
    %  PANEL (c): RR interval vs Age — bottom right
    %  ================================================================
    ax3 = axes('Position', [0.57 0.08 0.40 0.38]);
    hold on; box on;

    h_lines2 = gobjects(2, 1);

    for g = [2, 1]
        idx = masks{g};
        mdl = fitlm(all_data.Age(idx), all_data.RR(idx));
        xfit = linspace(age_lim(1), age_lim(2), 300);
        [yfit, ci] = predict(mdl, xfit', 'Alpha', 0.05);

        fill([xfit fliplr(xfit)], [ci(:,1)' fliplr(ci(:,2)')], colors(g,:), ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h_lines2(g) = plot(xfit, yfit, 'Color', colors(g,:) * 0.75, 'LineWidth', line_w);

        scatter(all_data.Age(idx), all_data.RR(idx), dot_size, colors(g,:), ...
                'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    end

    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('Median RR interval (s)', 'FontSize', lab_fs);

    leg3 = legend(h_lines2, { ...
        sprintf('Clinically Normal (n=%s)', format_comma(n_cn)), ...
        sprintf('Pathological (n=%s)', format_comma(n_path))}, ...
        'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
    set(leg3, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);

    xlim(age_lim);
    set(ax3, 'FontSize', ax_fs, 'LineWidth', 0.5, 'TickDir', 'out');

    % Panel label
    axes(ax3);
    xl3 = xlim; yl3 = ylim;
    text(xl3(1) - 0.15*diff(xl3), yl3(2) + 0.04*diff(yl3), '\bfc', ...
         'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    hold off;

    %% ================================================================
    %  SAVE
    %  ================================================================

    out_pdf = fullfile(paths.figures, 'SI_Fig6_code15.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig6_code15.png');
    out_fig = fullfile(paths.figures, 'SI_Fig6_code15.fig');

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf('  Saved: %s (raster, 300 dpi)\n', out_png);

    exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 300);
    fprintf('  Saved: %s (raster-in-PDF, 300 dpi)\n', out_pdf);

    fprintf('  Skipped: %s (too many graphic objects)\n', out_fig);

    fprintf('\nSupplementary Figure 6 saved.\n');
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