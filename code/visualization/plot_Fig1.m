function plot_Fig1()
% PLOT_FIG1 - Figure 1: The cardiac duty cycle converges on 1/e in
%   healthy hearts and deviates with pathological aging
%
%   Panel a: ΔCDC from the theoretical optimum (1/e) across the adult
%     lifespan, by clinical group. Individual subjects shown as
%     semi-transparent dots; linear regression with 95% CI shading
%     per group. 
%
%   Panel b: Thermodynamic diastole duration (ms) across the adult
%     lifespan, by clinical group. Derived from subject-level median
%     CDC and HR as: diastole = RR × (1 − CDC), where RR = 60000/HR.
%
% Nature Aging formatting: double-column (183 mm), bold lowercase
%   panel labels outside axes.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD DATA
    %  ================================================================

    S = load(fullfile(paths.results, 'hierarchical_results.mat'), ...
             'all_data', 'group_modes', 'group_mode_cis', 'inv_e');
    D = S.all_data;

    % Bootstrap mode results (order: HC, CN, Path)
    group_modes    = S.group_modes;
    group_mode_cis = S.group_mode_cis;

    % Derived quantities
    D.delta_CDC = D.CDC - inv_e;              
    D.RR_ms     = 60000 ./ D.HR;              
    D.Dias_ms   = D.RR_ms .* (1 - D.CDC);    

    % Group masks
    is_hc   = D.Group == 'HealthyControl';
    is_cn   = D.Group == 'ClinicallyNormal';
    is_path = D.Group == 'Pathological';

    n_hc   = sum(is_hc);
    n_cn   = sum(is_cn);
    n_path = sum(is_path);

    fprintf('Figure 1 — data loaded:\n');
    fprintf('  Healthy Control:   N = %s\n', format_comma(n_hc));
    fprintf('  Clinically Normal: N = %s\n', format_comma(n_cn));
    fprintf('  Pathological:      N = %s\n', format_comma(n_path));

    %% ================================================================
    %  COLOUR SCHEME & LABELS
    %  ================================================================

    col_hc   = [0.20 0.55 0.85];   % Blue
    col_cn   = [0.25 0.70 0.35];   % Green
    col_path = [0.85 0.25 0.20];   % Red 

    colors = [col_hc; col_cn; col_path];
    masks  = {is_hc, is_cn, is_path};
    ns     = [n_hc, n_cn, n_path];

    group_names = {'Healthy Control', 'Clinically Normal', 'Pathological'};
    labels_a = cell(3, 1);
    for g = 1:3
        m  = group_modes(g);
        lo = group_mode_cis(g, 1);
        hi = group_mode_cis(g, 2);
        in_ci = inv_e >= lo && inv_e <= hi;
        if in_ci
            ci_note = ', 1/\ite\rm \in CI';
        else
            ci_note = '';
        end
        % Splitting into two lines with \n to prevent horizontal stretching
        labels_a{g} = sprintf('%s (\\itn\\rm = %s)\nmode: %.3f [%.3f, %.3f]%s', ...
            group_names{g}, format_comma(ns(g)), m, lo, hi, ci_note);
    end

    labels_b = { ...
        sprintf('Healthy Control (\\itn\\rm = %s)', format_comma(n_hc)), ...
        sprintf('Clinically Normal (\\itn\\rm = %s)', format_comma(n_cn)), ...
        sprintf('Pathological (\\itn\\rm = %s)', format_comma(n_path))};

    %% ================================================================
    %  FIGURE SETUP (SCALED FOR LEGIBILITY)
    %  ================================================================
    
    % Draw the figure at 2x scale. When scaled to 50% in print, 
    % fonts will naturally hit the required Nature standard sizes.
    fig_scale = 2; 
    
    fig_w_cm = 18.3 * fig_scale;
    fig_h_cm = 9.0 * fig_scale; 

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 fig_w_cm fig_h_cm]);

    % Scaled font and line parameters
    ax_fs    = 8 * fig_scale;    % tick labels 
    lab_fs   = 9 * fig_scale;    % axis labels 
    panel_fs = 11 * fig_scale;   % panel labels 
    leg_fs   = 8 * fig_scale;    % legend labels
    
    alpha_dot = 0.08;
    dot_size  = 4 * fig_scale;
    line_w    = 1.5 * fig_scale;
    ax_line_w = 0.5 * fig_scale;

    age_lim = [15 95];
    xfit = linspace(age_lim(1), age_lim(2), 300)';

    %% ================================================================
    %  PANEL (a): ΔCDC from optimal (1/e) vs Age
    %  ================================================================

    ax1 = subplot(1, 2, 1);
    hold on; box on;

    scatter(D.Age(is_path), D.delta_CDC(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.delta_CDC(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.delta_CDC(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, 'HandleVisibility', 'off');

    yline(0, 'k--', 'LineWidth', line_w * 0.7, 'HandleVisibility', 'off');
    text(ax1, age_lim(2) - 1, -0.008, 'Optimal (1/\ite\rm)', ...
        'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    h_lines = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.delta_CDC(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h_lines(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, 'LineWidth', line_w);
    end

    xlim(age_lim);
    ylim([-0.08 0.12]);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', ax_line_w, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015]);

    % Panel label A 
    text(ax1, -0.12, 1.04, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'Clipping', 'off', 'VerticalAlignment', 'bottom');

    % Legend placed in the bottom right corner
    leg1 = legend(ax1, h_lines, labels_a, 'Location', 'southeast', 'Box', 'on');
    set(leg1, 'FontSize', leg_fs, 'Color', 'w');

    hold off;

    %% ================================================================
    %  PANEL (b): Thermodynamic diastole (ms) vs Age
    %  ================================================================

    ax2 = subplot(1, 2, 2);
    hold on; box on;

    scatter(D.Age(is_path), D.Dias_ms(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.Dias_ms(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, 'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.Dias_ms(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, 'HandleVisibility', 'off');

    h_lines2 = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.Dias_ms(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h_lines2(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, 'LineWidth', line_w);
    end

    xlim(age_lim);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('Thermodynamic diastole (ms)', 'FontSize', lab_fs);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', ax_line_w, ...
        'TickDir', 'out', 'TickLength', [0.015 0.015]);

    leg2 = legend(ax2, h_lines2, labels_b, 'Location', 'northeast', ...
        'FontSize', leg_fs, 'Box', 'off');
    leg2.ItemTokenSize = [15 15];

    % Panel label B 
    text(ax2, -0.12, 1.04, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'Clipping', 'off', 'VerticalAlignment', 'bottom');

    hold off;

    %% ================================================================
    %  EXPLICIT LAYOUT AND SAVE
    %  ================================================================

    % Hardcode bounds [left bottom width height] to guarantee no overlap
    ax1_pos = [0.10, 0.15, 0.38, 0.75];
    ax2_pos = [0.58, 0.15, 0.38, 0.75];

    set(ax1, 'Position', ax1_pos);
    set(ax2, 'Position', ax2_pos);

    out_pdf = fullfile(paths.figures, 'Fig1_cdc_aging.pdf');
    out_png = fullfile(paths.figures, 'Fig1_cdc_aging.png');
    out_fig = fullfile(paths.figures, 'Fig1_cdc_aging.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nFigure 1 saved.\n');
end

%% ========================================================================
%  HELPERS
%  ========================================================================

function save_large_figure(fig, out_pdf, out_png, out_fig, w_cm, h_cm)
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [w_cm h_cm]);
    set(fig, 'PaperPosition', [0 0 w_cm h_cm]);

    % PNG — raster
    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf('  Saved: %s (raster, 300 dpi)\n', out_png);

    % PDF — raster-in-PDF to prevent memory crash with 20k scatter points
    exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 300);
    fprintf('  Saved: %s (raster-in-PDF, 300 dpi)\n', out_pdf);

    fprintf('  Skipped: %s (too many graphic objects for .fig format)\n', out_fig);
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
