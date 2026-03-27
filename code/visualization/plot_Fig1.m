function plot_Fig1()
% PLOT_FIG1 - Figure 1: The cardiac duty cycle converges on 1/e in
%   healthy hearts and deviates with pathological aging
%
%   Panel a: ΔCDC from the theoretical optimum (1/e) across the adult
%     lifespan, by clinical group. Individual subjects shown as
%     semi-transparent dots; linear regression with 95% CI shading
%     per group. The near-zero slope in healthy controls vs the
%     positive slopes in clinical groups is the core finding.
%
%   Panel b: Thermodynamic diastole duration (ms) across the adult
%     lifespan, by clinical group. Derived from subject-level median
%     CDC and HR as: diastole = RR × (1 − CDC), where RR = 60000/HR.
%
% Linear regression was chosen over LOESS because the healthy control
% group (N = 1,165) is drawn from databases with coarse or bimodal age
% distributions (Fantasia: bimodal at ~26/~74 yr; Autonomic Aging:
% age-range midpoints), causing LOESS to overfit to cluster structure.
% OLS regression with 95% CI bands is robust to uneven age density and
% transparently communicates uncertainty (wider bands for smaller N).
%
% Data source: hierarchical_results.mat (via analyze_hierarchical_model.m)
%
% Nature Aging formatting: double-column (183 mm), bold lowercase
%   panel labels outside axes, 7-pt tick labels, no sgtitle.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD DATA
    %  ================================================================

    S = load(fullfile(paths.results, 'hierarchical_results.mat'), 'all_data');
    D = S.all_data;

    % Derived quantities
    D.delta_CDC = D.CDC - inv_e;              % deviation from optimum
    D.RR_ms     = 60000 ./ D.HR;              % RR interval (ms)
    D.Dias_ms   = D.RR_ms .* (1 - D.CDC);    % thermodynamic diastole (ms)

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
    fprintf('  Total:             N = %s\n', format_comma(height(D)));

    %% ================================================================
    %  COLOUR SCHEME (consistent across all main figures)
    %  ================================================================

    col_hc   = [0.20 0.55 0.85];   % Blue  — verified healthy volunteers
    col_cn   = [0.25 0.70 0.35];   % Green — clinically normal ECG
    col_path = [0.85 0.25 0.20];   % Red   — pathological

    colors = [col_hc; col_cn; col_path];
    masks  = {is_hc, is_cn, is_path};
    ns     = [n_hc, n_cn, n_path];

    % Legend labels with sample sizes (Nature verbal cues)
    labels = { ...
        sprintf('Healthy Control (\\itn\\rm = %s)', format_comma(n_hc)), ...
        sprintf('Clinically Normal (\\itn\\rm = %s)', format_comma(n_cn)), ...
        sprintf('Pathological (\\itn\\rm = %s)', format_comma(n_path))};

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================
    %  Nature Aging double-column width: 183 mm.
    %  Two panels, ~80 mm height.

    fig_w_cm = 18.3;
    fig_h_cm = 8.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 5 fig_w_cm fig_h_cm]);

    % Font sizes for final print dimensions (Nature minimum: 5 pt)
    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    panel_fs = 10;   % panel labels (a, b)

    % Scatter aesthetics
    alpha_dot = 0.08;
    dot_size  = 3;

    % Shared axis limits
    age_lim = [15 95];
    xfit = linspace(age_lim(1), age_lim(2), 300)';

    %% ================================================================
    %  PANEL (a): ΔCDC from optimal (1/e) vs Age
    %  ================================================================

    ax1 = subplot(1, 2, 1);
    hold on; box on;

    % Individual scatter — back-to-front so HC visible on top
    scatter(D.Age(is_path), D.delta_CDC(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.delta_CDC(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.delta_CDC(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, ...
        'HandleVisibility', 'off');

    % Reference line at ΔCDC = 0 (CDC = 1/e)
    yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    text(age_lim(2) - 1, -0.008, 'Optimal (1/\ite\rm)', ...
        'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

    % Linear regression with 95% CI per group
    h_lines = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.delta_CDC(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        h_lines(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, ...
            'LineWidth', 2.0);
    end

    xlim(age_lim);
    ylim([-0.08 0.12]);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    % Legend moved to southwest with black-rimmed white box
    leg1 = legend(h_lines, labels, 'Location', 'southwest', ...
        'FontSize', ax_fs - 0.5, 'Box', 'on');
    set(leg1, 'Color', 'w', 'EdgeColor', 'k', 'FontSize', ax_fs - 0.5);

    text(-0.14, 1.08, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): Thermodynamic diastole (ms) vs Age
    %  ================================================================

    ax2 = subplot(1, 2, 2);
    hold on; box on;

    % Individual scatter
    scatter(D.Age(is_path), D.Dias_ms(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.Dias_ms(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.Dias_ms(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, ...
        'HandleVisibility', 'off');

    % Linear regression with 95% CI per group
    h_lines2 = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.Dias_ms(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        h_lines2(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, ...
            'LineWidth', 2.0);

        beta = mdl.Coefficients.Estimate(2);
        p_val = mdl.Coefficients.pValue(2);
        fprintf('  Panel b — %s: slope = %+.2f ms/yr, p = %.2e\n', ...
            labels{g}, beta, p_val);
    end

    xlim(age_lim);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('Thermodynamic diastole (ms)', 'FontSize', lab_fs);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    leg2 = legend(h_lines2, labels, ...
        'Location', 'northeast', 'FontSize', ax_fs - 0.5, 'Box', 'off');
    leg2.ItemTokenSize = [12 8];

    % Panel label
    text(-0.14, 1.08, '\bfb', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================

    set(ax1, 'Position', [0.08  0.17  0.40  0.73]);
    set(ax2, 'Position', [0.56  0.17  0.40  0.73]);

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'Fig1_cdc_aging.pdf');
    out_png = fullfile(paths.figures, 'Fig1_cdc_aging.png');
    out_fig = fullfile(paths.figures, 'Fig1_cdc_aging.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nFigure 1 saved:\n');
    fprintf('  %s  (vector, for submission)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  SUMMARY STATISTICS (for legend / Methods)
    %  ================================================================

    fprintf('\n--- Summary for figure legend ---\n');
    fprintf('Total subjects: N = %s\n', format_comma(height(D)));
    fprintf('  HC:   N = %5s, age %d-%d yr\n', format_comma(n_hc), ...
        round(min(D.Age(is_hc))), round(max(D.Age(is_hc))));
    fprintf('  CN:   N = %5s, age %d-%d yr\n', format_comma(n_cn), ...
        round(min(D.Age(is_cn))), round(max(D.Age(is_cn))));
    fprintf('  Path: N = %5s, age %d-%d yr\n', format_comma(n_path), ...
        round(min(D.Age(is_path))), round(max(D.Age(is_path))));

    fprintf('\nDiastole summary (ms):\n');
    for g = 1:3
        m = masks{g};
        fprintf('  %-20s  %.0f +/- %.0f ms\n', ...
            char(regexp(labels{g}, '^[^(]+', 'match')), ...
            mean(D.Dias_ms(m)), std(D.Dias_ms(m)));
    end
end


%% ========================================================================
%  HELPERS
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
