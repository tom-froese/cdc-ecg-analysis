function plot_SI_Fig2()
% PLOT_SI_FIG2 - Supplementary Figure 2: Systole and RR interval across
%   the lifespan
%
%   Panel a: Mechanical systole duration (ms) vs Age by clinical group.
%     Shows the gradual age-dependent prolongation of systole that is
%     common to all groups but accelerated in pathological hearts.
%     Together with Fig. 1b (diastole), this reveals the compensatory
%     mechanism: healthy hearts maintain CDC near 1/e by adjusting
%     cycle length in proportion to systolic lengthening.
%
%   Panel b: RR interval (ms) vs Age by clinical group.
%     Confirms that healthy controls slow their heart rate with age,
%     providing the cycle-length headroom for diastolic compensation.
%
% Data source: hierarchical_results.mat (via analyze_hierarchical_model.m)
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();

    %% ================================================================
    %  LOAD DATA
    %  ================================================================

    S = load(fullfile(paths.results, 'hierarchical_results.mat'), 'all_data');
    D = S.all_data;

    % Derived quantities
    D.RR_ms  = 60000 ./ D.HR;              % RR interval (ms)
    D.RT_ms  = D.CDC .* D.RR_ms;           % mechanical systole (ms)

    % Group masks
    is_hc   = D.Group == 'HealthyControl';
    is_cn   = D.Group == 'ClinicallyNormal';
    is_path = D.Group == 'Pathological';

    n_hc   = sum(is_hc);
    n_cn   = sum(is_cn);
    n_path = sum(is_path);

    fprintf('SI Figure 2 — data loaded:\n');
    fprintf('  Healthy Control:   N = %s\n', format_comma(n_hc));
    fprintf('  Clinically Normal: N = %s\n', format_comma(n_cn));
    fprintf('  Pathological:      N = %s\n', format_comma(n_path));

    %% ================================================================
    %  COLOUR SCHEME
    %  ================================================================

    col_hc   = [0.20 0.55 0.85];
    col_cn   = [0.25 0.70 0.35];
    col_path = [0.85 0.25 0.20];

    colors = [col_hc; col_cn; col_path];
    masks  = {is_hc, is_cn, is_path};
    ns     = [n_hc, n_cn, n_path];

    labels = { ...
        sprintf('Healthy Control (\\itn\\rm = %s)', format_comma(n_hc)), ...
        sprintf('Clinically Normal (\\itn\\rm = %s)', format_comma(n_cn)), ...
        sprintf('Pathological (\\itn\\rm = %s)', format_comma(n_path))};

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================

    fig_w_cm = 18.3;
    fig_h_cm = 8.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 5 fig_w_cm fig_h_cm]);

    ax_fs    = 7;
    lab_fs   = 8;
    panel_fs = 10;

    alpha_dot = 0.08;
    dot_size  = 3;

    age_lim = [15 95];
    xfit = linspace(age_lim(1), age_lim(2), 300)';

    %% ================================================================
    %  PANEL (a): Mechanical systole (ms) vs Age
    %  ================================================================

    ax1 = subplot(1, 2, 1);
    hold on; box on;

    % Scatter
    scatter(D.Age(is_path), D.RT_ms(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.RT_ms(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.RT_ms(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, ...
        'HandleVisibility', 'off');

    % Linear regression with 95% CI
    h_lines = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.RT_ms(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h_lines(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, ...
            'LineWidth', 2.0);

        beta = mdl.Coefficients.Estimate(2);
        p_val = mdl.Coefficients.pValue(2);
        fprintf('  Panel a — %s: slope = %+.2f ms/yr, p = %.2e\n', ...
            char(regexp(labels{g}, '^[^(]+', 'match')), beta, p_val);
    end

    xlim(age_lim);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('Mechanical systole (ms)', 'FontSize', lab_fs);
    set(ax1, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    leg1 = legend(h_lines, labels, ...
        'Location', 'northwest', 'FontSize', ax_fs - 0.5, 'Box', 'off');
    leg1.ItemTokenSize = [12 8];

    text(-0.14, 1.08, '\bfa', 'Units', 'normalized', ...
        'FontSize', panel_fs, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');

    hold off;

    %% ================================================================
    %  PANEL (b): RR interval (ms) vs Age
    %  ================================================================

    ax2 = subplot(1, 2, 2);
    hold on; box on;

    % Scatter
    scatter(D.Age(is_path), D.RR_ms(is_path), dot_size, ...
        col_path, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_cn), D.RR_ms(is_cn), dot_size, ...
        col_cn, 'filled', 'MarkerFaceAlpha', alpha_dot, ...
        'HandleVisibility', 'off');
    scatter(D.Age(is_hc), D.RR_ms(is_hc), dot_size, ...
        col_hc, 'filled', 'MarkerFaceAlpha', alpha_dot * 2.5, ...
        'HandleVisibility', 'off');

    % Linear regression with 95% CI
    h_lines2 = gobjects(3, 1);
    for g = 1:3
        mask = masks{g};
        mdl = fitlm(D.Age(mask), D.RR_ms(mask));
        [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

        fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
            colors(g,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        h_lines2(g) = plot(xfit, yfit, '-', 'Color', colors(g,:) * 0.8, ...
            'LineWidth', 2.0);

        beta = mdl.Coefficients.Estimate(2);
        p_val = mdl.Coefficients.pValue(2);
        fprintf('  Panel b — %s: slope = %+.2f ms/yr, p = %.2e\n', ...
            char(regexp(labels{g}, '^[^(]+', 'match')), beta, p_val);
    end

    xlim(age_lim);
    xlabel('Age (years)', 'FontSize', lab_fs);
    ylabel('RR interval (ms)', 'FontSize', lab_fs);
    set(ax2, 'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);

    leg2 = legend(h_lines2, labels, ...
        'Location', 'southeast', 'FontSize', ax_fs - 0.5, 'Box', 'off');
    leg2.ItemTokenSize = [12 8];

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

    out_pdf = fullfile(paths.figures, 'SI_Fig2_systole_rr.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig2_systole_rr.png');
    out_fig = fullfile(paths.figures, 'SI_Fig2_systole_rr.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nSI Figure 2 saved:\n');
    fprintf('  %s  (vector)\n', out_pdf);
    fprintf('  %s  (300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  SUMMARY
    %  ================================================================

    fprintf('\n--- Systole summary (ms) ---\n');
    for g = 1:3
        m = masks{g};
        fprintf('  %-20s  %.0f +/- %.0f ms\n', ...
            char(regexp(labels{g}, '^[^(]+', 'match')), ...
            mean(D.RT_ms(m)), std(D.RT_ms(m)));
    end
    fprintf('\n--- RR interval summary (ms) ---\n');
    for g = 1:3
        m = masks{g};
        fprintf('  %-20s  %.0f +/- %.0f ms\n', ...
            char(regexp(labels{g}, '^[^(]+', 'match')), ...
            mean(D.RR_ms(m)), std(D.RR_ms(m)));
    end
end


%% ========================================================================
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
