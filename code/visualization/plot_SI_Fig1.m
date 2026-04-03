function plot_SI_Fig1()
% PLOT_SI_FIG1 - Supplementary Figure 1: CDC predicts healthy HR and
%   identifies hidden mortality risk beyond standard ECG screening
%
%   a) The 1/e ratio predicts the healthy resting heart rate.
%      Healthy controls from 4 databases (N = 1,177): the empirical
%      CDC-HR regression intersects CDC = 1/e at 66 bpm [65.8, 66.7],
%      within the textbook optimal resting HR range.
%
%   b) CDC deviation predicts mortality in clinically normal patients.
%      CODE-15% patients with normal ECG: within 5-bpm HR strata,
%      those whose CDC deviates further from 1/e show higher all-cause
%      mortality. Only bins where every CDC-deviation tertile contains
%      >= 20 mortality events are shown, ensuring stable rate estimates.
%
% Loads pre-computed results from analyze_cdc_vs_hr.m.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    S = load(fullfile(paths.results, 'cdc_vs_hr_results.mat'));
    inv_e = S.inv_e;

    fig = figure('Position', [60 120 1400 620], 'Color', 'w');

    col_near = [0.30 0.65 0.40];
    col_mod  = [0.85 0.65 0.13];
    col_far  = [0.80 0.30 0.25];
    tert_colors = [col_near; col_mod; col_far];
    col_cdc = [0.20 0.45 0.75];

    %% ================================================================
    %  LEFT PANEL (a): 1/e predicts the healthy resting HR
    %  ================================================================

    ax1 = subplot(1, 2, 1);
    hold on; grid on; box on;

    healthy = S.healthy.data;
    hc_bin_stats = S.healthy.bin_stats;
    valid_bins = find([hc_bin_stats.n] >= 10);

    hr_lo = min([hc_bin_stats(valid_bins).center]) - 5;
    hr_hi = max([hc_bin_stats(valid_bins).center]) + 5;

    % Individual scatter (faded)
    scatter(healthy.HR, healthy.CDC, 6, [0.70 0.70 0.70], 'filled', ...
        'MarkerFaceAlpha', 0.15, 'HandleVisibility', 'off');

    % Bin means with 95% CI
    x_c  = [hc_bin_stats(valid_bins).center];
    y_m  = [hc_bin_stats(valid_bins).mean_cdc];
    y_lo = [hc_bin_stats(valid_bins).ci_lo];
    y_hi = [hc_bin_stats(valid_bins).ci_hi];

    errorbar(x_c, y_m, y_m - y_lo, y_hi - y_m, ...
        'o-', 'Color', col_cdc * 0.7, 'MarkerFaceColor', col_cdc, ...
        'MarkerEdgeColor', col_cdc * 0.5, 'MarkerSize', 9, ...
        'LineWidth', 2.2, 'CapSize', 8);

    % Regression line
    xfit = linspace(hr_lo, hr_hi, 200);
    yfit = S.healthy.intercept + S.healthy.slope * xfit;
    plot(xfit, yfit, '-', 'Color', [col_cdc 0.4], 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    % 1/e reference line
    yline(inv_e, 'k--', 'LineWidth', 2.0, 'Label', '1/e', ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 10);

    % Intersection
    hr_opt    = S.healthy.hr_opt;
    hr_opt_ci = S.healthy.hr_opt_ci;

    patch([hr_opt_ci(1) hr_opt_ci(2) hr_opt_ci(2) hr_opt_ci(1)], ...
          [inv_e-0.003 inv_e-0.003 inv_e+0.003 inv_e+0.003], ...
          col_cdc, 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
          'HandleVisibility', 'off');
    plot(hr_opt, inv_e, 'p', 'MarkerSize', 16, ...
        'MarkerFaceColor', [0.85 0.25 0.15], 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.2);
    plot([hr_opt hr_opt], [0.26 inv_e], ':', 'Color', [0.5 0.2 0.2], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Intersection annotation — below the 1/e line
    text(hr_opt + 1.5, 0.31, ...
        sprintf('%.1f bpm\n[%.1f, %.1f]', hr_opt, hr_opt_ci(1), hr_opt_ci(2)), ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.5 0.1 0.1], ...
        'VerticalAlignment', 'top');

    % R^2
    text(0.97, 0.05, sprintf('R^2 = %.3f', S.healthy.r2_cdc_on_hr), ...
        'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);

    % Bin N labels
    for i = 1:length(valid_bins)
        bi = valid_bins(i);
        text(hc_bin_stats(bi).center, y_lo(i) - 0.007, ...
            sprintf('n=%s', format_comma(hc_bin_stats(bi).n)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, ...
            'Color', [0.35 0.35 0.35]);
    end

    xlabel('Heart Rate (bpm)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Cardiac Duty Cycle (RT/RR)', 'FontSize', 13, 'FontWeight', 'bold');
    title({'{\bfa}  1/e predicts the healthy resting heart rate', ...
           sprintf('Healthy controls from 4 databases (N = %s)', ...
                   format_comma(height(healthy)))}, ...
        'FontSize', 12, 'FontWeight', 'normal');
    xlim([hr_lo hr_hi]); ylim([0.26 0.48]);
    set(gca, 'FontSize', 11, 'LineWidth', 0.8);
    hold off;

  
    %% ================================================================
    %  RIGHT PANEL (b): Hidden risk in clinically normal patients
    %  ================================================================

    ax2 = subplot(1, 2, 2);
    hold on; grid on; box on;

    hr_labels   = S.code15_cn.hr_strata.hr_labels;
    hr_edges    = S.code15_cn.hr_strata.hr_edges;
    mort_rate   = S.code15_cn.hr_strata.mort_rate;
    mort_n      = S.code15_cn.hr_strata.mort_n;
    mort_dead   = S.code15_cn.hr_strata.mort_dead;
    tert_labels = S.code15_cn.hr_strata.tert_labels;
    n_cn        = S.code15_cn.n;

    MIN_DEATHS = 20;
    idx_r = find(all(mort_dead >= MIN_DEATHS, 2));
    n_valid = length(idx_r);

    if n_valid < 2
        text(0.5, 0.5, 'Insufficient data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        hold off; return;
    end

    % Grouped bar chart
    bar_data = mort_rate(idx_r, :) * 100;
    b = bar(1:n_valid, bar_data, 'grouped');
    for ti = 1:3
        b(ti).FaceColor = tert_colors(ti, :);
        b(ti).EdgeColor = 'k'; b(ti).LineWidth = 0.8;
    end

    % 95% Wald CI error bars
    for ti = 1:3
        xpos = b(ti).XEndPoints;
        for bi = 1:n_valid
            hi = idx_r(bi);
            p_hat = mort_rate(hi, ti);
            se = sqrt(p_hat * (1-p_hat) / mort_n(hi, ti));
            errorbar(xpos(bi), p_hat*100, 1.96*se*100, 'k', ...
                'LineWidth', 0.8, 'CapSize', 4, 'HandleVisibility', 'off');
        end
    end

    % Y-axis headroom
    max_bar = max(bar_data(:));
    y_max = max_bar * 1.80;
    ylim([0 y_max]);

    % Fold-change annotations
    for bi = 1:n_valid
        hi = idx_r(bi);
        if mort_rate(hi,1) > 0 && ~isnan(mort_rate(hi,3))
            fold = mort_rate(hi,3) / mort_rate(hi,1);
            max_rate = max(mort_rate(hi,:));
            max_ti = find(mort_rate(hi,:) == max_rate, 1);
            max_se = sqrt(max_rate*(1-max_rate)/mort_n(hi,max_ti));
            y_top = max_rate*100 + 1.96*max_se*100;
            text(bi, y_top + 0.12, sprintf('%.1f\\times', fold), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', [0.5 0.1 0.1]);
        end
    end

    % Sample-size info (rotated)
    for bi = 1:n_valid
        hi = idx_r(bi);
        n_total_bin = sum(mort_n(hi, :));
        n_dead_bin  = sum(mort_dead(hi, :));
        text(bi, y_max * 0.88, ...
            sprintf('n=%s (%d)', format_comma(n_total_bin), n_dead_bin), ...
            'HorizontalAlignment', 'left', 'FontSize', 8, ...
            'Color', [0.35 0.35 0.35], 'Rotation', 45);
    end

    % === LEGEND: bottom-left with white background + black rim ===
    leg = legend(tert_labels, 'Location', 'southwest', 'FontSize', 10, 'Box', 'on');
    set(leg, 'Color', 'w', 'EdgeColor', 'k', 'FontSize', 10);   % <-- this fixes the error

    % === CMH annotation: middle-right with boxed background ===
    cmh_p = S.code15_cn.hr_strata.cmh_p;
    if cmh_p < 0.001, cmh_str = 'CMH p < 0.001';
    else, cmh_str = sprintf('CMH p = %.3f', cmh_p); end
    text(0.97, 0.55, cmh_str, 'Units', 'normalized', ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'BackgroundColor', [1 1 1 0.85], 'EdgeColor', 'k', 'Margin', 3);

    % Inclusion criteria (already boxed at bottom-right)
    text(0.97, 0.04, sprintf('Bins with \\geq %d events per tertile', MIN_DEATHS), ...
        'Units', 'normalized', 'FontSize', 8, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
        'Color', [0.3 0.3 0.3], ...
        'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);

    % Clean x-ticks
    set(gca, 'XTick', 1:n_valid, 'XTickLabel', hr_labels(idx_r), ...
        'FontSize', 11, 'LineWidth', 0.8);

    xlabel('Heart Rate (bpm)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('All-Cause Mortality (%)', 'FontSize', 13, 'FontWeight', 'bold');
    title({'{\bfb}  CDC deviation predicts mortality in clinically normal patients', ...
           sprintf('CODE-15%% normal ECG subgroup (N = %s)', format_comma(n_cn))}, ...
        'FontSize', 12, 'FontWeight', 'normal');
    hold off;

    %% ================================================================
    %  SAVE
    %  ================================================================

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [35 16]);
    set(fig, 'PaperPosition', [0 0 35 16]);

    print(fig, fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.pdf'), '-dpdf', '-painters');
    print(fig, fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.png'), '-dpng', '-r300');
    savefig(fig, fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.fig'));

    fprintf('Saved: %s\n', fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.pdf'));
    fprintf('Saved: %s\n', fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.png'));
    fprintf('Saved: %s\n', fullfile(paths.figures, 'SI_Fig1_cdc_vs_hr.fig'));
end

function s = format_comma(n)
    s = num2str(n);
    if n >= 1000
        idx = length(s) - 2;
        while idx > 1
            s = [s(1:idx-1) ',' s(idx:end)]; idx = idx - 3;
        end
    end
end