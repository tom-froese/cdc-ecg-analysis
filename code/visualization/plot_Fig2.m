function plot_Fig2()
% PLOT_FIG2 - Figure 2: Proximity to 1/e predicts survival at every age
%
%   Grouped bar chart showing all-cause mortality (%) by age decade and
%   within-stratum CDC-deviation tertile (Near 1/e, Moderate, Far from
%   1/e). Fold-change annotations above each age group. CMH test
%   p-value annotated.
%
% Data source: age_stratified_mortality_results.mat
%   (via analyze_age_stratified_mortality.m)
%
% Nature Aging formatting: 1.5-column width (120 mm), bold lowercase
%   panel label, 7-pt tick labels, Wald 95% CI error bars.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();

    S = load(fullfile(paths.results, 'age_stratified_mortality_results.mat'));
    age_labels = S.age_labels;
    n_age = length(age_labels);

    mort_rate  = S.tertile.mort_rate;
    mort_n     = S.tertile.mort_n;
    tert_labels = S.tertile.labels;
    n_dev = size(mort_rate, 2);

    %% ================================================================
    %  COLOUR SCHEME
    %  ================================================================
    %  Green / amber / red for Near / Moderate / Far from 1/e.
    %  Slightly differentiated from the HC/CN/Path blues and greens in
    %  Figure 1 to avoid confusion across figures.

    col_near = [0.30 0.65 0.40];   % Green  — nearest to 1/e
    col_mod  = [0.85 0.65 0.13];   % Amber  — moderate deviation
    col_far  = [0.80 0.30 0.25];   % Red    — farthest from 1/e
    colors = [col_near; col_mod; col_far];

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================
    %  1.5-column width (120 mm) works well for 6 age groups × 3 bars.

    fig_w_cm = 12.0;
    fig_h_cm = 8.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [3 5 fig_w_cm fig_h_cm]);

    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    panel_fs = 10;   % not used (single-panel figure) but kept for consistency

    %% ================================================================
    %  GROUPED BAR CHART
    %  ================================================================

    ax = gca;
    hold on; box on;

    bar_data = mort_rate * 100;
    b = bar(1:n_age, bar_data, 'grouped');
    for di = 1:n_dev
        b(di).FaceColor = colors(di, :);
        b(di).EdgeColor = 'k';
        b(di).LineWidth = 0.5;
    end

    % --- 95% Wald CI error bars ---
    for di = 1:n_dev
        xpos = b(di).XEndPoints;
        for ai = 1:n_age
            if mort_n(ai, di) > 0
                p_hat = mort_rate(ai, di);
                se = sqrt(p_hat * (1 - p_hat) / mort_n(ai, di));
                errorbar(xpos(ai), p_hat * 100, 1.96 * se * 100, ...
                    'k', 'LineWidth', 0.6, 'CapSize', 3, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % --- Fold-change annotations (Far / Near) ---
    for ai = 1:n_age
        if mort_rate(ai, 1) > 0 && ~isnan(mort_rate(ai, n_dev))
            fold = mort_rate(ai, n_dev) / mort_rate(ai, 1);
            % Position above the tallest bar + error bar
            p_hat = mort_rate(ai, n_dev);
            se = sqrt(p_hat * (1 - p_hat) / max(mort_n(ai, n_dev), 1));
            y_top = p_hat * 100 + 1.96 * se * 100;
            y_pos = max(max(bar_data(ai, :)), y_top) + 0.4;
            text(ai, y_pos, sprintf('%.1f\\times', fold), ...
                'HorizontalAlignment', 'center', 'FontSize', ax_fs - 0.5, ...
                'FontWeight', 'bold', 'Color', [0.4 0.1 0.1]);
        end
    end

    % --- CMH p-value annotation ---
    if isfield(S, 'tertile') && isfield(S.tertile, 'cmh_p')
        cmh_p = S.tertile.cmh_p;
    else
        cmh_p = NaN;  % will be computed if re-running updated analysis
    end

    if ~isnan(cmh_p)
        if cmh_p < 0.001
            cmh_str = 'CMH \itp\rm < 0.001';
        else
            cmh_str = sprintf('CMH \\itp\\rm = %.3f', cmh_p);
        end
        text(0.97, 0.95, cmh_str, 'Units', 'normalized', ...
            'FontSize', ax_fs, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.5 0.5 0.5]);
    end

    % --- Sample size annotations (rotated, above bar groups) ---
    for ai = 1:n_age
        n_total_bin = sum(mort_n(ai, :));
        if n_total_bin > 0
            % Position at bottom of plot
            text(ai, -0.5, sprintf('N=%s', format_comma(n_total_bin)), ...
                'HorizontalAlignment', 'center', 'FontSize', ax_fs - 1, ...
                'Color', [0.4 0.4 0.4]);
        end
    end

    hold off;

    %% ================================================================
    %  AXES FORMATTING
    %  ================================================================

    set(ax, 'XTick', 1:n_age, 'XTickLabel', age_labels, ...
        'FontSize', ax_fs, 'LineWidth', 0.5, ...
        'TickDir', 'out', 'TickLength', [0.02 0.02]);
    xlabel('Age group (years)', 'FontSize', lab_fs);
    ylabel('All-cause mortality (%)', 'FontSize', lab_fs);

    % Adjust y-axis to accommodate fold-change labels and N labels
    yl = ylim;
    ylim([min(-1.0, yl(1)), yl(2) * 1.12]);

    leg = legend(b, tert_labels, ...
        'Location', 'northwest', 'FontSize', ax_fs - 0.5, 'Box', 'off');
    leg.ItemTokenSize = [10 8];

    %% ================================================================
    %  SAVE
    %  ================================================================

    set(ax, 'Position', [0.14  0.17  0.82  0.76]);

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    out_pdf = fullfile(paths.figures, 'Fig2_mortality_tertiles.pdf');
    out_png = fullfile(paths.figures, 'Fig2_mortality_tertiles.png');
    out_fig = fullfile(paths.figures, 'Fig2_mortality_tertiles.fig');

    print(fig, out_pdf, '-dpdf', '-painters');
    print(fig, out_png, '-dpng', '-r300');
    savefig(fig, out_fig);

    fprintf('\nFigure 2 saved:\n');
    fprintf('  %s  (vector, for submission)\n', out_pdf);
    fprintf('  %s  (raster, 300 dpi)\n', out_png);
    fprintf('  %s  (editable)\n', out_fig);

    %% ================================================================
    %  SUMMARY STATISTICS (for legend / Methods)
    %  ================================================================

    fprintf('\n--- Summary for figure legend ---\n');
    fprintf('Total N = %s\n', format_comma(S.n_total));
    fprintf('Deaths  = %s (%.1f%%)\n', format_comma(S.n_deceased), ...
        100 * S.n_deceased / S.n_total);
    fprintf('Age groups: %s\n', strjoin(age_labels, ', '));

    fprintf('\nMortality by tertile (all ages pooled):\n');
    for di = 1:n_dev
        n_all = sum(mort_n(:, di));
        d_all = sum(S.tertile.mort_d(:, di));
        fprintf('  %-20s  %.2f%%  (%d / %s)\n', tert_labels{di}, ...
            100 * d_all / n_all, d_all, format_comma(n_all));
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
