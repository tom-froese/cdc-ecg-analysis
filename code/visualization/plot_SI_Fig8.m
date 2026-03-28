function plot_SI_Fig8()
% PLOT_SI_FIG8 - Supplementary Figure 8: Sex-stratified CDC aging
%   trajectories confirm 1/e convergence in both sexes
%
% Four-panel figure:
%   a: ΔCDC vs Age — Female, all three clinical groups
%   b: ΔCDC vs Age — Male, all three clinical groups
%   c: Thermodynamic diastole (ms) vs Age — Female, all groups
%   d: Thermodynamic diastole (ms) vs Age — Male, all groups
%
% Design rationale:
%   The panels mirror the structure of Figure 1 (which shows ΔCDC and
%   diastole vs Age pooled across sexes), split into female and male
%   subpanels to demonstrate that the 1/e convergence pattern — and the
%   diastolic compensation mechanism — hold in both sexes. Within each
%   panel, all three clinical groups are overlaid, preserving the
%   group-comparison structure of Fig. 1. This layout makes it visually
%   immediate that the aging trajectories are consistent across sexes,
%   which is the key claim for the reviewer.
%
% Colour scheme matches Figure 1 (plot_Fig1.m):
%   Blue    [0.20 0.55 0.85]  = Healthy Control
%   Green   [0.25 0.70 0.35]  = Clinically Normal
%   Red     [0.85 0.25 0.20]  = Pathological
%
% Data source: sex_stratified_results.mat (via analyze_sex_stratified.m)
%
% Nature Aging formatting: double-column width (183 mm), bold lowercase
%   panel labels outside axes, 7-pt tick labels.
%
% Tom Froese, OIST Embodied Cognitive Science Unit

    paths = config();
    inv_e = 1 / exp(1);

    %% ================================================================
    %  LOAD PRECOMPUTED RESULTS
    %  ================================================================

    S = load(fullfile(paths.results, 'sex_stratified_results.mat'));

    D            = S.all_data;
    groups       = S.groups;
    group_labels = S.group_labels;
    sexes        = S.sexes;
    sex_labels   = S.sex_labels;
    ols_dcdc     = S.ols_dcdc;
    ols_dias     = S.ols_dias;
    interaction  = S.interaction;
    desc         = S.desc;

    %% ================================================================
    %  COLOUR SCHEME (matches Figure 1)
    %  ================================================================

    col_hc   = [0.20 0.55 0.85];   % Blue  — Healthy Control
    col_cn   = [0.25 0.70 0.35];   % Green — Clinically Normal
    col_path = [0.85 0.25 0.20];   % Red   — Pathological
    colors = [col_hc; col_cn; col_path];

    %% ================================================================
    %  FIGURE SETUP — Nature Aging formatting
    %  ================================================================

    fig_w_cm = 18.3;   % double-column width
    fig_h_cm = 16.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 1 fig_w_cm fig_h_cm]);

    ax_fs    = 7;    % tick labels
    lab_fs   = 8;    % axis labels
    title_fs = 8;    % panel titles
    panel_fs = 10;   % panel letters
    leg_fs   = 6;    % legend

    % Scatter aesthetics
    dot_size  = 2;
    alpha_dot = 0.06;

    % Shared axis limits
    age_lim   = [15 95];
    dcdc_lim  = [-0.08 0.12];

    panel_labels = {'a', 'b', 'c', 'd'};
    ax = gobjects(4, 1);

    %% ================================================================
    %  PANELS (a,b): ΔCDC vs Age — Female, then Male
    %  ================================================================

    for si = 1:2
        ax(si) = subplot(2, 2, si);
        hold on; box on;

        sex_mask = (D.Sex == sexes{si});

        % Reference line at ΔCDC = 0 (optimal)
        yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

        % Scatter and OLS per group (back-to-front: Path, CN, HC)
        h_lines = gobjects(3, 1);
        xfit = linspace(age_lim(1), age_lim(2), 300)';

        for gi = [3, 2, 1]
            mask = sex_mask & (D.Group == groups{gi});

            % Scatter
            alpha_g = alpha_dot;
            if gi == 1, alpha_g = alpha_dot * 2.5; end  % HC more visible
            scatter(D.Age(mask), D.delta_CDC(mask), dot_size, ...
                colors(gi, :), 'filled', 'MarkerFaceAlpha', alpha_g, ...
                'HandleVisibility', 'off');

            % OLS trend + CI
            slope = ols_dcdc(gi, si).slope;
            intercept = ols_dcdc(gi, si).intercept;
            yfit = intercept + slope * xfit;

            % Refit for CI band (need the model object)
            mdl = fitlm(D.Age(mask), D.delta_CDC(mask));
            [~, ci] = predict(mdl, xfit, 'Alpha', 0.05);

            fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
                colors(gi, :), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');

            h_lines(gi) = plot(xfit, yfit, '-', ...
                'Color', colors(gi, :), 'LineWidth', 1.8);
        end

        % Optimal label
        if si == 1
            text(age_lim(2) - 1, 0.004, 'Optimal (1/\ite\rm)', ...
                'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end

        xlim(age_lim);
        ylim(dcdc_lim);
        xlabel('Age (years)', 'FontSize', lab_fs);
        ylabel('\DeltaCDC from 1/\ite', 'FontSize', lab_fs);
        title(sex_labels{si}, 'FontSize', title_fs, 'FontWeight', 'bold');

        set(ax(si), 'FontSize', ax_fs, 'LineWidth', 0.5, ...
            'TickDir', 'out', 'TickLength', [0.02 0.02]);

        % Legend (panel a only)
        if si == 1
            leg_strs = cell(3, 1);
            for gi = 1:3
                leg_strs{gi} = sprintf('%s (n=%s)', ...
                    group_labels{gi}, format_comma(desc(gi, si).n));
            end
            leg = legend(h_lines, leg_strs, ...
                'Location', 'northwest', 'FontSize', leg_fs, 'Box', 'off');
            leg.ItemTokenSize = [12 8];
        end

        % Slope annotations (top-right, compact)
        for gi = 1:3
            slope = ols_dcdc(gi, si).slope;
            p_val = ols_dcdc(gi, si).p_slope;

            if p_val < 0.001
                p_str = '<0.001';
            else
                p_str = sprintf('%.3f', p_val);
            end

            y_pos = 0.97 - (gi - 1) * 0.09;
            text(0.97, y_pos, ...
                sprintf('%+.5f/yr (\\itp\\rm=%s)', slope, p_str), ...
                'Units', 'normalized', 'FontSize', ax_fs - 0.5, ...
                'Color', colors(gi, :) * 0.7, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end

        % Panel label
        text(-0.14, 1.08, ['\bf' panel_labels{si}], 'Units', 'normalized', ...
            'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

        hold off;
    end

    %% ================================================================
    %  PANELS (c,d): DIASTOLE vs Age — Female, then Male
    %  ================================================================

    for si = 1:2
        ax(si + 2) = subplot(2, 2, si + 2);
        hold on; box on;

        sex_mask = (D.Sex == sexes{si});

        h_lines2 = gobjects(3, 1);
        xfit = linspace(age_lim(1), age_lim(2), 300)';

        for gi = [3, 2, 1]
            mask = sex_mask & (D.Group == groups{gi});

            alpha_g = alpha_dot;
            if gi == 1, alpha_g = alpha_dot * 2.5; end
            scatter(D.Age(mask), D.Dias_ms(mask), dot_size, ...
                colors(gi, :), 'filled', 'MarkerFaceAlpha', alpha_g, ...
                'HandleVisibility', 'off');

            mdl = fitlm(D.Age(mask), D.Dias_ms(mask));
            [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

            fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
                colors(gi, :), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');

            h_lines2(gi) = plot(xfit, yfit, '-', ...
                'Color', colors(gi, :), 'LineWidth', 1.8);
        end

        xlim(age_lim);
        xlabel('Age (years)', 'FontSize', lab_fs);
        ylabel('Thermodynamic diastole (ms)', 'FontSize', lab_fs);
        title(sex_labels{si}, 'FontSize', title_fs, 'FontWeight', 'bold');

        set(ax(si + 2), 'FontSize', ax_fs, 'LineWidth', 0.5, ...
            'TickDir', 'out', 'TickLength', [0.02 0.02]);

        % Legend (panel c only)
        if si == 1
            leg_strs2 = cell(3, 1);
            for gi = 1:3
                leg_strs2{gi} = sprintf('%s (n=%s)', ...
                    group_labels{gi}, format_comma(desc(gi, si).n));
            end
            leg2 = legend(h_lines2, leg_strs2, ...
                'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'off');
            leg2.ItemTokenSize = [12 8];
        end

        % Slope annotations
        for gi = 1:3
            slope = ols_dias(gi, si).slope;
            p_val = ols_dias(gi, si).p_slope;

            if p_val < 0.001
                p_str = '<0.001';
            else
                p_str = sprintf('%.3f', p_val);
            end

            y_pos = 0.97 - (gi - 1) * 0.09;
            text(0.97, y_pos, ...
                sprintf('%+.2f ms/yr (\\itp\\rm=%s)', slope, p_str), ...
                'Units', 'normalized', 'FontSize', ax_fs - 0.5, ...
                'Color', colors(gi, :) * 0.7, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end

        % Panel label
        text(-0.14, 1.08, ['\bf' panel_labels{si + 2}], 'Units', 'normalized', ...
            'FontSize', panel_fs, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

        hold off;
    end

    %% ================================================================
    %  LAYOUT AND SAVE
    %  ================================================================

    set(ax(1), 'Position', [0.08  0.58  0.40  0.36]);
    set(ax(2), 'Position', [0.56  0.58  0.40  0.36]);
    set(ax(3), 'Position', [0.08  0.08  0.40  0.36]);
    set(ax(4), 'Position', [0.56  0.08  0.40  0.36]);

    out_pdf = fullfile(paths.figures, 'SI_Fig8_sex_stratified.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig8_sex_stratified.png');
    out_fig = fullfile(paths.figures, 'SI_Fig8_sex_stratified.fig');

    save_large_figure(fig, out_pdf, out_png, out_fig, fig_w_cm, fig_h_cm);

    fprintf('\nSupplementary Figure 8 saved.\n');

    %% ================================================================
    %  CONSOLE SUMMARY (for SI figure legend)
    %  ================================================================

    fprintf('\n--- Summary for SI Fig 8 legend ---\n');

    n_total = height(D);
    n_f = sum(D.Sex == 'F');
    n_m = sum(D.Sex == 'M');
    fprintf('N = %s subjects (%s F, %s M)\n', ...
        format_comma(n_total), format_comma(n_f), format_comma(n_m));

    fprintf('\nGroup x Sex breakdown:\n');
    for gi = 1:3
        fprintf('  %s:  F=%s, M=%s\n', group_labels{gi}, ...
            format_comma(desc(gi, 1).n), format_comma(desc(gi, 2).n));
    end

    fprintf('\ndCDC slopes (per year):\n');
    for gi = 1:3
        for si = 1:2
            fprintf('  %-22s %-8s  %+.6f/yr (p=%.2e)\n', ...
                group_labels{gi}, sex_labels{si}, ...
                ols_dcdc(gi, si).slope, ols_dcdc(gi, si).p_slope);
        end
    end

    fprintf('\nDiastole slopes (ms per year):\n');
    for gi = 1:3
        for si = 1:2
            fprintf('  %-22s %-8s  %+.4f ms/yr (p=%.2e)\n', ...
                group_labels{gi}, sex_labels{si}, ...
                ols_dias(gi, si).slope, ols_dias(gi, si).p_slope);
        end
    end

    fprintf('\nSex x Age interaction (within group):\n');
    for gi = 1:3
        fprintf('  %-22s  beta=%+.6f, p=%.2e\n', ...
            interaction(gi).group, interaction(gi).beta, interaction(gi).p);
    end

    fprintf('\nCohen''s d (M-F):\n');
    for gi = 1:3
        fprintf('  %-22s  CDC: d=%+.3f,  HR: d=%+.3f\n', ...
            S.cohens_d(gi).group, S.cohens_d(gi).d_cdc, S.cohens_d(gi).d_hr);
    end
end


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function save_large_figure(fig, out_pdf, out_png, out_fig, w_cm, h_cm)
% SAVE_LARGE_FIGURE - Save figure without crashing on large scatter data
%
% For figures with many graphic objects, MATLAB's painters renderer and
% savefig serialize every element, consuming gigabytes of memory.
%
% Strategy:
%   PNG: exportgraphics (raster, always safe)
%   PDF: exportgraphics with ContentType 'image' (raster-in-PDF wrapper,
%        avoids the painters memory explosion while producing a PDF that
%        embeds at 300 dpi)
%   FIG: skipped (the .fig format stores all graphic objects and can
%        itself become multi-GB)

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

    close(fig);
end


function s = format_comma(n)
% FORMAT_COMMA - Format integer with thousands separators
    s = num2str(n);
    if n >= 1000
        idx = length(s) - 2;
        while idx > 1
            s = [s(1:idx-1) ',' s(idx:end)];
            idx = idx - 3;
        end
    end
end
