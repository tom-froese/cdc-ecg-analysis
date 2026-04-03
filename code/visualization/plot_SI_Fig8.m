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
% Colour scheme matches Figure 1:
%   Blue    [0.20 0.55 0.85]  = Healthy Control
%   Green   [0.25 0.70 0.35]  = Clinically Normal
%   Red     [0.85 0.25 0.20]  = Pathological
%
% Uses manual axes positioning and data-coordinate annotations to
% avoid exportgraphics displacement bugs.
%
% Data source: sex_stratified_results.mat (via analyze_sex_stratified.m)
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
    desc         = S.desc;

    %% ================================================================
    %  COLOUR SCHEME (matches Figure 1)
    %  ================================================================

    col_hc   = [0.20 0.55 0.85];
    col_cn   = [0.25 0.70 0.35];
    col_path = [0.85 0.25 0.20];
    colors = [col_hc; col_cn; col_path];

    %% ================================================================
    %  FIGURE SETUP
    %  ================================================================

    fig_w_cm = 18.3;
    fig_h_cm = 18.0;

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 1 fig_w_cm fig_h_cm]);

    % Font sizes — enlarged for readability
    ax_fs    = 9;
    lab_fs   = 10;
    title_fs = 11;
    panel_fs = 13;
    leg_fs   = 7.5;
    slope_fs = 7.5;

    % Scatter aesthetics
    dot_size  = 2;
    alpha_dot = 0.06;
    line_w    = 2.0;

    % Shared axis limits
    age_lim   = [15 95];
    dcdc_lim  = [-0.08 0.12];

    % Manual panel positions: [left, bottom, width, height]
    pw = 0.38;
    ph = 0.35;
    positions = {
        [0.10, 0.58, pw, ph],   % a: top-left (Female ΔCDC)
        [0.57, 0.58, pw, ph],   % b: top-right (Male ΔCDC)
        [0.10, 0.08, pw, ph],   % c: bottom-left (Female diastole)
        [0.57, 0.08, pw, ph],   % d: bottom-right (Male diastole)
    };

    panel_labels = {'a', 'b', 'c', 'd'};
    ax = gobjects(4, 1);

    %% ================================================================
    %  PANELS (a,b): ΔCDC vs Age — Female, then Male
    %  ================================================================

    for si = 1:2
        ax(si) = axes('Position', positions{si});
        hold on; box on;

        sex_mask = (D.Sex == sexes{si});

        % Reference line at ΔCDC = 0
        yline(0, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

        % Scatter and OLS per group (back-to-front: Path, CN, HC)
        h_lines = gobjects(3, 1);
        xfit = linspace(age_lim(1), age_lim(2), 300)';

        for gi = [3, 2, 1]
            mask = sex_mask & (D.Group == groups{gi});

            alpha_g = alpha_dot;
            if gi == 1, alpha_g = alpha_dot * 2.5; end
            scatter(D.Age(mask), D.delta_CDC(mask), dot_size, ...
                colors(gi, :), 'filled', 'MarkerFaceAlpha', alpha_g, ...
                'HandleVisibility', 'off');

            mdl = fitlm(D.Age(mask), D.delta_CDC(mask));
            [yfit, ci] = predict(mdl, xfit, 'Alpha', 0.05);

            fill([xfit; flipud(xfit)], [ci(:,1); flipud(ci(:,2))], ...
                colors(gi, :), 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');

            h_lines(gi) = plot(xfit, yfit, '-', ...
                'Color', colors(gi, :), 'LineWidth', line_w);
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
            leg1 = legend(h_lines, leg_strs, ...
                'Location', 'southwest', 'FontSize', leg_fs, 'Box', 'on');
            set(leg1, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
            leg1.ItemTokenSize = [12 8];
        end

        % Slope annotations (data coordinates, bottom-right)
        axes(ax(si));
        xl = xlim; yl = ylim;
        x_ann = xl(2) - 0.03 * diff(xl);

        for gi = 1:3
            slope = ols_dcdc(gi, si).slope;
            p_val = ols_dcdc(gi, si).p_slope;

            if p_val < 0.001
                p_str = '<0.001';
            else
                p_str = sprintf('%.3f', p_val);
            end

            y_ann = yl(1) + (0.18 - (gi-1)*0.06) * diff(yl);
            text(x_ann, y_ann, ...
                sprintf('%+.5f/yr (\\itp\\rm=%s)', slope, p_str), ...
                'FontSize', slope_fs, ...
                'Color', colors(gi, :) * 0.7, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
        end

        % Optimal label (panel a only)
        if si == 1
            text(xl(2) - 0.02*diff(xl), 0.004, 'Optimal (1/\ite\rm)', ...
                'FontSize', ax_fs, 'Color', [0.3 0.3 0.3], ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end

        % Panel label (data coordinates)
        text(xl(1) - 0.14*diff(xl), yl(2) + 0.04*diff(yl), ...
            ['\bf' panel_labels{si}], ...
            'FontSize', panel_fs, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'bottom');

        hold off;
    end

    %% ================================================================
    %  PANELS (c,d): DIASTOLE vs Age — Female, then Male
    %  ================================================================

    for si = 1:2
        ax(si+2) = axes('Position', positions{si+2});
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
                'Color', colors(gi, :), 'LineWidth', line_w);
        end

        xlim(age_lim);
        xlabel('Age (years)', 'FontSize', lab_fs);
        ylabel('Thermodynamic diastole (ms)', 'FontSize', lab_fs);
        title(sex_labels{si}, 'FontSize', title_fs, 'FontWeight', 'bold');

        set(ax(si+2), 'FontSize', ax_fs, 'LineWidth', 0.5, ...
            'TickDir', 'out', 'TickLength', [0.02 0.02]);

        % Legend (panel c only)
        if si == 1
            leg_strs2 = cell(3, 1);
            for gi = 1:3
                leg_strs2{gi} = sprintf('%s (n=%s)', ...
                    group_labels{gi}, format_comma(desc(gi, si).n));
            end
            leg2 = legend(h_lines2, leg_strs2, ...
                'Location', 'northeast', 'FontSize', leg_fs, 'Box', 'on');
            set(leg2, 'Color', 'w', 'EdgeColor', [0.7 0.7 0.7]);
            leg2.ItemTokenSize = [12 8];
        end

        % Slope annotations (data coordinates, bottom-left)
        axes(ax(si+2));
        xl = xlim; yl = ylim;
        x_ann = xl(1) + 0.03 * diff(xl);

        for gi = 1:3
            slope = ols_dias(gi, si).slope;
            p_val = ols_dias(gi, si).p_slope;

            if p_val < 0.001
                p_str = '<0.001';
            else
                p_str = sprintf('%.3f', p_val);
            end

            y_ann = yl(1) + (0.18 - (gi-1)*0.06) * diff(yl);
            text(x_ann, y_ann, ...
                sprintf('%+.2f ms/yr (\\itp\\rm=%s)', slope, p_str), ...
                'FontSize', slope_fs, ...
                'Color', colors(gi, :) * 0.7, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        end

        % Panel label (data coordinates)
        text(xl(1) - 0.14*diff(xl), yl(2) + 0.04*diff(yl), ...
            ['\bf' panel_labels{si+2}], ...
            'FontSize', panel_fs, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'bottom');

        hold off;
    end

    %% ================================================================
    %  SAVE
    %  ================================================================

    out_pdf = fullfile(paths.figures, 'SI_Fig8_sex_stratified.pdf');
    out_png = fullfile(paths.figures, 'SI_Fig8_sex_stratified.png');
    out_fig = fullfile(paths.figures, 'SI_Fig8_sex_stratified.fig');

    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [fig_w_cm fig_h_cm]);
    set(fig, 'PaperPosition', [0 0 fig_w_cm fig_h_cm]);

    exportgraphics(fig, out_png, 'Resolution', 300);
    fprintf('  Saved: %s (raster, 300 dpi)\n', out_png);

    exportgraphics(fig, out_pdf, 'ContentType', 'image', 'Resolution', 300);
    fprintf('  Saved: %s (raster-in-PDF, 300 dpi)\n', out_pdf);

    fprintf('  Skipped: %s (too many graphic objects)\n', out_fig);
    fprintf('\nSupplementary Figure 8 saved.\n');

    %% ================================================================
    %  CONSOLE SUMMARY
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
            S.interaction(gi).group, S.interaction(gi).beta, S.interaction(gi).p);
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