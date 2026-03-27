function [mode_est, ci] = bootstrap_mode(data, n_boot)
% BOOTSTRAP_MODE - Estimate the mode with bootstrap confidence intervals
%
% Usage:
%   [mode_est, ci] = bootstrap_mode(data, n_boot)
%
% Estimates the mode of a distribution using kernel density estimation,
% with 95% bootstrap confidence intervals.
%
% Inputs:
%   data   - Vector of observations
%   n_boot - Number of bootstrap resamples (default: 5000)
%
% Outputs:
%   mode_est - Point estimate of the mode
%   ci       - [lower, upper] 95% confidence interval

    if nargin < 2
        n_boot = 5000;
    end

    data = data(:);
    data = data(~isnan(data));

    if length(data) < 3
        mode_est = NaN;
        ci = [NaN, NaN];
        return;
    end

    % Point estimate via kernel density
    [f, x] = ksdensity(data, 'NumPoints', 300);
    [~, idx] = max(f);
    mode_est = x(idx);

    % Bootstrap confidence interval
    n = length(data);
    boot_modes = zeros(n_boot, 1);

    rng(42);  % Reproducibility
    for b = 1:n_boot
        boot_sample = data(randi(n, n, 1));
        [bf, bx] = ksdensity(boot_sample, 'NumPoints', 200);
        [~, bidx] = max(bf);
        boot_modes(b) = bx(bidx);
    end

    ci = [prctile(boot_modes, 2.5), prctile(boot_modes, 97.5)];
end
