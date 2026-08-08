function [slope_results] = analyze_single_arm_slopes(run_data, options)
    % analyze_single_arm_slopes - Analyzes slopes of dF/F traces for single food arm maze
    % Compares sess0 (no access) vs sess1 (food available) for towards and away runs
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_single_arm_runs_with_speed_check
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to analyze (default: {'sess0', 'sess1'})
    %       .directions - Cell array of directions (default: {'towards', 'away'})
    %       .max_distance - Maximum distance for analysis (default: 200)
    %       .slope_alpha - Alpha level for significance (default: 0.05)
    %       .min_data_points - Minimum data points for regression (default: 10)
    %       .slope_range - Distance range around food for testing (default: 50 units)
    %       .smoothing - Window size for smoothing before slope calculation (default: 5)
    %       .plot_results - Whether to create plots (default: true)
    %       .expected_direction - 'auto', 'positive', 'negative', or 'none' (default: 'auto')
    %       .bin_width - Width of distance bins for averaging (default: 1)

    if nargin < 2, options = struct(); end

    if ~isfield(options, 'sessions'),           options.sessions           = {'sess0', 'sess1'}; end
    if ~isfield(options, 'directions'),         options.directions         = {'towards', 'away'}; end
    if ~isfield(options, 'max_distance'),       options.max_distance       = 200; end
    if ~isfield(options, 'slope_alpha'),        options.slope_alpha        = 0.05; end
    if ~isfield(options, 'min_data_points'),    options.min_data_points    = 10; end
    if ~isfield(options, 'slope_range'),        options.slope_range        = 50; end
    if ~isfield(options, 'smoothing'),          options.smoothing          = 5; end
    if ~isfield(options, 'plot_results'),       options.plot_results       = true; end
    if ~isfield(options, 'expected_direction'), options.expected_direction = 'auto'; end
    if ~isfield(options, 'bin_width'),          options.bin_width          = 1; end

    % Session mapping
    session_numbers = [];
    for i = 1:length(options.sessions)
        session_numbers(end+1) = str2double(options.sessions{i}(5:end));
    end
    [session_numbers, sort_idx] = sort(session_numbers);
    options.sessions = options.sessions(sort_idx);

    session_labels = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        switch session_numbers(i)
            case 0,   session_labels{i} = 'No Access';
            case 1,   session_labels{i} = 'Food Available';
            otherwise, session_labels{i} = sprintf('Session %d', session_numbers(i));
        end
    end

    fprintf('Analyzing single-arm slopes for sessions (bin width: %.1f): ', options.bin_width);
    for i = 1:length(options.sessions)
        fprintf('%s (%s) ', options.sessions{i}, session_labels{i});
    end
    fprintf('\n');

    bin_width   = options.bin_width;
    dist_bins   = 0:bin_width:options.max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;

    % Initialise per-mouse data store
    mouse_data = struct();
    for s = 1:length(session_numbers)
        for d = 1:length(options.directions)
            key = sprintf('sess%d_%s', session_numbers(s), options.directions{d});
            mouse_data.(key) = cell(length(dist_bins)-1, 1);
            for i = 1:length(dist_bins)-1
                mouse_data.(key){i} = containers.Map();
            end
        end
    end

    % Process each mouse
    for m = 1:length(run_data)
        mouse_id   = run_data(m).mouse_id;
        mouse_sess = run_data(m).session;
        sess_str   = ['sess' num2str(mouse_sess)];
        if ~ismember(sess_str, options.sessions), continue; end

        runs = run_data(m).runs;

        mouse_temp_bins = struct();
        for d = 1:length(options.directions)
            key = sprintf('sess%d_%s', mouse_sess, options.directions{d});
            mouse_temp_bins.(key) = cell(length(dist_bins)-1, 1);
            for i = 1:length(dist_bins)-1, mouse_temp_bins.(key){i} = []; end
        end

        all_dff = [];
        for r = 1:length(runs), all_dff = [all_dff; runs(r).dff]; end
        dff_mean = mean(all_dff);
        dff_std  = std(all_dff);
        if dff_std < 1e-10
            warning('Mouse %s has near-zero dF/F std, skipping.', mouse_id); continue;
        end

        for r = 1:length(runs)
            run = runs(r);
            if ~ismember(run.type, options.directions), continue; end
            valid_idx = run.distance <= options.max_distance;
            if ~any(valid_idx), continue; end
            z_dff = (run.dff(valid_idx) - dff_mean) / dff_std;
            dist  = run.distance(valid_idx);
            key   = sprintf('sess%d_%s', mouse_sess, run.type);
            for i = 1:length(dist_bins)-1
                idx = dist >= dist_bins(i) & dist < dist_bins(i+1);
                if any(idx)
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; z_dff(idx)];
                end
            end
        end

        for d = 1:length(options.directions)
            key = sprintf('sess%d_%s', mouse_sess, options.directions{d});
            for i = 1:length(dist_bins)-1
                if ~isempty(mouse_temp_bins.(key){i})
                    mouse_data.(key){i}(mouse_id) = mean(mouse_temp_bins.(key){i});
                end
            end
        end
    end

    % Analyse slopes
    slope_results = struct();
    slope_results.summary  = struct();
    slope_results.detailed = struct();

    for s = 1:length(session_numbers)
        sess       = session_numbers(s);
        sess_label = session_labels{s};

        for d = 1:length(options.directions)
            direction = options.directions{d};
            fprintf('\n===== Slope Analysis: %s, %s (±%d units from food) =====\n', ...
                    sess_label, direction, options.slope_range);

            key = sprintf('sess%d_%s', sess, direction);
            if ~isfield(mouse_data, key)
                fprintf('No data available.\n'); continue;
            end

            means = nan(length(bin_centers), 1);
            n_mice_per_bin = zeros(length(bin_centers), 1);
            for i = 1:length(dist_bins)-1
                bmd = mouse_data.(key){i};
                if bmd.Count > 0
                    vals = cell2mat(values(bmd));
                    means(i)         = mean(vals);
                    n_mice_per_bin(i) = length(vals);
                end
            end

            valid = ~isnan(means);
            if sum(valid) > options.smoothing
                x_valid      = bin_centers(valid);
                y_smoothed   = movmean(means(valid), options.smoothing);
                n_mice_valid = n_mice_per_bin(valid);
            else
                x_valid      = bin_centers(valid);
                y_smoothed   = means(valid);
                n_mice_valid = n_mice_per_bin(valid);
            end

            if strcmp(direction, 'towards'), x_valid = -x_valid; end

            if strcmp(direction, 'towards')
                range_min = -options.slope_range; range_max = 0;
                expected_direction = ternary(strcmp(options.expected_direction,'auto'), ...
                                            'negative', options.expected_direction);
            else
                range_min = 0; range_max = options.slope_range;
                expected_direction = ternary(strcmp(options.expected_direction,'auto'), ...
                                            'positive', options.expected_direction);
            end

            ri      = x_valid >= range_min & x_valid <= range_max;
            x_slope = x_valid(ri);
            y_slope = y_smoothed(ri);
            n_mice_slope = n_mice_valid(ri);

            if length(x_slope) >= options.min_data_points
                slope_stats = performMixedEffectsRegression(x_slope, y_slope, key, ...
                                                             mouse_data, bin_centers, options);

                correct_direction = strcmp(expected_direction,'positive') && slope_stats.slope > 0 || ...
                                    strcmp(expected_direction,'negative') && slope_stats.slope < 0 || ...
                                    strcmp(expected_direction,'none');

                sig_str = '';
                if slope_stats.is_significant
                    if     slope_stats.p_value < 0.001, sig_str = ' ***';
                    elseif slope_stats.p_value < 0.01,  sig_str = ' **';
                    else,                                sig_str = ' *';
                    end
                end

                fprintf('Slope = %.4f ± %.4f, t = %.2f, p = %.4f, R² = %.3f (n=%d mice, avg %.1f mice/bin)%s\n', ...
                    slope_stats.slope, slope_stats.se_slope, slope_stats.t_stat, ...
                    slope_stats.p_value, slope_stats.r_squared, slope_stats.n_mice, ...
                    mean(n_mice_slope), sig_str);

                result_key = sprintf('sess%d_%s', sess, direction);
                slope_results.detailed.(result_key)                   = slope_stats;
                slope_results.detailed.(result_key).correct_direction = correct_direction;
                slope_results.detailed.(result_key).expected_direction = expected_direction;
                slope_results.detailed.(result_key).x_data            = x_slope;
                slope_results.detailed.(result_key).y_data            = y_slope;
                slope_results.detailed.(result_key).x_all             = x_valid;
                slope_results.detailed.(result_key).y_all             = y_smoothed;
                slope_results.detailed.(result_key).n_mice_per_bin    = n_mice_valid;
                slope_results.detailed.(result_key).session_label     = sess_label;
            else
                fprintf('Insufficient data points (%d) for slope analysis\n', length(x_slope));
            end
        end
    end

    createSlopeSummary(slope_results, options);
    if options.plot_results
        plotSingleArmSlopeResults(slope_results, options, session_labels);
    end
end

% ==========================================================================

function slope_stats = performMixedEffectsRegression(x_slope, y_slope, data_key, ...
                                                      mouse_data, bin_centers, options)
    mouse_x_data   = [];
    mouse_y_data   = [];
    mouse_ids_data = {};

    for i = 1:length(bin_centers)
        bin_center = bin_centers(i);
        x_coord = ternary(contains(data_key,'towards'), -bin_center, bin_center);

        if any(abs(x_slope - x_coord) < 0.01)
            if isfield(mouse_data, data_key) && i <= length(mouse_data.(data_key))
                bmd = mouse_data.(data_key){i};
                if bmd.Count > 0
                    mids  = keys(bmd);
                    mvals = cell2mat(values(bmd));
                    n     = length(mvals);
                    mouse_x_data   = [mouse_x_data;   repmat(x_coord, n, 1)];
                    mouse_y_data   = [mouse_y_data;   mvals(:)];
                    mouse_ids_data = [mouse_ids_data; mids(:)];
                end
            end
        end
    end

    unique_mice = unique(mouse_ids_data);
    n_mice      = length(unique_mice);

    if length(mouse_x_data) < options.min_data_points || n_mice < 2
        slope_stats        = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = max(1, n_mice);
        slope_stats.method = 'simple_regression';
        slope_stats.mouse_slopes = [];
        slope_stats.mouse_ids    = {};
        return;
    end

    mouse_slopes     = [];
    mouse_intercepts = [];
    valid_mice       = {};   % tracks IDs that actually produced a slope

    for m = 1:length(unique_mice)
        mid     = unique_mice{m};
        idx     = strcmp(mouse_ids_data, mid);
        x_mouse = mouse_x_data(idx);
        y_mouse = mouse_y_data(idx);

        if sum(idx) >= 3
            try
                beta = [ones(length(x_mouse),1), x_mouse] \ y_mouse;
                mouse_intercepts(end+1) = beta(1);  %#ok<AGROW>
                mouse_slopes(end+1)     = beta(2);  %#ok<AGROW>
                valid_mice{end+1}       = mid;      %#ok<AGROW>
            catch
                % skip mice where regression fails
            end
        end
    end

    if length(mouse_slopes) >= 2
        [~, p_value, ~, st] = ttest(mouse_slopes(:), 0);
        mean_slope = mean(mouse_slopes);
        se_slope   = std(mouse_slopes) / sqrt(length(mouse_slopes));

        X_all     = [ones(length(mouse_x_data),1), mouse_x_data];
        beta_all  = X_all \ mouse_y_data;
        y_pred    = X_all * beta_all;
        ss_tot    = sum((mouse_y_data - mean(mouse_y_data)).^2);
        r_squared = 1 - sum((mouse_y_data - y_pred).^2) / ss_tot;

        slope_stats                  = struct();
        slope_stats.slope            = mean_slope;
        slope_stats.intercept        = mean(mouse_intercepts);
        slope_stats.se_slope         = se_slope;
        slope_stats.t_stat           = st.tstat;
        slope_stats.p_value          = p_value;
        slope_stats.r_squared        = r_squared;
        slope_stats.is_significant   = p_value < options.slope_alpha;
        slope_stats.n_points         = length(mouse_x_data);
        slope_stats.n_mice           = length(mouse_slopes);
        slope_stats.method           = 'mixed_effects';
        slope_stats.mouse_slopes     = mouse_slopes(:);
        slope_stats.mouse_intercepts = mouse_intercepts(:);
        slope_stats.mouse_ids        = valid_mice(:);  % 1-to-1 with mouse_slopes
    else
        slope_stats        = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = length(unique_mice);
        slope_stats.method = 'simple_regression';
        slope_stats.mouse_slopes = [];
        slope_stats.mouse_ids    = {};
    end
end

% ==========================================================================

function slope_stats = performLinearRegression(x_data, y_data, alpha)
    n   = length(x_data);
    X   = [ones(n,1), x_data(:)];
    beta = X \ y_data(:);
    y_pred    = X * beta;
    residuals = y_data(:) - y_pred;
    mse       = sum(residuals.^2) / (n - 2);
    cov_mat   = mse * inv(X'*X);
    se_slope  = sqrt(cov_mat(2,2));
    t_stat    = beta(2) / se_slope;
    p_value   = 2 * (1 - tcdf(abs(t_stat), n-2));
    ss_tot    = sum((y_data(:) - mean(y_data)).^2);
    r_squared = 1 - sum(residuals.^2) / ss_tot;

    slope_stats              = struct();
    slope_stats.slope        = beta(2);
    slope_stats.intercept    = beta(1);
    slope_stats.se_slope     = se_slope;
    slope_stats.t_stat       = t_stat;
    slope_stats.p_value      = p_value;
    slope_stats.r_squared    = r_squared;
    slope_stats.is_significant = p_value < alpha;
    slope_stats.n_points     = n;
    slope_stats.method       = 'simple_regression';
end

% ==========================================================================

function createSlopeSummary(slope_results, options)
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        slope_results.summary.message = 'No slope analyses completed'; return;
    end

    all_slopes = []; all_p = []; all_r2 = [];
    n_sig = 0; n_correct = 0;
    for i = 1:length(detailed_fields)
        st = slope_results.detailed.(detailed_fields{i});
        all_slopes(end+1) = st.slope;
        all_p(end+1)      = st.p_value;
        all_r2(end+1)     = st.r_squared;
        if st.is_significant,  n_sig     = n_sig + 1;     end
        if st.correct_direction, n_correct = n_correct + 1; end
    end

    n = length(detailed_fields);
    slope_results.summary.total_analyses    = n;
    slope_results.summary.significant_slopes = n_sig;
    slope_results.summary.correct_direction = n_correct;
    slope_results.summary.mean_slope        = mean(all_slopes);
    slope_results.summary.std_slope         = std(all_slopes);
    slope_results.summary.mean_r_squared    = mean(all_r2);
    slope_results.summary.median_p_value    = median(all_p);

    fprintf('\n===== SINGLE-ARM SLOPE ANALYSIS SUMMARY =====\n');
    fprintf('Total analyses: %d\n', n);
    fprintf('Significant (p<%.3f): %d (%.1f%%)\n', options.slope_alpha, n_sig, 100*n_sig/n);
    fprintf('Correct direction: %d (%.1f%%)\n', n_correct, 100*n_correct/n);
    fprintf('Mean slope: %.4f ± %.4f\n', mean(all_slopes), std(all_slopes));
    fprintf('Mean R²: %.3f\n', mean(all_r2));
    fprintf('Median p: %.4f\n', median(all_p));
end

% ==========================================================================

function plotSingleArmSlopeResults(slope_results, options, session_labels)
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields), fprintf('No data to plot.\n'); return; end

    sessions   = [0, 1];
    directions = options.directions;
    n_sess = length(sessions); n_dir = length(directions);
    session_colors = {[0.2,0.4,0.8], [0.8,0.2,0.2]};

    figure('Name','Single-Arm Slope Analysis Results', ...
           'Position',[100,100,400*n_dir,300*n_sess]);
    plot_idx = 1;
    for s = 1:n_sess
        for d = 1:n_dir
            subplot(n_sess, n_dir, plot_idx);
            hold on;
            fname = sprintf('sess%d_%s', sessions(s), directions{d});
            if isfield(slope_results.detailed, fname)
                st = slope_results.detailed.(fname);
                plot(st.x_all, st.y_all, 'o', 'Color', session_colors{s}, ...
                     'MarkerSize', 4, 'DisplayName', ...
                     sprintf('%s (R²=%.3f)', st.session_label, st.r_squared));
                if st.is_significant && st.correct_direction
                    xf = linspace(min(st.x_all), max(st.x_all), 100);
                    plot(xf, st.intercept + st.slope*xf, '-', ...
                         'Color', session_colors{s}, 'LineWidth', 2, ...
                         'HandleVisibility','off');
                end
                if st.is_significant
                    if     st.p_value < 0.001, ps = '***';
                    elseif st.p_value < 0.01,  ps = '**';
                    else,                       ps = '*';
                    end
                    text(0.1, 0.9, ['p ' ps], 'Units','normalized', ...
                         'FontSize',12,'FontWeight','bold','Color',session_colors{s});
                end
            end
            title(sprintf('%s - %s', session_labels{s}, directions{d}));
            xlabel(ternary(strcmp(directions{d},'towards'), ...
                           'Distance to Food', 'Distance from Food'));
            ylabel('Z-scored dF/F');
            legend('Location','best'); grid on;
            plot_idx = plot_idx + 1;
        end
    end
    sgtitle('Single-Arm Slope Analysis: dF/F vs Distance to Food','FontSize',14);
end

% ==========================================================================

function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end