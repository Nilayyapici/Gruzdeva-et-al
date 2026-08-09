function [slope_results] = analyze_combined_nonfood_slopes(run_data, options)
    % analyze_combined_nonfood_slopes - Analyzes slopes of dF/F traces with combined non-food arms
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_mouse_runs
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to analyze (default: all sessions)
    %       .arms - Cell array of arms to analyze (default: {'food', 'nonfood'})
    %       .directions - Cell array of directions (default: {'towards', 'away'})
    %       .max_distance - Maximum distance for analysis (default: 200)
    %       .slope_alpha - Alpha level for significance (default: 0.05)
    %       .min_data_points - Minimum data points for regression (default: 10)
    %       .slope_range - Distance range around food for testing (default: 50 units)
    %                     For 'towards': tests from -slope_range to 0
    %                     For 'away': tests from 0 to +slope_range
    %       .smoothing - Window size for smoothing before slope calculation (default: 5)
    %       .plot_results - Whether to create plots (default: true)
    %       .expected_direction - 'auto', 'positive', 'negative', or 'none' (default: 'auto')
    
    % Set default options
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'sessions')
        % Find all unique session numbers in run_data
        all_sessions = [];
        for m = 1:length(run_data)
            if ~ismember(run_data(m).session, all_sessions)
                all_sessions = [all_sessions, run_data(m).session];
            end
        end
        sessions_to_plot_nums = sort(all_sessions);
        
        % Convert numbers to strings
        options.sessions = cell(1, length(sessions_to_plot_nums));
        for i = 1:length(sessions_to_plot_nums)
            options.sessions{i} = ['sess' num2str(sessions_to_plot_nums(i))];
        end
    end
    
    if ~isfield(options, 'arms'), options.arms = {'food', 'nonfood'}; end % Combined non-food
    if ~isfield(options, 'directions'), options.directions = {'towards', 'away'}; end
    if ~isfield(options, 'max_distance'), options.max_distance = 200; end
    if ~isfield(options, 'slope_alpha'), options.slope_alpha = 0.05; end
    if ~isfield(options, 'min_data_points'), options.min_data_points = 10; end
    if ~isfield(options, 'slope_range'), options.slope_range = 50; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'plot_results'), options.plot_results = true; end
    if ~isfield(options, 'expected_direction'), options.expected_direction = 'auto'; end
    
    % Create session mapping
    session_numbers = [];
    for i = 1:length(options.sessions)
        session_str = options.sessions{i};
        session_num = str2double(session_str(5:end));
        session_numbers = [session_numbers, session_num];
    end
    [session_numbers, sort_idx] = sort(session_numbers);
    options.sessions = options.sessions(sort_idx);
    
    % Create session labels
    session_labels = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        sess_num = session_numbers(i);
        switch sess_num
            case 0
                session_labels{i} = 'Before';
            case 1
                session_labels{i} = 'Learning';
            case 2
                session_labels{i} = 'Test';
            otherwise
                session_labels{i} = sprintf('Session %d', sess_num);
        end
    end
    
    fprintf('Analyzing slopes with combined non-food arms for sessions: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf('\n');
    
    % Create distance bins
    bin_width = 0.5;
    dist_bins = 0:bin_width:options.max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialize data collection structure - store per-mouse data properly
    mouse_data = struct(); % Store individual mouse averages for proper statistics
    
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for a = 1:length(options.arms)
                arm = options.arms{a};
                key = sprintf('sess%d_%s_%s', sess, direction, arm);
                mouse_data.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_data.(key){i} = containers.Map(); % Mouse ID -> mean value for this bin
                end
            end
        end
    end
    
    % Process each mouse
    for m = 1:length(run_data)
        mouse_id = run_data(m).mouse_id;
        mouse_sess = run_data(m).session;
        
        % Skip if session not in sessions to analyze
        sess_str = ['sess' num2str(mouse_sess)];
        if ~ismember(sess_str, options.sessions)
            continue;
        end
        
        runs = run_data(m).runs;
        
        % Initialize temporary storage for this mouse's binned data
        mouse_temp_bins = struct();
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for a = 1:length(options.arms)
                arm = options.arms{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, direction, arm);
                mouse_temp_bins.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_temp_bins.(key){i} = [];
                end
            end
        end
        
        % Collect all dF/F values for z-scoring
        all_dff_values = [];
        for r = 1:length(runs)
            all_dff_values = [all_dff_values; runs(r).dff];
        end
        
        % Calculate z-score parameters
        dff_mean = mean(all_dff_values);
        dff_std = std(all_dff_values);
        
        if dff_std < 1e-10
            warning('Mouse %s has near-zero dF/F standard deviation, skipping.', mouse_id);
            continue;
        end
        
        % Process runs
        for r = 1:length(runs)
            run = runs(r);
            
            % Determine if this is food or combined non-food arm
            if strcmp(run.arm, 'food')
                combined_arm = 'food';
            elseif ismember(run.arm, {'nonfood1', 'nonfood2'})
                combined_arm = 'nonfood'; % Combine both non-food arms
            else
                continue; % Skip unrecognized arms
            end
            
            % Skip if not a recognized direction
            if ~ismember(run.type, options.directions)
                continue;
            end
            
            % Apply distance limit
            valid_distance_idx = run.distance <= options.max_distance;
            if ~any(valid_distance_idx)
                continue;
            end
            
            % Filter data
            filtered_distance = run.distance(valid_distance_idx);
            filtered_dff = run.dff(valid_distance_idx);
            
            % Z-score the dF/F
            z_scored_dff = (filtered_dff - dff_mean) / dff_std;
            
            % Bin the data
            for i = 1:length(dist_bins)-1
                indices = filtered_distance >= dist_bins(i) & filtered_distance < dist_bins(i+1);
                if any(indices)
                    key = sprintf('sess%d_%s_%s', mouse_sess, run.type, combined_arm);
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; z_scored_dff(indices)];
                end
            end
        end
        
        % Calculate per-mouse averages for each bin
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for a = 1:length(options.arms)
                arm = options.arms{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, direction, arm);
                
                for i = 1:length(dist_bins)-1
                    if ~isempty(mouse_temp_bins.(key){i})
                        mouse_bin_avg = mean(mouse_temp_bins.(key){i});
                        
                        % If this mouse already has data in this bin (from multiple non-food arms), 
                        % average with existing
                        if mouse_data.(key){i}.isKey(mouse_id)
                            existing_avg = mouse_data.(key){i}(mouse_id);
                            mouse_data.(key){i}(mouse_id) = (existing_avg + mouse_bin_avg) / 2;
                        else
                            mouse_data.(key){i}(mouse_id) = mouse_bin_avg;
                        end
                    end
                end
            end
        end
    end
    
    % Initialize results structure
    slope_results = struct();
    slope_results.summary = struct();
    slope_results.detailed = struct();
    
    % Analyze slopes for each combination
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        sess_label = session_labels{s};
        
        for d = 1:length(options.directions)
            direction = options.directions{d};
            
            for a = 1:length(options.arms)
                arm = options.arms{a};
                
                fprintf('\n===== Slope Analysis: %s, %s, %s (±%d units from food) =====\n', ...
                    sess_label, arm, direction, options.slope_range);
                
                key = sprintf('sess%d_%s_%s', sess, direction, arm);
                
                if ~isfield(mouse_data, key)
                    fprintf('No data available for this combination.\n');
                    continue;
                end
                
                % Calculate mean dF/F for each bin using proper mouse-level statistics
                means = nan(length(bin_centers), 1);
                n_mice_per_bin = zeros(length(bin_centers), 1);
                
                for i = 1:length(dist_bins)-1
                    bin_mouse_data = mouse_data.(key){i};
                    if bin_mouse_data.Count > 0
                        % Get per-mouse means for this bin
                        mouse_means_in_bin = cell2mat(values(bin_mouse_data));
                        
                        % Average across mice (proper statistical unit)
                        means(i) = mean(mouse_means_in_bin);
                        n_mice_per_bin(i) = length(mouse_means_in_bin);
                    end
                end
                
                % Apply smoothing
                valid = ~isnan(means);
                if sum(valid) > options.smoothing
                    x_valid = bin_centers(valid);
                    y_valid = means(valid);
                    y_smoothed = movmean(y_valid, options.smoothing);
                    n_mice_valid = n_mice_per_bin(valid);
                else
                    x_valid = bin_centers(valid);
                    y_smoothed = means(valid);
                    n_mice_valid = n_mice_per_bin(valid);
                end
                
                % Transform coordinates for 'towards' direction
                if strcmp(direction, 'towards')
                    x_valid = -x_valid; % Make negative for towards runs
                end
                
                % Define slope testing range - limited distance around food (x=0)
                if strcmp(direction, 'towards')
                    % For towards: test from -slope_range to 0 (approaching food)
                    range_min = -options.slope_range;
                    range_max = 0;
                    if strcmp(options.expected_direction, 'auto')
                        expected_direction = 'negative'; % Activity decreases approaching food
                    else
                        expected_direction = options.expected_direction;
                    end
                else
                    % For away: test from 0 to +slope_range (leaving food)
                    range_min = 0;
                    range_max = options.slope_range;
                    if strcmp(options.expected_direction, 'auto')
                        expected_direction = 'positive'; % Activity increases leaving food
                    else
                        expected_direction = options.expected_direction;
                    end
                end
                
                % Filter data to slope testing range
                range_indices = x_valid >= range_min & x_valid <= range_max;
                x_slope = x_valid(range_indices);
                y_slope = y_smoothed(range_indices);
                n_mice_slope = n_mice_valid(range_indices);
                
                % Perform regression accounting for between-mouse variability
                if length(x_slope) >= options.min_data_points
                    [slope_stats] = performMixedEffectsRegression(x_slope, y_slope, key, mouse_data, bin_centers, options);
                    
                    % Check direction
                    if strcmp(expected_direction, 'positive')
                        correct_direction = slope_stats.slope > 0;
                        direction_text = sprintf(' (%s direction)', expected_direction);
                    elseif strcmp(expected_direction, 'negative')
                        correct_direction = slope_stats.slope < 0;
                        direction_text = sprintf(' (%s direction)', expected_direction);
                    else
                        correct_direction = true; % No direction expectation
                        direction_text = '';
                    end
                    
                    if ~correct_direction
                        direction_text = sprintf(' (opposite to %s)', expected_direction);
                    end
                    
                    % Significance stars
                    significance_text = '';
                    if slope_stats.is_significant
                        if slope_stats.p_value < 0.001
                            significance_text = ' ***';
                        elseif slope_stats.p_value < 0.01
                            significance_text = ' **';
                        elseif slope_stats.p_value < 0.05
                            significance_text = ' *';
                        end
                    end
                    
                    % Print results with sample size info
                    fprintf('Slope = %.4f ± %.4f, t = %.2f, p = %.4f, R² = %.3f (n=%d mice, avg %.1f mice/bin)%s%s\n', ...
                        slope_stats.slope, slope_stats.se_slope, slope_stats.t_stat, ...
                        slope_stats.p_value, slope_stats.r_squared, slope_stats.n_mice, ...
                        mean(n_mice_slope), direction_text, significance_text);
                    
                    % Store results
                    result_key = sprintf('sess%d_%s_%s', sess, direction, arm);
                    slope_results.detailed.(result_key) = slope_stats;
                    slope_results.detailed.(result_key).correct_direction = correct_direction;
                    slope_results.detailed.(result_key).expected_direction = expected_direction;
                    slope_results.detailed.(result_key).x_data = x_slope;
                    slope_results.detailed.(result_key).y_data = y_slope;
                    slope_results.detailed.(result_key).x_all = x_valid;
                    slope_results.detailed.(result_key).y_all = y_smoothed;
                    slope_results.detailed.(result_key).n_mice_per_bin = n_mice_valid;
                    
                else
                    fprintf('Insufficient data points (%d) for slope analysis\n', length(x_slope));
                end
            end
        end
    end
    
    % Create summary statistics
    createSlopeSummary(slope_results, options);
    
    % Create plots if requested
    if options.plot_results
        plotCombinedSlopeResults(slope_results, options, session_labels);
    end
end

function [slope_stats] = performMixedEffectsRegression(x_slope, y_slope, data_key, mouse_data, bin_centers, options)
    % Perform regression accounting for between-mouse variability
    % Uses individual mouse data points rather than just the mean trace
    
    % Extract individual mouse data points for the slope range
    mouse_x_data = [];
    mouse_y_data = [];
    mouse_ids_data = [];
    
    % Find which bins correspond to our slope range
    for i = 1:length(bin_centers)
        bin_center = bin_centers(i);
        
        % Transform bin center based on direction (same as in main function)
        if contains(data_key, 'towards')
            x_coord = -bin_center;
        else
            x_coord = bin_center;
        end
        
        % Check if this bin is in our slope range (with small tolerance for floating point)
        if any(abs(x_slope - x_coord) < 0.01)
            if isfield(mouse_data, data_key) && i <= length(mouse_data.(data_key))
                bin_mouse_data = mouse_data.(data_key){i};
                if bin_mouse_data.Count > 0
                    mouse_ids = keys(bin_mouse_data);
                    mouse_values = cell2mat(values(bin_mouse_data));
                    
                    % Add data points
                    mouse_x_data = [mouse_x_data; repmat(x_coord, length(mouse_values), 1)];
                    mouse_y_data = [mouse_y_data; mouse_values(:)];
                    mouse_ids_data = [mouse_ids_data; mouse_ids(:)];
                end
            end
        end
    end
    
    % Check if we have enough data
    unique_mice = unique(mouse_ids_data);
    n_mice = length(unique_mice);
    
    if length(mouse_x_data) < options.min_data_points || n_mice < 2
        % Fall back to simple regression if not enough mouse data
        slope_stats = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = max(1, n_mice); % Use actual number or 1 as minimum
        slope_stats.method = 'simple_regression';
        return;
    end
    
    % Calculate per-mouse slopes
    mouse_slopes = [];
    mouse_intercepts = [];
    
    for m = 1:length(unique_mice)
        mouse_id = unique_mice{m};
        mouse_indices = strcmp(mouse_ids_data, mouse_id);
        
        if sum(mouse_indices) >= 3 % Need at least 3 points per mouse
            x_mouse = mouse_x_data(mouse_indices);
            y_mouse = mouse_y_data(mouse_indices);
            
            % Simple regression for this mouse
            try
                X_mouse = [ones(length(x_mouse), 1), x_mouse];
                beta_mouse = X_mouse \ y_mouse;
                mouse_intercepts = [mouse_intercepts; beta_mouse(1)];
                mouse_slopes = [mouse_slopes; beta_mouse(2)];
            catch
                % Skip this mouse if regression fails
                continue;
            end
        end
    end
    
    % Test if the average slope across mice is significantly different from zero
    if length(mouse_slopes) >= 2
        % One-sample t-test on the mouse slopes
        [~, p_value, ~, stats_struct] = ttest(mouse_slopes, 0);
        mean_slope = mean(mouse_slopes);
        se_slope = std(mouse_slopes) / sqrt(length(mouse_slopes));
        t_stat = stats_struct.tstat;
        
        % Calculate R-squared using the overall fit
        X_all = [ones(length(mouse_x_data), 1), mouse_x_data];
        beta_overall = X_all \ mouse_y_data;
        y_pred = X_all * beta_overall;
        ss_tot = sum((mouse_y_data - mean(mouse_y_data)).^2);
        ss_res = sum((mouse_y_data - y_pred).^2);
        r_squared = 1 - (ss_res / ss_tot);
        
        % Store results
        slope_stats = struct();
        slope_stats.slope = mean_slope;
        slope_stats.intercept = mean(mouse_intercepts);
        slope_stats.se_slope = se_slope;
        slope_stats.t_stat = t_stat;
        slope_stats.p_value = p_value;
        slope_stats.r_squared = r_squared;
        slope_stats.is_significant = p_value < options.slope_alpha;
        slope_stats.n_points = length(mouse_x_data);
        slope_stats.n_mice = length(mouse_slopes);
        slope_stats.method = 'mixed_effects';
        slope_stats.mouse_slopes = mouse_slopes;
        slope_stats.mouse_intercepts = mouse_intercepts;
        
    else
        % Fall back to simple regression
        slope_stats = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = length(unique_mice);
        slope_stats.method = 'simple_regression';
    end
end

function [slope_stats] = performLinearRegression(x_data, y_data, alpha)
    % Perform simple linear regression (fallback method)
    n = length(x_data);
    X = [ones(n, 1), x_data(:)]; % Design matrix with intercept
    
    % Calculate regression coefficients
    beta = X \ y_data(:);
    intercept = beta(1);
    slope = beta(2);
    
    % Calculate statistics
    y_pred = X * beta;
    residuals = y_data(:) - y_pred;
    p = 2; % Number of parameters
    
    % Standard error of slope
    mse = sum(residuals.^2) / (n - p);
    cov_matrix = mse * inv(X' * X);
    se_slope = sqrt(cov_matrix(2, 2));
    
    % T-test for slope
    t_stat = slope / se_slope;
    df = n - p;
    p_value = 2 * (1 - tcdf(abs(t_stat), df));
    
    % R-squared
    ss_tot = sum((y_data(:) - mean(y_data)).^2);
    ss_res = sum(residuals.^2);
    r_squared = 1 - (ss_res / ss_tot);
    
    % Store results
    slope_stats = struct();
    slope_stats.slope = slope;
    slope_stats.intercept = intercept;
    slope_stats.se_slope = se_slope;
    slope_stats.t_stat = t_stat;
    slope_stats.p_value = p_value;
    slope_stats.r_squared = r_squared;
    slope_stats.is_significant = p_value < alpha;
    slope_stats.n_points = n;
    slope_stats.method = 'simple_regression';
end

function createSlopeSummary(slope_results, options)
    % Create summary statistics across all conditions
    
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        slope_results.summary.message = 'No slope analyses completed';
        return;
    end
    
    % Collect all slopes and statistics
    all_slopes = [];
    all_p_values = [];
    all_r_squared = [];
    significant_count = 0;
    correct_direction_count = 0;
    
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        
        all_slopes = [all_slopes; stats.slope];
        all_p_values = [all_p_values; stats.p_value];
        all_r_squared = [all_r_squared; stats.r_squared];
        
        if stats.is_significant
            significant_count = significant_count + 1;
        end
        
        if stats.correct_direction
            correct_direction_count = correct_direction_count + 1;
        end
    end
    
    % Summary statistics
    slope_results.summary.total_analyses = length(detailed_fields);
    slope_results.summary.significant_slopes = significant_count;
    slope_results.summary.correct_direction = correct_direction_count;
    slope_results.summary.mean_slope = mean(all_slopes);
    slope_results.summary.std_slope = std(all_slopes);
    slope_results.summary.mean_r_squared = mean(all_r_squared);
    slope_results.summary.median_p_value = median(all_p_values);
    
    % Print summary
    fprintf('\n===== SLOPE ANALYSIS SUMMARY (Combined Non-Food) =====\n');
    fprintf('Total analyses performed: %d\n', slope_results.summary.total_analyses);
    fprintf('Significant slopes (p < %.3f): %d (%.1f%%)\n', options.slope_alpha, ...
        significant_count, 100*significant_count/length(detailed_fields));
    fprintf('Slopes in expected direction: %d (%.1f%%)\n', ...
        correct_direction_count, 100*correct_direction_count/length(detailed_fields));
    fprintf('Mean slope: %.4f ± %.4f\n', slope_results.summary.mean_slope, slope_results.summary.std_slope);
    fprintf('Mean R²: %.3f\n', slope_results.summary.mean_r_squared);
    fprintf('Median p-value: %.4f\n', slope_results.summary.median_p_value);
end

function plotCombinedSlopeResults(slope_results, options, session_labels)
    % Create plots showing slope analysis results with combined non-food arms
    
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No data to plot.\n');
        return;
    end
    
    % Organize data by session, direction, and arm
    sessions = unique(cellfun(@(x) str2double(x(5)), detailed_fields, 'UniformOutput', true));
    directions = options.directions;
    arms = options.arms; % Should be {'food', 'nonfood'}
    
    % Create figure
    n_sessions = length(sessions);
    n_directions = length(directions);
    
    figure('Name', 'Slope Analysis Results (Combined Non-Food)', 'Position', [100, 100, 400*n_directions, 300*n_sessions]);
    
    plot_idx = 1;
    for s = 1:n_sessions
        sess_num = sessions(s);
        
        for d = 1:n_directions
            direction = directions{d};
            
            subplot(n_sessions, n_directions, plot_idx);
            hold on;
            
            % Define colors for each arm type
            arm_colors = {[0, 0, 0.8], [0.8, 0, 0]}; % Blue for food, Red for non-food
            
            for a = 1:length(arms)
                arm = arms{a};
                field_name = sprintf('sess%d_%s_%s', sess_num, direction, arm);
                
                if isfield(slope_results.detailed, field_name)
                    stats = slope_results.detailed.(field_name);
                    
                    % Plot the data and regression line
                    x_data = stats.x_all;
                    y_data = stats.y_all;
                    
                    % Create proper arm label
                    if strcmp(arm, 'food')
                        arm_label = 'Food';
                    else
                        arm_label = 'Non-Food';
                    end
                    
                    % Plot data points
                    plot(x_data, y_data, 'o', 'Color', arm_colors{a}, 'MarkerSize', 4, ...
                         'DisplayName', sprintf('%s (R²=%.3f)', arm_label, stats.r_squared));
                    
                    % Plot regression line if significant
                    if stats.is_significant && stats.correct_direction
                        x_fit = linspace(min(x_data), max(x_data), 100);
                        y_fit = stats.intercept + stats.slope * x_fit;
                        plot(x_fit, y_fit, '-', 'Color', arm_colors{a}, 'LineWidth', 2, ...
                             'HandleVisibility', 'off');
                    end
                end
            end
            
            % Formatting
            sess_idx = find(sessions == sess_num);
            title(sprintf('%s - %s', session_labels{sess_idx}, direction));
            if strcmp(direction, 'towards')
                xlabel('Distance to Food');
            else
                xlabel('Distance from Food');
            end
            ylabel('Z-scored dF/F');
            legend('Location', 'best');
            legend('boxoff')
            grid off;
            
            plot_idx = plot_idx + 1;
        end
    end
    
    sgtitle('Slope Analysis: Linear Trends Toward Food Location (Combined Non-Food)', 'FontSize', 14);
end