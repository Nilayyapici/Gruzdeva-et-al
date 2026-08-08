function [slope_results] = analyze_behavior_slopes(run_data, options)
    % analyze_behavior_slopes - Analyzes slopes of dF/F traces for behavior-classified runs
    % Compares eating vs visit runs separately for towards and away directions
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_runs_with_behavior_classification
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to analyze (default: all available)
    %       .directions - Cell array of directions (default: {'towards', 'away'})
    %       .behaviors - Cell array of behaviors (default: {'eating', 'visit'})
    %       .max_distance - Maximum distance for analysis (default: 200)
    %       .slope_alpha - Alpha level for significance (default: 0.05)
    %       .min_data_points - Minimum data points for regression (default: 10)
    %       .slope_range - Distance range around food for testing (default: 50 units)
    %       .smoothing - Window size for smoothing (default: 5)
    %       .plot_results - Whether to create plots (default: true)
    %       .bin_width - Width of distance bins (default: 1)
    %       .min_runs_per_mouse - Minimum runs per mouse to include (default: 3)
    %                            Mouse must have at least this many runs of a specific
    %                            behavior/direction to be included in that analysis
    
    % Set default options
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'sessions'), options.sessions = []; end % Auto-detect
    if ~isfield(options, 'directions'), options.directions = {'towards', 'away'}; end
    if ~isfield(options, 'behaviors'), options.behaviors = {'eating', 'visit'}; end
    if ~isfield(options, 'max_distance'), options.max_distance = 200; end
    if ~isfield(options, 'slope_alpha'), options.slope_alpha = 0.05; end
    if ~isfield(options, 'min_data_points'), options.min_data_points = 10; end
    if ~isfield(options, 'slope_range'), options.slope_range = 50; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'plot_results'), options.plot_results = true; end
    if ~isfield(options, 'bin_width'), options.bin_width = 1; end
    if ~isfield(options, 'min_runs_per_mouse'), options.min_runs_per_mouse = 3; end
    
    % Auto-detect sessions if not provided
    if isempty(options.sessions)
        all_sessions = [];
        for m = 1:length(run_data)
            if ~ismember(run_data(m).session, all_sessions)
                all_sessions = [all_sessions, run_data(m).session];
            end
        end
        session_numbers = sort(all_sessions);
        options.sessions = cell(1, length(session_numbers));
        for i = 1:length(session_numbers)
            options.sessions{i} = ['sess' num2str(session_numbers(i))];
        end
    else
        session_numbers = [];
        for i = 1:length(options.sessions)
            session_numbers = [session_numbers, str2double(options.sessions{i}(5:end))];
        end
    end
    
    fprintf('Analyzing behavior-specific slopes (bin width: %.1f)\n', options.bin_width);
    fprintf('Minimum runs per mouse: %d\n', options.min_runs_per_mouse);
    fprintf('Sessions: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf('\n');
    
    % Create distance bins
    bin_width = options.bin_width;
    dist_bins = 0:bin_width:options.max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialize data collection structure
    mouse_data = struct();
    
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                key = sprintf('sess%d_%s_%s', sess, direction, behavior);
                mouse_data.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_data.(key){i} = containers.Map();
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
        
        % Count runs for each direction/behavior combination for this mouse
        run_counts = struct();
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                key = sprintf('%s_%s', direction, behavior);
                
                % Count runs matching this direction and behavior
                count = 0;
                for r = 1:length(runs)
                    if strcmp(runs(r).type, direction) && strcmp(runs(r).behavior, behavior)
                        count = count + 1;
                    end
                end
                run_counts.(key) = count;
            end
        end
        
        % Initialize temporary storage for this mouse's binned data
        mouse_temp_bins = struct();
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                key = sprintf('sess%d_%s_%s', mouse_sess, direction, behavior);
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
            
            % Skip if not a recognized direction or behavior
            if ~ismember(run.type, options.directions) || ~ismember(run.behavior, options.behaviors)
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
                    key = sprintf('sess%d_%s_%s', mouse_sess, run.type, run.behavior);
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; z_scored_dff(indices)];
                end
            end
        end
        
        % Calculate per-mouse averages for each bin (only if mouse has enough runs)
        for d = 1:length(options.directions)
            direction = options.directions{d};
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                key = sprintf('sess%d_%s_%s', mouse_sess, direction, behavior);
                count_key = sprintf('%s_%s', direction, behavior);
                
                % Only include this mouse if it has enough runs for this combination
                if run_counts.(count_key) >= options.min_runs_per_mouse
                    for i = 1:length(dist_bins)-1
                        if ~isempty(mouse_temp_bins.(key){i})
                            mouse_bin_avg = mean(mouse_temp_bins.(key){i});
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
        
        for d = 1:length(options.directions)
            direction = options.directions{d};
            
            fprintf('\n===== Session %d, %s Direction =====\n', sess, direction);
            
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                
                key = sprintf('sess%d_%s_%s', sess, direction, behavior);
                
                if ~isfield(mouse_data, key)
                    fprintf('  Behavior: %s - No data available.\n', behavior);
                    continue;
                end
                
                % Calculate mean dF/F for each bin
                means = nan(length(bin_centers), 1);
                n_mice_per_bin = zeros(length(bin_centers), 1);
                
                for i = 1:length(dist_bins)-1
                    bin_mouse_data = mouse_data.(key){i};
                    if bin_mouse_data.Count > 0
                        mouse_means_in_bin = cell2mat(values(bin_mouse_data));
                        means(i) = mean(mouse_means_in_bin);
                        n_mice_per_bin(i) = length(mouse_means_in_bin);
                    end
                end
                
                % Count total unique mice for this condition
                all_mice_this_condition = containers.Map();
                for i = 1:length(dist_bins)-1
                    bin_mouse_data = mouse_data.(key){i};
                    if bin_mouse_data.Count > 0
                        mouse_ids_in_bin = keys(bin_mouse_data);
                        for m_id = 1:length(mouse_ids_in_bin)
                            all_mice_this_condition(mouse_ids_in_bin{m_id}) = 1;
                        end
                    end
                end
                n_mice_total = all_mice_this_condition.Count;
                
                fprintf('  Behavior: %s (±%d units, n=%d mice with ≥%d runs)\n', ...
                    behavior, options.slope_range, n_mice_total, options.min_runs_per_mouse);
                
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
                    x_valid = -x_valid;
                end
                
                % Define slope testing range
                if strcmp(direction, 'towards')
                    range_min = -options.slope_range;
                    range_max = 0;
                else
                    range_min = 0;
                    range_max = options.slope_range;
                end
                
                % Filter data to slope testing range
                range_indices = x_valid >= range_min & x_valid <= range_max;
                x_slope = x_valid(range_indices);
                y_slope = y_smoothed(range_indices);
                n_mice_slope = n_mice_valid(range_indices);
                
                % Perform regression
                if length(x_slope) >= options.min_data_points
                    [slope_stats] = performMixedEffectsRegression(x_slope, y_slope, key, mouse_data, bin_centers, options);
                    
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
                    
                    fprintf('    Slope = %.4f ± %.4f, t = %.2f, p = %.4f, R² = %.3f (n=%d mice)%s\n', ...
                        slope_stats.slope, slope_stats.se_slope, slope_stats.t_stat, ...
                        slope_stats.p_value, slope_stats.r_squared, slope_stats.n_mice, significance_text);
                    
                    % Store results
                    slope_results.detailed.(key) = slope_stats;
                    slope_results.detailed.(key).x_data = x_slope;
                    slope_results.detailed.(key).y_data = y_slope;
                    slope_results.detailed.(key).x_all = x_valid;
                    slope_results.detailed.(key).y_all = y_smoothed;
                    slope_results.detailed.(key).n_mice_per_bin = n_mice_valid;
                    slope_results.detailed.(key).session = sess;
                    slope_results.detailed.(key).direction = direction;
                    slope_results.detailed.(key).behavior = behavior;
                else
                    fprintf('    Insufficient data points (%d)\n', length(x_slope));
                end
            end
        end
    end
    
    % Create summary
    createBehaviorSlopeSummary(slope_results, options);
    
    % Create plots if requested
    if options.plot_results
        plotBehaviorSlopeResults(slope_results, options);
    end
end

function [slope_stats] = performMixedEffectsRegression(x_slope, y_slope, data_key, mouse_data, bin_centers, options)
    % Perform regression accounting for between-mouse variability
    
    mouse_x_data = [];
    mouse_y_data = [];
    mouse_ids_data = [];
    
    % Extract individual mouse data points
    for i = 1:length(bin_centers)
        bin_center = bin_centers(i);
        
        % Transform bin center based on direction
        if contains(data_key, 'towards')
            x_coord = -bin_center;
        else
            x_coord = bin_center;
        end
        
        % Check if this bin is in our slope range
        if any(abs(x_slope - x_coord) < 0.01)
            if isfield(mouse_data, data_key) && i <= length(mouse_data.(data_key))
                bin_mouse_data = mouse_data.(data_key){i};
                if bin_mouse_data.Count > 0
                    mouse_ids = keys(bin_mouse_data);
                    mouse_values = cell2mat(values(bin_mouse_data));
                    
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
        slope_stats = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = max(1, n_mice);
        slope_stats.method = 'simple_regression';
        return;
    end
    
    % Calculate per-mouse slopes
    mouse_slopes = [];
    mouse_intercepts = [];
    
    for m = 1:length(unique_mice)
        mouse_id = unique_mice{m};
        mouse_indices = strcmp(mouse_ids_data, mouse_id);
        
        if sum(mouse_indices) >= 3
            x_mouse = mouse_x_data(mouse_indices);
            y_mouse = mouse_y_data(mouse_indices);
            
            try
                X_mouse = [ones(length(x_mouse), 1), x_mouse];
                beta_mouse = X_mouse \ y_mouse;
                mouse_intercepts = [mouse_intercepts; beta_mouse(1)];
                mouse_slopes = [mouse_slopes; beta_mouse(2)];
            catch
                continue;
            end
        end
    end
    
    % Test if average slope is significantly different from zero
    if length(mouse_slopes) >= 2
        [~, p_value, ~, stats_struct] = ttest(mouse_slopes, 0);
        mean_slope = mean(mouse_slopes);
        se_slope = std(mouse_slopes) / sqrt(length(mouse_slopes));
        t_stat = stats_struct.tstat;
        
        % Calculate R-squared
        X_all = [ones(length(mouse_x_data), 1), mouse_x_data];
        beta_overall = X_all \ mouse_y_data;
        y_pred = X_all * beta_overall;
        ss_tot = sum((mouse_y_data - mean(mouse_y_data)).^2);
        ss_res = sum((mouse_y_data - y_pred).^2);
        r_squared = 1 - (ss_res / ss_tot);
        
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
        slope_stats = performLinearRegression(x_slope, y_slope, options.slope_alpha);
        slope_stats.n_mice = length(unique_mice);
        slope_stats.method = 'simple_regression';
    end
end

function [slope_stats] = performLinearRegression(x_data, y_data, alpha)
    % Perform simple linear regression
    n = length(x_data);
    X = [ones(n, 1), x_data(:)];
    
    beta = X \ y_data(:);
    intercept = beta(1);
    slope = beta(2);
    
    y_pred = X * beta;
    residuals = y_data(:) - y_pred;
    p = 2;
    
    mse = sum(residuals.^2) / (n - p);
    cov_matrix = mse * inv(X' * X);
    se_slope = sqrt(cov_matrix(2, 2));
    
    t_stat = slope / se_slope;
    df = n - p;
    p_value = 2 * (1 - tcdf(abs(t_stat), df));
    
    ss_tot = sum((y_data(:) - mean(y_data)).^2);
    ss_res = sum(residuals.^2);
    r_squared = 1 - (ss_res / ss_tot);
    
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

function createBehaviorSlopeSummary(slope_results, options)
    % Create summary statistics
    
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        slope_results.summary.message = 'No slope analyses completed';
        return;
    end
    
    all_slopes = [];
    all_p_values = [];
    all_r_squared = [];
    significant_count = 0;
    
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        
        all_slopes = [all_slopes; stats.slope];
        all_p_values = [all_p_values; stats.p_value];
        all_r_squared = [all_r_squared; stats.r_squared];
        
        if stats.is_significant
            significant_count = significant_count + 1;
        end
    end
    
    slope_results.summary.total_analyses = length(detailed_fields);
    slope_results.summary.significant_slopes = significant_count;
    slope_results.summary.mean_slope = mean(all_slopes);
    slope_results.summary.std_slope = std(all_slopes);
    slope_results.summary.mean_r_squared = mean(all_r_squared);
    slope_results.summary.median_p_value = median(all_p_values);
    
    fprintf('\n===== BEHAVIOR SLOPE ANALYSIS SUMMARY =====\n');
    fprintf('Total analyses: %d\n', slope_results.summary.total_analyses);
    fprintf('Significant slopes: %d (%.1f%%)\n', ...
        significant_count, 100*significant_count/length(detailed_fields));
    fprintf('Mean slope: %.4f ± %.4f\n', ...
        slope_results.summary.mean_slope, slope_results.summary.std_slope);
    fprintf('Mean R²: %.3f\n', slope_results.summary.mean_r_squared);
end

function plotBehaviorSlopeResults(slope_results, options)
    % Create plots showing behavior-specific slopes
    
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No data to plot.\n');
        return;
    end
    
    % Organize by direction
    directions = options.directions;
    
    for d = 1:length(directions)
        direction = directions{d};
        
        figure('Name', sprintf('Behavior Slopes: %s', direction), ...
               'Position', [100 + (d-1)*800, 100, 1200, 600]);
        
        % Find available sessions for this direction
        sessions = [];
        for i = 1:length(detailed_fields)
            if contains(detailed_fields{i}, direction)
                parts = split(detailed_fields{i}, '_');
                sess_num = str2double(parts{1}(5:end));
                if ~ismember(sess_num, sessions)
                    sessions = [sessions, sess_num];
                end
            end
        end
        sessions = sort(sessions);
        
        % Colors for behaviors
        eating_color = [0.8, 0.2, 0.2];  % Red
        visit_color = [0.2, 0.4, 0.8];   % Blue
        
        for s_idx = 1:length(sessions)
            sess = sessions(s_idx);
            subplot(1, length(sessions), s_idx);
            hold on;
            
            for b = 1:length(options.behaviors)
                behavior = options.behaviors{b};
                key = sprintf('sess%d_%s_%s', sess, direction, behavior);
                
                if isfield(slope_results.detailed, key)
                    stats = slope_results.detailed.(key);
                    
                    % Choose color
                    if strcmp(behavior, 'eating')
                        color = eating_color;
                    else
                        color = visit_color;
                    end
                    
                    % Plot data
                    x_data = stats.x_all;
                    y_data = stats.y_all;
                    
                    plot(x_data, y_data, 'o-', 'Color', color, 'MarkerSize', 4, ...
                         'LineWidth', 1.5, 'DisplayName', sprintf('%s (R²=%.3f)', behavior, stats.r_squared));
                    
                    % Plot regression line if significant
                    if stats.is_significant
                        x_fit = linspace(min(x_data), max(x_data), 100);
                        y_fit = stats.intercept + stats.slope * x_fit;
                        plot(x_fit, y_fit, '--', 'Color', color, 'LineWidth', 2, ...
                             'HandleVisibility', 'off');
                        
                        % Add significance indicator
                        if stats.p_value < 0.001
                            sig_text = '***';
                        elseif stats.p_value < 0.01
                            sig_text = '**';
                        else
                            sig_text = '*';
                        end
                        
                        text(0.9, 0.9 - (b-1)*0.1, sig_text, 'Units', 'normalized', ...
                             'FontSize', 12, 'Color', color, 'FontWeight', 'bold');
                    end
                end
            end
            
            title(sprintf('Session %d', sess), 'FontSize', 12);
            if strcmp(direction, 'towards')
                xlabel('Distance to Food', 'FontSize', 11);
            else
                xlabel('Distance from Food', 'FontSize', 11);
            end
            if s_idx == 1
                ylabel('Z-scored dF/F', 'FontSize', 11);
            end
            legend('Location', 'best', 'FontSize', 9);
            grid on;
            box off;
        end
        
        sgtitle(sprintf('%s Direction: Eating vs Visit', upper(direction)), ...
                'FontSize', 14, 'FontWeight', 'bold');
    end
end