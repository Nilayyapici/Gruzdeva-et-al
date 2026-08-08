function plot_fast_slow_comparison(run_data, options)
    % PLOT_FAST_SLOW_COMPARISON_BINNED Plots dF/F vs distance for fast, middle and slow runs
    % using a binning approach with proper SEM bands and customizable speed thresholds
    %
    % Options:
    %   .fast_threshold - Speed threshold for fast runs (default: median speed)
    %   .slow_threshold - Speed threshold for slow runs (default: same as fast_threshold)
    %   .normalize_method - Method for normalizing dF/F (default: 'zscore')
    %   .ylim - Y-axis limits [min max] (default: [-1.5 1.5])
    %   .xlim - X-axis limits [min max] (default: [0 210])
    %   .show_individual_runs - Whether to show individual runs (default: true)
    %   .show_middle_runs - Whether to plot middle runs (default: false)
    %   .bin_width - Width of distance bins (default: 5)
    %   .smoothing_window - Smoothing window size (default: 3, set to 0 to disable)
    %   .min_points_per_bin - Minimum points per bin for mean/SEM calc (default: 3)
    %   .session - Specific session to analyze (optional)
    %   .run_category - Category of runs to include: 'food', 'not_food', 'all' (default: 'food')
    
    % Default options
    if nargin < 2
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'fast_threshold') && ~isfield(options, 'speed_threshold')
        options.fast_threshold = 'median';
    elseif ~isfield(options, 'fast_threshold') && isfield(options, 'speed_threshold')
        options.fast_threshold = options.speed_threshold; % Backwards compatibility
    end
    
    if ~isfield(options, 'slow_threshold')
        if isfield(options, 'fast_threshold')
            options.slow_threshold = options.fast_threshold; % Default to same threshold
        else
            options.slow_threshold = 'median';
        end
    end
    
    if ~isfield(options, 'normalize_method'), options.normalize_method = 'zscore'; end
    if ~isfield(options, 'ylim'), options.ylim = [-1.5 1.5]; end
    if ~isfield(options, 'xlim'), options.xlim = [0 210]; end
    if ~isfield(options, 'show_individual_runs'), options.show_individual_runs = true; end
    if ~isfield(options, 'show_middle_runs'), options.show_middle_runs = false; end % Default: don't show middle runs
    if ~isfield(options, 'bin_width'), options.bin_width = 5; end
    if ~isfield(options, 'smoothing_window'), options.smoothing_window = 3; end
    if ~isfield(options, 'min_points_per_bin'), options.min_points_per_bin = 3; end
    if ~isfield(options, 'run_category'), options.run_category = 'food'; end % Options: 'food', 'not_food', 'all'
    
    % Print data structure info
    fprintf('Input run_data contains %d mice\n', length(run_data));
    
    % Filter by session if specified
    if isfield(options, 'session')
        fprintf('Filtering for session %d...\n', options.session);
        run_data = filter_by_session(run_data, options.session);
        session_str = sprintf('Session %d: ', options.session);
    else
        session_str = 'All Sessions: ';
    end
    
    fprintf('After filtering: %d mice\n', length(run_data));
    
    % Determine which run types to include based on run_category option
    switch options.run_category
        case 'food'
            % Only include runs that reached/started from food
            valid_towards_types = {'towards'};
            valid_away_types = {'away'};
            category_str = ' (Food Runs)';
        case 'not_food'
            % Only include runs that didn't reach/start from food
            valid_towards_types = {'not_food_towards'};
            valid_away_types = {'not_food_away'};
            category_str = ' (Non-Food Runs)';
        case 'all'
            % Include all run types
            valid_towards_types = {'towards', 'not_food_towards'};
            valid_away_types = {'away', 'not_food_away'};
            category_str = ' (All Runs)';
        otherwise
            error('Invalid run_category option. Must be ''food'', ''not_food'', or ''all''');
    end
    
    % Check data structure for debugging
    total_towards = 0;
    total_away = 0;
    for m = 1:length(run_data)
        for r = 1:length(run_data(m).runs)
            run_type = run_data(m).runs(r).type;
            if any(strcmp(run_type, valid_towards_types))
                total_towards = total_towards + 1;
            elseif any(strcmp(run_type, valid_away_types))
                total_away = total_away + 1;
            end
        end
    end
    fprintf('Total runs after filtering%s: %d towards, %d away\n', category_str, total_towards, total_away);
    
    % Calculate all speeds for thresholding (only from valid run categories)
    all_speeds = [];
    for m = 1:length(run_data)
        for r = 1:length(run_data(m).runs)
            run = run_data(m).runs(r);
            
            % Check if this run type should be included
            is_valid_run = any(strcmp(run.type, valid_towards_types)) || ...
                          any(strcmp(run.type, valid_away_types));
            
            if ~is_valid_run
                continue;
            end
            
            if isfield(run, 'calculated_speed') && ~isempty(run.calculated_speed)
                all_speeds = [all_speeds, mean(run.calculated_speed)];
            elseif ~isempty(run.speed)
                all_speeds = [all_speeds, mean(run.speed)];
            end
        end
    end
    
    % Determine fast threshold
    if ischar(options.fast_threshold)
        if strcmp(options.fast_threshold, 'median')
            fast_threshold = median(all_speeds);
        elseif strcmp(options.fast_threshold, 'mean')
            fast_threshold = mean(all_speeds);
        elseif strcmp(options.fast_threshold, 'percentile')
            if isfield(options, 'fast_percentile')
                fast_threshold = prctile(all_speeds, options.fast_percentile);
            else
                fast_threshold = prctile(all_speeds, 75); % Default to 75th percentile
            end
        else
            error('Unknown fast_threshold method');
        end
    else
        fast_threshold = options.fast_threshold;
    end
    
    % Determine slow threshold
    if ischar(options.slow_threshold)
        if strcmp(options.slow_threshold, 'median')
            slow_threshold = median(all_speeds);
        elseif strcmp(options.slow_threshold, 'mean')
            slow_threshold = mean(all_speeds);
        elseif strcmp(options.slow_threshold, 'percentile')
            if isfield(options, 'slow_percentile')
                slow_threshold = prctile(all_speeds, options.slow_percentile);
            else
                slow_threshold = prctile(all_speeds, 25); % Default to 25th percentile
            end
        else
            error('Unknown slow_threshold method');
        end
    else
        slow_threshold = options.slow_threshold;
    end
    
    % Check if thresholds are the same
    if fast_threshold == slow_threshold
        % Using the same threshold in both directions (traditional fast vs slow)
        fprintf('Speed threshold for both fast (>) and slow (<): %.2f cm/s\n', fast_threshold);
        title_str = sprintf('Fast (>%.2f) vs Slow (<%.2f) cm/s%s', fast_threshold, slow_threshold, category_str);
    else
        % Using different thresholds (e.g., very fast vs very slow)
        fprintf('Fast threshold (>=): %.2f cm/s\n', fast_threshold);
        fprintf('Slow threshold (<=): %.2f cm/s\n', slow_threshold);
        
        if options.show_middle_runs
            title_str = sprintf('Fast (≥%.2f) vs Middle vs Slow (≤%.2f) cm/s%s', fast_threshold, slow_threshold, category_str);
        else
            title_str = sprintf('Fast (≥%.2f) vs Slow (≤%.2f) cm/s%s', fast_threshold, slow_threshold, category_str);
        end
    end
    
    % Create figure
    figure('Position', [100, 100, 1000, 500], 'Name', [session_str 'Fast vs Slow Runs: dF/F vs Distance' category_str]);
    
    % Process towards runs
    subplot(1, 2, 1);
    process_and_plot_binned(run_data, valid_towards_types, fast_threshold, slow_threshold, options);
    title('Towards Runs: dF/F vs Distance');
    set(gca, 'XDir', 'reverse'); % Reverse x-axis for towards runs
    
    % Process away runs
    subplot(1, 2, 2);
    process_and_plot_binned(run_data, valid_away_types, fast_threshold, slow_threshold, options);
    title('Away Runs: dF/F vs Distance');
    
    % Add overall title
    if isfield(options, 'session')
        sgtitle(sprintf('Session %d: %s', options.session, title_str), 'FontSize', 14);
    else
        sgtitle(sprintf('All Sessions: %s', title_str), 'FontSize', 14);
    end
end

function filtered_data = filter_by_session(run_data, session_number)
    % Filter run_data to include only the specified session
    filtered_data = struct('mouse_id', {}, 'session', {}, 'runs', {});
    
    for m = 1:length(run_data)
        if run_data(m).session == session_number
            filtered_data = [filtered_data; run_data(m)];
        end
    end
    
    % Check if any data was found for the session
    if isempty(filtered_data)
        warning('No data found for session %d', session_number);
    else
        fprintf('Found %d mice for session %d\n', length(filtered_data), session_number);
    end
end

function process_and_plot_binned(run_data, valid_run_types, fast_threshold, slow_threshold, options)
    % Using a binning approach with proper SEM bands
    % valid_run_types is now a cell array of run types to include
    
    % Create bins for distance
    bin_edges = options.xlim(1):options.bin_width:options.xlim(2);
    bin_centers = bin_edges(1:end-1) + options.bin_width/2;
    num_bins = length(bin_centers);
    
    % Initialize arrays to collect dF/F values for each bin
    fast_bin_dff = cell(num_bins, 1);
    middle_bin_dff = cell(num_bins, 1);  % New array for middle runs
    slow_bin_dff = cell(num_bins, 1);
    
    for i = 1:num_bins
        fast_bin_dff{i} = [];
        middle_bin_dff{i} = [];  % Initialize middle bins
        slow_bin_dff{i} = [];
    end
    
    % Arrays to store all run data for individual plotting
    fast_runs_dist = {};
    fast_runs_dff = {};
    middle_runs_dist = {};  % New arrays for middle runs
    middle_runs_dff = {};
    slow_runs_dist = {};
    slow_runs_dff = {};
    
    % Process each run
    fast_count = 0;
    middle_count = 0;  % Count middle runs
    slow_count = 0;
    
    for m = 1:length(run_data)
        for r = 1:length(run_data(m).runs)
            run = run_data(m).runs(r);
            
            % Skip if not one of the valid run types
            if ~any(strcmp(run.type, valid_run_types))
                continue;
            end
            
            % Check if we have both distance and dF/F data
            if isempty(run.distance) || isempty(run.dff)
                continue;
            end
            
            % Determine which speed field to use (calculated_speed if available, otherwise speed)
            if isfield(run, 'calculated_speed') && ~isempty(run.calculated_speed)
                run_speed = run.calculated_speed;
            else
                run_speed = run.speed;
            end
            
            if isempty(run_speed)
                continue;
            end
            
            % Check for NaN or Inf values
            if any(isnan(run.distance)) || any(isnan(run.dff)) || any(isnan(run_speed)) || ...
               any(isinf(run.distance)) || any(isinf(run.dff)) || any(isinf(run_speed))
                continue;
            end
            
            % Calculate average speed for this run
            avg_speed = mean(run_speed);
            
            % Normalize dF/F
            dff = run.dff;
            switch options.normalize_method
                case 'zscore'
                    if std(dff) > 0
                        dff_norm = (dff - mean(dff)) / std(dff);
                    else
                        continue;
                    end
                case 'minmax'
                    range_dff = max(dff) - min(dff);
                    if range_dff > 0
                        dff_norm = (dff - min(dff)) / range_dff;
                    else
                        continue;
                    end
                case 'none'
                    dff_norm = dff;
                otherwise
                    dff_norm = dff;
            end
            
            % Classify runs based on separate thresholds for fast and slow
            if avg_speed >= fast_threshold
                % Fast run (greater than or equal to fast threshold)
                if options.show_individual_runs
                    fast_runs_dist{end+1} = run.distance;
                    fast_runs_dff{end+1} = dff_norm;
                end
                
                % Bin the data
                for i = 1:length(bin_edges)-1
                    bin_indices = run.distance >= bin_edges(i) & run.distance < bin_edges(i+1);
                    if any(bin_indices)
                        fast_bin_dff{i} = [fast_bin_dff{i}; dff_norm(bin_indices)];
                    end
                end
                
                fast_count = fast_count + 1;
            elseif avg_speed <= slow_threshold
                % Slow run (less than or equal to slow threshold)
                if options.show_individual_runs
                    slow_runs_dist{end+1} = run.distance;
                    slow_runs_dff{end+1} = dff_norm;
                end
                
                % Bin the data
                for i = 1:length(bin_edges)-1
                    bin_indices = run.distance >= bin_edges(i) & run.distance < bin_edges(i+1);
                    if any(bin_indices)
                        slow_bin_dff{i} = [slow_bin_dff{i}; dff_norm(bin_indices)];
                    end
                end
                
                slow_count = slow_count + 1;
            else
                % Middle speed run (between slow and fast thresholds)
                if options.show_individual_runs && options.show_middle_runs
                    middle_runs_dist{end+1} = run.distance;
                    middle_runs_dff{end+1} = dff_norm;
                end
                
                % Bin the data
                for i = 1:length(bin_edges)-1
                    bin_indices = run.distance >= bin_edges(i) & run.distance < bin_edges(i+1);
                    if any(bin_indices)
                        middle_bin_dff{i} = [middle_bin_dff{i}; dff_norm(bin_indices)];
                    end
                end
                
                middle_count = middle_count + 1;
            end
        end
    end
    
    % Determine run type label for printing
    if length(valid_run_types) == 1
        if strcmp(valid_run_types{1}, 'towards') || strcmp(valid_run_types{1}, 'not_food_towards')
            run_type_label = 'Towards';
        else
            run_type_label = 'Away';
        end
    else
        run_type_label = 'Combined';
    end
    
    % Print summary
    if fast_threshold == slow_threshold
        fprintf('%s runs: %d fast (>%.2f), %d slow (<%.2f)\n', ...
                run_type_label, fast_count, fast_threshold, slow_count, slow_threshold);
    else
        fprintf('%s runs: %d fast (≥%.2f), %d middle (%.2f-%.2f), %d slow (≤%.2f)\n', ...
                run_type_label, fast_count, fast_threshold, middle_count, slow_threshold, fast_threshold, slow_count, slow_threshold);
    end
    
    % Calculate mean and SEM for each bin
    fast_means = nan(num_bins, 1);
    fast_sems = nan(num_bins, 1);
    fast_counts = zeros(num_bins, 1);
    
    middle_means = nan(num_bins, 1);  % Calculate middle means and SEMs
    middle_sems = nan(num_bins, 1);
    middle_counts = zeros(num_bins, 1);
    
    slow_means = nan(num_bins, 1);
    slow_sems = nan(num_bins, 1);
    slow_counts = zeros(num_bins, 1);
    
    for i = 1:num_bins
        if length(fast_bin_dff{i}) >= options.min_points_per_bin
            fast_means(i) = mean(fast_bin_dff{i});
            fast_sems(i) = std(fast_bin_dff{i}) / sqrt(length(fast_bin_dff{i}));
            fast_counts(i) = length(fast_bin_dff{i});
        end
        
        if length(middle_bin_dff{i}) >= options.min_points_per_bin
            middle_means(i) = mean(middle_bin_dff{i});
            middle_sems(i) = std(middle_bin_dff{i}) / sqrt(length(middle_bin_dff{i}));
            middle_counts(i) = length(middle_bin_dff{i});
        end
        
        if length(slow_bin_dff{i}) >= options.min_points_per_bin
            slow_means(i) = mean(slow_bin_dff{i});
            slow_sems(i) = std(slow_bin_dff{i}) / sqrt(length(slow_bin_dff{i}));
            slow_counts(i) = length(slow_bin_dff{i});
        end
    end
    
    % Apply optional smoothing
    if options.smoothing_window > 0
        valid_fast = ~isnan(fast_means);
        if sum(valid_fast) > options.smoothing_window
            fast_means(valid_fast) = movmean(fast_means(valid_fast), options.smoothing_window);
            fast_sems(valid_fast) = movmean(fast_sems(valid_fast), options.smoothing_window);
        end
        
        valid_middle = ~isnan(middle_means);  % Smooth middle data
        if sum(valid_middle) > options.smoothing_window
            middle_means(valid_middle) = movmean(middle_means(valid_middle), options.smoothing_window);
            middle_sems(valid_middle) = movmean(middle_sems(valid_middle), options.smoothing_window);
        end
        
        valid_slow = ~isnan(slow_means);
        if sum(valid_slow) > options.smoothing_window
            slow_means(valid_slow) = movmean(slow_means(valid_slow), options.smoothing_window);
            slow_sems(valid_slow) = movmean(slow_sems(valid_slow), options.smoothing_window);
        end
    end
    
    % Plot data
    hold on;
    
    % Plot individual runs if requested
    if options.show_individual_runs
        for i = 1:length(fast_runs_dist)
            h = plot(fast_runs_dist{i}, fast_runs_dff{i}, 'Color', [0.2, 0.6, 0.8, 0.1], 'LineWidth', 0.5);
            set(h, 'HandleVisibility', 'off');
        end
        
        if options.show_middle_runs
            for i = 1:length(middle_runs_dist)  % Plot middle individual runs
                h = plot(middle_runs_dist{i}, middle_runs_dff{i}, 'Color', [0.5, 0.5, 0.5, 0.1], 'LineWidth', 0.5);
                set(h, 'HandleVisibility', 'off');
            end
        end
        
        for i = 1:length(slow_runs_dist)
            h = plot(slow_runs_dist{i}, slow_runs_dff{i}, 'Color', [0.8, 0.4, 0.2, 0.1], 'LineWidth', 0.5);
            set(h, 'HandleVisibility', 'off');
        end
    end
    
    % Only include bins with valid data
    valid_fast = ~isnan(fast_means) & fast_counts >= options.min_points_per_bin;
    valid_middle = ~isnan(middle_means) & middle_counts >= options.min_points_per_bin;
    valid_slow = ~isnan(slow_means) & slow_counts >= options.min_points_per_bin;
    
    % Plot middle run means with shaded SEM if enabled
    if options.show_middle_runs && any(valid_middle)
        % Extract valid data
        valid_x = bin_centers(valid_middle);
        valid_means = middle_means(valid_middle);
        valid_sems = middle_sems(valid_middle);
        
        % Plot the mean line
        h_middle = plot(valid_x, valid_means, '-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, ...
                       'DisplayName', sprintf('Middle (%.2f-%.2f cm/s)', slow_threshold, fast_threshold));
        
        % Create SEM shading
        lower_bound = valid_means - valid_sems;
        upper_bound = valid_means + valid_sems;
        
        % Create shaded area using patch
        x_patch = [valid_x, fliplr(valid_x)];
        y_patch = [upper_bound', fliplr(lower_bound')];
        h_patch = patch(x_patch, y_patch, [0.5, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        set(h_patch, 'HandleVisibility', 'off');  % Hide from legend
    end
    
    % Plot slow run means with shaded SEM
    if any(valid_slow)
        % Extract valid data
        valid_x = bin_centers(valid_slow);
        valid_means = slow_means(valid_slow);
        valid_sems = slow_sems(valid_slow);
        
        % Plot the mean line
        h_slow = plot(valid_x, valid_means, '-', 'Color', [0.8, 0.4, 0.2], 'LineWidth', 2, ...
                     'DisplayName', sprintf('Slow (≤%.2f cm/s)', slow_threshold));
        
        % Create SEM shading
        lower_bound = valid_means - valid_sems;
        upper_bound = valid_means + valid_sems;
        
        % Create shaded area using patch
        x_patch = [valid_x, fliplr(valid_x)];
        y_patch = [upper_bound', fliplr(lower_bound')];
        h_patch = patch(x_patch, y_patch, [0.8, 0.4, 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        set(h_patch, 'HandleVisibility', 'off');  % Hide from legend
    end
    
    % Plot fast run means with shaded SEM
    if any(valid_fast)
        % Extract valid data
        valid_x = bin_centers(valid_fast);
        valid_means = fast_means(valid_fast);
        valid_sems = fast_sems(valid_fast);
        
        % Plot the mean line
        h_fast = plot(valid_x, valid_means, '-', 'Color', [0.2, 0.6, 0.8], 'LineWidth', 2, ...
                     'DisplayName', sprintf('Fast (≥%.2f cm/s)', fast_threshold));
        
        % Create SEM shading
        lower_bound = valid_means - valid_sems;
        upper_bound = valid_means + valid_sems;
        
        % Create shaded area using patch
        x_patch = [valid_x, fliplr(valid_x)];
        y_patch = [upper_bound', fliplr(lower_bound')];
        h_patch = patch(x_patch, y_patch, [0.2, 0.6, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        set(h_patch, 'HandleVisibility', 'off');  % Hide from legend
    end
    
    % Add reference line at y=0
    h_ref = plot(get(gca, 'XLim'), [0 0], 'k--', 'LineWidth', 1);
    set(h_ref, 'HandleVisibility', 'off');
    
    % Set axis limits and labels
    xlim(options.xlim);
    ylim(options.ylim);
    xlabel('Distance from Food (cm)');
    ylabel('Normalized dF/F');
    
    % Add counts
    text(0.05, 0.95, sprintf('Fast: n=%d runs', fast_count), ...
         'Units', 'normalized', 'Color', [0.2, 0.6, 0.8], ...
         'FontWeight', 'bold', 'FontSize', 8);
    text(0.05, 0.9, sprintf('Slow: n=%d runs', slow_count), ...
         'Units', 'normalized', 'Color', [0.8, 0.4, 0.2], ...
         'FontWeight', 'bold', 'FontSize', 8);
         
    % Show middle run count (whether plotting them or not)
    if fast_threshold ~= slow_threshold
        if options.show_middle_runs
            text(0.05, 0.85, sprintf('Middle: n=%d runs', middle_count), ...
                 'Units', 'normalized', 'Color', [0.5, 0.5, 0.5], ...
                 'FontWeight', 'bold', 'FontSize', 8);
        else
            text(0.05, 0.85, sprintf('Middle: n=%d runs (not shown)', middle_count), ...
                 'Units', 'normalized', 'Color', [0.5, 0.5, 0.5], ...
                 'FontSize', 8);
        end
    end
    
    % Add legend
    legend('show', 'Location', 'best');
    legend('boxoff');
    
    % Add grid
    grid on;
    box off;
    
    % Add label for empty plot if needed
    if fast_count == 0 && slow_count == 0 && (middle_count == 0 || ~options.show_middle_runs)
        text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
end