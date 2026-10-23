function [run_data] = analyze_cheeseboard_runs(mice, options)
    % Analyze mouse runs to and from food hole in cheeseboard maze
    %
    % Inputs:
    %   mice - Cell array containing cheeseboard mouse data
    %   options - Structure with analysis parameters:
    %     .time_limit - Time limit in minutes (optional, analyzes full session if empty)
    %     .food_area - Distance threshold for food area in cm (default: 2)
    %     .group_filter - Group to analyze ('saline', 'CNO', or 'all') (default: 'all')
    %     .memory - Memory strength to analyze ('strong', 'weak', or 'all') (default: 'all')
    %               Requires memory classifications in column 7 of mice array
    %     .session_filter - Session to analyze ('pre', 'test', or 'both') (default: 'both')
    %     .exclude_grooming - Remove grooming periods (default: true)
    %     .validate_plots - Boolean to control validation plotting for all mice (default: false)
    %     .plot_specific_mouse - Mouse name to show validation plots for (optional)
    %                            Example: 'MDRE8_saline' or 'F13_CNO'
    %                            If specified, validation plots will be shown for this mouse
    %                            regardless of the validate_plots setting
    %     .threshold - Velocity threshold for run detection (default: 0.1)
    %     .colormap - Colormap for dF/F visualization (default: 'parula')
    %                 Standard MATLAB: 'parula', 'jet', 'hot', 'cool', 'viridis', 'plasma', 'inferno'
    %                 Custom options: 'bluewhitered', 'blueWhiteRed', 'redblue'
    %     .caxis - Color axis limits for dF/F z-score plots (default: [] for auto)
    %              Example: [-2, 2] to set limits from -2 to +2 z-scores
    %              If empty, uses 5th-95th percentile of z-scored dF/F
    %
    % Outputs:
    %   run_data - Structure containing run information with dF/F vs distance slopes
    %              Each run includes both raw dF/F and z-scored dF/F
    %              Slopes are calculated using z-scored dF/F for standardization
    
    % Default parameters
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'time_limit')
        options.time_limit = []; % No time limit by default
    end
    
    if ~isfield(options, 'food_area')
        options.food_area = 2; % 2 cm default
    end
    
    if ~isfield(options, 'group_filter')
        options.group_filter = 'all';
    end
    
    if ~isfield(options, 'memory')
        options.memory = 'all'; % All memory types by default
    end
    
    if ~isfield(options, 'session_filter')
        options.session_filter = 'both';
    end
    
    if ~isfield(options, 'exclude_grooming')
        options.exclude_grooming = true;
    end
    
    if ~isfield(options, 'validate_plots')
        options.validate_plots = false;
    end
    
    if ~isfield(options, 'plot_specific_mouse')
        options.plot_specific_mouse = ''; % No specific mouse by default
    end
    
    if ~isfield(options, 'threshold')
        options.threshold = 0.1; % Distance derivative threshold
    end
    
    if ~isfield(options, 'colormap')
        options.colormap = 'parula'; % Default colormap
    end
    
    if ~isfield(options, 'caxis')
        options.caxis = []; % Auto-scale by default (empty means use percentiles)
    end
    
    % Convert time limit to seconds if specified
    if ~isempty(options.time_limit)
        time_limit_sec = options.time_limit * 60;
        use_time_limit = true;
    else
        use_time_limit = false;
    end
    
    fprintf('=== CHEESEBOARD RUN ANALYSIS ===\n');
    fprintf('Food area threshold: %.1f cm\n', options.food_area);
    fprintf('Group filter: %s\n', options.group_filter);
    fprintf('Memory filter: %s\n', options.memory);
    fprintf('Session filter: %s\n', options.session_filter);
    fprintf('Exclude grooming: %s\n', mat2str(options.exclude_grooming));
    if use_time_limit
        fprintf('Time limit: %.1f minutes\n', options.time_limit);
    else
        fprintf('Time limit: Full session\n');
    end
    fprintf('Distance derivative threshold: %.3f\n', options.threshold);
    fprintf('Colormap for validation: %s\n', options.colormap);
    if ~isempty(options.caxis)
        fprintf('Color axis limits: [%.2f, %.2f]\n', options.caxis(1), options.caxis(2));
    else
        fprintf('Color axis limits: Auto (5th-95th percentile)\n');
    end
    fprintf('Validate plots for all mice: %s\n', mat2str(options.validate_plots));
    if ~isempty(options.plot_specific_mouse)
        fprintf('Show validation plots for specific mouse: %s\n', options.plot_specific_mouse);
    end
    
    % Initialize output structure
    run_data = struct('mouse_id', {}, 'session', {}, 'runs', {}, 'food_position', {}, 'memory_strength', {});
    
    % Get mouse information
    mouse_names = mice(:,1);
    n_mice = size(mice, 1);
    
    % Check if memory classification exists (column 7)
    has_memory_classification = size(mice, 2) >= 7 && ~all(cellfun(@isempty, mice(:,7)));
    if has_memory_classification
        memory_classifications = mice(:,7);
        fprintf('Memory classifications found in column 7\n');
    else
        if ~strcmp(options.memory, 'all')
            error('Memory classifications not found in mice array (column 7), but memory filter "%s" was specified', options.memory);
        end
        fprintf('No memory classifications found - ignoring memory filter\n');
    end
    
    % Filter mice by group and memory
    selected_indices = [];
    for i = 1:n_mice
        mouse_name = mouse_names{i};
        
        % Apply group filter
        group_passed = false;
        if strcmp(options.group_filter, 'all')
            group_passed = true;
        elseif strcmp(options.group_filter, 'saline') && contains(mouse_name, 'saline')
            group_passed = true;
        elseif strcmp(options.group_filter, 'CNO') && contains(mouse_name, 'CNO')
            group_passed = true;
        end
        
        % Apply memory filter
        memory_passed = false;
        if strcmp(options.memory, 'all')
            memory_passed = true;
        elseif has_memory_classification
            mouse_memory = memory_classifications{i};
            if strcmp(options.memory, 'strong') && strcmp(mouse_memory, 'strong memory')
                memory_passed = true;
            elseif strcmp(options.memory, 'weak') && strcmp(mouse_memory, 'weak memory')
                memory_passed = true;
            end
        end
        
        % Include mouse if both filters pass
        if group_passed && memory_passed
            selected_indices = [selected_indices; i];
        end
    end
    
    fprintf('Analyzing %d mice\n', length(selected_indices));
    
    % Print breakdown of selected mice
    if has_memory_classification
        saline_strong = 0; saline_weak = 0; cno_strong = 0; cno_weak = 0;
        for i = 1:length(selected_indices)
            idx = selected_indices(i);
            mouse_name = mouse_names{idx};
            memory_type = memory_classifications{idx};
            
            if contains(mouse_name, 'saline')
                if strcmp(memory_type, 'strong memory')
                    saline_strong = saline_strong + 1;
                else
                    saline_weak = saline_weak + 1;
                end
            elseif contains(mouse_name, 'CNO')
                if strcmp(memory_type, 'strong memory')
                    cno_strong = cno_strong + 1;
                else
                    cno_weak = cno_weak + 1;
                end
            end
        end
        fprintf('  - Saline strong memory: %d mice\n', saline_strong);
        fprintf('  - Saline weak memory: %d mice\n', saline_weak);
        fprintf('  - CNO strong memory: %d mice\n', cno_strong);
        fprintf('  - CNO weak memory: %d mice\n', cno_weak);
    end
    
    % Check if specific mouse exists in the dataset
    specific_mouse_found = false;
    if ~isempty(options.plot_specific_mouse)
        % Check if the specified mouse exists in the selected mice
        for i = 1:length(selected_indices)
            mouse_idx = selected_indices(i);
            mouse_name = mouse_names{mouse_idx};
            if contains(mouse_name, options.plot_specific_mouse) || strcmp(mouse_name, options.plot_specific_mouse)
                specific_mouse_found = true;
                fprintf('Found specified mouse for validation plots: %s\n', mouse_name);
                break;
            end
        end
        
        if ~specific_mouse_found
            fprintf('Warning: Specified mouse "%s" not found in the dataset\n', options.plot_specific_mouse);
            fprintf('Available mice:\n');
            for i = 1:length(selected_indices)
                mouse_idx = selected_indices(i);
                mouse_name = mouse_names{mouse_idx};
                if has_memory_classification
                    memory_type = memory_classifications{mouse_idx};
                    fprintf('  %s (%s)\n', mouse_name, memory_type);
                else
                    fprintf('  %s\n', mouse_name);
                end
            end
        end
    end
    
    % Process each selected mouse
    for idx = 1:length(selected_indices)
        mouse_idx = selected_indices(idx);
        mouse_name = mouse_names{mouse_idx};
        food_pos = mice{mouse_idx, 2};
        
        % Get memory strength if available
        memory_strength = '';
        if has_memory_classification && ~isempty(mice{mouse_idx, 7})
            memory_strength = mice{mouse_idx, 7};
        end
        
        fprintf('\nProcessing mouse %d/%d: %s', idx, length(selected_indices), mouse_name);
        if ~isempty(memory_strength)
            fprintf(' (%s)', memory_strength);
        end
        fprintf('\n');
        fprintf('  Food position: [%.1f, %.1f] cm\n', food_pos(1), food_pos(2));
        
        % Determine if we should show validation plots for this mouse
        show_plots_for_this_mouse = false;
        
        if options.validate_plots
            % Show plots for all mice if validate_plots is true
            show_plots_for_this_mouse = true;
        elseif ~isempty(options.plot_specific_mouse)
            % Show plots only for the specified mouse
            if contains(mouse_name, options.plot_specific_mouse) || strcmp(mouse_name, options.plot_specific_mouse)
                show_plots_for_this_mouse = true;
                fprintf('  -> This mouse selected for validation plots\n');
            end
        end
        
        % Determine which sessions to analyze
        sessions_to_analyze = {};
        session_data = {};
        
        if strcmp(options.session_filter, 'pre') || strcmp(options.session_filter, 'both')
            sessions_to_analyze{end+1} = 'pre';
            session_data{end+1} = mice{mouse_idx, 3};
        end
        
        if strcmp(options.session_filter, 'test') || strcmp(options.session_filter, 'both')
            sessions_to_analyze{end+1} = 'test';
            session_data{end+1} = mice{mouse_idx, 4};
        end
        
        % Process each session
        for sess_idx = 1:length(sessions_to_analyze)
            session_name = sessions_to_analyze{sess_idx};
            data = session_data{sess_idx};
            
            if isempty(data) || size(data, 1) < 10
                fprintf('  %s session: insufficient data\n', session_name);
                continue;
            end
            
            fprintf('  Processing %s session (%d data points)\n', session_name, size(data, 1));
            
            % Apply time limit if specified
            if use_time_limit
                time_mask = (data(:, 1) - data(1, 1)) <= time_limit_sec;
                data = data(time_mask, :);
                fprintf('    After time limit: %d data points\n', size(data, 1));
            end
            
            % Remove grooming periods if requested
            if options.exclude_grooming
                non_grooming_mask = data(:, 10) == 0;
                data = data(non_grooming_mask, :);
                fprintf('    After grooming removal: %d data points\n', size(data, 1));
            end
            
            % Check if we have enough data for filtering
            if size(data, 1) < 10
                fprintf('    Insufficient data after filtering\n');
                continue;
            end
            
            % Apply smoothing filters
            data_filtered = apply_smoothing_filters(data);
            
            % Calculate z-scored dF/F for the session
            dff_mean = mean(data_filtered(:, 6), 'omitnan');
            dff_std = std(data_filtered(:, 6), 'omitnan');
            
            % Add z-scored dF/F as column 11
            if dff_std > 0
                data_filtered(:, 11) = (data_filtered(:, 6) - dff_mean) / dff_std;
            else
                data_filtered(:, 11) = zeros(size(data_filtered, 1), 1);
            end
            
            fprintf('    dF/F statistics: mean=%.3f, std=%.3f\n', dff_mean, dff_std);
            
            % Calculate distance derivative
            distance_derivative = calculate_distance_derivative(data_filtered);
            
            % Detect runs
            detected_runs = detect_runs(data_filtered, distance_derivative, options);
            
            if ~isempty(detected_runs)
                % Create mouse entry
                mouse_entry = struct();
                mouse_entry.mouse_id = mouse_name;
                mouse_entry.session = session_name;
                mouse_entry.runs = detected_runs;
                mouse_entry.food_position = food_pos;
                mouse_entry.memory_strength = memory_strength;
                
                run_data = [run_data; mouse_entry];
                
                fprintf('    Found %d valid runs (%d towards, %d away)\n', ...
                        length(detected_runs), ...
                        sum(strcmp({detected_runs.type}, 'towards')), ...
                        sum(strcmp({detected_runs.type}, 'away')));
                
                % Generate validation plots if requested for this mouse
                if show_plots_for_this_mouse
                    fprintf('    -> Generating validation plots for %s %s session\n', mouse_name, session_name);
                    plot_cheeseboard_run_validation(data_filtered, distance_derivative, detected_runs, options, mouse_name, session_name);
                end
            else
                fprintf('    No valid runs detected\n');
            end
        end
    end
    
    % Print summary
    print_run_summary(run_data, options);
end

function data_filtered = apply_smoothing_filters(data)
    % Apply smoothing filters to the data
    
    data_filtered = data;
    
    % Check if we have enough data points for filtering
    if size(data, 1) < 4
        return;
    end
    
    % Design filters
    [b1, a1] = butter(1, 0.1, 'low');  % Low-pass filter for signals
    
    try
        % Filter signals
        data_filtered(:, 4) = filtfilt(b1, a1, data(:, 4));  % 465nm signal
        data_filtered(:, 5) = filtfilt(b1, a1, data(:, 5));  % 405nm signal
        data_filtered(:, 6) = filtfilt(b1, a1, data(:, 6));  % dF/F
        data_filtered(:, 7) = filtfilt(b1, a1, data(:, 7));  % Speed
        data_filtered(:, 9) = filtfilt(b1, a1, data(:, 9));  % Distance
    catch
        % If filtering fails, use original data
        warning('Filtering failed, using original data');
    end
end

function distance_derivative = calculate_distance_derivative(data)
    % Calculate the derivative of distance to food
    
    distance = data(:, 9);  % Distance column
    
    % Calculate derivative
    derivative = [0; diff(distance)];  % Prepend zero to maintain size
    
    % Apply additional smoothing to derivative if we have enough points
    if length(derivative) >= 4
        try
            [b, a] = butter(1, 0.5, 'low');
            distance_derivative = filtfilt(b, a, derivative);
        catch
            distance_derivative = derivative;
        end
    else
        distance_derivative = derivative;
    end
end

function runs = detect_runs(data, distance_derivative, options)
    % Detect runs towards and away from food
    
    runs = [];
    run_counter = 1;
    
    % Parameters
    threshold = options.threshold;
    food_area = options.food_area;
    
    % Track run state
    in_towards_run = false;
    in_away_run = false;
    current_run_start = 0;
    
    ii = 2;
    while ii <= size(data, 1) - 1
        if ii > length(distance_derivative)
            break;
        end
        
        % Check for TOWARDS run (distance decreasing, derivative negative)
        if ~in_towards_run && distance_derivative(ii) < -threshold && data(ii, 9) <= food_area
            % Start of a new towards run
            in_towards_run = true;
            current_run_start = find(distance_derivative(1:ii) > -threshold, 1, 'last') + 1;
            if isempty(current_run_start)
                current_run_start = 1;
            end
            
        elseif in_towards_run
            % Already in a towards run, check if it's ending
            if distance_derivative(ii) >= -threshold || data(ii, 9) > food_area
                % End of towards run
                in_towards_run = false;
                % Only save if mouse got close to food
                if min(data(current_run_start:ii, 9)) <= food_area
                    run_indices = (current_run_start:ii)';
                    
                    % Create run entry
                    run_entry = create_run_entry(run_counter, 'towards', run_indices, data, food_area);
                    
                    if ~isempty(run_entry)
                        runs = [runs; run_entry];
                        run_counter = run_counter + 1;
                    end
                end
            end
        end
        
        % Check for AWAY run (distance increasing, derivative positive)
        if ~in_away_run && distance_derivative(ii) > threshold && data(ii, 9) <= food_area
            % Start of a new away run
            in_away_run = true;
            current_run_start = find(distance_derivative(1:ii) < threshold, 1, 'last') + 1;
            if isempty(current_run_start)
                current_run_start = 1;
            end
            
        elseif in_away_run
            % Already in an away run, check if it's ending
            if distance_derivative(ii) <= threshold
                % End of away run
                in_away_run = false;
                % Only save if mouse started close to food
                if any(data(current_run_start:current_run_start+min(10,ii-current_run_start), 9) <= food_area)
                    run_indices = (current_run_start:ii)';
                    
                    % Create run entry
                    run_entry = create_run_entry(run_counter, 'away', run_indices, data, food_area);
                    
                    if ~isempty(run_entry)
                        runs = [runs; run_entry];
                        run_counter = run_counter + 1;
                    end
                end
            end
        end
        
        ii = ii + 1;
    end
    
    % Handle any ongoing runs at the end of the data
    if in_towards_run && min(data(current_run_start:end, 9)) <= food_area
        run_indices = (current_run_start:size(data, 1))';
        run_entry = create_run_entry(run_counter, 'towards', run_indices, data, food_area);
        if ~isempty(run_entry)
            runs = [runs; run_entry];
        end
    end
    
    if in_away_run && any(data(current_run_start:current_run_start+min(10,size(data,1)-current_run_start), 9) <= food_area)
        run_indices = (current_run_start:size(data, 1))';
        run_entry = create_run_entry(run_counter, 'away', run_indices, data, food_area);
        if ~isempty(run_entry)
            runs = [runs; run_entry];
        end
    end
end

function run_entry = create_run_entry(run_id, run_type, run_indices, data, food_area)
    % Create a run entry structure
    
    % Minimum run duration (in seconds)
    min_duration = 0.5;
    
    % Check if run is long enough
    duration = data(run_indices(end), 1) - data(run_indices(1), 1);
    if duration < min_duration
        run_entry = [];
        return;
    end
    
    % Extract run data
    distance_vec = data(run_indices, 9);
    dff_vec = data(run_indices, 6);
    dff_z_vec = data(run_indices, 11);  % Z-scored dF/F
    
    % Calculate dF/F vs distance slope (using z-scored dF/F)
    dff_slope = calculate_dff_distance_slope(distance_vec, dff_z_vec, run_type);
    
    % Calculate speed from x,y coordinates and compare with original
    [calculated_speed, speed_correlation] = calculate_speed_from_position(...
        data(run_indices, 1),...  % time
        data(run_indices, 2),...  % x
        data(run_indices, 3),...  % y
        data(run_indices, 7));    % original speed
    
    % Create run structure
    run_entry = struct();
    run_entry.id = run_id;
    run_entry.type = run_type;
    run_entry.indices = run_indices;
    run_entry.time = data(run_indices, 1);
    run_entry.x_coordinate = data(run_indices, 2);
    run_entry.y_coordinate = data(run_indices, 3);
    run_entry.distance = distance_vec;
    run_entry.dff = dff_vec;
    run_entry.dff_z = dff_z_vec;  % Z-scored dF/F
    run_entry.dff_slope = dff_slope;
    run_entry.speed = data(run_indices, 7);  % Original speed
    run_entry.calculated_speed = calculated_speed;  % Calculated speed
    run_entry.speed_correlation = speed_correlation;  % Correlation between speeds
    run_entry.signal_465 = data(run_indices, 4);
    run_entry.signal_405 = data(run_indices, 5);
    run_entry.zones = data(run_indices, 8);
    run_entry.grooming = data(run_indices, 10);
    run_entry.start_time = data(run_indices(1), 1);
    run_entry.end_time = data(run_indices(end), 1);
    run_entry.duration = duration;
    run_entry.min_distance = min(distance_vec);
    run_entry.max_distance = max(distance_vec);
    run_entry.distance_change = distance_vec(end) - distance_vec(1);
end

function dff_slope = calculate_dff_distance_slope(distance_vec, dff_vec, run_type)
    % Calculate the slope of dF/F relative to distance from food
    
    % Check if we have enough data points
    if length(distance_vec) < 3 || length(dff_vec) < 3
        dff_slope = NaN;
        return;
    end
    
    % Remove any NaN or infinite values
    valid_idx = ~isnan(distance_vec) & ~isnan(dff_vec) & isfinite(distance_vec) & isfinite(dff_vec);
    
    if sum(valid_idx) < 3
        dff_slope = NaN;
        return;
    end
    
    distance_clean = distance_vec(valid_idx);
    dff_clean = dff_vec(valid_idx);
    
    % Transform distance based on run type
    if strcmp(run_type, 'towards')
        % For towards runs: use negative distance so that approaching food
        % becomes increasing x-axis values. Positive slope = dF/F increases approaching food
        x_values = -distance_clean;
    elseif strcmp(run_type, 'away')
        % For away runs: use distance as-is. Negative slope = dF/F decreases moving away from food
        x_values = distance_clean;
    else
        error('run_type must be either "towards" or "away"');
    end
    
    % Perform linear regression: dF/F = slope * x_values + intercept
    try
        p = polyfit(x_values, dff_clean, 1);
        dff_slope = p(1); % The slope coefficient
    catch
        % If polyfit fails, use manual calculation
        try
            n = length(x_values);
            sum_x = sum(x_values);
            sum_dff = sum(dff_clean);
            sum_x2 = sum(x_values.^2);
            sum_x_dff = sum(x_values .* dff_clean);
            
            % Calculate slope using least squares formula
            denominator = n * sum_x2 - sum_x^2;
            if abs(denominator) < 1e-12
                dff_slope = NaN;
            else
                dff_slope = (n * sum_x_dff - sum_x * sum_dff) / denominator;
            end
        catch
            dff_slope = NaN;
        end
    end
end

function plot_cheeseboard_run_validation(data, distance_derivative, runs, options, mouse_name, session_name)
    % Create validation plot for run detection with dF/F scatter overlay
    
    figure('Position', [100, 100, 1400, 800], 'Name', sprintf('Run Validation: %s - %s', mouse_name, session_name));
    
    % Set the chosen colormap with error handling
    try
        if strcmp(options.colormap, 'bluewhitered') || strcmp(options.colormap, 'blueWhiteRed')
            % Create custom blue-white-red colormap
            cmap = create_blue_white_red_colormap();
            colormap(cmap);
        elseif strcmp(options.colormap, 'redblue')
            % Create custom red-blue colormap
            cmap = create_red_blue_colormap();
            colormap(cmap);
        else
            % Use standard MATLAB colormap
            colormap(options.colormap);
        end
    catch
        % Fallback to parula if colormap fails
        warning('Colormap "%s" not found, using parula instead', options.colormap);
        colormap('parula');
    end
    
    % Top subplot: Distance and derivative
    subplot(2, 2, [1, 2]);
    
    % Plot distance to food
    plot(data(:, 1), data(:, 9), 'k-', 'LineWidth', 1.5);
    hold on;
    
    % Plot food area threshold
    yline(options.food_area, 'k--', 'Food Area Threshold', 'LineWidth', 1);
    
    % Highlight runs with different colors for towards/away
    for i = 1:length(runs)
        run_indices = runs(i).indices;
        if strcmp(runs(i).type, 'towards')
            scatter(data(run_indices, 1), data(run_indices, 9), 15, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        else
            scatter(data(run_indices, 1), data(run_indices, 9), 15, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        end
    end
    
    title('Distance to Food and Run Detection');
    xlabel('Time (s)');
    ylabel('Distance (cm)');
    grid off;
    box off;
    
    % Calculate global dF/F z-score range for consistent color scaling
    all_dff_z_in_runs = [];
    for i = 1:length(runs)
        all_dff_z_in_runs = [all_dff_z_in_runs; runs(i).dff_z];
    end
    
    % Determine color axis limits
    if ~isempty(options.caxis)
        % Use user-specified caxis
        dff_min = options.caxis(1);
        dff_max = options.caxis(2);
    elseif ~isempty(all_dff_z_in_runs)
        % Auto-scale using percentiles
        dff_min = prctile(all_dff_z_in_runs, 5);   % 5th percentile
        dff_max = prctile(all_dff_z_in_runs, 95);  % 95th percentile
    else
        % Default range for z-scores
        dff_min = -2;
        dff_max = 2;
    end
    
    % Bottom left subplot: Towards runs with dF/F scatter
    subplot(2, 2, 3);
    
    % Plot full trajectory in light gray
    plot(data(:, 2), data(:, 3), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    hold on;
    
    
    % Plot towards runs with dF/F coloring
    towards_count = 0;
    for i = 1:length(runs)
        if strcmp(runs(i).type, 'towards')
            run_indices = runs(i).indices;
            towards_count = towards_count + 1;
            
            % Plot trajectory as thin black line
            plot(data(run_indices, 2), data(run_indices, 3), 'k-', 'LineWidth', 1);
            
            % Scatter z-scored dF/F values
            scatter(data(run_indices, 2), data(run_indices, 3), 10, runs(i).dff_z, 'filled', 'MarkerFaceAlpha', 0.8);
            
            % % Add run number at start of run
            % text(data(run_indices(1), 2), data(run_indices(1), 3), num2str(towards_count), ...
            %     'FontSize', 5, 'FontWeight', 'bold', 'Color', 'black', ...
            %     'BackgroundColor', 'white', 'EdgeColor', 'white', 'Margin', 1);
        end
    end
    
    % Set consistent color limits
    if ~isempty(all_dff_z_in_runs) || ~isempty(options.caxis)
        caxis([dff_min, dff_max]);
    end
    % Plot food location
    [~, food_idx] = min(data(:, 9));
    plot(data(food_idx, 2), data(food_idx, 3), 'ko', 'MarkerSize', 7, 'LineWidth', 1, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k');
    title(sprintf('Towards Runs (n=%d) - Z-scored dF/F', towards_count));
    xlabel('X Position');
    ylabel('Y Position');
    axis equal;
    axis off;
    
    % Add colorbar
    c1 = colorbar('Location', 'eastoutside');
    c1.Label.String = 'dF/F (z-score)';
    c1.Label.FontSize = 10;
    
    % Bottom right subplot: Away runs with dF/F scatter
    subplot(2, 2, 4);
    
    % Plot full trajectory in light gray
    plot(data(:, 2), data(:, 3), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    hold on;
    
    % Plot away runs with dF/F coloring
    away_count = 0;
    for i = 1:length(runs)
        if strcmp(runs(i).type, 'away')
            run_indices = runs(i).indices;
            away_count = away_count + 1;
            
            % Plot trajectory as thin black line
            plot(data(run_indices, 2), data(run_indices, 3), 'k-', 'LineWidth', 1);
            
            % Scatter z-scored dF/F values
            scatter(data(run_indices, 2), data(run_indices, 3), 10, runs(i).dff_z, 'filled', 'MarkerFaceAlpha', 0.8);
            
            % % Add run number at start of run
            %  text(data(run_indices(1), 2), data(run_indices(1), 3), num2str(away_count), ...
            %     'FontSize', 5, 'FontWeight', 'bold', 'Color', 'black', ...
            %     'BackgroundColor', 'white', 'EdgeColor', 'white', 'Margin', 1);
        end
    end
    
    % Set consistent color limits
    if ~isempty(all_dff_z_in_runs) || ~isempty(options.caxis)
        caxis([dff_min, dff_max]);
    end
    % Plot food location
    plot(data(food_idx, 2), data(food_idx, 3), 'ko', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'none');
    title(sprintf('Away Runs (n=%d) - Z-scored dF/F', away_count));
    xlabel('X Position');
    ylabel('Y Position');
    axis equal;
    axis off;
    
    % Add colorbar
    c2 = colorbar('Location', 'eastoutside');
    c2.Label.String = 'dF/F (z-score)';
    c2.Label.FontSize = 10;
    
    % Add overall title with run summary
    sgtitle(sprintf('%s - %s Session: %d towards, %d away runs detected (Colormap: %s)', ...
            mouse_name, session_name, towards_count, away_count, options.colormap));
end

function print_run_summary(run_data, options)
    % Print summary statistics
    
    if isempty(run_data)
        fprintf('\nNo valid runs found with current parameters\n');
        return;
    end
    
    % Count runs by type and session
    total_runs = 0;
    towards_runs = 0;
    away_runs = 0;
    pre_runs = 0;
    test_runs = 0;
    all_slopes = [];
    all_durations = [];
    all_speed_correlations = [];
    
    for i = 1:length(run_data)
        mouse_runs = run_data(i).runs;
        total_runs = total_runs + length(mouse_runs);
        
        if strcmp(run_data(i).session, 'pre')
            pre_runs = pre_runs + length(mouse_runs);
        else
            test_runs = test_runs + length(mouse_runs);
        end
        
        for j = 1:length(mouse_runs)
            if strcmp(mouse_runs(j).type, 'towards')
                towards_runs = towards_runs + 1;
            else
                away_runs = away_runs + 1;
            end
            
            if ~isnan(mouse_runs(j).dff_slope)
                all_slopes = [all_slopes; mouse_runs(j).dff_slope];
            end
            
            if ~isnan(mouse_runs(j).speed_correlation)
                all_speed_correlations = [all_speed_correlations; mouse_runs(j).speed_correlation];
            end
            
            all_durations = [all_durations; mouse_runs(j).duration];
        end
    end
    
    fprintf('\n=== RUN ANALYSIS SUMMARY ===\n');
    fprintf('Total mice analyzed: %d\n', length(run_data));
    fprintf('Total runs detected: %d\n', total_runs);
    fprintf('  Towards runs: %d (%.1f%%)\n', towards_runs, towards_runs/total_runs*100);
    fprintf('  Away runs: %d (%.1f%%)\n', away_runs, away_runs/total_runs*100);
    
    if strcmp(options.session_filter, 'both')
        fprintf('  Pre-test runs: %d\n', pre_runs);
        fprintf('  Test runs: %d\n', test_runs);
    end
    
    fprintf('\nRun duration statistics:\n');
    fprintf('  Mean: %.2f ± %.2f seconds\n', mean(all_durations), std(all_durations));
    fprintf('  Range: %.2f to %.2f seconds\n', min(all_durations), max(all_durations));
    
    if ~isempty(all_slopes)
        fprintf('\ndF/F vs distance slope statistics:\n');
        fprintf('  Mean: %.6f ± %.6f per distance unit\n', mean(all_slopes), std(all_slopes));
        fprintf('  Range: %.6f to %.6f per distance unit\n', min(all_slopes), max(all_slopes));
        fprintf('  Valid slopes: %d/%d runs (%.1f%%)\n', length(all_slopes), total_runs, length(all_slopes)/total_runs*100);
    end
    
    if ~isempty(all_speed_correlations)
        fprintf('\nSpeed calculation validation:\n');
        fprintf('  Mean correlation (calculated vs original): %.4f ± %.4f\n', mean(all_speed_correlations), std(all_speed_correlations));
        fprintf('  Range: %.4f to %.4f\n', min(all_speed_correlations), max(all_speed_correlations));
        fprintf('  Valid correlations: %d/%d runs (%.1f%%)\n', length(all_speed_correlations), total_runs, length(all_speed_correlations)/total_runs*100);
    end
    
    fprintf('\nAnalysis parameters used:\n');
    fprintf('  Food area threshold: %.1f cm\n', options.food_area);
    fprintf('  Distance derivative threshold: %.3f\n', options.threshold);
    fprintf('  Group filter: %s\n', options.group_filter);
    fprintf('  Memory filter: %s\n', options.memory);
    fprintf('  Session filter: %s\n', options.session_filter);
    fprintf('  Grooming excluded: %s\n', mat2str(options.exclude_grooming));
    fprintf('  Colormap used: %s\n', options.colormap);
    if ~isempty(options.caxis)
        fprintf('  Color axis limits: [%.2f, %.2f]\n', options.caxis(1), options.caxis(2));
    else
        fprintf('  Color axis limits: Auto (percentile-based)\n');
    end
    if ~isempty(options.plot_specific_mouse)
        fprintf('  Validation plots shown for: %s\n', options.plot_specific_mouse);
    end
end

function [calculated_speed, correlation] = calculate_speed_from_position(time, x, y, original_speed)
    % Calculate speed from position data (x,y) and time, with conversion to cm/s
    % Inputs:
    %   time - Vector of time points
    %   x - Vector of x coordinates (in pixels)
    %   y - Vector of y coordinates (in pixels)
    %   original_speed - Original speed from data (for correlation)
    %
    % Outputs:
    %   calculated_speed - Vector of calculated speeds (in cm/s)
    %   correlation - Correlation coefficient between original and calculated speed
    
    % Conversion factor from pixels to cm (adjust if needed for cheeseboard)
    pixels_to_cm = 0.17; % 0.17 cm per pixel
    
    % Check if we have enough points to calculate speed
    if length(time) < 2
        calculated_speed = zeros(size(time));
        correlation = NaN;
        return;
    end
    
    % Calculate displacements
    dx = diff(x);
    dy = diff(y);
    dt = diff(time);
    
    % Calculate speed as distance/time
    displacement_pixels = sqrt(dx.^2 + dy.^2);
    displacement_cm = displacement_pixels * pixels_to_cm; % Convert to cm
    instant_speed = displacement_cm ./ dt; % Speed in cm/s
    
    % Handle potential division by zero or very small dt
    instant_speed(dt < 1e-6) = 0;
    instant_speed(isnan(instant_speed)) = 0;
    instant_speed(isinf(instant_speed)) = 0;
  
    % Smooth the calculated speed (using a simple moving average)
    window_size = min(5, length(instant_speed));
    if window_size > 1
        smoothed_speed = movmean(instant_speed, window_size);
    else
        smoothed_speed = instant_speed;
    end
    
    % Match lengths (append first value at the beginning)
    calculated_speed = [smoothed_speed(1); smoothed_speed];
    
    % Calculate correlation with original speed (if available)
    if nargin >= 4 && ~isempty(original_speed)
        try
            correlation = corr(calculated_speed, original_speed);
        catch
            % Handle any errors in correlation calculation
            correlation = NaN;
        end
    else
        correlation = NaN;
    end
end

function cmap = create_blue_white_red_colormap()
    % Create a blue-white-red colormap (64 colors)
    n = 64;
    half = round(n/2);
    
    % Blue to white
    blue_to_white = [linspace(0, 1, half)', linspace(0, 1, half)', ones(half, 1)];
    
    % White to red
    white_to_red = [ones(n-half, 1), linspace(1, 0, n-half)', linspace(1, 0, n-half)'];
    
    cmap = [blue_to_white; white_to_red];
end

function cmap = create_red_blue_colormap()
    % Create a red-blue colormap (64 colors)
    n = 64;
    cmap = [linspace(1, 0, n)', zeros(n, 1), linspace(0, 1, n)'];
end