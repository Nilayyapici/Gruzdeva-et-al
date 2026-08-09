function [run_data] = analyze_food_runs_by_visits(mice_all, group_filter, source_filter, validate_plots)
    % Analyze mouse runs towards and away from food based on eating/food visit events
    % Runs are defined by eating/food visit events:
    %   - AWAY run: starts at end of eating/food visit, ends when distance starts decreasing
    %   - TOWARDS run: starts when distance decreases, ends at next eating/food visit
    % Excludes grooming periods and eating/food visit points from runs
    %
    % Inputs:
    %   mice_all - Cell array with mouse data
    %   group_filter - 'fasted', 'fed', or 'all'
    %   source_filter - 'food', 'gel', or 'all'
    %   validate_plots - true/false to generate validation plots
    %
    % Output:
    %   run_data - Structure array with detected runs for each mouse/session
    
    % Default parameters
    if nargin < 2
        group_filter = 'all';
    end
    if nargin < 3
        source_filter = 'all';
    end
    if nargin < 4
        validate_plots = false;
    end
    
    % Constants for data columns
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    COL_SPEED = 4;
    COL_DIST = 5;
    COL_PATH = 6;
    COL_DOOR = 7;
    COL_FOOD_VISIT = 8;
    COL_EATING = 9;
    COL_GROOMING = 10;
    COL_DFF = 11;
    
    % Initialize output structure
    run_data = struct('mouse_id', {}, 'session', {}, 'group', {}, 'source', {}, 'runs', {});
    
    % Loop through each mouse
    for i = 1:size(mice_all, 1)
        % Extract mouse information
        session_info = mice_all{i, 1};
        group = mice_all{i, 2};
        food_type = mice_all{i, 3};
        data = mice_all{i, 4};
        
        % Apply filters
        if ~strcmp(group_filter, 'all') && ~strcmp(group, group_filter)
            continue;
        end
        if ~strcmp(source_filter, 'all') && ~strcmp(food_type, source_filter)
            continue;
        end
        
        % Extract session number
        if contains(session_info, 'sess0')
            session_number = 0;
        elseif contains(session_info, 'sess1')
            session_number = 1;
        elseif contains(session_info, 'sess2')
            session_number = 2;
        elseif contains(session_info, 'sess3')
            session_number = 3;
        else
            continue; % Skip if no valid session info
        end
        
        % Extract mouse ID
        mouse_id = session_info;
        if contains(session_info, '_sess')
            mouse_id = extractBefore(session_info, '_sess');
        end
        
        % Check if we have enough data
        if isempty(data) || size(data, 1) < 10
            continue;
        end
        
        % Remove grooming periods from data
        valid_idx = data(:, COL_GROOMING) == 0;
        data_clean = data(valid_idx, :);
        
        if isempty(data_clean) || size(data_clean, 1) < 10
            continue;
        end
        
        % Apply Butterworth filter to distance to reduce noise
        % Check if we have enough points for filtering
        if size(data_clean, 1) >= 12  % filtfilt requires at least 3*padlen where padlen=3*(filter_order)
            [b_dist, a_dist] = butter(1, 0.05, 'low');  % Low-pass filter for distance
            data_clean(:, COL_DIST) = filtfilt(b_dist, a_dist, data_clean(:, COL_DIST));
        end
        
        % Detect eating/food visit events
        eating_or_visit = (data_clean(:, COL_EATING) == 1) | (data_clean(:, COL_FOOD_VISIT) == 1);
        
        % Find transitions: start and end of eating/food visit bouts
        eating_diff = diff([0; eating_or_visit; 0]);
        eating_starts = find(eating_diff == 1);   % Start of eating/visit
        eating_ends = find(eating_diff == -1) - 1; % End of eating/visit
        
        if isempty(eating_ends)
            continue; % No eating/food visits detected
        end
        
        % Initialize runs array
        runs = [];
        run_counter = 1;
        
        % For validation plots
        if validate_plots
            towards_periods = [];
            away_periods = [];
        end
        
        % Process runs between eating/food visit events
        for j = 1:length(eating_ends)
            % === AWAY RUN: from end of current eating to where distance starts decreasing ===
            if eating_ends(j) < size(data_clean, 1) - 5
                % Start of away run (right after eating ends)
                away_start = eating_ends(j) + 1;
                
                % Find where distance starts decreasing (local maximum in distance)
                % Look ahead to find when mouse starts moving back towards food
                search_end = min(away_start + 500, size(data_clean, 1));
                
                % Calculate smoothed distance derivative
                dist_segment = data_clean(away_start:search_end, COL_DIST);
                
                if length(dist_segment) > 10
                    % Distance is already filtered, just calculate derivative
                    dist_deriv = diff(dist_segment);
                    
                    % Find first point where distance starts consistently decreasing
                    away_end_local = find_run_transition(dist_deriv, 'positive_to_negative');
                    
                    if ~isempty(away_end_local)
                        away_end = away_start + away_end_local - 1;
                    else
                        % If no clear transition, use the maximum distance point
                        [~, max_idx] = max(dist_segment);
                        away_end = away_start + max_idx - 1;
                    end
                    
                    % Make sure away run doesn't overlap with next eating event
                    if j < length(eating_starts)
                        away_end = min(away_end, eating_starts(j+1) - 1);
                    end
                    
                    % Extract away run data
                    if away_end > away_start + 2 % At least 3 points
                        run_indices = away_start:away_end;
                        
                        % Exclude any eating/food visit points that might be included
                        run_valid = ~eating_or_visit(run_indices);
                        run_indices = run_indices(run_valid);
                        
                        if length(run_indices) >= 3
                            % Create run entry
                            run_entry = create_run_entry(run_counter, 'away', run_indices, data_clean);
                            runs = [runs; run_entry];
                            run_counter = run_counter + 1;
                            
                            if validate_plots
                                away_periods = [away_periods; run_indices(1), run_indices(end)];
                            end
                        end
                    end
                end
            end
            
            % === TOWARDS RUN: from where distance was increasing to beginning of next eating event ===
            if j < length(eating_starts)
                % End of towards run (right before next eating starts)
                towards_end = eating_starts(j+1) - 1;
                
                % Find start of towards run by looking backward from eating start
                % to find where distance was last increasing (positive derivative)
                % or where it transitioned from increasing to decreasing
                
                % Determine search window start
                if exist('away_end', 'var') && away_end < towards_end
                    search_start = away_end + 1;
                else
                    search_start = eating_ends(j) + 1;
                end
                
                if towards_end > search_start + 5
                    % Calculate distance derivative in this window
                    search_window = search_start:towards_end;
                    dist_segment = data_clean(search_window, COL_DIST);
                    
                    % Distance is already filtered, just calculate derivative
                    dist_deriv = diff(dist_segment);
                    
                    % Find where distance transitioned from increasing to decreasing
                    % Look backward from the end to find last point where distance was increasing
                    towards_start_local = [];
                    
                    % Find the last point where derivative was positive (increasing distance)
                    % Then the towards run starts when it becomes negative (decreasing)
                    increasing_points = find(dist_deriv > 0);
                    
                    if ~isempty(increasing_points)
                        % Start of towards run is after the last increasing point
                        towards_start_local = increasing_points(end) + 1;
                    else
                        % If distance was always decreasing, start from beginning of window
                        towards_start_local = 1;
                    end
                    
                    % Convert back to full data indices
                    towards_start = search_start + towards_start_local - 1;
                    
                    % Extract towards run data
                    if towards_end > towards_start + 2 % At least 3 points
                        run_indices = towards_start:towards_end;
                        
                        % Exclude any eating/food visit points
                        run_valid = ~eating_or_visit(run_indices);
                        run_indices = run_indices(run_valid);
                        
                        if length(run_indices) >= 3
                            % Create run entry
                            run_entry = create_run_entry(run_counter, 'towards', run_indices, data_clean);
                            runs = [runs; run_entry];
                            run_counter = run_counter + 1;
                            
                            if validate_plots
                                towards_periods = [towards_periods; run_indices(1), run_indices(end)];
                            end
                        end
                    end
                end
            end
        end
        
        % Store runs for this mouse/session
        if ~isempty(runs)
            mouse_entry = struct(...
                'mouse_id', mouse_id, ...
                'session', session_number, ...
                'group', group, ...
                'source', food_type, ...
                'runs', runs);
            
            run_data = [run_data; mouse_entry];
            
            % Generate validation plot
            if validate_plots
                plot_run_validation(data_clean, towards_periods, away_periods, ...
                    eating_or_visit, session_info, group, food_type);
            end
        end
    end
    
    % Print summary
    print_summary(run_data, group_filter, source_filter);
end

function transition_idx = find_run_transition(derivative, transition_type)
    % Find where derivative transitions from positive to negative or vice versa
    % Uses a more robust approach with smoothing and consistency checks
    
    if isempty(derivative) || length(derivative) < 3
        transition_idx = [];
        return;
    end
    
    % Smooth derivative to reduce noise
    if length(derivative) > 5
        deriv_smooth = movmean(derivative, 3);
    else
        deriv_smooth = derivative;
    end
    
    if strcmp(transition_type, 'positive_to_negative')
        % Find where derivative changes from positive to negative
        % Look for sustained negative values
        negative_runs = find(deriv_smooth < 0);
        
        if ~isempty(negative_runs)
            % Find first sustained negative period (at least 3 consecutive points)
            for i = 1:length(negative_runs)-2
                if negative_runs(i+1) == negative_runs(i)+1 && ...
                   negative_runs(i+2) == negative_runs(i)+2
                    transition_idx = negative_runs(i);
                    return;
                end
            end
            % If no sustained period, use first negative point
            transition_idx = negative_runs(1);
        else
            transition_idx = [];
        end
        
    elseif strcmp(transition_type, 'negative_to_positive')
        % Find where derivative changes from negative to positive
        positive_runs = find(deriv_smooth > 0);
        
        if ~isempty(positive_runs)
            % Find first sustained positive period
            for i = 1:length(positive_runs)-2
                if positive_runs(i+1) == positive_runs(i)+1 && ...
                   positive_runs(i+2) == positive_runs(i)+2
                    transition_idx = positive_runs(i);
                    return;
                end
            end
            % If no sustained period, use first positive point
            transition_idx = positive_runs(1);
        else
            transition_idx = [];
        end
    else
        transition_idx = [];
    end
end

function run_entry = create_run_entry(run_id, run_type, indices, data)
    % Create a run entry structure with all relevant data
    
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    COL_SPEED = 4;
    COL_DIST = 5;
    COL_PATH = 6;
    COL_DFF = 11;
    
    run_entry = struct(...
        'id', run_id, ...
        'type', run_type, ...
        'indices', indices, ...
        'time', data(indices, COL_TIME), ...
        'x', data(indices, COL_X), ...
        'y', data(indices, COL_Y), ...
        'speed', data(indices, COL_SPEED), ...
        'distance', data(indices, COL_DIST), ...
        'path_length', data(indices, COL_PATH), ...
        'dff', data(indices, COL_DFF), ...
        'start_time', data(indices(1), COL_TIME), ...
        'end_time', data(indices(end), COL_TIME), ...
        'duration', data(indices(end), COL_TIME) - data(indices(1), COL_TIME), ...
        'start_dist', data(indices(1), COL_DIST), ...
        'end_dist', data(indices(end), COL_DIST), ...
        'dist_change', data(indices(end), COL_DIST) - data(indices(1), COL_DIST), ...
        'mean_dff', mean(data(indices, COL_DFF), 'omitnan'), ...
        'mean_speed', mean(data(indices, COL_SPEED), 'omitnan'));
end

function plot_run_validation(data, towards_periods, away_periods, eating_or_visit, ...
                             session_info, group, food_type)
    % Create validation plot showing detected runs
    
    COL_TIME = 1;
    COL_DIST = 5;
    COL_DFF = 11;
    
    % Create figure
    figure('Position', [100, 100, 1400, 600], 'Name', ['Run Validation: ', session_info]);
    
    % Subplot 1: Distance to food with runs highlighted
    subplot(2, 1, 1);
    
    % Plot distance
    plot(data(:, COL_TIME), data(:, COL_DIST), 'k-', 'LineWidth', 1);
    hold on;
    
    % Highlight eating/food visit periods in yellow
    eating_idx = find(eating_or_visit);
    if ~isempty(eating_idx)
        scatter(data(eating_idx, COL_TIME), data(eating_idx, COL_DIST), ...
            20, 'y', 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    % Highlight towards runs in blue
    for i = 1:size(towards_periods, 1)
        idx_start = towards_periods(i, 1);
        idx_end = towards_periods(i, 2);
        scatter(data(idx_start:idx_end, COL_TIME), data(idx_start:idx_end, COL_DIST), ...
            15, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % Highlight away runs in red
    for i = 1:size(away_periods, 1)
        idx_start = away_periods(i, 1);
        idx_end = away_periods(i, 2);
        scatter(data(idx_start:idx_end, COL_TIME), data(idx_start:idx_end, COL_DIST), ...
            15, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Distance to Food', 'FontSize', 12);
    title(sprintf('Distance to Food - %s (%s, %s)', session_info, group, food_type), ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    legend('Distance', 'Eating/Food Visit', 'Towards Runs', 'Away Runs', ...
        'Location', 'best');
    grid on;
    box off;
    hold off;
    
    % Subplot 2: dF/F signal with runs highlighted
    subplot(2, 1, 2);
    
    % Plot dF/F
    plot(data(:, COL_TIME), data(:, COL_DFF), 'k-', 'LineWidth', 0.5);
    hold on;
    
    % Highlight eating/food visit periods
    if ~isempty(eating_idx)
        scatter(data(eating_idx, COL_TIME), data(eating_idx, COL_DFF), ...
            20, 'y', 'filled', 'MarkerFaceAlpha', 0.6);
    end
    
    % Highlight towards runs
    for i = 1:size(towards_periods, 1)
        idx_start = towards_periods(i, 1);
        idx_end = towards_periods(i, 2);
        scatter(data(idx_start:idx_end, COL_TIME), data(idx_start:idx_end, COL_DFF), ...
            15, 'b', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % Highlight away runs
    for i = 1:size(away_periods, 1)
        idx_start = away_periods(i, 1);
        idx_end = away_periods(i, 2);
        scatter(data(idx_start:idx_end, COL_TIME), data(idx_start:idx_end, COL_DFF), ...
            15, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('dF/F', 'FontSize', 12);
    title(sprintf('dF/F Signal - %s', session_info), 'FontWeight', 'bold', 'Interpreter', 'none');
    grid on;
    box off;
    hold off;
end

function print_summary(run_data, group_filter, source_filter)
    % Print summary statistics
    
    if isempty(run_data)
        fprintf('No runs detected for group="%s", source="%s"\n', group_filter, source_filter);
        return;
    end
    
    % Count runs by type and session
    total_runs = 0;
    towards_runs = 0;
    away_runs = 0;
    session_counts = zeros(4, 1); % Sessions 0-3
    
    for i = 1:length(run_data)
        n_runs = length(run_data(i).runs);
        total_runs = total_runs + n_runs;
        session_counts(run_data(i).session + 1) = session_counts(run_data(i).session + 1) + n_runs;
        
        for j = 1:n_runs
            if strcmp(run_data(i).runs(j).type, 'towards')
                towards_runs = towards_runs + 1;
            else
                away_runs = away_runs + 1;
            end
        end
    end
    
    fprintf('\n=== Run Detection Summary ===\n');
    fprintf('Filters: group="%s", source="%s"\n', group_filter, source_filter);
    fprintf('Total mice/sessions: %d\n', length(run_data));
    fprintf('Total runs detected: %d\n', total_runs);
    fprintf('  Towards runs: %d (%.1f%%)\n', towards_runs, 100*towards_runs/total_runs);
    fprintf('  Away runs: %d (%.1f%%)\n', away_runs, 100*away_runs/total_runs);
    fprintf('\nRuns by session:\n');
    for s = 0:3
        if session_counts(s+1) > 0
            fprintf('  Session %d: %d runs\n', s, session_counts(s+1));
        end
    end
    fprintf('=============================\n\n');
end