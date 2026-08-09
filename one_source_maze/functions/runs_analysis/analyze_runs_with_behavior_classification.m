function [run_data] = analyze_runs_with_behavior_classification(mice_all, food_area, group_filter, source_filter, validate_plots, max_time_gap, max_spatial_dist, lookback_window)
    % Analyze mouse runs and classify them by what they end/start with
    % 
    % Key difference from analyze_single_arm_runs_with_speed_check:
    % - TOWARDS runs: Classified by what happens at/after the run END in ORIGINAL data
    % - AWAY runs: Classified by what happens at/before the run START in ORIGINAL data
    %
    % Classification:
    %   'eating' - Run ends/starts at eating
    %   'visit' - Run ends/starts at food interaction (no eating)
    %   'neither' - Run ends/starts without eating or visiting
    
    % Default parameters
    if nargin < 3, group_filter = 'fasted'; end
    if nargin < 4, source_filter = 'all'; end
    if nargin < 5, validate_plots = false; end
    if nargin < 6, max_time_gap = 1.0; end
    if nargin < 7, max_spatial_dist = 5.0; end
    if nargin < 8, lookback_window = 10; end  % seconds to look at endpoints
    
    % Constants for data columns
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    COL_SPEED = 4;
    COL_DIST = 5;
    COL_PATH = 6;
    COL_DOOR = 7;
    COL_INTERACTION = 8;
    COL_EATING = 9;
    COL_GROOMING = 10;
    COL_DFF = 11;
    
    threshold = 0.3;
    MAX_TIME_GAP = max_time_gap;
    MAX_SPATIAL_DIST = max_spatial_dist;
    pixels_to_cm = 0.17;
    
    run_data = struct('mouse_id', {}, 'session', {}, 'runs', {});
    
    % Loop through each entry in mice_all
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        group = mice_all{i, 2};
        food_type = mice_all{i, 3};
        data_original = mice_all{i, 4};  % Keep original data
        
        % Process mice based on filters
        if (strcmp(group_filter, 'all') || strcmp(group, group_filter)) && ...
           (strcmp(source_filter, 'all') || strcmp(food_type, source_filter))
            
            % Session identification
            session_number = contains(session_info, 'sess0') * 0 + ...
                             contains(session_info, 'sess1') * 1 + ...
                             contains(session_info, 'sess2') * 2 + ...
                             contains(session_info, 'sess3') * 3;
            
            % Apply door filtering to original data
            if session_number == 2
                data_original = data_original(data_original(:,COL_DOOR) == 0, :);
            elseif session_number == 1
                data_original = data_original(data_original(:,COL_DOOR) == 1, :);
                if ~isempty(data_original)
                    data_original(:,COL_TIME) = data_original(:,COL_TIME) - data_original(1,COL_TIME);
                end
            elseif session_number == 0
                data_original = data_original(data_original(:,COL_DOOR) == 0, :);
            elseif session_number == 3
                data_original = data_original(data_original(:,COL_DOOR) == 1, :);
                if ~isempty(data_original)
                    data_original(:,COL_TIME) = data_original(:,COL_TIME) - data_original(1,COL_TIME);
                end
            else
                continue;
            end
            
            if isempty(data_original)
                continue;
            end
            
            % Keep original data for classification
            % For run detection, we'll use original data but mark eating/grooming points
            data_for_detection = data_original;
            
            if size(data_for_detection, 1) < 4
                continue;
            end
            
            % Apply filters for distance and dF/F
            [b, a] = butter(1, 0.05, 'low');
            [b1, a1] = butter(1, 0.05, 'low');
            
            filtered_dist = filtfilt(b1, a1, data_for_detection(:, COL_DIST));
            dist_derivative = [diff(filtered_dist); 0];
            
            if size(dist_derivative, 1) >= 4
                dist_derivative = filtfilt(b, a, dist_derivative);
            end
            
            data_for_detection(:, COL_DFF) = filtfilt(b1, a1, data_for_detection(:, COL_DFF));
            
            % Create mask for valid points (not eating/grooming) for run data extraction
            valid_points = data_for_detection(:,COL_GROOMING) == 0 & ...
                          data_for_detection(:,COL_EATING) == 0;
            
            % Extract mouse ID
            mouse_id = session_info;
            if contains(session_info, '_sess')
                mouse_id = extractBefore(session_info, '_sess');
            end
            
            mouse_entry = struct('mouse_id', mouse_id, 'session', session_number, 'runs', []);
            run_counter = 1;
            
            % Track validation plot data
            if validate_plots
                towards_periods = zeros(0, 2);
                away_periods = zeros(0, 2);
            end
            
            in_towards_run = false;
            in_away_run = false;
            current_run_start = 0;
            
            % Analyze runs in ORIGINAL data (detecting based on distance changes)
            ii = 2;
            while ii <= size(data_for_detection, 1) - 1
                % TOWARDS run detection
                if ~in_towards_run && dist_derivative(ii) < -threshold && data_for_detection(ii, COL_DIST) <= food_area
                    in_towards_run = true;
                    current_run_start = find(dist_derivative(1:ii) > -threshold, 1, 'last') + 1;
                    if isempty(current_run_start), current_run_start = 1; end
                    
                elseif in_towards_run
                    if dist_derivative(ii) >= -threshold || data_for_detection(ii, COL_DIST) > food_area
                        in_towards_run = false;
                        if min(data_for_detection(current_run_start:ii, COL_DIST)) <= food_area
                            % Get run indices
                            run_indices = (current_run_start:ii)';
                            run_indices = filter_run_artifacts(data_for_detection, run_indices, MAX_TIME_GAP, MAX_SPATIAL_DIST, pixels_to_cm);
                            
                            if length(run_indices) >= 3
                                % Classify by looking at ORIGINAL data at run END
                                run_end_idx = run_indices(end);
                                run_end_time = data_for_detection(run_end_idx, COL_TIME);
                                
                                % Check if run END is during eating/visit, OR look forward in time
                                lookahead_mask = data_for_detection(:, COL_TIME) >= run_end_time & ...
                                               data_for_detection(:, COL_TIME) <= run_end_time + lookback_window;
                                
                                has_eating = any(data_for_detection(lookahead_mask, COL_EATING) == 1);
                                has_interaction = any(data_for_detection(lookahead_mask, COL_INTERACTION) == 1);
                                
                                if has_eating
                                    behavior_class = 'eating';
                                elseif has_interaction
                                    behavior_class = 'visit';
                                else
                                    behavior_class = 'neither';
                                end
                                
                                % For run data, only use valid (non-eating/grooming) points
                                run_valid_mask = valid_points(run_indices);
                                run_indices_valid = run_indices(run_valid_mask);
                                
                                if length(run_indices_valid) >= 3
                                    % Calculate speed
                                    [calculated_speed, correlation] = calculate_speed_from_position(...
                                        data_for_detection(run_indices_valid, COL_TIME),...
                                        data_for_detection(run_indices_valid, COL_X),...
                                        data_for_detection(run_indices_valid, COL_Y),...
                                        data_for_detection(run_indices_valid, COL_SPEED));
                                    
                                    % Create run entry
                                    run_entry = struct('id', run_counter, ...
                                                      'type', 'towards', ...
                                                      'behavior', behavior_class, ...
                                                      'indices', run_indices_valid, ...
                                                      'time', data_for_detection(run_indices_valid, COL_TIME), ...
                                                      'xcoordinate', data_for_detection(run_indices_valid, COL_X), ...
                                                      'ycoordinate', data_for_detection(run_indices_valid, COL_Y), ...
                                                      'distance', data_for_detection(run_indices_valid, COL_DIST), ...
                                                      'dff', data_for_detection(run_indices_valid, COL_DFF), ...
                                                      'speed', data_for_detection(run_indices_valid, COL_SPEED), ...
                                                      'calculated_speed', calculated_speed, ...
                                                      'speed_correlation', correlation, ...
                                                      'start_time', data_for_detection(run_indices_valid(1), COL_TIME), ...
                                                      'end_time', data_for_detection(run_indices_valid(end), COL_TIME), ...
                                                      'duration', data_for_detection(run_indices_valid(end), COL_TIME) - data_for_detection(run_indices_valid(1), COL_TIME));
                                    
                                    mouse_entry.runs = [mouse_entry.runs; run_entry];
                                    run_counter = run_counter + 1;
                                    
                                    if validate_plots
                                        towards_periods(end+1,:) = [run_indices(1), run_indices(end)];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % AWAY run detection
                if ~in_away_run && dist_derivative(ii) > threshold && data_for_detection(ii, COL_DIST) <= food_area
                    in_away_run = true;
                    current_run_start = find(dist_derivative(1:ii) < threshold, 1, 'last') + 1;
                    if isempty(current_run_start), current_run_start = 1; end
                    
                elseif in_away_run
                    if dist_derivative(ii) <= threshold
                        in_away_run = false;
                        if any(data_for_detection(current_run_start:current_run_start+min(10,ii-current_run_start), COL_DIST) <= food_area)
                            % Get run indices
                            run_indices = (current_run_start:ii)';
                            run_indices = filter_run_artifacts(data_for_detection, run_indices, MAX_TIME_GAP, MAX_SPATIAL_DIST, pixels_to_cm);
                            
                            if length(run_indices) >= 3
                                % Classify by looking at ORIGINAL data at run START
                                run_start_idx = run_indices(1);
                                run_start_time = data_for_detection(run_start_idx, COL_TIME);
                                
                                % Check if run START is during eating/visit, OR look backward in time
                                lookback_mask = data_for_detection(:, COL_TIME) >= run_start_time - lookback_window & ...
                                              data_for_detection(:, COL_TIME) <= run_start_time;
                                
                                has_eating = any(data_for_detection(lookback_mask, COL_EATING) == 1);
                                has_interaction = any(data_for_detection(lookback_mask, COL_INTERACTION) == 1);
                                
                                if has_eating
                                    behavior_class = 'eating';
                                elseif has_interaction
                                    behavior_class = 'visit';
                                else
                                    behavior_class = 'neither';
                                end
                                
                                % For run data, only use valid (non-eating/grooming) points
                                run_valid_mask = valid_points(run_indices);
                                run_indices_valid = run_indices(run_valid_mask);
                                
                                if length(run_indices_valid) >= 3
                                    % Calculate speed
                                    [calculated_speed, correlation] = calculate_speed_from_position(...
                                        data_for_detection(run_indices_valid, COL_TIME),...
                                        data_for_detection(run_indices_valid, COL_X),...
                                        data_for_detection(run_indices_valid, COL_Y),...
                                        data_for_detection(run_indices_valid, COL_SPEED));
                                    
                                    % Create run entry
                                    run_entry = struct('id', run_counter, ...
                                                      'type', 'away', ...
                                                      'behavior', behavior_class, ...
                                                      'indices', run_indices_valid, ...
                                                      'time', data_for_detection(run_indices_valid, COL_TIME), ...
                                                      'xcoordinate', data_for_detection(run_indices_valid, COL_X), ...
                                                      'ycoordinate', data_for_detection(run_indices_valid, COL_Y), ...
                                                      'distance', data_for_detection(run_indices_valid, COL_DIST), ...
                                                      'dff', data_for_detection(run_indices_valid, COL_DFF), ...
                                                      'speed', data_for_detection(run_indices_valid, COL_SPEED), ...
                                                      'calculated_speed', calculated_speed, ...
                                                      'speed_correlation', correlation, ...
                                                      'start_time', data_for_detection(run_indices_valid(1), COL_TIME), ...
                                                      'end_time', data_for_detection(run_indices_valid(end), COL_TIME), ...
                                                      'duration', data_for_detection(run_indices_valid(end), COL_TIME) - data_for_detection(run_indices_valid(1), COL_TIME));
                                    
                                    mouse_entry.runs = [mouse_entry.runs; run_entry];
                                    run_counter = run_counter + 1;
                                    
                                    if validate_plots
                                        away_periods(end+1,:) = [run_indices(1), run_indices(end)];
                                    end
                                end
                            end
                        end
                    end
                end
                
                ii = ii + 1;
            end
            
            % Add mouse data if it has runs
            if ~isempty(mouse_entry.runs)
                run_data = [run_data; mouse_entry];
                
                if validate_plots
                    plot_classification_validation(data_for_detection, mouse_entry.runs, ...
                        towards_periods, away_periods, food_area, session_info, lookback_window);
                end
            end
        end
    end
    
    % Print summary
    if ~isempty(run_data)
        total_runs = 0;
        towards_eating = 0;
        towards_visit = 0;
        towards_neither = 0;
        away_eating = 0;
        away_visit = 0;
        away_neither = 0;
        
        for i = 1:length(run_data)
            for j = 1:length(run_data(i).runs)
                run = run_data(i).runs(j);
                total_runs = total_runs + 1;
                
                if strcmp(run.type, 'towards')
                    if strcmp(run.behavior, 'eating')
                        towards_eating = towards_eating + 1;
                    elseif strcmp(run.behavior, 'visit')
                        towards_visit = towards_visit + 1;
                    else
                        towards_neither = towards_neither + 1;
                    end
                elseif strcmp(run.type, 'away')
                    if strcmp(run.behavior, 'eating')
                        away_eating = away_eating + 1;
                    elseif strcmp(run.behavior, 'visit')
                        away_visit = away_visit + 1;
                    else
                        away_neither = away_neither + 1;
                    end
                end
            end
        end
        
        fprintf('\n=== Run Classification Summary ===\n');
        fprintf('Total: %d mice with %d runs\n', length(run_data), total_runs);
        fprintf('Filters: Group="%s", Source="%s"\n', group_filter, source_filter);
        fprintf('Time window: %.1fs\n\n', lookback_window);
        
        fprintf('TOWARDS runs (classified by what follows):\n');
        fprintf('  Eating:  %d (%.1f%%)\n', towards_eating, 100*towards_eating/total_runs);
        fprintf('  Visit:   %d (%.1f%%)\n', towards_visit, 100*towards_visit/total_runs);
        fprintf('  Neither: %d (%.1f%%)\n\n', towards_neither, 100*towards_neither/total_runs);
        
        fprintf('AWAY runs (classified by what preceded):\n');
        fprintf('  Eating:  %d (%.1f%%)\n', away_eating, 100*away_eating/total_runs);
        fprintf('  Visit:   %d (%.1f%%)\n', away_visit, 100*away_visit/total_runs);
        fprintf('  Neither: %d (%.1f%%)\n', away_neither, 100*away_neither/total_runs);
    else
        fprintf('No valid runs found for group "%s" and source "%s"\n', group_filter, source_filter);
    end
end

% Include all the helper functions from the original
function filtered_indices = filter_run_artifacts(data, indices, max_time_gap, max_spatial_dist, pixels_to_cm)
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    
    if length(indices) < 2
        filtered_indices = indices;
        return;
    end
    
    segments = {};
    current_segment = indices(1);
    
    for i = 2:length(indices)
        curr_idx = indices(i);
        prev_idx = indices(i-1);
        
        time_gap = data(curr_idx, COL_TIME) - data(prev_idx, COL_TIME);
        dx = data(curr_idx, COL_X) - data(prev_idx, COL_X);
        dy = data(curr_idx, COL_Y) - data(prev_idx, COL_Y);
        spatial_dist = sqrt(dx^2 + dy^2) * pixels_to_cm;
        
        if time_gap <= max_time_gap && spatial_dist <= max_spatial_dist
            current_segment = [current_segment; curr_idx];
        else
            if length(current_segment) >= 3
                segments{end+1} = current_segment;
            end
            current_segment = curr_idx;
        end
    end
    
    if length(current_segment) >= 3
        segments{end+1} = current_segment;
    end
    
    if isempty(segments)
        filtered_indices = [];
    else
        [~, longest_idx] = max(cellfun(@length, segments));
        filtered_indices = segments{longest_idx};
    end
end

function [calculated_speed, correlation] = calculate_speed_from_position(time, x, y, original_speed)
    pixels_to_cm = 0.17;
    
    if length(time) < 2
        calculated_speed = zeros(size(time));
        correlation = NaN;
        return;
    end
    
    dx = diff(x);
    dy = diff(y);
    dt = diff(time);
    
    displacement_pixels = sqrt(dx.^2 + dy.^2);
    displacement_cm = displacement_pixels * pixels_to_cm;
    instant_speed = displacement_cm ./ dt;
    
    instant_speed(dt < 1e-6) = 0;
    instant_speed(isnan(instant_speed)) = 0;
    instant_speed(isinf(instant_speed)) = 0;
    
    window_size = min(5, length(instant_speed));
    if window_size > 1
        smoothed_speed = movmean(instant_speed, window_size);
    else
        smoothed_speed = instant_speed;
    end
    
    calculated_speed = [smoothed_speed(1); smoothed_speed];
    
    if nargin >= 4 && ~isempty(original_speed)
        try
            correlation = corr(calculated_speed, original_speed);
        catch
            correlation = NaN;
        end
    else
        correlation = NaN;
    end
end

function plot_classification_validation(data_original, runs, towards_periods, away_periods, food_area, session_info, lookback_window)
    % Create validation plot showing eating/visit patches and classified runs
    
    COL_TIME = 1;
    COL_DIST = 5;
    COL_INTERACTION = 8;
    COL_EATING = 9;
    
    % Create figure
    figure('Position', [100, 100, 1400, 700], 'Name', ['Classification Validation: ', session_info]);
    
    % Plot distance to food from original data
    plot(data_original(:, COL_TIME), data_original(:, COL_DIST), 'k-', 'LineWidth', 1);
    hold on;
    
    % Plot food area threshold
    yline(food_area, 'k--', 'Food Area', 'LineWidth', 2);
    
    % Find eating and visit periods in original data
    eating_starts = find(diff([0; data_original(:, COL_EATING)]) == 1);
    eating_ends = find(diff([data_original(:, COL_EATING); 0]) == -1);
    
    visit_starts = find(diff([0; data_original(:, COL_INTERACTION)]) == 1);
    visit_ends = find(diff([data_original(:, COL_INTERACTION); 0]) == -1);
    
    % Remove visits that overlap with eating
    visit_mask = true(length(visit_starts), 1);
    for v = 1:length(visit_starts)
        visit_time_start = data_original(visit_starts(v), COL_TIME);
        visit_time_end = data_original(visit_ends(v), COL_TIME);
        
        % Check if any eating period overlaps
        for e = 1:length(eating_starts)
            eating_time_start = data_original(eating_starts(e), COL_TIME);
            eating_time_end = data_original(eating_ends(e), COL_TIME);
            
            % Check for overlap
            if (visit_time_start <= eating_time_end && visit_time_end >= eating_time_start)
                visit_mask(v) = false;
                break;
            end
        end
    end
    visit_starts = visit_starts(visit_mask);
    visit_ends = visit_ends(visit_mask);
    
    % Get y-axis limits for patches
    y_lim = ylim;
    
    % Plot eating patches (red with transparency)
    for e = 1:length(eating_starts)
        if eating_starts(e) <= size(data_original, 1) && eating_ends(e) <= size(data_original, 1)
            t_start = data_original(eating_starts(e), COL_TIME);
            t_end = data_original(eating_ends(e), COL_TIME);
            patch([t_start, t_end, t_end, t_start], ...
                  [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                  [1, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    
    % Plot visit patches (yellow with transparency)
    for v = 1:length(visit_starts)
        if visit_starts(v) <= size(data_original, 1) && visit_ends(v) <= size(data_original, 1)
            t_start = data_original(visit_starts(v), COL_TIME);
            t_end = data_original(visit_ends(v), COL_TIME);
            patch([t_start, t_end, t_end, t_start], ...
                  [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                  [1, 1, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    
    % Define colors for run classifications
    % Towards runs: blue shades
    towards_eating_color = [0, 0, 1];      % Dark blue
    towards_visit_color = [0.3, 0.6, 1];   % Light blue
    towards_neither_color = [0.7, 0.7, 0.7]; % Gray
    
    % Away runs: red shades
    away_eating_color = [1, 0, 0];         % Dark red
    away_visit_color = [1, 0.4, 0.4];      % Light red/pink
    away_neither_color = [0.7, 0.7, 0.7];  % Gray
    
    % Plot runs colored by their classification
    for r = 1:length(runs)
        run = runs(r);
        
        % Get color based on type and behavior
        if strcmp(run.type, 'towards')
            if strcmp(run.behavior, 'eating')
                run_color = towards_eating_color;
            elseif strcmp(run.behavior, 'visit')
                run_color = towards_visit_color;
            else
                run_color = towards_neither_color;
            end
        else % away
            if strcmp(run.behavior, 'eating')
                run_color = away_eating_color;
            elseif strcmp(run.behavior, 'visit')
                run_color = away_visit_color;
            else
                run_color = away_neither_color;
            end
        end
        
        % Plot the run
        scatter(run.time, run.distance, 12, run_color, 'filled', 'MarkerFaceAlpha', 0.8);
    end
    
    % Set labels and title
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Distance to Food', 'FontSize', 14);
    title(sprintf('Run Classification: %s (lookback window = %.1fs)', session_info, lookback_window), 'FontSize', 12);
    
    % Create custom legend
    legend_entries = {};
    legend_handles = [];
    
    % Add eating/visit patches to legend
    h1 = patch(NaN, NaN, [1, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend_handles = [legend_handles, h1];
    legend_entries{end+1} = 'Eating';
    
    h2 = patch(NaN, NaN, [1, 1, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend_handles = [legend_handles, h2];
    legend_entries{end+1} = 'Visit (no eating)';
    
    % Add towards runs to legend
    h3 = scatter(NaN, NaN, 20, towards_eating_color, 'filled');
    legend_handles = [legend_handles, h3];
    legend_entries{end+1} = 'Towards → Eating';
    
    h4 = scatter(NaN, NaN, 20, towards_visit_color, 'filled');
    legend_handles = [legend_handles, h4];
    legend_entries{end+1} = 'Towards → Visit';
    
    h5 = scatter(NaN, NaN, 20, towards_neither_color, 'filled');
    legend_handles = [legend_handles, h5];
    legend_entries{end+1} = 'Towards → Neither';
    
    % Add away runs to legend
    h6 = scatter(NaN, NaN, 20, away_eating_color, 'filled');
    legend_handles = [legend_handles, h6];
    legend_entries{end+1} = 'Eating → Away';
    
    h7 = scatter(NaN, NaN, 20, away_visit_color, 'filled');
    legend_handles = [legend_handles, h7];
    legend_entries{end+1} = 'Visit → Away';
    
    h8 = scatter(NaN, NaN, 20, away_neither_color, 'filled');
    legend_handles = [legend_handles, h8];
    legend_entries{end+1} = 'Neither → Away';
    
    legend(legend_handles, legend_entries, 'Location', 'best', 'FontSize', 9);
    
    % Set reasonable xlim to see detail
    xlim([0 1000]);
    
    grid off;
    box off;
    hold off;
end