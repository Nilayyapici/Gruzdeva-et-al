function validateDistanceLimit(mice, options)
    % VALIDATEDISTANCELIMIT Visualizes the distance limit threshold on mouse trajectory
    %   mice: cell array with mouse data
    %   options: struct with field dist_limit
    %
    % This function helps validate the dist_limit parameter by:
    %   - Plotting the full trajectory of a random mouse
    %   - Highlighting points within the distance limit in red
    %   - Showing the estimated food location and distance circle
    %   - Providing statistics about data inclusion/exclusion
    
    % Set default distance limit if not provided
    if ~isfield(options, 'dist_limit')
        options.dist_limit = 5;
        fprintf('Using default dist_limit = 5\n');
    end
    
    % Define constants for data columns
    COL_X = 2;       % X coordinate
    COL_Y = 3;       % Y coordinate  
    COL_DIST = 5;    % Distance to food
    COL_GROOM = 10;  % Grooming
    COL_EATING = 9;  % Eating
    
    % Filter mice based on options if available
    valid_mice = [];
    for i = 1:size(mice, 1)
        % Check if this mouse matches the state/source filters
        include_mouse = true;
        
        if isfield(options, 'state') && ~isempty(options.state)
            if ischar(options.state)
                states_to_check = {options.state};
            else
                states_to_check = options.state;
            end
            if ~any(strcmp(mice{i, 2}, states_to_check))
                include_mouse = false;
            end
        end
        
        if isfield(options, 'source') && ~isempty(options.source)
            if ischar(options.source)
                sources_to_check = {options.source};
            else
                sources_to_check = options.source;
            end
            if ~any(strcmp(mice{i, 3}, sources_to_check))
                include_mouse = false;
            end
        end
        
        if include_mouse
            valid_mice = [valid_mice; i];
        end
    end
    
    if isempty(valid_mice)
        error('No mice match the specified filters');
    end
    
    % Select a random mouse from valid mice
    random_idx = valid_mice(randi(length(valid_mice)));
    
    % Extract data for the selected mouse
    mouse_info = mice{random_idx, 1};
    group = mice{random_idx, 2};
    stimulus = mice{random_idx, 3};
    data = mice{random_idx, 4};
    
    % Check if discovery frame is available
    discovery_frame = [];
    if size(mice, 2) >= 6
        discovery_frame = mice{random_idx, 6};
    end
    
    % Create figure
    figure('Position', [100, 100, 1000, 700], 'Name', ['Distance Limit Validation: ', mouse_info]);
    
    % Create subplot for trajectory
    subplot(2, 2, [1, 2]);
    
    % Plot the full trajectory in gray
    plot(data(:, COL_X), data(:, COL_Y), '-', 'LineWidth', 1, 'Color', [0.8 0.8 0.8]);
    hold on;
    
    % Find points within distance limit
    within_limit = data(:, COL_DIST) < options.dist_limit;
    
    % Apply grooming and eating filters if specified in options
    if isfield(options, 'remove_grooming') && options.remove_grooming
        within_limit = within_limit & data(:, COL_GROOM) == 0 & data(:, COL_EATING) == 0;
    end
    
    % Scatter red points for data within distance limit (after all filters)
    if any(within_limit)
        scatter(data(within_limit, COL_X), data(within_limit, COL_Y), 15, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    end
    
    % Estimate food location (point with minimum distance)
    [min_dist, min_idx] = min(data(:, COL_DIST));
    food_x = data(min_idx, COL_X);
    food_y = data(min_idx, COL_Y);
    
    % Plot food location
    scatter(food_x, food_y, 150, 'g', 'filled', 'Marker', 's', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Draw circle showing distance limit around food (approximate)
    % Note: This assumes distance is calculated in the same units as x,y coordinates
    theta = linspace(0, 2*pi, 100);
    circle_x = food_x + options.dist_limit * cos(theta);
    circle_y = food_y + options.dist_limit * sin(theta);
    plot(circle_x, circle_y, 'r--', 'LineWidth', 2);
    
    % Mark discovery frame if available
    if ~isempty(discovery_frame) && discovery_frame <= size(data, 1)
        scatter(data(discovery_frame, COL_X), data(discovery_frame, COL_Y), 100, 'b', 'filled', 'Marker', '^');
    end
    
    % Format plot
    title(sprintf('Distance Limit Validation (limit = %.1f)\nMouse: %s (%s, %s)', ...
          options.dist_limit, mouse_info, group, stimulus), 'FontSize', 12);
    xlabel('X Coordinate', 'FontSize', 11);
    ylabel('Y Coordinate', 'FontSize', 11);
    
    % Create legend
    legend_items = {'Full Trajectory', 'Points < Distance Limit', 'Food Location', 'Distance Limit Circle'};
    if ~isempty(discovery_frame)
        legend_items{end+1} = 'Discovery Frame';
    end
    legend(legend_items, 'Location', 'best', 'FontSize', 9);
    
    axis equal;
    grid on;
    box off;
    
    % Subplot for distance over time
    subplot(2, 2, 3);
    plot(1:length(data), data(:, COL_DIST), 'k-', 'LineWidth', 1);
    hold on;
    yline(options.dist_limit, 'r--', 'Distance Limit', 'LineWidth', 2);
    
    % Highlight points within limit
    if any(within_limit)
        scatter(find(within_limit), data(within_limit, COL_DIST), 10, 'r', 'filled');
    end
    
    % Plot eating events
    eating_frames = find(data(:, COL_EATING) == 1);
    if ~isempty(eating_frames)
        % Plot vertical patches for eating events (semi-transparent)
        y_limits = ylim;
        patch_width = 0.5; % Width of the vertical patch
        for i = 1:length(eating_frames)
            % Create a thin vertical rectangle using patch
            x_patch = [eating_frames(i) - patch_width/2, eating_frames(i) + patch_width/2, ...
                      eating_frames(i) + patch_width/2, eating_frames(i) - patch_width/2];
            y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
            patch(x_patch, y_patch, [1, 0.5, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        % Also add scatter points at the actual distance during eating
        scatter(eating_frames, data(eating_frames, COL_DIST), 25, [1, 0.5, 0], 'filled', 'Marker', 'd');
    end
    
    % Mark discovery frame if available
    if ~isempty(discovery_frame)
        xline(discovery_frame, 'b-', 'Discovery', 'LineWidth', 2);
    end
    
    title('Distance to Food Over Time', 'FontSize', 11);
    xlabel('Frame Number', 'FontSize', 10);
    ylabel('Distance to Food', 'FontSize', 10);
    
    % Add legend for distance plot
    legend_items_dist = {'Distance to Food', 'Distance Limit', 'Points < Limit'};
    if ~isempty(eating_frames)
        legend_items_dist{end+1} = 'Eating Events';
    end
    if ~isempty(discovery_frame)
        legend_items_dist{end+1} = 'Discovery';
    end
    legend(legend_items_dist, 'Location', 'best', 'FontSize', 8);
    
    grid on;
    box off;
    
    % Subplot for statistics
    subplot(2, 2, 4);
    axis off;
    
    % Calculate statistics
    n_total = size(data, 1);
    n_within_dist = sum(data(:, COL_DIST) < options.dist_limit);
    n_within_all_filters = sum(within_limit);
    
    percentage_dist = (n_within_dist / n_total) * 100;
    percentage_all = (n_within_all_filters / n_total) * 100;
    
    % Additional filter statistics
    n_grooming = sum(data(:, COL_GROOM) == 1);
    n_eating = sum(data(:, COL_EATING) == 1);
    
    % Eating event statistics
    eating_frames = find(data(:, COL_EATING) == 1);
    n_eating_within_limit = sum(data(eating_frames, COL_DIST) < options.dist_limit);
    if ~isempty(eating_frames)
        mean_eating_distance = mean(data(eating_frames, COL_DIST));
        min_eating_distance = min(data(eating_frames, COL_DIST));
    else
        mean_eating_distance = NaN;
        min_eating_distance = NaN;
    end
    
    % Create statistics text
    stats_text = sprintf([...
        'DISTANCE LIMIT VALIDATION\n\n' ...
        'Mouse: %s\n' ...
        'Group: %s\n' ...
        'Stimulus: %s\n\n' ...
        'Distance Limit: %.1f\n' ...
        'Min Distance: %.2f\n\n' ...
        'FILTERING RESULTS:\n' ...
        'Total frames: %d\n' ...
        'Within distance limit: %d (%.1f%%)\n' ...
        'After all filters: %d (%.1f%%)\n\n' ...
        'EATING EVENTS:\n' ...
        'Total eating frames: %d (%.1f%%)\n' ...
        'Eating within limit: %d (%.1f%%)\n' ...
        'Mean eating distance: %.2f\n' ...
        'Min eating distance: %.2f\n\n' ...
        'OTHER EXCLUDED:\n' ...
        'Grooming frames: %d (%.1f%%)\n'], ...
        mouse_info, group, stimulus, ...
        options.dist_limit, min_dist, ...
        n_total, n_within_dist, percentage_dist, ...
        n_within_all_filters, percentage_all, ...
        n_eating, (n_eating/n_total)*100, ...
        n_eating_within_limit, (n_eating_within_limit/max(1,n_eating))*100, ...
        mean_eating_distance, min_eating_distance, ...
        n_grooming, (n_grooming/n_total)*100);
    
    text(0.05, 0.95, stats_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'FontSize', 10, 'FontName', 'FixedWidth', ...
         'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', 'black', 'Margin', 5);
    
    % Print summary to command window
    fprintf('\n=== DISTANCE LIMIT VALIDATION ===\n');
    fprintf('Mouse: %s (%s, %s)\n', mouse_info, group, stimulus);
    fprintf('Distance limit: %.1f\n', options.dist_limit);
    fprintf('Points within distance limit: %d/%d (%.1f%%)\n', ...
            n_within_dist, n_total, percentage_dist);
    fprintf('Points after all filters: %d/%d (%.1f%%)\n', ...
            n_within_all_filters, n_total, percentage_all);
    fprintf('Data reduction: %.1f%% → %.1f%%\n', percentage_dist, percentage_all);
    fprintf('\nEating events: %d total, %d within limit (%.1f%% of eating)\n', ...
            n_eating, n_eating_within_limit, (n_eating_within_limit/max(1,n_eating))*100);
    if ~isempty(eating_frames)
        fprintf('Eating distance stats: mean=%.2f, min=%.2f\n', mean_eating_distance, min_eating_distance);
    end
    fprintf('================================\n\n');
    
    hold off;
end