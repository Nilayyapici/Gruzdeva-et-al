function trajectoryCheeseboard(mice, condition, options)
% TRAJECTORYCHEESEBOARD - Plots x,y trajectories for each mouse (Pre and Test sessions)
%
% Usage:
%   trajectoryCheeseboard(mice, 'saline')
%   trajectoryCheeseboard(mice, 'CNO')
%   trajectoryCheeseboard(mice, 'strong')
%   trajectoryCheeseboard(mice, 'weak')
%   
%   % With options:
%   options.time_limit = 10;            % Analyze only first 10 minutes (optional)
%   options.exclude_grooming = false;   % Exclude grooming periods (default: false)
%   options.speed_threshold = 2;        % Minimum speed to include (cm/s, optional)
%   options.figure_size = [100, 100, 1000, 500]; % Figure position and size
%   options.x_col = 2;                  % Column index for x-coordinate (default: 2)
%   options.y_col = 3;                  % Column index for y-coordinate (default: 3)

%% Input validation and defaults
if nargin < 2
    error('Both mice data and condition are required');
end

if nargin < 3
    options = struct();
end

% Debug: Show what options were passed
fprintf('=== OPTIONS DEBUG ===\n');
if isempty(fieldnames(options))
    fprintf('No options provided, using defaults\n');
else
    disp(options);
end

% Set defaults
if ~isfield(options, 'exclude_grooming')
    options.exclude_grooming = false;
end

if ~isfield(options, 'figure_size') || isempty(options.figure_size) || length(options.figure_size) ~= 4
    options.figure_size = [100, 100, 1000, 500];
    fprintf('Using default figure size: [100, 100, 1000, 500]\n');
else
    fprintf('Using provided figure size: [%.0f, %.0f, %.0f, %.0f]\n', options.figure_size);
end

if ~isfield(options, 'x_col')
    options.x_col = 2;
end

if ~isfield(options, 'y_col')
    options.y_col = 3;
end

if isfield(options, 'time_limit') && ~isempty(options.time_limit)
    use_time_limit = true;
    time_limit_min = options.time_limit;
    time_limit_sec = time_limit_min * 60;
    fprintf('*** TIME LIMIT ENABLED: %.1f minutes (%.0f seconds) ***\n', time_limit_min, time_limit_sec);
else
    use_time_limit = false;
    fprintf('*** NO TIME LIMIT - analyzing full sessions ***\n');
end

if isfield(options, 'speed_threshold') && ~isempty(options.speed_threshold)
    use_speed_threshold = true;
    speed_threshold = options.speed_threshold;
    fprintf('*** SPEED THRESHOLD ENABLED: %.2f cm/s ***\n', speed_threshold);
else
    use_speed_threshold = false;
    fprintf('*** NO SPEED THRESHOLD - including all speeds ***\n');
end

%% Extract mouse information and filter by condition
mouse_names = mice(:,1);
saline_idx = contains(mouse_names, 'saline');
cno_idx = contains(mouse_names, 'CNO');

% Check if memory classification exists (column 7)
has_memory_classification = size(mice, 2) >= 7 && ~all(cellfun(@isempty, mice(:,7)));
if has_memory_classification
    memory_classifications = mice(:,7);
    strong_memory_idx = strcmp(memory_classifications, 'strong memory');
    weak_memory_idx = strcmp(memory_classifications, 'weak memory');
end

% Select mice based on condition
switch lower(condition)
    case 'saline'
        selected_idx = saline_idx;
        condition_name = 'Saline';
    case 'cno'
        selected_idx = cno_idx;
        condition_name = 'CNO';
    case 'strong'
        if ~has_memory_classification
            error('Memory classifications not found. Run classify_memory_strength first.');
        end
        selected_idx = strong_memory_idx;
        condition_name = 'Strong Memory';
    case 'weak'
        if ~has_memory_classification
            error('Memory classifications not found. Run classify_memory_strength first.');
        end
        selected_idx = weak_memory_idx;
        condition_name = 'Weak Memory';
    otherwise
        error('Invalid condition. Use: saline, CNO, strong, or weak');
end

selected_mice = find(selected_idx);
n_mice = length(selected_mice);

if n_mice == 0
    error('No mice found for condition: %s', condition);
end

fprintf('\nPlotting trajectories for %d mice in %s condition\n', n_mice, condition_name);

%% Plot trajectory for each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice{mouse_idx, 1};
    
    % Get pre and test data
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};
    food_pos = mice{mouse_idx, 2};  % Food position [x, y]
    
    fprintf('\n=== Mouse %d/%d: %s ===\n', i, n_mice, mouse_name);
    fprintf('  Original data - Pre: %d points, Test: %d points\n', size(pre_data,1), size(test_data,1));
    
    % Apply filters
    pre_data_filtered = apply_filters(pre_data, use_time_limit, time_limit_sec, ...
                                      use_speed_threshold, speed_threshold, ...
                                      options.exclude_grooming);
    test_data_filtered = apply_filters(test_data, use_time_limit, time_limit_sec, ...
                                       use_speed_threshold, speed_threshold, ...
                                       options.exclude_grooming);
    
    fprintf('  After filters - Pre: %d points, Test: %d points\n', ...
            size(pre_data_filtered,1), size(test_data_filtered,1));
    
    % Extract x, y coordinates
    x_pre = pre_data_filtered(:, options.x_col);
    y_pre = pre_data_filtered(:, options.y_col);
    x_test = test_data_filtered(:, options.x_col);
    y_test = test_data_filtered(:, options.y_col);
    
    % Create figure for this mouse
    try
        fig = figure('Name', sprintf('Trajectory: %s', mouse_name), ...
                     'Position', options.figure_size, 'Color', 'white');
    catch
        fig = figure('Name', sprintf('Trajectory: %s', mouse_name), 'Color', 'white');
        fprintf('Warning: Could not set figure position, using default\n');
    end
    
    % Plot PRE session (left subplot)
    subplot(1, 2, 1);
    plot(x_pre, y_pre, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    hold on;
    
    % Mark start and end points
    if ~isempty(x_pre)
        scatter(x_pre(1), y_pre(1), 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        scatter(x_pre(end), y_pre(end), 100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    % Mark food position
    if ~isempty(food_pos) && length(food_pos) >= 2
        scatter(food_pos(1), food_pos(2), 150, [1 0.6 0], 'p', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end
    
    axis equal;
    grid off;
    box off;
    axis off;
    xlabel('X Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Before Session', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'Trajectory', 'Start', 'End', 'Food'}, 'Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
    hold off;
    
    % Plot TEST session (right subplot)
    subplot(1, 2, 2);
    plot(x_test, y_test, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1);
    hold on;
    
    % Mark start and end points
    if ~isempty(x_test)
        scatter(x_test(1), y_test(1), 100, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        scatter(x_test(end), y_test(end), 100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    % Mark food position
    if ~isempty(food_pos) && length(food_pos) >= 2
        scatter(food_pos(1), food_pos(2), 150, [1 0.6 0], 'p', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 2);
    end
    
    axis equal;
    grid off;
    box off;
    axis off;
    xlabel('X Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y Position (cm)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Test Session', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'Trajectory', 'Start', 'End', 'Food'}, 'Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
    hold off;
    
    % Add overall title with filter information
    title_text = sprintf('%s - %s', mouse_name, condition_name);
    filter_parts = {};
    if use_time_limit
        filter_parts{end+1} = sprintf('%.1f min', time_limit_min);
    end
    if use_speed_threshold
        filter_parts{end+1} = sprintf('Speed >= %.1f cm/s', speed_threshold);
    end
    if options.exclude_grooming
        filter_parts{end+1} = 'No grooming';
    end
    if ~isempty(filter_parts)
        title_text = sprintf('%s (%s)', title_text, strjoin(filter_parts, ', '));
    end
    
    sgtitle(title_text, 'FontSize', 16, 'FontWeight', 'bold');
end

fprintf('\n=== Trajectory plotting complete ===\n');
fprintf('Generated %d figures (one per mouse)\n', n_mice);

end

%% Helper function to apply filters
function filtered_data = apply_filters(data, use_time_limit, time_limit_sec, ...
                                       use_speed_threshold, speed_threshold, ...
                                       exclude_grooming)
    filtered_data = data;
    
    % Apply time limit if specified
    if use_time_limit
        time = filtered_data(:, 1);
        start_time = min(time);
        time_mask = (time - start_time) <= time_limit_sec;
        filtered_data = filtered_data(time_mask, :);
    end
    
    % Apply speed threshold if specified
    if use_speed_threshold
        speed = filtered_data(:, 7);
        speed_mask = speed >= speed_threshold;
        filtered_data = filtered_data(speed_mask, :);
    end
    
    % Apply grooming exclusion if specified
    if exclude_grooming
        grooming = filtered_data(:, 10);
        non_grooming_mask = grooming == 0;
        filtered_data = filtered_data(non_grooming_mask, :);
    end
end