function plotMouseCheeseboardData(mice, mouse_name, options)
    % Plot df/f vs distance for a single mouse in cheeseboard maze (pre-test and test sessions)
    % Each row contains: heatmap, time series, and correlation plot
    %
    % Parameters:
    %   mice: cell array with mouse data structure:
    %     mice{i,1} - mouse name
    %     mice{i,2} - [x,y] food position
    %     mice{i,3} - pre-test data (Nx10 double)
    %     mice{i,4} - test data (Nx10 double)
    %     mice{i,5} - [x0,y0] starting position
    %     mice{i,6} - latencies [pre;test]
    %     mice{i,7} - 'strong memory' or 'weak memory'
    %     mice{i,8} - food zone change (test - pre) as percentage
    %
    %   mouse_name: string with mouse ID (e.g., 'FDRE14_saline', 'FDRE15_CNO')
    %   options: struct with options for filtering
    %     - dist_limit: minimum distance threshold (default: 0)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - smooth_window: window size for time series smoothing (default: 20)
    %     - colormap: colormap for heatmap (default: parula)
    %     - sigma: sigma for Gaussian smoothing (default: 2)
    %     - resolution: grid resolution for heatmap (default: 100)
    %     - time_xlim: x-axis limits for time series plots [min max] (optional)
    %     - group_filter: filter by experimental group ('saline', 'CNO', or 'all') (default: 'all')
    %     - dff_ylim: manual y-axis limits for dF/F [min max] (optional)
    %     - dist_ylim: manual y-axis limits for distance [min max] (optional)
    
    % Define column indices for cheeseboard data
    COL_TIME = 1;     % Time
    COL_X = 2;        % X position
    COL_Y = 3;        % Y position
    COL_465 = 4;      % 465nm signal
    COL_405 = 5;      % 405nm signal
    COL_DFF = 6;      % dF/F
    COL_SPEED = 7;    % Speed
    COL_ZONES = 8;    % Zones (0-outside, 1-food, 2-area2, 3-area3)
    COL_DISTANCE = 9; % Distance to food
    COL_GROOM = 10;   % Grooming
    
    % Set default options if not provided
    if nargin < 3
        options = struct();
    end
    
    if ~isfield(options, 'dist_limit')
        options.dist_limit = 0;
    end
    
    if ~isfield(options, 'remove_grooming')
        options.remove_grooming = true;
    end
    
    if ~isfield(options, 'smooth_window')
        options.smooth_window = 20;
    end
    
    if ~isfield(options, 'colormap')
        options.colormap = parula;
    end
    
    if ~isfield(options, 'sigma')
        options.sigma = 2;
    end
    
    if ~isfield(options, 'resolution')
        options.resolution = 100;
    end
    
    if ~isfield(options, 'group_filter')
        options.group_filter = 'all';
    end
    
    % Manual y-axis limits for time series plots
    if ~isfield(options, 'dff_ylim')
        options.dff_ylim = []; % Auto if empty, or [min, max] for manual
    end

    if ~isfield(options, 'dist_ylim')
        options.dist_ylim = []; % Auto if empty, or [min, max] for manual
    end
    
    % Validate group option
    valid_groups = {'all', 'saline', 'CNO'};
    if ~ismember(options.group_filter, valid_groups)
        error('Invalid group_filter. Must be one of: ''all'', ''saline'', or ''CNO''');
    end
    
    % Find the mouse in the data
    mouse_idx = [];
    for i = 1:size(mice, 1)
        if strcmp(mice{i,1}, mouse_name)
            % Check group filter
            if ~strcmp(options.group_filter, 'all')
                if strcmp(options.group_filter, 'saline') && ~contains(mouse_name, 'saline')
                    continue;
                elseif strcmp(options.group_filter, 'CNO') && ~contains(mouse_name, 'CNO')
                    continue;
                end
            end
            mouse_idx = i;
            break;
        end
    end
    
    if isempty(mouse_idx)
        if strcmp(options.group_filter, 'all')
            error('Mouse not found: %s', mouse_name);
        else
            error('Mouse not found: %s in group: %s', mouse_name, options.group_filter);
        end
    end
    
    % Extract mouse data
    food_pos = mice{mouse_idx, 2};
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};
    start_pos = mice{mouse_idx, 5};
    
    % Calculate latencies directly from data (discovery distance = 0.5 cm)
    discovery_distance = 0.5; % cm
    latencies = calculate_latencies(pre_data, test_data, discovery_distance);
    
    % Get additional info if available
    memory_strength = '';
    food_zone_change = '';
    if size(mice, 2) >= 7 && ~isempty(mice{mouse_idx, 7})
        memory_strength = mice{mouse_idx, 7};
    end
    if size(mice, 2) >= 8 && ~isempty(mice{mouse_idx, 8})
        food_zone_change = sprintf('%.1f%%', mice{mouse_idx, 8});
    end
    
    % Determine treatment group from mouse name
    if contains(mouse_name, 'saline')
        treatment_group = 'Saline';
    elseif contains(mouse_name, 'CNO')
        treatment_group = 'CNO';
    else
        treatment_group = 'Unknown';
    end
    
    % Store session data
    session_data = {pre_data, test_data};
    session_names = {'Pre-Test', 'Test'};
    num_sessions = 2;
    
    % Process data - remove grooming if requested
    processed_data = cell(1, num_sessions);
    for s = 1:num_sessions
        data = session_data{s};
        
        % Remove grooming periods if requested
        if options.remove_grooming && size(data, 2) >= COL_GROOM
            data = data(data(:, COL_GROOM) == 0, :);
        end
        
        % Check if we have enough data after filtering
        if size(data, 1) < 10
            warning('Very few data points remaining for %s session after filtering', session_names{s});
        end
        
        processed_data{s} = data;
    end
    
    fprintf('Plotting mouse %s (%s)\n', mouse_name, treatment_group);
    fprintf('Food position: [%.1f, %.1f] cm\n', food_pos(1), food_pos(2));
    
    % Print latencies with error handling
    if ~isnan(latencies(1)) && ~isnan(latencies(2))
        fprintf('Pre-test latency: %.1f s, Test latency: %.1f s\n', latencies(1), latencies(2));
    elseif ~isnan(latencies(1))
        fprintf('Pre-test latency: %.1f s, Test latency: N/A\n', latencies(1));
    elseif ~isnan(latencies(2))
        fprintf('Pre-test latency: N/A, Test latency: %.1f s\n', latencies(2));
    else
        fprintf('Pre-test latency: N/A, Test latency: N/A\n');
    end
    
    if options.remove_grooming
        fprintf('Grooming periods: REMOVED from all analyses\n');
    else
        fprintf('Grooming periods: INCLUDED in analyses\n');
    end
    if ~isempty(memory_strength)
        fprintf('Memory strength: %s\n', memory_strength);
    end
    if ~isempty(food_zone_change)
        fprintf('Food zone change: %s\n', food_zone_change);
    end
    
    % Create figure with a row for each session
    row_height = 250; % pixels per row
    title_height = 80; % pixels for title
    bottom_margin = 10; % Additional pixels for bottom margin
    figure('Position', [50, 50, 1400, row_height*num_sessions+title_height+bottom_margin]);

    % Calculate global color limits for heatmaps for consistency
    global_heatmap_min = Inf;
    global_heatmap_max = -Inf;

    % Calculate global y-axis limits for time series consistency
    global_dff_min = Inf;
    global_dff_max = -Inf;
    global_dist_min = Inf;
    global_dist_max = -Inf;

    % Find global limits across both sessions
    for s = 1:num_sessions
        data = processed_data{s};
        if isempty(data), continue; end

        dff_vals = data(:, COL_DFF);
        dist_vals = data(:, COL_DISTANCE);

        global_dff_min = min(global_dff_min, min(dff_vals));
        global_dff_max = max(global_dff_max, max(dff_vals));
        global_dist_min = min(global_dist_min, min(dist_vals));
        global_dist_max = max(global_dist_max, max(dist_vals));
    end

    % Calculate global heatmap limits
    for s = 1:num_sessions
        data = processed_data{s};
        if isempty(data), continue; end

        % Extract position and dff data
        x_pos = data(:, COL_X);
        y_pos = data(:, COL_Y);
        dff_all = data(:, COL_DFF);
        
        % Create a grid for visualization
        resolution = options.resolution;
        x_min = min(x_pos);
        x_max = max(x_pos);
        y_min = min(y_pos);
        y_max = max(y_pos);
        
        % Initialize grid for averaged values and visit count
        dff_grid = zeros(resolution, resolution);
        visit_count = zeros(resolution, resolution);
        
        % Calculate which grid cell each data point belongs to
        x_idx = round((x_pos - x_min) / (x_max - x_min) * (resolution-1)) + 1;
        y_idx = round((y_pos - y_min) / (y_max - y_min) * (resolution-1)) + 1;
        
        % Enforce valid indices
        x_idx = max(1, min(resolution, x_idx));
        y_idx = max(1, min(resolution, y_idx));
        
        % Sum dff values and count visits for each grid cell
        for i = 1:length(x_pos)
            xi = x_idx(i);
            yi = y_idx(i);
            dff_grid(yi, xi) = dff_grid(yi, xi) + dff_all(i);
            visit_count(yi, xi) = visit_count(yi, xi) + 1;
        end
        
        % Calculate mean dff for each cell
        visited_mask = visit_count > 0;
        mean_dff_grid = zeros(size(dff_grid));
        mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
        
        % Create a masked array with NaNs for unvisited locations
        masked_dff = mean_dff_grid;
        masked_dff(~visited_mask) = NaN;
        
        % Apply Gaussian smoothing to the grid
        smoothed_dff = smoothHeatmap(masked_dff, options.sigma);
        
        % Update global limits for consistent color scale
        valid_values = smoothed_dff(~isnan(smoothed_dff));
        if ~isempty(valid_values)
            global_heatmap_min = min(global_heatmap_min, min(valid_values));
            global_heatmap_max = max(global_heatmap_max, max(valid_values));
        end
    end
    
    % Ensure symmetric color scale for heatmaps
    if ~isinf(global_heatmap_min) && ~isinf(global_heatmap_max)
        heatmap_limit = max(abs(global_heatmap_min), abs(global_heatmap_max));
        global_clim = [-heatmap_limit, heatmap_limit];
    else
        global_clim = [-1, 1]; % Default range
    end
    
    % Layout parameters
    heatmap_width = 0.25;
    time_width = 0.45;
    corr_width = 0.15;
    padding = 0.03;
    
    % Calculate positioning
    total_height = row_height*num_sessions+title_height+bottom_margin;
    bottom_margin_frac = bottom_margin / total_height;
    available_height = 1 - bottom_margin_frac - 0.05;
    row_frac = available_height / num_sessions;
    row_height_frac = row_frac * 0.8;
    
    % Process each session and create plots
    for s = 1:num_sessions
        session_name = session_names{s};
        data = processed_data{s};
        
        if isempty(data)
            fprintf('Warning: No data available for %s session\n', session_name);
            continue;
        end
        
        % Define subplot positions for this row
        row_top = 1 - 0.05 - (s-1) * row_frac;
        row_bottom = row_top - row_height_frac;
        
        pos_heatmap = [0.02, row_bottom, heatmap_width, row_height_frac];
        pos_time = [0.05+heatmap_width+padding-0.032, row_bottom, time_width, row_height_frac];
        pos_corr = [0.05+heatmap_width+padding+time_width+padding, row_bottom, corr_width, row_height_frac];
        
        %% HEATMAP PLOT
        heatmap_axes = subplot('Position', pos_heatmap);
        
        % Extract position and dff data
        x_pos = data(:, COL_X);
        y_pos = data(:, COL_Y);
        dff_all = data(:, COL_DFF);
        
        % Create heatmap (same process as above)
        resolution = options.resolution;
        x_min = min(x_pos);
        x_max = max(x_pos);
        y_min = min(y_pos);
        y_max = max(y_pos);
        x_edges = linspace(x_min, x_max, resolution);
        y_edges = linspace(y_min, y_max, resolution);
        
        % Initialize grid
        dff_grid = zeros(resolution, resolution);
        visit_count = zeros(resolution, resolution);
        
        % Calculate grid indices
        x_idx = round((x_pos - x_min) / (x_max - x_min) * (resolution-1)) + 1;
        y_idx = round((y_pos - y_min) / (y_max - y_min) * (resolution-1)) + 1;
        x_idx = max(1, min(resolution, x_idx));
        y_idx = max(1, min(resolution, y_idx));
        
        % Fill grid
        for i = 1:length(x_pos)
            xi = x_idx(i);
            yi = y_idx(i);
            dff_grid(yi, xi) = dff_grid(yi, xi) + dff_all(i);
            visit_count(yi, xi) = visit_count(yi, xi) + 1;
        end
        
        % Calculate mean and smooth
        visited_mask = visit_count > 0;
        mean_dff_grid = zeros(size(dff_grid));
        mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
        masked_dff = mean_dff_grid;
        masked_dff(~visited_mask) = NaN;
        smoothed_dff = smoothHeatmap(masked_dff, options.sigma);
        
        % Plot heatmap
        imagesc(x_edges, y_edges, smoothed_dff, 'AlphaData', ~isnan(smoothed_dff));
        hold on;
        
        % Add food location
        plot(food_pos(1), food_pos(2), 'mo', 'MarkerSize', 12, 'LineWidth', 3, 'MarkerFaceColor', 'm');
        
        % % Add start position
        % plot(start_pos(1), start_pos(2), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'green');
        
        % Format heatmap
        axis equal tight
        colormap(heatmap_axes, options.colormap);
        c = colorbar;
        caxis(global_clim);
        c.Label.String = 'Mean \DeltaF/F';
        c.Label.FontSize = 10;
        xlabel('X Position (cm)', 'FontSize', 10);
        ylabel('Y Position (cm)', 'FontSize', 10);
        title([session_name, ': Mean dF/F Heatmap'], 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
        box off
        axis off
        
        % % Add legend for markers
        % if s == 1 % Only on first plot to avoid clutter
        %     legend_entries = {'Food', 'Start'};
        %     legend('boxoff', 'Location', 'northeast', 'FontSize', 8);
        % end
        
        %% TIME SERIES PLOT
        time_axes = subplot('Position', pos_time);
        
        % Filter data based on distance threshold
        distances = data(:, COL_DISTANCE);
        valid_indices = distances > options.dist_limit;
        
        % Calculate time vector using filtered data (grooming already removed from processed_data)
        % Since grooming periods are removed, we need to create a continuous time vector
        time_data = data(:, COL_TIME);
        if options.remove_grooming
            % Create continuous time vector by calculating cumulative time differences
            time_diffs = [0; diff(time_data)];
            % Cap large time gaps (where grooming was removed) to median interval
            median_interval = median(time_diffs(time_diffs > 0));
            time_diffs(time_diffs > 5 * median_interval) = median_interval;
            time_vector = cumsum(time_diffs);
        else
            time_vector = time_data - time_data(1); % Start from 0
        end
        
        % Smooth data for plotting (using the already filtered data)
        smooth_dff = movmean(data(:, COL_DFF), options.smooth_window);
        smooth_dist = movmean(distances, options.smooth_window);
        
        % Plot dF/F vs time (left y-axis)
        yyaxis left;
        plot(time_vector, smooth_dff, 'LineWidth', 1.5, 'Color', [0, 0.4, 0.8]);
        ylabel('\DeltaF/F', 'FontSize', 10, 'Color', [0, 0.4, 0.8]);
        ax = gca;
        ax.YAxis(1).Color = [0, 0.4, 0.8];
        
        % Plot distance vs time (right y-axis)
        yyaxis right;
        plot(time_vector, smooth_dist, 'LineWidth', 1.5, 'Color', [0.9, 0.3, 0.5]);
        ylabel('Distance to Food (cm)', 'FontSize', 10, 'Color', [0.9, 0.3, 0.5]);
        ax.YAxis(2).Color = [0.9, 0.3, 0.5];
        
        % Format time series plot
        xlabel('Time (s)', 'FontSize', 10);
        
        % Set y-axis limits (manual override or consistent auto limits)
        yyaxis left;
        if ~isempty(options.dff_ylim) && length(options.dff_ylim) == 2
            ylim(options.dff_ylim); % Use manual limits
        elseif ~isinf(global_dff_min) && ~isinf(global_dff_max)
            ylim([global_dff_min, global_dff_max]); % Use auto limits
        end

        yyaxis right;
        if ~isempty(options.dist_ylim) && length(options.dist_ylim) == 2
            ylim(options.dist_ylim); % Use manual limits
        elseif ~isinf(global_dist_min) && ~isinf(global_dist_max)
            ylim([global_dist_min, global_dist_max]); % Use auto limits
        end

        if options.remove_grooming
            if ~isnan(latencies(s))
                title_str = sprintf('%s: dF/F and Distance vs Time (Latency: %.1fs, Grooming Removed)', session_name, latencies(s));
            else
                title_str = sprintf('%s: dF/F and Distance vs Time (Latency: N/A, Grooming Removed)', session_name);
            end
        else
            if ~isnan(latencies(s))
                title_str = sprintf('%s: dF/F and Distance vs Time (Latency: %.1fs)', session_name, latencies(s));
            else
                title_str = sprintf('%s: dF/F and Distance vs Time (Latency: N/A)', session_name);
            end
        end
        title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box off
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
        
        % Apply custom time limits if provided
        if isfield(options, 'time_xlim') && ~isempty(options.time_xlim) && length(options.time_xlim) == 2
            xlim(options.time_xlim);
        end
        
        %% CORRELATION SCATTER PLOT
        corr_axes = subplot('Position', pos_corr);
        
        % Extract valid data for correlation
        dff_valid = data(valid_indices, COL_DFF);
        dist_valid = distances(valid_indices);
        
        % Calculate correlation
        if length(dff_valid) > 1 && length(dist_valid) > 1
            [rho, pval] = corr(dist_valid, dff_valid, 'Type', 'Pearson');
        else
            rho = NaN;
            pval = NaN;
        end
        
        % Plot scatter
        scatter(dist_valid, dff_valid, 25, 'filled', 'MarkerFaceColor', [0.5, 0.7, 0.9], 'MarkerFaceAlpha', 0.7);
        hold on;
        
        % Add regression line
        if length(dist_valid) > 1 && ~isnan(rho)
            p = polyfit(dist_valid, dff_valid, 1);
            x_range = linspace(min(dist_valid), max(dist_valid), 100);
            y_fit = polyval(p, x_range);
            plot(x_range, y_fit, 'LineWidth', 2, 'Color', [0.3, 0.5, 0.7]);
        end
        
        % Add correlation coefficient
        if ~isnan(rho)
            text(0.05, 0.95, sprintf('r = %.3f\np = %.3f', rho, pval), ...
                 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        else
            text(0.05, 0.95, 'r = N/A', ...
                 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
        
        % Format correlation plot
        xlabel('Distance to Food (cm)', 'FontSize', 10);
        ylabel('\DeltaF/F', 'FontSize', 10);
        title('Distance-dF/F Correlation', 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box off
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
    end
    
    % Add main title for the whole figure
    axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
    
    % Create comprehensive title - fix underscore issue
    clean_mouse_name = strrep(mouse_name, '_', '\_'); % Escape underscores for LaTeX
    title_parts = {sprintf('Mouse %s', clean_mouse_name)};

    if ~isempty(memory_strength)
        title_parts{end+1} = memory_strength;
    end
    
    if ~isempty(food_zone_change)
        title_parts{end+1} = sprintf('Time increase: %s', food_zone_change);
    end
    
    title_str = strjoin(title_parts, ' | ');
    
    text(0.5, 0, title_str, 'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
end

% Helper function to calculate latencies
function latencies = calculate_latencies(pre_data, test_data, discovery_distance)
    % Calculate latencies for pre-test and test sessions
    
    % Initialize default latencies
    latencies = [NaN, NaN];
    
    % Check if data exists and has required columns
    if isempty(pre_data) || size(pre_data, 2) < 9
        fprintf('Warning: Pre-test data missing or insufficient columns\n');
        latencies(1) = NaN;
    else
        % Pre-test latency
        pre_distance = pre_data(:, 9);  % Distance column
        pre_time = pre_data(:, 1);      % Time column
        
        if ~isempty(pre_distance) && ~isempty(pre_time)
            pre_discovery_idx = find(pre_distance <= discovery_distance, 1, 'first');
            if ~isempty(pre_discovery_idx)
                latencies(1) = pre_time(pre_discovery_idx) - pre_time(1);
            else
                latencies(1) = max(pre_time) - min(pre_time); % Session duration if not discovered
            end
        else
            fprintf('Warning: Pre-test distance or time data is empty\n');
            latencies(1) = NaN;
        end
    end
    
    % Check if test data exists and has required columns
    if isempty(test_data) || size(test_data, 2) < 9
        fprintf('Warning: Test data missing or insufficient columns\n');
        latencies(2) = NaN;
    else
        % Test latency
        test_distance = test_data(:, 9);  % Distance column
        test_time = test_data(:, 1);      % Time column
        
        if ~isempty(test_distance) && ~isempty(test_time)
            test_discovery_idx = find(test_distance <= discovery_distance, 1, 'first');
            if ~isempty(test_discovery_idx)
                latencies(2) = test_time(test_discovery_idx) - test_time(1);
            else
                latencies(2) = max(test_time) - min(test_time); % Session duration if not discovered
            end
        else
            fprintf('Warning: Test distance or time data is empty\n');
            latencies(2) = NaN;
        end
    end
end

% Helper function to smooth a heatmap while preserving NaN values
function smoothed_data = smoothHeatmap(data, sigma)
    % Create a Gaussian kernel
    kernel_size = 7; % Must be odd
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*sigma^2));
    kernel = kernel / sum(kernel(:));
    
    % Apply the smoothing while preserving the mask
    smoothed_data = NaN(size(data));
    valid_mask = ~isnan(data);
    
    % Pad the arrays to handle edges
    padded_data = padarray(data, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);
    
    % Convolve with the kernel
    for i = 1+floor(kernel_size/2):size(padded_data,1)-floor(kernel_size/2)
        for j = 1+floor(kernel_size/2):size(padded_data,2)-floor(kernel_size/2)
            if padded_mask(i, j)
                window = padded_data(i-floor(kernel_size/2):i+floor(kernel_size/2), j-floor(kernel_size/2):j+floor(kernel_size/2));
                window_mask = ~isnan(window);
                if sum(window_mask(:)) > 0
                    weighted_values = window .* kernel;
                    weighted_values(~window_mask) = 0;
                    % Normalize by the valid weights
                    weight_sum = sum(kernel(window_mask));
                    if weight_sum > 0
                        smoothed_data(i-floor(kernel_size/2), j-floor(kernel_size/2)) = sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end
end