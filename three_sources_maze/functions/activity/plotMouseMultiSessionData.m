function plotMouseMultiSessionData(mice_all, mouse_name, arm_type, options)
    % Plot df/f vs distance for a single mouse in all sessions (arranged in rows)
    % Each row contains: heatmap, time series, and correlation plot
    %
    % Parameters:
    %   mice_all: cell array with mouse data
    %   mouse_name: string with mouse ID (e.g., 'F13', 'F14')
    %   arm_type: 'food', 'nonfood1', or 'nonfood2'
    %   options: struct with options for filtering
    %     - dist_limit: minimum distance threshold (default: 0)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - smooth_window: window size for time series smoothing (default: 20)
    %     - colormap: colormap for heatmap (default: parula)
    %     - sigma: sigma for Gaussian smoothing (default: 2)
    %     - resolution: grid resolution for heatmap (default: 100)
    %     - time_xlim: x-axis limits for time series plots [min max] (optional)
    %     - group: filter by experimental group ('control', 'CNO', 'chrimson', or 'all') (default: 'all')
    
    % Define column indices for data
    COL_TIME = 1;     % Time
    COL_X = 2;        % X position
    COL_Y = 3;        % Y position
    COL_DFF = 6;      % dF/F
    COL_SPEED = 7;    % Speed
    COL_ZONES = 8;    % Zones
    COL_DIST_FOOD = 9;       % Distance to food arm
    COL_DIST_NONFOOD2 = 10;  % Distance to non-food arm 2
    COL_DIST_NONFOOD1 = 11;  % Distance to non-food arm 1
    COL_DOOR = 12;    % Door status (0=closed, 1=open, 2=closed again)
    COL_GROOM = 13;   % Grooming
    
    % Set default options if not provided
    if nargin < 4
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
    
    if ~isfield(options, 'group')
        options.group = 'all';
    end
    
    % Validate arm type
    valid_arms = {'food', 'nonfood1', 'nonfood2'};
    if ~ismember(arm_type, valid_arms)
        error('Invalid arm_type. Must be one of: ''food'', ''nonfood1'', or ''nonfood2''');
    end
    
    % Validate group option
    valid_groups = {'all', 'control', 'CNO', 'chrimson'};
    if ~ismember(options.group, valid_groups)
        error('Invalid group. Must be one of: ''all'', ''control'', ''CNO'', or ''chrimson''');
    end
    
    % Determine distance column based on arm_type
    if strcmp(arm_type, 'food')
        dist_col = COL_DIST_FOOD;
        arm_label = 'Food Arm';
    elseif strcmp(arm_type, 'nonfood1')
        dist_col = COL_DIST_NONFOOD1;
        arm_label = 'Non-Food Arm 1';
    else % nonfood2
        dist_col = COL_DIST_NONFOOD2;
        arm_label = 'Non-Food Arm 2';
    end
    
    % Find all sessions matching this mouse name and group
    session_indices = {};
    session_data = {};
    mouse_id = mouse_name;
    group = '';
    food_arm_location = '';
    memory_strength = '';
    
    % First, collect all sessions for this mouse that match the group filter
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        current_group = mice_all{i, 2};
        
        % Check if this row belongs to our mouse
        if contains(session_info, [mouse_name, '_sess'])
            % Skip if group filter is active and doesn't match
            if ~strcmp(options.group, 'all') && ~strcmp(current_group, options.group)
                continue;
            end
            
            % Extract session type (sess0, sess1, sess2)
            if contains(session_info, 'sess0')
                session_type = 'sess0';
            elseif contains(session_info, 'sess1')
                session_type = 'sess1';
            elseif contains(session_info, 'sess2')
                session_type = 'sess2';
            else
                continue; % Skip if not a recognized session
            end
            
            % Store group and food arm info (same for all sessions of this mouse)
            if isempty(group)
                group = current_group;
                food_arm_location = mice_all{i, 3};
                if size(mice_all, 2) >= 6
                    memory_strength = mice_all{i, 6};
                end
            end
            
            % Process the data for this session
            raw_data = mice_all{i, 4};
            
            % Process session data based on session type
            if strcmp(session_type, 'sess2')
                % For session 2, keep only closed door condition
                filtered_data = raw_data(raw_data(:,COL_DOOR) == 0, :);
            elseif strcmp(session_type, 'sess1')
                % For session 1, make sure we're only using the door=1 (open) condition
                filtered_data = raw_data(raw_data(:,COL_DOOR) == 1, :);
                
                % Fix time by making it relative to the start
                if ~isempty(filtered_data)
                    filtered_data(:,COL_TIME) = filtered_data(:,COL_TIME) - filtered_data(1,COL_TIME);
                end
            elseif strcmp(session_type, 'sess0')
                % For session 0, make sure we're only using the door=0 (closed) condition
                filtered_data = raw_data(raw_data(:,COL_DOOR) == 0, :);
            end
            
            % Remove grooming periods if requested
            if options.remove_grooming && ~isempty(filtered_data)
                filtered_data = filtered_data(filtered_data(:,COL_GROOM) == 0, :);
            end
            
            % Only keep the session if we have data after filtering
            if ~isempty(filtered_data) && size(filtered_data, 1) > 10
                session_indices{end+1} = session_type;
                session_data{end+1} = filtered_data;
            end
        end
    end

    % Check if we found any sessions
    if isempty(session_indices)
        if strcmp(options.group, 'all')
            error('No valid sessions found for mouse: %s', mouse_name);
        else
            error('No valid sessions found for mouse: %s in group: %s', mouse_name, options.group);
        end
    end

    % Sort sessions in order (sess0, sess1, sess2)
    [sorted_sessions, sort_idx] = sort(session_indices);
    sorted_data = session_data(sort_idx);

    % Debug info
    if strcmp(options.group, 'all')
        fprintf('Found %d sessions for mouse %s: %s\n', length(sorted_sessions), mouse_name, strjoin(sorted_sessions, ', '));
    else
        fprintf('Found %d sessions for mouse %s (group: %s): %s\n', length(sorted_sessions), mouse_name, options.group, strjoin(sorted_sessions, ', '));
    end

    % Create figure with a row for each session
    num_sessions = length(sorted_sessions);

    % Increase row height and add extra space at the bottom
    row_height = 250; % pixels per row
    title_height = 80; % pixels for title
    bottom_margin = 10; % Additional pixels for bottom margin
    figure('Position', [50, 50, 1400, row_height*num_sessions+title_height+bottom_margin]);

    % Define global color limits for heatmaps for consistency
    global_heatmap_min = Inf;
    global_heatmap_max = -Inf;

    % Calculate global heatmap limits - no need to filter for grooming again here
    % as filtering is already applied to sorted_data
    for s = 1:num_sessions
        data = sorted_data{s};

        % Extract x,y position and dff data
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

        % Calculate mean dff for each cell and create a mask for visited locations
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
    heatmap_limit = max(abs(global_heatmap_min), abs(global_heatmap_max));
    global_clim = [-heatmap_limit, heatmap_limit];

    % Adjust the layout - make the time series plot wider and correlation plot narrower
    heatmap_width = 0.25;
    time_width = 0.45; % Increased from 0.3
    corr_width = 0.10; % Decreased from 0.25

    % More space between plots
    padding = 0.05;

    % Adjust the vertical positioning to account for bottom margin
    total_height = row_height*num_sessions+title_height+bottom_margin;
    bottom_margin_frac = bottom_margin / total_height; % Bottom margin as fraction of figure height

    % Process each session and create plots
    for s = 1:num_sessions
        session_type = sorted_sessions{s};
        data = sorted_data{s};
        
        % Get formatted session name
        formatted_session = strrep(session_type, 'sess', 'Session ');

        % Define subplot positions for this row with increased spacing between rows
        % Use absolute positioning based on normalized coordinates
        % Define subplot positions for this row with proper spacing
        available_height = 1 - bottom_margin_frac - 0.05; % Available height (minus bottom margin and title space)
        row_frac = available_height / num_sessions; % Height of each row as fraction of available space
        row_height_frac = row_frac * 0.8; % 80% of row height for the plots
        row_top = 1 - 0.05 - (s-1) * row_frac; % Position from top
        row_bottom = row_top - row_height_frac; % Bottom of the plot

        % Position the plots in this row
        pos_heatmap = [0.05, row_bottom, heatmap_width, row_height_frac];
        pos_time = [0.05+heatmap_width, row_bottom, time_width, row_height_frac];
        pos_corr = [0.05+heatmap_width+padding+time_width+padding, row_bottom, corr_width, row_height_frac];
        
        %% HEATMAP PLOT
        heatmap_axes = subplot('Position', pos_heatmap);
        
        % Extract x,y position and dff data
        x_pos = data(:, COL_X);
        y_pos = data(:, COL_Y);
        dff_all = data(:, COL_DFF);
        
        % Create a grid for visualization
        resolution = options.resolution;
        x_min = min(x_pos);
        x_max = max(x_pos);
        y_min = min(y_pos);
        y_max = max(y_pos);
        x_edges = linspace(x_min, x_max, resolution);
        y_edges = linspace(y_min, y_max, resolution);
        
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
        
        % Calculate mean dff for each cell and create a mask for visited locations
        visited_mask = visit_count > 0;
        mean_dff_grid = zeros(size(dff_grid));
        mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
        
        % Create a masked array with NaNs for unvisited locations
        masked_dff = mean_dff_grid;
        masked_dff(~visited_mask) = NaN;
        
        % Apply Gaussian smoothing to the grid
        smoothed_dff = smoothHeatmap(masked_dff, options.sigma);
        
        % Plot the smoothed heatmap
        imagesc(x_edges, y_edges, smoothed_dff, 'AlphaData', ~isnan(smoothed_dff));
        
        % Format heatmap
        axis equal tight
        colormap(heatmap_axes, options.colormap);
        c = colorbar;
        caxis(global_clim); % Use consistent color limits
        c.Label.String = 'Mean \Delta F/F';
        c.Label.FontSize = 10;
        xlabel('X Position (cm)', 'FontSize', 10);
        ylabel('Y Position (cm)', 'FontSize', 10);
        title([formatted_session, ': Mean DFF Heatmap'], 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
        box off
        axis off
        
        %% TIME SERIES PLOT
        time_axes = subplot('Position', pos_time);
        
        % Filter data based on options
        if options.remove_grooming
            % Remove rows with grooming - not needed, already filtered in the session data
            valid_indices = data(:, dist_col) > options.dist_limit;
        else
            % Only filter by distance threshold
            valid_indices = data(:, dist_col) > options.dist_limit;
        end
        
        % Calculate frame rate from time data
        frame_rate = 1/mean(diff(data(:, COL_TIME))); % Hz
        time_vector = (0:size(data, 1)-1) / frame_rate;
        
        % Smooth dF/F and distance data for clarity
        smooth_dff = movmean(data(:, COL_DFF), options.smooth_window);
        smooth_dist = movmean(data(:, dist_col), options.smooth_window);
        
        % Plot data vs time (left y-axis)
        yyaxis left;
        plot(time_vector, smooth_dff, 'LineWidth', 1.5, 'Color', [0, 0.4, 0.8]);
        ylabel('\Delta F/F', 'FontSize', 10, 'Color', [0, 0.4, 0.8]);
        ylim([-2,3]);
        % Set y-axis color
        ax = gca;
        ax.YAxis(1).Color = [0, 0.4, 0.8];
        
        % Plot data vs time (right y-axis)
        yyaxis right;
        plot(time_vector, smooth_dist, 'LineWidth', 1.5, 'Color', [0.9, 0.3, 0.5]);
        ylabel(['Distance (cm)'], 'FontSize', 10, 'Color', [0.9, 0.3, 0.5]);
        % Set y-axis color
        ax.YAxis(2).Color = [0.9, 0.3, 0.5];
        
        % Format time series plot
        xlabel('Time (s)', 'FontSize', 10);
        title([formatted_session, ': dF/F and Distance vs Time'], 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box off;
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
        
        % Apply custom time limits if provided
        if isfield(options, 'time_xlim') && ~isempty(options.time_xlim) && length(options.time_xlim) == 2
            xlim(options.time_xlim);
        end
        
        %% CORRELATION SCATTER PLOT
        corr_axes = subplot('Position', pos_corr);
        
        % Extract valid data for correlation
        dff_valid = data(valid_indices, COL_DFF);
        dist_valid = data(valid_indices, dist_col);
        
        % Calculate correlation
        [rho, pval] = corr(dist_valid, dff_valid, 'Type', 'Pearson');
        
        % Plot scatter
        scatter(dist_valid, dff_valid, 25, 'filled', 'MarkerFaceColor', [0.5, 0.7, 0.9], 'MarkerFaceAlpha', 0.7);
        hold on;
        
        % Add regression line
        if length(dist_valid) > 1
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
        end
        
        % Format correlation plot
        xlabel(['Distance (cm)'], 'FontSize', 10);
        ylabel('\Delta F/F', 'FontSize', 10);
        title('Correlation', 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box off;
        set(gca, 'FontSize', 9, 'LineWidth', 1.5);
        
    end
    
    % Add main title for the whole figure with proper position
    % Create a title at the top of the figure with some padding
    axes('Position', [0, 0.95, 1, 0.05], 'Visible', 'off');
    
    % Add memory strength if available
    if ~isempty(memory_strength)
        if strcmp(options.group, 'all')
            title_str = sprintf('Mouse %s: Analysis for %s (Group: %s, Food Arm: %s, %s)', ...
                mouse_name, arm_label, group, food_arm_location, memory_strength);
        else
            title_str = sprintf('Mouse %s: Analysis for %s (Group: %s, Food Arm: %s, %s)', ...
                mouse_name, arm_label, options.group, food_arm_location, memory_strength);
        end
    else
        if strcmp(options.group, 'all')
            title_str = sprintf('Mouse %s: Analysis for %s (Group: %s, Food Arm: %s)', ...
                mouse_name, arm_label, group, food_arm_location);
        else
            title_str = sprintf('Mouse %s: Analysis for %s (Group: %s, Food Arm: %s)', ...
                mouse_name, arm_label, options.group, food_arm_location);
        end
    end
    
    text(0.5, 0, title_str, 'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
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