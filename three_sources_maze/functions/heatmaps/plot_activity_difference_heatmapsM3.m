function plot_activity_difference_heatmapsM3(mice_all, condition, colormap_name, options)
    % PLOT_ACTIVITY_DIFFERENCE_HEATMAPS Creates z-scored dF/F activity difference heatmaps (Learning-Before, Test-Before)
    % Similar to plot_spatial_difference_heatmaps but uses z-scored dF/F activity instead of occupancy
    % Uses the subtraction approach from plotMouseDffDifferenceHeatmap (z-score per period then subtract)
    %
    % Inputs:
    %   mice_all - Cell array containing mouse data with structure:
    %              Column 1: time, Column 2: x, Column 3: y, Column 6: dF/F, 
    %              Column 7: speed, Column 9: distance to food, Column 12: door, Column 13: grooming
    %   condition - String specifying the condition to analyze (e.g., 'control', 'CNO', 'chrimson')
    %   colormap_name - String specifying the colormap (e.g., 'RdBu', 'coolwarm', 'seismic')
    %   options - Structure with optional parameters
    %
    % Options:
    %   .bin_size - Spatial bin size for heatmap (default: 4)
    %   .x_range - X coordinate range [min max] (default: auto-detected)
    %   .y_range - Y coordinate range [min max] (default: auto-detected)
    %   .min_time_threshold - Minimum time spent to include data (default: 1 second)
    %   .diff_clim - Custom difference color limits [min max] (default: auto-calculated)
    %   .exclude_grooming - Exclude data points during grooming behavior (default: false)
    %   .speed_threshold - Minimum speed threshold to include data points (default: 0, no threshold)
    %   .resolution - Grid resolution for spatial maps (default: 100)
    %   .sigma - Gaussian smoothing sigma (default: 4 - increased for smoother maps)
    
    % Column definitions for data structure
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    COL_DFF = 6;
    COL_SPEED = 7;
    COL_DIST_FOOD = 9;
    COL_DOOR = 12;
    COL_GROOM = 13;
    
    % Set default options
    if nargin < 4
        options = struct();
    end

    if ~isfield(options, 'min_time_threshold'), options.min_time_threshold = 1; end
    if ~isfield(options, 'exclude_grooming'), options.exclude_grooming = false; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 0; end  % 0 means no threshold
    if ~isfield(options, 'resolution'), options.resolution = 100; end
    if ~isfield(options, 'sigma'), options.sigma = 4; end  % Increased to fill gaps with sparse data
    
    % Set default colormap for difference plots
    if nargin < 3 || isempty(colormap_name)
        colormap_name = 'RdBu'; % Diverging colormap good for differences
    end
    
    % Session information
    sessions = {'sess0', 'sess1', 'sess2'};
    session_names = {'Before', 'Learning', 'Test'};
    
    % Get unique food arm locations for this condition
    condition_indices = strcmp(mice_all(:,2), condition) & strcmp(mice_all(:,5), 'ok_signal');
    if sum(condition_indices) == 0
        error('No data found for condition "%s" with good signal quality', condition);
    end
    
    food_arms = unique(mice_all(condition_indices, 3));
    fprintf('Found food arm locations for condition "%s": %s\n', condition, strjoin(food_arms, ', '));
    
    % Initialize data collection structures for each food arm location
    arm_data = struct();
    all_coordinates = [];
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        arm_data.(food_arm_loc) = struct();
        arm_data.(food_arm_loc).mouse_data = {}; % Store individual mouse difference maps
        arm_data.(food_arm_loc).mouse_ids = {};
    end
    
    % Track filtering statistics
    total_points_before_filtering = 0;
    total_points_after_filtering = 0;
    speed_filtered_points = 0;
    grooming_filtered_points = 0;
    
    % First pass: collect all coordinates to determine common ranges
    mice_count_by_food_arm = containers.Map();
    for fa = 1:length(food_arms)
        mice_count_by_food_arm(food_arms{fa}) = 0;
    end
    
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        mouse_condition = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        signal_quality = mice_all{i, 5};
        
        % Skip if not the desired condition or bad signal
        if ~strcmp(mouse_condition, condition) || strcmp(signal_quality, 'bad_signal')
            continue;
        end
        
        % Track original data size
        original_size = size(data, 1);
        total_points_before_filtering = total_points_before_filtering + original_size;
        
        % Apply filtering and collect coordinates
        filtered_data = data;
        
        % Apply grooming filter if requested
        if options.exclude_grooming && size(filtered_data, 2) >= COL_GROOM
            grooming_indices = filtered_data(:, COL_GROOM) == 0;
            grooming_filtered_points = grooming_filtered_points + sum(~grooming_indices);
            filtered_data = filtered_data(grooming_indices, :);
        end
        
        % Apply speed threshold if requested
        if options.speed_threshold > 0 && size(filtered_data, 2) >= COL_SPEED
            speed_indices = filtered_data(:, COL_SPEED) >= options.speed_threshold;
            speed_filtered_points = speed_filtered_points + sum(~speed_indices);
            filtered_data = filtered_data(speed_indices, :);
        end
        
        % Track final data size
        total_points_after_filtering = total_points_after_filtering + size(filtered_data, 1);
        
        % Collect all x,y coordinates for range determination
        if ~isempty(filtered_data)
            all_coordinates = [all_coordinates; filtered_data(:, [COL_X, COL_Y])];
        end
    end
    
    % Set coordinate ranges (ensure consistency across all plots)
    if ~isfield(options, 'x_range')
        options.x_range = [min(all_coordinates(:,1)) - 10, max(all_coordinates(:,1)) + 10];
    end
    if ~isfield(options, 'y_range')
        options.y_range = [min(all_coordinates(:,2)) - 10, max(all_coordinates(:,2)) + 10];
    end
    
    fprintf('Activity grid: %d x %d resolution\n', options.resolution, options.resolution);
    fprintf('Coordinate ranges: X [%.1f, %.1f], Y [%.1f, %.1f]\n', ...
            options.x_range(1), options.x_range(2), options.y_range(1), options.y_range(2));
    fprintf('Smoothing parameter (sigma): %.1f\n', options.sigma);
    
    % Print filtering statistics
    fprintf('\n=== DATA FILTERING STATISTICS ===\n');
    fprintf('Total data points before filtering: %d\n', total_points_before_filtering);
    if options.exclude_grooming
        fprintf('Points removed by grooming filter: %d (%.1f%%)\n', ...
                grooming_filtered_points, grooming_filtered_points/total_points_before_filtering*100);
    end
    if options.speed_threshold > 0
        fprintf('Points removed by speed threshold (< %.2f): %d (%.1f%%)\n', ...
                options.speed_threshold, speed_filtered_points, speed_filtered_points/total_points_before_filtering*100);
    end
    fprintf('Total data points after filtering: %d (%.1f%% retained)\n', ...
            total_points_after_filtering, total_points_after_filtering/total_points_before_filtering*100);
    
    % Create common bounds structure
    bounds = struct('x_min', options.x_range(1), 'x_max', options.x_range(2), ...
                   'y_min', options.y_range(1), 'y_max', options.y_range(2));
    
    % Second pass: process each mouse and create individual difference maps
    unique_mice_by_food_arm = containers.Map();
    for fa = 1:length(food_arms)
        unique_mice_by_food_arm(food_arms{fa}) = {};
    end
    
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        mouse_condition = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        signal_quality = mice_all{i, 5};
        
        % Skip if not the desired condition or bad signal
        if ~strcmp(mouse_condition, condition) || strcmp(signal_quality, 'bad_signal')
            continue;
        end
        
        % Determine session number
        sess_idx = 0;
        for s = 1:length(sessions)
            if contains(session_info, sessions{s})
                sess_idx = s;
                break;
            end
        end
        
        if sess_idx == 0
            continue;
        end
        
        % Extract mouse ID for tracking
        mouse_id = '';
        if contains(session_info, '_')
            parts = strsplit(session_info, '_');
            mouse_id = parts{1};
        else
            mouse_id = session_info;
        end
        
        % Track unique mice
        current_mice = unique_mice_by_food_arm(food_arm);
        if ~any(strcmp(current_mice, mouse_id))
            current_mice{end+1} = mouse_id;
            unique_mice_by_food_arm(food_arm) = current_mice;
        end
        
        % Process individual mouse data and create difference maps
        if sess_idx == 1 % Before session - we'll need this for difference calculation
            % Store before session data for this mouse
            before_key = sprintf('%s_before', mouse_id);
            arm_data.(food_arm).(before_key) = data;
        elseif sess_idx > 1 % Learning or Test session
            % Look for corresponding before session data
            before_key = sprintf('%s_before', mouse_id);
            
            if isfield(arm_data.(food_arm), before_key)
                before_data = arm_data.(food_arm).(before_key);
                current_data = data;
                
                % Create difference map for this mouse (current session - before session)
                mouse_diff_map = createMouseActivityDifferenceMap(before_data, current_data, ...
                                                                bounds, options, sess_idx);
                
                if ~isempty(mouse_diff_map) && ~all(isnan(mouse_diff_map(:)))
                    % Store the difference map
                    session_name = session_names{sess_idx};
                    if ~isfield(arm_data.(food_arm), session_name)
                        arm_data.(food_arm).(session_name) = {};
                    end
                    arm_data.(food_arm).(session_name){end+1} = mouse_diff_map;
                    
                    % Store mouse ID for this difference map
                    mouse_key = sprintf('%s_mouse_ids', session_name);
                    if ~isfield(arm_data.(food_arm), mouse_key)
                        arm_data.(food_arm).(mouse_key) = {};
                    end
                    arm_data.(food_arm).(mouse_key){end+1} = mouse_id;
                end
            end
        end
    end
    
    % Update mice counts
    for fa = 1:length(food_arms)
        mice_count_by_food_arm(food_arms{fa}) = length(unique_mice_by_food_arm(food_arms{fa}));
    end
    
    % Average difference maps across mice for each food arm and session
    averaged_difference_data = struct();
    all_differences = [];
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        averaged_difference_data.(food_arm_loc) = struct();
        
        % Process Learning - Before and Test - Before
        for sess_idx = 2:3 % Learning and Test sessions
            session_name = session_names{sess_idx};
            
            if isfield(arm_data.(food_arm_loc), session_name) && ...
               ~isempty(arm_data.(food_arm_loc).(session_name))
                
                % Stack all difference maps for this condition
                diff_maps = arm_data.(food_arm_loc).(session_name);
                
                if ~isempty(diff_maps)
                    % Convert cell array to 3D array
                    map_stack = cat(3, diff_maps{:});
                    
                    % Average across mice (3rd dimension)
                    averaged_map = nanmean(map_stack, 3);
                    
                    averaged_difference_data.(food_arm_loc).(session_name) = averaged_map;
                    all_differences = [all_differences; averaged_map(:)];
                    
                    fprintf('Averaged %d difference maps for %s - %s\n', ...
                            length(diff_maps), food_arm_loc, session_name);
                end
            end
        end
    end
    
    % Calculate global difference color limits (symmetric around 0)
    if isfield(options, 'diff_clim') && ~isempty(options.diff_clim)
        diff_clim_max = max(abs(options.diff_clim));
        diff_clim_min = -diff_clim_max;
        fprintf('Using custom difference color limits: [%.2f, %.2f] z-scored dF/F change\n', diff_clim_min, diff_clim_max);
    else
        if ~isempty(all_differences)
            all_differences_clean = all_differences(~isnan(all_differences));
            if ~isempty(all_differences_clean)
                diff_clim_max = max(abs(prctile(all_differences_clean, [5, 95]))); % Use 5th-95th percentile range
                diff_clim_min = -diff_clim_max;
            else
                diff_clim_max = 1;
                diff_clim_min = -1;
            end
        else
            diff_clim_max = 1;
            diff_clim_min = -1;
        end
        fprintf('Auto-calculated difference color limits: [%.2f, %.2f] z-scored dF/F change\n', diff_clim_min, diff_clim_max);
    end
    
    % Create coordinate arrays for plotting
    x_coords = linspace(bounds.x_min, bounds.x_max, options.resolution);
    y_coords = linspace(bounds.y_min, bounds.y_max, options.resolution);
    
    % Create plots for each food arm location
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        
        % Add filtering status to figure name
        filtering_status = '';
        if options.exclude_grooming
            filtering_status = [filtering_status, ' - No Grooming'];
        end
        if options.speed_threshold > 0
            filtering_status = [filtering_status, sprintf(' - Speed ≥ %.1f', options.speed_threshold)];
        end
        
        figure('Name', sprintf('Activity Difference Heatmaps - %s - Food Arm: %s - %s%s', ...
                              condition, food_arm_loc, colormap_name, filtering_status), ...
               'Position', [100 + (fa-1)*50, 100 + (fa-1)*50, 800, 400]);
        
        % Plot difference heatmaps
        diff_titles = {'Learning - Before', 'Test - Before'};
        diff_sessions = {'Learning', 'Test'};
        
        for diff_idx = 1:2
            subplot(1, 2, diff_idx);
            session_name = diff_sessions{diff_idx};
            
            if isfield(averaged_difference_data.(food_arm_loc), session_name)
                diff_map = averaged_difference_data.(food_arm_loc).(session_name);
                
                if ~isempty(diff_map) && ~all(isnan(diff_map(:)))
                    % Display the difference map
                    imagesc(x_coords, y_coords, diff_map, 'AlphaData', ~isnan(diff_map));
                    axis xy equal tight;
                    
                    % Apply colormap and set limits
                    colormap(gca, colormap_name);
                    set(gca, 'CLim', [diff_clim_min, diff_clim_max]);
                    
                    % Add colorbar
                    cb = colorbar;
                    cb.Label.String = '\Delta(Z-scored \Delta F/F)';
                    
                    title(sprintf('%s\n(Δ Z-scored ΔF/F)', diff_titles{diff_idx}), 'FontSize', 14, 'FontWeight', 'bold');
                    xlabel('X Position', 'FontSize', 12);
                    ylabel('Y Position', 'FontSize', 12);
                    
                else
                    title(sprintf('%s\n(No Data)', diff_titles{diff_idx}), 'FontSize', 14);
                    set(gca, 'XLim', options.x_range, 'YLim', options.y_range);
                end
            else
                title(sprintf('%s\n(No Data)', diff_titles{diff_idx}), 'FontSize', 14);
                set(gca, 'XLim', options.x_range, 'YLim', options.y_range);
            end
            
            axis off;
            box off;
        end
        
        % Add overall title
        main_title = sprintf('Activity Difference Analysis - %s - Food Arm: %s (n=%d mice)', ...
                           condition, food_arm_loc, mice_count_by_food_arm(food_arm_loc));
        if options.exclude_grooming
            main_title = [main_title, ' - Grooming Excluded'];
        end
        if options.speed_threshold > 0
            main_title = [main_title, sprintf(' - Speed ≥ %.1f', options.speed_threshold)];
        end
        
        sgtitle(main_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    % Print summary statistics
    fprintf('\n=== ACTIVITY DIFFERENCE ANALYSIS SUMMARY ===\n');
    fprintf('Condition: %s\n', condition);
    fprintf('Analysis: Z-scored dF/F difference heatmaps (Learning-Before, Test-Before)\n');
    
    % Print filtering settings
    fprintf('\nFILTERING SETTINGS:\n');
    if options.exclude_grooming
        fprintf('Grooming exclusion: ENABLED\n');
    else
        fprintf('Grooming exclusion: DISABLED\n');
    end
    if options.speed_threshold > 0
        fprintf('Speed threshold: %.2f (excluding low-speed data points)\n', options.speed_threshold);
    else
        fprintf('Speed threshold: DISABLED (all speeds included)\n');
    end
    
    fprintf('\nANALYSIS PARAMETERS:\n');
    fprintf('Difference color limits: [%.2f, %.2f] z-scored dF/F change\n', diff_clim_min, diff_clim_max);
    fprintf('Resolution: %d x %d grid\n', options.resolution, options.resolution);
    fprintf('Smoothing sigma: %.1f\n', options.sigma);
    fprintf('Coordinate ranges: X [%.1f, %.1f], Y [%.1f, %.1f]\n', ...
            options.x_range(1), options.x_range(2), options.y_range(1), options.y_range(2));
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        fprintf('\nFood Arm Location: %s (%d mice)\n', food_arm_loc, mice_count_by_food_arm(food_arm_loc));
        
        % Calculate and report change statistics for each session comparison
        diff_sessions = {'Learning', 'Test'};
        diff_titles = {'Learning - Before', 'Test - Before'};
        
        for sess_idx = 1:2
            session_name = diff_sessions{sess_idx};
            diff_title = diff_titles{sess_idx};
            
            if isfield(averaged_difference_data.(food_arm_loc), session_name)
                diff_map = averaged_difference_data.(food_arm_loc).(session_name);
                
                if ~isempty(diff_map) && ~all(isnan(diff_map(:)))
                    valid_pixels = diff_map(~isnan(diff_map));
                    
                    pos_change = sum(valid_pixels > 0.1); % Areas with >0.1 z-score increase
                    neg_change = sum(valid_pixels < -0.1); % Areas with >0.1 z-score decrease
                    total_pixels = length(valid_pixels);
                    
                    fprintf('  %s: %d pixels increased (%.1f%%), %d pixels decreased (%.1f%%)\n', ...
                            diff_title, pos_change, pos_change/total_pixels*100, ...
                            neg_change, neg_change/total_pixels*100);
                    
                    % Report extreme changes
                    max_increase = max(valid_pixels);
                    max_decrease = min(valid_pixels);
                    mean_change = mean(valid_pixels);
                    fprintf('    Max increase: %.2f z-score, Max decrease: %.2f z-score, Mean: %.2f z-score\n', ...
                            max_increase, max_decrease, mean_change);
                else
                    fprintf('  %s: No valid data available\n', diff_title);
                end
            else
                fprintf('  %s: No data available\n', diff_title);
            end
        end
    end
end

function difference_map = createMouseActivityDifferenceMap(before_data, current_data, bounds, options, sess_idx)
    % Create difference map for a single mouse (current session - before session)
    % Uses z-scoring approach from plotMouseDffDifferenceHeatmap
    
    % Column definitions
    COL_X = 2;
    COL_Y = 3;
    COL_DFF = 6;
    COL_SPEED = 7;
    COL_DIST_FOOD = 9;
    COL_GROOM = 13;
    
    % Filter before period data
    valid_before = true(size(before_data, 1), 1);
    
    % Apply grooming filter
    if options.exclude_grooming && size(before_data, 2) >= COL_GROOM
        valid_before = valid_before & (before_data(:, COL_GROOM) == 0);
    end
    
    % Apply speed threshold
    if options.speed_threshold > 0 && size(before_data, 2) >= COL_SPEED
        valid_before = valid_before & (before_data(:, COL_SPEED) >= options.speed_threshold);
    end
    
    % Filter current period data
    valid_current = true(size(current_data, 1), 1);
    
    % Apply grooming filter
    if options.exclude_grooming && size(current_data, 2) >= COL_GROOM
        valid_current = valid_current & (current_data(:, COL_GROOM) == 0);
    end
    
    % Apply speed threshold
    if options.speed_threshold > 0 && size(current_data, 2) >= COL_SPEED
        valid_current = valid_current & (current_data(:, COL_SPEED) >= options.speed_threshold);
    end
    
    % Extract filtered data
    x_before = before_data(valid_before, COL_X);
    y_before = before_data(valid_before, COL_Y);
    dff_before = before_data(valid_before, COL_DFF);
    
    x_current = current_data(valid_current, COL_X);
    y_current = current_data(valid_current, COL_Y);
    dff_current = current_data(valid_current, COL_DFF);
    
    % Check if we have enough data
    sampling_rate = 10; % Assuming 10 Hz
    min_points = options.min_time_threshold * sampling_rate;
    
    if length(x_before) < min_points || length(x_current) < min_points
        difference_map = [];
        return;
    end
    
    % Z-score dF/F values for each period separately
    if ~isempty(dff_before) && std(dff_before) > 0
        dff_before_zscore = (dff_before - mean(dff_before)) / std(dff_before);
    else
        difference_map = [];
        return;
    end
    
    if ~isempty(dff_current) && std(dff_current) > 0
        dff_current_zscore = (dff_current - mean(dff_current)) / std(dff_current);
    else
        difference_map = [];
        return;
    end
    
    % Create spatial grids for both periods using common bounds
    before_grid = createSpatialGridWithCommonBounds(x_before, y_before, dff_before_zscore, ...
                                                   bounds, options.resolution, options.sigma);
    current_grid = createSpatialGridWithCommonBounds(x_current, y_current, dff_current_zscore, ...
                                                    bounds, options.resolution, options.sigma);
    
    % Compute difference (current - before)
    valid_both = ~isnan(before_grid) & ~isnan(current_grid);
    difference_map = NaN(size(before_grid));
    difference_map(valid_both) = current_grid(valid_both) - before_grid(valid_both);
end

function spatial_grid = createSpatialGridWithCommonBounds(x_coords, y_coords, dff_values, bounds, resolution, sigma)
    % Create a smoothed spatial grid using specified bounds for coordinate alignment
    % Uses dynamic kernel size that scales with sigma (same as plotMouseDffDifferenceHeatmap)
    
    % Initialize grid for values and visit count
    dff_grid = zeros(resolution, resolution);
    visit_count = zeros(resolution, resolution);
    
    % Use provided bounds for coordinate mapping
    x_min = bounds.x_min;
    x_max = bounds.x_max;
    y_min = bounds.y_min;
    y_max = bounds.y_max;
    
    % Check for valid bounds
    if isempty(x_coords) || x_max <= x_min || y_max <= y_min
        spatial_grid = NaN(resolution, resolution);
        return;
    end
    
    % Map coordinates to grid indices using common bounds
    x_idx = round((x_coords - x_min) / (x_max - x_min) * (resolution-1)) + 1;
    y_idx = round((y_coords - y_min) / (y_max - y_min) * (resolution-1)) + 1;
    
    % Enforce valid indices
    x_idx = max(1, min(resolution, x_idx));
    y_idx = max(1, min(resolution, y_idx));
    
    % Accumulate dF/F values and counts for each grid cell
    for i = 1:length(x_coords)
        xi = x_idx(i);
        yi = y_idx(i);
        dff_grid(yi, xi) = dff_grid(yi, xi) + dff_values(i);
        visit_count(yi, xi) = visit_count(yi, xi) + 1;
    end
    
    % Calculate mean dF/F for each cell
    visited_mask = visit_count > 0;
    mean_dff_grid = zeros(size(dff_grid));
    mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
    
    % Create masked array with NaNs for unvisited locations
    masked_dff = mean_dff_grid;
    masked_dff(~visited_mask) = NaN;
    
    % Apply Gaussian smoothing while preserving NaN values
    % Dynamic kernel size that scales with sigma (same as plotMouseDffDifferenceHeatmap)
    kernel_size = 2*ceil(3*sigma) + 1; % Scale with sigma
    kernel_size = min(kernel_size, 21); % Cap at reasonable size
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), ...
                      -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*sigma^2));
    kernel = kernel / sum(kernel(:));
    
    % Initialize smoothed grid
    smoothed_dff = NaN(size(masked_dff));
    valid_mask = ~isnan(masked_dff);
    
    % Pad arrays to handle edges
    padded_dff = padarray(masked_dff, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);
    
    % Apply smoothing
    for i = 1+floor(kernel_size/2):size(padded_dff,1)-floor(kernel_size/2)
        for j = 1+floor(kernel_size/2):size(padded_dff,2)-floor(kernel_size/2)
            if padded_mask(i, j)
                window = padded_dff(i-floor(kernel_size/2):i+floor(kernel_size/2), ...
                                   j-floor(kernel_size/2):j+floor(kernel_size/2));
                window_mask = ~isnan(window);
                if sum(window_mask(:)) > 0
                    weighted_values = window .* kernel;
                    weighted_values(~window_mask) = 0;
                    weight_sum = sum(kernel(window_mask));
                    if weight_sum > 0
                        smoothed_dff(i-floor(kernel_size/2), j-floor(kernel_size/2)) = ...
                            sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end
    
    spatial_grid = smoothed_dff;
end