function plotAveragedWaterObjectAnalysis(mice_all, mouse_prefixes, options)
    % Plot averaged heatmaps across all mice for water and object sessions
    % Creates two figures with averaged z-scored dF/F heatmaps:
    % 1. Water analysis: before, water, difference (water-before)
    % 2. Object analysis: before, object, difference (object-before)
    %
    % Parameters:
    %   mice_all: cell array with all mouse data
    %   mouse_prefixes: cell array of mouse ID prefixes (e.g., {'F29', 'M28', ...})
    %   options: struct with options for filtering and plotting
    
    % Define column indices
    COL_TIME = 1; COL_X = 2; COL_Y = 3; COL_DFF = 6;
    COL_DIST_WATER = 9; COL_DIST_OBJECT = 10; COL_GROOM = 13;
    
    % Set default options
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    if ~isfield(options, 'dist_limit'), options.dist_limit = 5; end
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'colormap'), options.colormap = cool; end
    if ~isfield(options, 'sigma'), options.sigma = 4; end
    if ~isfield(options, 'resolution'), options.resolution = 100; end
    
    % Initialize storage for individual mouse grids
    water_grids = struct();
    water_grids.before = {};
    water_grids.water = {};
    water_grids.condition = '';
    water_grids.stimulus = '';
    
    object_grids = struct();
    object_grids.before = {};
    object_grids.object = {};
    object_grids.condition = '';
    object_grids.stimulus = '';
    
    % Collect all coordinates to determine common bounds
    all_x_water = []; all_y_water = [];
    all_x_object = []; all_y_object = [];
    
    fprintf('Collecting data from %d mice...\n', length(mouse_prefixes));
    
    % First pass: collect coordinates for common bounds
    valid_mice_water = {};
    valid_mice_object = {};
    
    for i = 1:length(mouse_prefixes)
        mouse_prefix = mouse_prefixes{i};
        
        % Find sessions for this mouse
        mouse_indices = find(contains(mice_all(:,1), mouse_prefix));
        
        before_idx = []; water_idx = []; object_idx = [];
        
        for j = mouse_indices'
            session_name = mice_all{j,1};
            if contains(session_name, 'before')
                before_idx = j;
            elseif contains(session_name, 'water')
                water_idx = j;
            elseif contains(session_name, 'object')
                object_idx = j;
            end
        end
        
        % Check if we have all water sessions
        if ~isempty(before_idx) && ~isempty(water_idx)
            before_data = mice_all{before_idx, 4};
            water_data = mice_all{water_idx, 4};
            
            % Get condition and stimulus
            if isempty(water_grids.condition)
                water_grids.condition = mice_all{before_idx, 2};
                water_grids.stimulus = mice_all{before_idx, 3};
            end
            
            % Collect coordinates
            [x_before, y_before] = getFilteredCoordinates(before_data, COL_X, COL_Y, COL_DIST_WATER, COL_GROOM, options);
            [x_water, y_water] = getFilteredCoordinates(water_data, COL_X, COL_Y, COL_DIST_WATER, COL_GROOM, options);
            
            all_x_water = [all_x_water; x_before; x_water];
            all_y_water = [all_y_water; y_before; y_water];
            
            valid_mice_water{end+1} = mouse_prefix;
        end
        
        % Check if we have all object sessions
        if ~isempty(before_idx) && ~isempty(object_idx)
            before_data = mice_all{before_idx, 4};
            object_data = mice_all{object_idx, 4};
            
            % Get condition and stimulus
            if isempty(object_grids.condition)
                object_grids.condition = mice_all{before_idx, 2};
                object_grids.stimulus = mice_all{before_idx, 3};
            end
            
            % Collect coordinates
            [x_before, y_before] = getFilteredCoordinates(before_data, COL_X, COL_Y, COL_DIST_OBJECT, COL_GROOM, options);
            [x_object, y_object] = getFilteredCoordinates(object_data, COL_X, COL_Y, COL_DIST_OBJECT, COL_GROOM, options);
            
            all_x_object = [all_x_object; x_before; x_object];
            all_y_object = [all_y_object; y_before; y_object];
            
            valid_mice_object{end+1} = mouse_prefix;
        end
    end
    
    fprintf('Found %d mice with water sessions\n', length(valid_mice_water));
    fprintf('Found %d mice with object sessions\n', length(valid_mice_object));
    
    % Calculate common bounds for water analysis
    if ~isempty(all_x_water)
        water_bounds = struct();
        water_bounds.x_min = min(all_x_water);
        water_bounds.x_max = max(all_x_water);
        water_bounds.y_min = min(all_y_water);
        water_bounds.y_max = max(all_y_water);
        
        fprintf('Water bounds - X: [%.1f, %.1f], Y: [%.1f, %.1f]\n', ...
                water_bounds.x_min, water_bounds.x_max, ...
                water_bounds.y_min, water_bounds.y_max);
    else
        error('No valid water data found');
    end
    
    % Calculate common bounds for object analysis
    if ~isempty(all_x_object)
        object_bounds = struct();
        object_bounds.x_min = min(all_x_object);
        object_bounds.x_max = max(all_x_object);
        object_bounds.y_min = min(all_y_object);
        object_bounds.y_max = max(all_y_object);
        
        fprintf('Object bounds - X: [%.1f, %.1f], Y: [%.1f, %.1f]\n', ...
                object_bounds.x_min, object_bounds.x_max, ...
                object_bounds.y_min, object_bounds.y_max);
    else
        error('No valid object data found');
    end
    
    % Second pass: create grids for each mouse
    fprintf('\nCreating individual mouse grids...\n');
    
    for i = 1:length(mouse_prefixes)
        mouse_prefix = mouse_prefixes{i};
        
        % Find sessions for this mouse
        mouse_indices = find(contains(mice_all(:,1), mouse_prefix));
        
        before_idx = []; water_idx = []; object_idx = [];
        
        for j = mouse_indices'
            session_name = mice_all{j,1};
            if contains(session_name, 'before')
                before_idx = j;
            elseif contains(session_name, 'water')
                water_idx = j;
            elseif contains(session_name, 'object')
                object_idx = j;
            end
        end
        
        % Process water sessions
        if ~isempty(before_idx) && ~isempty(water_idx)
            fprintf('  Processing water sessions for %s...\n', mouse_prefix);
            
            before_data = mice_all{before_idx, 4};
            water_data = mice_all{water_idx, 4};
            
            % Create grids with common bounds
            before_grid = createSpatialGridWithBounds(before_data, COL_X, COL_Y, COL_DFF, ...
                                                     COL_DIST_WATER, COL_GROOM, options, water_bounds);
            water_grid = createSpatialGridWithBounds(water_data, COL_X, COL_Y, COL_DFF, ...
                                                    COL_DIST_WATER, COL_GROOM, options, water_bounds);
            
            % Store grids
            water_grids.before{end+1} = before_grid;
            water_grids.water{end+1} = water_grid;
        end
        
        % Process object sessions
        if ~isempty(before_idx) && ~isempty(object_idx)
            fprintf('  Processing object sessions for %s...\n', mouse_prefix);
            
            before_data = mice_all{before_idx, 4};
            object_data = mice_all{object_idx, 4};
            
            % Create grids with common bounds
            before_grid = createSpatialGridWithBounds(before_data, COL_X, COL_Y, COL_DFF, ...
                                                     COL_DIST_OBJECT, COL_GROOM, options, object_bounds);
            object_grid = createSpatialGridWithBounds(object_data, COL_X, COL_Y, COL_DFF, ...
                                                    COL_DIST_OBJECT, COL_GROOM, options, object_bounds);
            
            % Store grids
            object_grids.before{end+1} = before_grid;
            object_grids.object{end+1} = object_grid;
        end
    end
    
    % Average grids across mice for water analysis
    if ~isempty(water_grids.before)
        fprintf('\nAveraging water grids across %d mice...\n', length(water_grids.before));
        
        avg_before_water = averageGrids(water_grids.before);
        avg_water = averageGrids(water_grids.water);
        avg_diff_water = computeAveragedDifference(avg_water, avg_before_water);
        
        % Plot water analysis
        figure('Position', [100, 100, 1600, 600]);
        sgtitle(sprintf('Averaged Water Analysis (N=%d mice, %s / %s)', ...
                       length(water_grids.before), upper(water_grids.condition), water_grids.stimulus), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % Before water heatmap
        subplot(1, 3, 1);
        plotHeatmap(avg_before_water, water_bounds, options.colormap);
        title('Before Water Discovery', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Water session heatmap
        subplot(1, 3, 2);
        plotHeatmap(avg_water, water_bounds, options.colormap);
        title('Water Session', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Difference (Water - Before)
        subplot(1, 3, 3);
        plotDifferenceHeatmap(avg_diff_water, water_bounds, options.colormap);
        title('Difference (Water - Before)', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Sync color limits for before and water plots
        syncColorLimitsForSubplots([1, 2]);
    end
    
    % Average grids across mice for object analysis
    if ~isempty(object_grids.before)
        fprintf('\nAveraging object grids across %d mice...\n', length(object_grids.before));
        
        avg_before_object = averageGrids(object_grids.before);
        avg_object = averageGrids(object_grids.object);
        avg_diff_object = computeAveragedDifference(avg_object, avg_before_object);
        
        % Plot object analysis
        figure('Position', [200, 200, 1600, 600]);
        sgtitle(sprintf('Averaged Object Analysis (N=%d mice, %s / %s)', ...
                       length(object_grids.before), upper(object_grids.condition), object_grids.stimulus), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        % Before object heatmap
        subplot(1, 3, 1);
        plotHeatmap(avg_before_object, object_bounds, options.colormap);
        title('Before Object Discovery', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Object session heatmap
        subplot(1, 3, 2);
        plotHeatmap(avg_object, object_bounds, options.colormap);
        title('Object Session', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Difference (Object - Before)
        subplot(1, 3, 3);
        plotDifferenceHeatmap(avg_diff_object, object_bounds, options.colormap);
        title('Difference (Object - Before)', 'FontSize', 14, 'FontWeight', 'bold');
        
        % Sync color limits for before and object plots
        syncColorLimitsForSubplots([1, 2]);
    end
end

%% Helper Functions

function [x_coords, y_coords] = getFilteredCoordinates(data, col_x, col_y, col_dist, col_groom, options)
    % Extract and filter coordinates based on distance and grooming
    
    x_coords = data(:, col_x);
    y_coords = data(:, col_y);
    
    % Apply filtering
    if options.remove_grooming && size(data, 2) >= col_groom
        valid_mask = data(:, col_dist) > options.dist_limit & data(:, col_groom) == 0;
    else
        valid_mask = data(:, col_dist) > options.dist_limit;
    end
    
    x_coords = x_coords(valid_mask);
    y_coords = y_coords(valid_mask);
end

function spatial_grid = createSpatialGridWithBounds(data, col_x, col_y, col_dff, col_dist, col_groom, options, bounds)
    % Create spatial grid using specified bounds for coordinate alignment
    % Same as in plotMouseWaterObjectAnalysis
    
    % Extract coordinates and values
    x_coords = data(:, col_x);
    y_coords = data(:, col_y);
    dff_values = data(:, col_dff);
    
    % Apply filtering
    if options.remove_grooming && size(data, 2) >= col_groom
        valid_mask = data(:, col_dist) > options.dist_limit & data(:, col_groom) == 0;
    else
        valid_mask = data(:, col_dist) > options.dist_limit;
    end
    
    % Filter data
    x_filt = x_coords(valid_mask);
    y_filt = y_coords(valid_mask);
    dff_filt = dff_values(valid_mask);
    
    % Z-score the dF/F values
    if ~isempty(dff_filt)
        dff_zscore = (dff_filt - mean(dff_filt)) / std(dff_filt);
    else
        dff_zscore = dff_filt;
    end
    
    % Use provided bounds
    x_min = bounds.x_min; x_max = bounds.x_max;
    y_min = bounds.y_min; y_max = bounds.y_max;
    
    % Create spatial bounds
    if isempty(x_filt) || x_max <= x_min || y_max <= y_min
        spatial_grid = NaN(options.resolution, options.resolution);
        return;
    end
    
    % Initialize grids
    resolution = options.resolution;
    dff_grid = zeros(resolution, resolution);
    visit_count = zeros(resolution, resolution);
    
    % Map coordinates to grid indices using common bounds
    x_idx = round((x_filt - x_min) / (x_max - x_min) * (resolution-1)) + 1;
    y_idx = round((y_filt - y_min) / (y_max - y_min) * (resolution-1)) + 1;
    
    % Enforce valid indices
    x_idx = max(1, min(resolution, x_idx));
    y_idx = max(1, min(resolution, y_idx));
    
    % Accumulate values
    for i = 1:length(x_filt)
        xi = x_idx(i); yi = y_idx(i);
        dff_grid(yi, xi) = dff_grid(yi, xi) + dff_zscore(i);
        visit_count(yi, xi) = visit_count(yi, xi) + 1;
    end
    
    % Calculate means
    visited_mask = visit_count > 0;
    mean_dff_grid = zeros(size(dff_grid));
    mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
    mean_dff_grid(~visited_mask) = NaN;
    
    % Apply Gaussian smoothing
    spatial_grid = applySmoothing(mean_dff_grid, options.sigma);
end

function smoothed_grid = applySmoothing(grid, sigma)
    % Apply Gaussian smoothing while preserving NaN values
    
    kernel_size = 11; % Must be odd
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), ...
                      -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*sigma^2));
    kernel = kernel / sum(kernel(:));
    
    % Initialize smoothed grid
    smoothed_grid = NaN(size(grid));
    valid_mask = ~isnan(grid);
    
    % Pad arrays to handle edges
    padded_grid = padarray(grid, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);
    
    % Apply smoothing
    for i = 1+floor(kernel_size/2):size(padded_grid,1)-floor(kernel_size/2)
        for j = 1+floor(kernel_size/2):size(padded_grid,2)-floor(kernel_size/2)
            if padded_mask(i, j)
                window = padded_grid(i-floor(kernel_size/2):i+floor(kernel_size/2), ...
                                   j-floor(kernel_size/2):j+floor(kernel_size/2));
                window_mask = ~isnan(window);
                if sum(window_mask(:)) > 0
                    weighted_values = window .* kernel;
                    weighted_values(~window_mask) = 0;
                    weight_sum = sum(kernel(window_mask));
                    if weight_sum > 0
                        smoothed_grid(i-floor(kernel_size/2), j-floor(kernel_size/2)) = ...
                            sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end
end

function avg_grid = averageGrids(grids)
    % Average multiple grids across mice, handling NaN values properly
    
    if isempty(grids)
        avg_grid = [];
        return;
    end
    
    % Stack all grids
    resolution = size(grids{1}, 1);
    grid_stack = NaN(resolution, resolution, length(grids));
    
    for i = 1:length(grids)
        grid_stack(:,:,i) = grids{i};
    end
    
    % Average across mice, ignoring NaN values
    avg_grid = nanmean(grid_stack, 3);
end

function diff_grid = computeAveragedDifference(grid1, grid2)
    % Compute difference between two averaged grids
    
    valid_both = ~isnan(grid1) & ~isnan(grid2);
    diff_grid = NaN(size(grid1));
    diff_grid(valid_both) = grid1(valid_both) - grid2(valid_both);
end

function plotHeatmap(grid, bounds, cmap)
    % Plot heatmap with proper axis scaling and NaN handling
    
    h = imagesc('XData', [bounds.x_min, bounds.x_max], ...
                'YData', [bounds.y_min, bounds.y_max], ...
                'CData', grid);
    
    % Make NaN values appear white (transparent)
    set(h, 'AlphaData', ~isnan(grid));
    
    axis xy; axis image;
    colormap(gca, cmap);
    colorbar;
    
    % Set white background for NaN areas
    set(gca, 'Color', 'white');
    
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
end

function plotDifferenceHeatmap(diff_grid, bounds, cmap)
    % Plot difference heatmap with diverging colormap and NaN handling
    
    h = imagesc('XData', [bounds.x_min, bounds.x_max], ...
                'YData', [bounds.y_min, bounds.y_max], ...
                'CData', diff_grid);
    
    % Make NaN values appear white (transparent)
    set(h, 'AlphaData', ~isnan(diff_grid));
    
    axis xy; axis image;
    
    % Use diverging colormap for differences
    colormap(gca, redblue(256));
    colorbar;
    
    % Set white background for NaN areas
    set(gca, 'Color', 'white');
    
    % Center colorbar at zero
    valid_vals = diff_grid(~isnan(diff_grid));
    if ~isempty(valid_vals)
        max_abs = max(abs(valid_vals));
        caxis([-max_abs, max_abs]);
    end
    
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
end

function syncColorLimitsForSubplots(subplot_indices)
    % Synchronize color limits across specified subplots
    
    all_values = [];
    
    for idx = subplot_indices
        subplot(1, 3, idx);
        h = findobj(gca, 'Type', 'Image');
        if ~isempty(h)
            data = get(h, 'CData');
            valid_data = data(~isnan(data));
            all_values = [all_values; valid_data(:)];
        end
    end
    
    if ~isempty(all_values)
        clim_range = [min(all_values), max(all_values)];
        for idx = subplot_indices
            subplot(1, 3, idx);
            caxis(clim_range);
        end
    end
end

function cmap = redblue(n)
    % Create red-blue diverging colormap
    
    if nargin < 1
        n = 256;
    end
    
    % Create diverging colormap: blue -> white -> red
    half = ceil(n/2);
    
    % Blue to white
    r1 = linspace(0, 1, half)';
    g1 = linspace(0, 1, half)';
    b1 = ones(half, 1);
    
    % White to red
    r2 = ones(n-half, 1);
    g2 = linspace(1, 0, n-half)';
    b2 = linspace(1, 0, n-half)';
    
    cmap = [r1, g1, b1; r2, g2, b2];
end