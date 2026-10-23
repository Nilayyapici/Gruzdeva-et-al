function plotMouseDffDifferenceHeatmap(mice, mouse_index, options)
    % Plot difference heatmap (after - before) for a single mouse
    % Shows how spatial activity patterns changed after food discovery
    % Uses z-scored dF/F values for each period with proper coordinate alignment
    %
    % Parameters:
    %   mice: cell array with mouse data
    %   mouse_index: index of the mouse to plot
    %   options: struct with options for filtering
    %     - dist_limit: minimum distance threshold (default: 5)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - colormap: colormap for the difference plot (default: redblue)
    %     - sigma: Gaussian smoothing sigma (default: 4 - higher to fill gaps with sparse data)
    %     - resolution: grid resolution (default: 100)
    %     - difference_limits: [min max] color limits for difference heatmap (optional)
    %     - zscore_difference: boolean, whether to z-score the difference map (default: false)
    
    % Define constants
    COL_DIST = 5;     % Distance to food
    COL_DOOR = 7;     % Door status
    COL_GROOM = 10;   % Grooming
    COL_EATING = 9;   % Eating
    COL_DFF = 11;     % DFF data
    COL_TIME = 1;     % Time
    COL_X = 2;        % X position
    COL_Y = 3;        % Y position
    
    % Set default options if not provided
    if ~isfield(options, 'dist_limit')
        options.dist_limit = 5;
    end
    
    if ~isfield(options, 'remove_grooming')
        options.remove_grooming = true;
    end
    
    if ~isfield(options, 'colormap')
        options.colormap = cool;
    end

    if ~isfield(options, 'sigma')
        options.sigma = 4;  % Increased to fill gaps with sparse data
    end

    if ~isfield(options, 'resolution')
        options.resolution = 100;
    end
    
    if ~isfield(options, 'zscore_difference')
        options.zscore_difference = false;
    end

    % Get mouse data
    mouse_id = mice{mouse_index, 1};
    condition = mice{mouse_index, 2};
    stimulus = mice{mouse_index, 3};
    data = mice{mouse_index, 4};
    discovery = mice{mouse_index, 6};
    
    % Find end frame (second closed door or end of data)
    if length(data) > 11000
        closed_indices = find(data(discovery:end, COL_DOOR) < 1);
        if ~isempty(closed_indices)
            end_frame = discovery + closed_indices(1) - 1;
        else
            end_frame = length(data);
        end
    else
        end_frame = length(data);
    end
    
    %% Calculate common bounds from filtered data only
    bounds = calculateCommonBoundsFromFilteredData(data, discovery, end_frame, ...
                                                  COL_X, COL_Y, COL_DIST, COL_GROOM, COL_EATING, options);
    
    % Create common grid edges for both periods
    resolution = options.resolution;
    x_edges = linspace(bounds.x_min, bounds.x_max, resolution);
    y_edges = linspace(bounds.y_min, bounds.y_max, resolution);
    
    %% BEFORE PERIOD - Create spatial activity map
    x_before = data(1:discovery, COL_X);
    y_before = data(1:discovery, COL_Y);
    dff_before_all = data(1:discovery, COL_DFF);
    
    % Filter data before discovery if specified
    if options.remove_grooming
        valid_before = data(1:discovery, COL_DIST) > options.dist_limit & ...
                      data(1:discovery, COL_GROOM) == 0 & ...
                      data(1:discovery, COL_EATING) == 0;
    else
        valid_before = data(1:discovery, COL_DIST) > options.dist_limit & ...
                      data(1:discovery, COL_EATING) == 0;
    end
    
    % Apply filtering
    x_before_filt = x_before(valid_before);
    y_before_filt = y_before(valid_before);
    dff_before_filt = dff_before_all(valid_before);
    
    % Z-score the dF/F values for the before period
    if ~isempty(dff_before_filt)
        dff_before_zscore = (dff_before_filt - mean(dff_before_filt)) / std(dff_before_filt);
    else
        dff_before_zscore = dff_before_filt;
    end
    
    % Create spatial grid for before period using common bounds
    before_grid = createSpatialGridWithCommonBounds(x_before_filt, y_before_filt, dff_before_zscore, ...
                                                   bounds, options.resolution, options.sigma);
    
    %% AFTER PERIOD - Create spatial activity map
    x_after = data(discovery:end_frame, COL_X);
    y_after = data(discovery:end_frame, COL_Y);
    dff_after_all = data(discovery:end_frame, COL_DFF);
    
    % Filter data after discovery if specified
    if options.remove_grooming
        valid_after = data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                     data(discovery:end_frame, COL_GROOM) == 0 & ...
                     data(discovery:end_frame, COL_EATING) == 0;
    else
        valid_after = data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                     data(discovery:end_frame, COL_EATING) == 0;
    end
    
    % Apply filtering
    x_after_filt = x_after(valid_after);
    y_after_filt = y_after(valid_after);
    dff_after_filt = dff_after_all(valid_after);
    
    % Z-score the dF/F values for the after period
    if ~isempty(dff_after_filt)
        dff_after_zscore = (dff_after_filt - mean(dff_after_filt)) / std(dff_after_filt);
    else
        dff_after_zscore = dff_after_filt;
    end
    
    % Create spatial grid for after period using common bounds
    after_grid = createSpatialGridWithCommonBounds(x_after_filt, y_after_filt, dff_after_zscore, ...
                                                  bounds, options.resolution, options.sigma);
    
    %% COMPUTE DIFFERENCE (After - Before)
    % Only compute difference where both periods have data
    valid_both = ~isnan(before_grid) & ~isnan(after_grid);
    difference_grid = NaN(size(before_grid));
    difference_grid(valid_both) = after_grid(valid_both) - before_grid(valid_both);
    
    % Optionally z-score the difference
    if options.zscore_difference
        valid_diff = difference_grid(~isnan(difference_grid));
        if ~isempty(valid_diff)
            mean_diff = mean(valid_diff);
            std_diff = std(valid_diff);
            difference_grid(~isnan(difference_grid)) = (difference_grid(~isnan(difference_grid)) - mean_diff) / std_diff;
        end
    end
    
    %% CREATE FIGURE WITH THREE SUBPLOTS
    figure('Position', [100, 100, 1600, 600]);
    
    % Subplot 1: Before period
    subplot(1, 3, 1);
    h = imagesc(x_edges, y_edges, before_grid, 'AlphaData', ~isnan(before_grid));
    h.Interpolation = 'nearest';
    axis equal tight
    colormap(gca, options.colormap);
    c1 = colorbar;
    c1.Label.String = 'Z-scored \Delta F/F';
    c1.Label.FontSize = 12;
    title('Before Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11);
    box off
    axis off
    
    % Subplot 2: After period
    subplot(1, 3, 2);
    h = imagesc(x_edges, y_edges, after_grid, 'AlphaData', ~isnan(after_grid));
    h.Interpolation = 'nearest';
    axis equal tight
    colormap(gca, options.colormap);
    c2 = colorbar;
    c2.Label.String = 'Z-scored \Delta F/F';
    c2.Label.FontSize = 12;
    title('After Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11);
    box off
    axis off
    
    % Set same color limits for before and after plots
    all_values = [before_grid(:); after_grid(:)];
    all_values = all_values(~isnan(all_values));
    if ~isempty(all_values)
        clim_range = [min(all_values), max(all_values)];
        subplot(1, 3, 1); caxis(clim_range);
        subplot(1, 3, 2); caxis(clim_range);
    end
    
    % Subplot 3: Difference (After - Before)
    subplot(1, 3, 3);
    h = imagesc(x_edges, y_edges, difference_grid, 'AlphaData', ~isnan(difference_grid));
    h.Interpolation = 'nearest';
    axis equal tight
    colormap(gca, options.colormap);
    c = colorbar;
    
    if options.zscore_difference
        c.Label.String = 'Z-scored \Delta(Z-scored \Delta F/F)';
        title('Difference (After - Before): Z-scored', 'FontSize', 14, 'FontWeight', 'bold');
    else
        c.Label.String = '\Delta(Z-scored \Delta F/F) (After - Before)';
        title('Difference (After - Before)', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    c.Label.FontSize = 12;
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11);
    box off
    axis off
    
    % Set color limits for difference plot
    if isfield(options, 'difference_limits') && length(options.difference_limits) == 2
        caxis(options.difference_limits);
    else
        % Use symmetric limits around zero
        diff_values = difference_grid(~isnan(difference_grid));
        if ~isempty(diff_values)
            max_abs = max(abs(diff_values));
            caxis([-max_abs, max_abs]);
        end
    end
    
    % Add overall title
    sgtitle(sprintf('Z-scored Spatial Activity Analysis: %s (%s / %s)', ...
            mouse_id(1:3), upper(condition), stimulus), ...
            'FontSize', 16, 'FontWeight', 'bold');
   
% ── Export spatial grids to Excel ──────────────────────────────────────
    excel_filename = sprintf('spatial_grids_%s_%s_%s.xlsx', ...
        mouse_id(1:3), condition, datestr(now, 'yyyymmdd'));
    
    % Sheet 1-3: the three grids (rows = y, cols = x)
    writematrix(before_grid,     excel_filename, 'Sheet', 'before_grid');
    writematrix(after_grid,      excel_filename, 'Sheet', 'after_grid');
    writematrix(difference_grid, excel_filename, 'Sheet', 'difference_grid');
    
    % Sheet 4: axis coordinates so grids can be reconstructed
    axes_table = table(x_edges(:), y_edges(:), ...
        'VariableNames', {'x_coords_cm', 'y_coords_cm'});
    writetable(axes_table, excel_filename, 'Sheet', 'axes_coords');
    
    fprintf('Spatial grids exported to: %s\n', excel_filename);

end

function bounds = calculateCommonBoundsFromFilteredData(data, discovery, end_frame, ...
                                                       col_x, col_y, col_dist, col_groom, col_eating, options)
    % Calculate common spatial bounds from filtered data of both periods
    
    % Filter before period data
    if options.remove_grooming
        valid_before = data(1:discovery, col_dist) > options.dist_limit & ...
                      data(1:discovery, col_groom) == 0 & ...
                      data(1:discovery, col_eating) == 0;
    else
        valid_before = data(1:discovery, col_dist) > options.dist_limit & ...
                      data(1:discovery, col_eating) == 0;
    end
    
    % Filter after period data  
    if options.remove_grooming
        valid_after = data(discovery:end_frame, col_dist) > options.dist_limit & ...
                     data(discovery:end_frame, col_groom) == 0 & ...
                     data(discovery:end_frame, col_eating) == 0;
    else
        valid_after = data(discovery:end_frame, col_dist) > options.dist_limit & ...
                     data(discovery:end_frame, col_eating) == 0;
    end
    
    % Get filtered coordinates from both periods
    x_before = data(1:discovery, col_x);
    y_before = data(1:discovery, col_y);
    x_after = data(discovery:end_frame, col_x);
    y_after = data(discovery:end_frame, col_y);
    
    x_before_filt = x_before(valid_before);
    y_before_filt = y_before(valid_before);
    x_after_filt = x_after(valid_after);
    y_after_filt = y_after(valid_after);
    
    % Combine all filtered coordinates to find common bounds
    all_x = [x_before_filt; x_after_filt];
    all_y = [y_before_filt; y_after_filt];
    
    if isempty(all_x)
        bounds = struct('x_min', 0, 'x_max', 1, 'y_min', 0, 'y_max', 1);
    else
        bounds = struct('x_min', min(all_x), 'x_max', max(all_x), ...
                       'y_min', min(all_y), 'y_max', max(all_y));
    end
end

function spatial_grid = createSpatialGridWithCommonBounds(x_coords, y_coords, dff_values, bounds, resolution, sigma)
    % Create a smoothed spatial grid using specified bounds for coordinate alignment
    % Inputs:
    %   x_coords, y_coords: position coordinates
    %   dff_values: corresponding dF/F values (already z-scored)
    %   bounds: struct with x_min, x_max, y_min, y_max
    %   resolution: grid resolution
    %   sigma: Gaussian smoothing parameter (higher = more spreading to fill gaps)
    % Output:
    %   spatial_grid: smoothed 2D grid of mean dF/F values
    
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
    % Larger kernel to match larger sigma
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