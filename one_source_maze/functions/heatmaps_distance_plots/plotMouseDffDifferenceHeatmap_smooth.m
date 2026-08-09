function plotMouseDffDifferenceHeatmap_smooth(mice, mouse_index, options)
    % Plot ultra-smooth difference heatmap (after - before) for a single mouse
    % Uses histcounts2 + imgaussfilt + optional interpolation for smooth results
    % Shows how spatial activity patterns changed after food discovery
    %
    % Parameters:
    %   mice: cell array with mouse data
    %   mouse_index: index of the mouse to plot
    %   options: struct with options for filtering
    %     - dist_limit: minimum distance threshold (default: 5)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - colormap: colormap for the difference plot (default: cool)
    %     - bin_size: spatial bin size in cm (default: 4)
    %     - smooth_factor: Gaussian smoothing factor for imgaussfilt (default: 2.0)
    %     - difference_limits: [min max] color limits for difference heatmap (optional)
    %     - zscore_difference: boolean, whether to z-score the difference map (default: false)
    %     - interpolate_display: interpolate to finer grid for display (default: true)
    %     - interp_factor: interpolation factor (default: 4)
    
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

    if ~isfield(options, 'bin_size')
        options.bin_size = 4;  % 4 cm bins
    end
    
    if ~isfield(options, 'smooth_factor')
        options.smooth_factor = 2.0;  % imgaussfilt smoothing
    end
    
    if ~isfield(options, 'zscore_difference')
        options.zscore_difference = false;
    end
    
    if ~isfield(options, 'interpolate_display')
        options.interpolate_display = true;
    end
    
    if ~isfield(options, 'interp_factor')
        options.interp_factor = 4;  % 4x finer display grid
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
    
    % Create spatial bins
    x_edges = bounds.x_min:options.bin_size:bounds.x_max;
    y_edges = bounds.y_min:options.bin_size:bounds.y_max;
    x_centers = x_edges(1:end-1) + options.bin_size/2;
    y_centers = y_edges(1:end-1) + options.bin_size/2;
    
    fprintf('Spatial grid: %d x %d bins (%.1f cm bins)\n', length(x_centers), length(y_centers), options.bin_size);
    fprintf('Smoothing factor: %.1f\n', options.smooth_factor);
    if options.interpolate_display
        fprintf('Display interpolation: %dx finer grid\n', options.interp_factor);
    end
    
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
    
    % Create smooth spatial grid for before period
    before_grid = createSmoothHeatmap(x_before_filt, y_before_filt, dff_before_zscore, ...
                                     x_edges, y_edges, x_centers, y_centers, options);
    
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
    
    % Create smooth spatial grid for after period
    after_grid = createSmoothHeatmap(x_after_filt, y_after_filt, dff_after_zscore, ...
                                    x_edges, y_edges, x_centers, y_centers, options);
    
    %% COMPUTE DIFFERENCE (After - Before)
    % Compute difference
    difference_grid = after_grid - before_grid;
    
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
    
    % Prepare interpolated grids if requested
    if options.interpolate_display
        [X, Y] = meshgrid(x_centers, y_centers);
        x_fine = linspace(x_centers(1), x_centers(end), length(x_centers) * options.interp_factor);
        y_fine = linspace(y_centers(1), y_centers(end), length(y_centers) * options.interp_factor);
        [X_fine, Y_fine] = meshgrid(x_fine, y_fine);
        
        % Interpolate all grids
        before_grid_fine = interp2(X, Y, before_grid', X_fine, Y_fine, 'cubic');
        after_grid_fine = interp2(X, Y, after_grid', X_fine, Y_fine, 'cubic');
        difference_grid_fine = interp2(X, Y, difference_grid', X_fine, Y_fine, 'cubic');
        
        x_plot = x_fine;
        y_plot = y_fine;
    else
        before_grid_fine = before_grid';
        after_grid_fine = after_grid';
        difference_grid_fine = difference_grid';
        
        x_plot = x_centers;
        y_plot = y_centers;
    end
    
    % Subplot 1: Before period
    subplot(1, 3, 1);
    imagesc(x_plot, y_plot, before_grid_fine, 'AlphaData', ~isnan(before_grid_fine));
    axis xy equal tight;
    colormap(gca, options.colormap);
    c1 = colorbar;
    c1.Label.String = 'Z-scored \Delta F/F';
    c1.Label.FontSize = 12;
    title('Before Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11);
    set(gca, 'Color', 'white');
    box off
    axis off
    
    % Subplot 2: After period
    subplot(1, 3, 2);
    imagesc(x_plot, y_plot, after_grid_fine, 'AlphaData', ~isnan(after_grid_fine));
    axis xy equal tight;
    colormap(gca, options.colormap);
    c2 = colorbar;
    c2.Label.String = 'Z-scored \Delta F/F';
    c2.Label.FontSize = 12;
    title('After Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('X Position (cm)', 'FontSize', 12);
    ylabel('Y Position (cm)', 'FontSize', 12);
    set(gca, 'FontSize', 11);
    set(gca, 'Color', 'white');
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
    imagesc(x_plot, y_plot, difference_grid_fine, 'AlphaData', ~isnan(difference_grid_fine));
    axis xy equal tight;
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
    set(gca, 'Color', 'white');
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
    sgtitle(sprintf('Smooth Z-scored Spatial Activity Analysis: %s (%s / %s)', ...
            mouse_id(1:3), upper(condition), stimulus), ...
            'FontSize', 16, 'FontWeight', 'bold');
   
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

function heatmap = createSmoothHeatmap(x_data, y_data, dff_data, x_edges, y_edges, x_centers, y_centers, options)
    % Create smooth heatmap using histcounts2 + imgaussfilt
    % Inputs:
    %   x_data, y_data: position coordinates
    %   dff_data: corresponding dF/F values (already z-scored)
    %   x_edges, y_edges: bin edges
    %   x_centers, y_centers: bin centers
    %   options: struct with smooth_factor
    % Output:
    %   heatmap: smoothed 2D grid of mean dF/F values
    
    % Initialize grids
    dff_grid = zeros(length(x_centers), length(y_centers));
    count_grid = zeros(length(x_centers), length(y_centers));
    
    % Bin the data
    x_bin = discretize(x_data, x_edges);
    y_bin = discretize(y_data, y_edges);
    
    % Accumulate dF/F values in each bin
    for i = 1:length(x_data)
        if ~isnan(x_bin(i)) && ~isnan(y_bin(i))
            dff_grid(x_bin(i), y_bin(i)) = dff_grid(x_bin(i), y_bin(i)) + dff_data(i);
            count_grid(x_bin(i), y_bin(i)) = count_grid(x_bin(i), y_bin(i)) + 1;
        end
    end
    
    % Calculate mean dF/F per bin
    mean_dff = zeros(size(dff_grid));
    valid_bins = count_grid > 0;
    mean_dff(valid_bins) = dff_grid(valid_bins) ./ count_grid(valid_bins);
    
    % Set unvisited bins to NaN
    mean_dff(~valid_bins) = NaN;
    
    % Apply Gaussian smoothing using imgaussfilt
    if options.smooth_factor > 0
        % Replace NaN with 0 for smoothing
        smooth_input = mean_dff;
        smooth_input(isnan(smooth_input)) = 0;
        
        % Smooth the data
        smoothed = imgaussfilt(smooth_input, options.smooth_factor);
        
        % Smooth the mask to know where we have data
        mask_smooth = imgaussfilt(double(valid_bins), options.smooth_factor);
        
        % Only keep smoothed values where we have enough data
        smoothed(mask_smooth < 0.1) = NaN;
        
        heatmap = smoothed;
    else
        heatmap = mean_dff;
    end
end