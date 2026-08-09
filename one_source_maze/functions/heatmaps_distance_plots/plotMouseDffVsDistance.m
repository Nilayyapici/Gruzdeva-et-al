function plotMouseDffVsDistance(mice, mouse_index, options)
% Plot df/f vs distance for a single mouse, before and after food discovery
% Also plots time series data of df/f and distance
% Heatmaps now show z-scored dF/F values
%
% Parameters:
%   mice: cell array with mouse data
%   mouse_index: index of the mouse to plot
%   options: struct with options for filtering
%     - dist_limit: minimum distance threshold (default: 5)
%     - remove_grooming: boolean, whether to remove grooming periods (default: true)
%     - remove_eating: boolean, whether to remove eating periods (default: true)
%     - time_before_limits: [min max] time limits in seconds for before discovery plot (optional)
%     - time_after_limits: [min max] time limits in seconds for after discovery plot (optional)
%     - frame_rate: acquisition rate in Hz (default: 10)
%     - colormap_limits: [min max] color limits for DFF heatmaps (optional)

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

% Option to remove grooming periods (default: true)
if ~isfield(options, 'remove_grooming')
    options.remove_grooming = true;
end

% Option to remove eating periods (default: true)
if ~isfield(options, 'remove_eating')
    options.remove_eating = true;
end

if ~isfield(options, 'colormap')
    options.colormap = cool;
end

if ~isfield(options, 'sigma')
    options.sigma = 2;
end

if ~isfield(options, 'resolution')
    options.resolution = 100;
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

% Create filtering masks for before discovery
valid_before = data(1:discovery, COL_DIST) > options.dist_limit;

if options.remove_grooming
    valid_before = valid_before & data(1:discovery, COL_GROOM) == 0;
end

if options.remove_eating
    valid_before = valid_before & data(1:discovery, COL_EATING) == 0;
end

% Get df/f and distance values before discovery
dff_before = data(1:discovery, COL_DFF);
dist_before = data(1:discovery, COL_DIST);

% Filter valid points
dff_before_valid = dff_before(valid_before);
dist_before_valid = dist_before(valid_before);

% Calculate correlation for before
[rho_before, pval_before] = corr(dff_before_valid, dist_before_valid, 'Type', 'Pearson');

% Create filtering masks for after discovery
valid_after = data(discovery:end_frame, COL_DIST) > options.dist_limit;

if options.remove_grooming
    valid_after = valid_after & data(discovery:end_frame, COL_GROOM) == 0;
end

if options.remove_eating
    valid_after = valid_after & data(discovery:end_frame, COL_EATING) == 0;
end

% Get df/f and distance values after discovery
dff_after = data(discovery:end_frame, COL_DFF);
dist_after = data(discovery:end_frame, COL_DIST);

% Filter valid points
dff_after_valid = dff_after(valid_after);
dist_after_valid = dist_after(valid_after);

% Calculate correlation for after
[rho_after, pval_after] = corr(dff_after_valid, dist_after_valid, 'Type', 'Pearson');

% Create time vectors using specified frame rate
frame_rate = 1/mean(diff(data(:, COL_TIME))); % Hz
time_before = (1:discovery) / frame_rate;
time_after = (1:length(data(discovery:end_frame, 1))) / frame_rate;

% Create a figure with 2 rows of subplots with custom widths
figure('Position', [100, 100, 1500, 800]);

% Create custom subplot positions for the layout with heatmaps on left
% Format: [left bottom width height]
pos_heatmap_before = [0.05, 0.55, 0.20, 0.35];   % Top left - heatmap before
pos_time_before = [0.30, 0.55, 0.45, 0.35];      % Top middle - time series before
pos_corr_before = [0.83, 0.55, 0.15, 0.35];      % Top right - correlation before

pos_heatmap_after = [0.05, 0.10, 0.20, 0.35];    % Bottom left - heatmap after
pos_time_after = [0.30, 0.10, 0.45, 0.35];       % Bottom middle - time series after
pos_corr_after = [0.83, 0.10, 0.15, 0.35];       % Bottom right - correlation after

%% TOP ROW: Before discovery plots

% Heatmap of x,y activity (before) - TOP LEFT
heatmap_before_axes = subplot('Position', pos_heatmap_before);

% Extract x,y position and dff data before discovery (FILTERED)
x_before_all = data(1:discovery, COL_X);
y_before_all = data(1:discovery, COL_Y);
dff_before_all = data(1:discovery, COL_DFF);

% Apply the same filtering to heatmap data
x_before = x_before_all(valid_before);
y_before = y_before_all(valid_before);
dff_before_filtered = dff_before_all(valid_before);

% Z-score the filtered dF/F values for the before period
if ~isempty(dff_before_filtered)
    dff_before_zscore = (dff_before_filtered - mean(dff_before_filtered)) / std(dff_before_filtered);

    % Create a grid for visualization
    resolution = options.resolution; % Grid size
    x_min = min(x_before);
    x_max = max(x_before);
    y_min = min(y_before);
    y_max = max(y_before);
    x_edges = linspace(x_min, x_max, resolution);
    y_edges = linspace(y_min, y_max, resolution);

    % Initialize grid for averaged values and visit count
    dff_grid = zeros(resolution, resolution);
    visit_count = zeros(resolution, resolution);

    % Calculate which grid cell each data point belongs to
    x_idx = round((x_before - x_min) / (x_max - x_min) * (resolution-1)) + 1;
    y_idx = round((y_before - y_min) / (y_max - y_min) * (resolution-1)) + 1;

    % Enforce valid indices
    x_idx = max(1, min(resolution, x_idx));
    y_idx = max(1, min(resolution, y_idx));

    % Sum z-scored dff values and count visits for each grid cell
    for i = 1:length(x_before)
        xi = x_idx(i);
        yi = y_idx(i);
        dff_grid(yi, xi) = dff_grid(yi, xi) + dff_before_zscore(i);
        visit_count(yi, xi) = visit_count(yi, xi) + 1;
    end

    % Calculate mean z-scored dff for each cell and create a mask for visited locations
    visited_mask = visit_count > 0;
    mean_dff_grid = zeros(size(dff_grid));
    mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);

    % Create a masked array with NaNs for unvisited locations
    masked_dff = mean_dff_grid;
    masked_dff(~visited_mask) = NaN;

    % Apply Gaussian smoothing
    masked_dff_before = masked_dff;

    % Create a Gaussian kernel
    kernel_size = 7; % Must be odd
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*options.sigma^2));
    kernel = kernel / sum(kernel(:));

    % Apply the smoothing while preserving the mask
    smoothed_dff = NaN(size(masked_dff));
    valid_mask = ~isnan(masked_dff);

    % Pad the arrays to handle edges
    padded_dff = padarray(masked_dff, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);

    % Convolve with the kernel
    for i = 1+floor(kernel_size/2):size(padded_dff,1)-floor(kernel_size/2)
        for j = 1+floor(kernel_size/2):size(padded_dff,2)-floor(kernel_size/2)
            if padded_mask(i, j)
                window = padded_dff(i-floor(kernel_size/2):i+floor(kernel_size/2), j-floor(kernel_size/2):j+floor(kernel_size/2));
                window_mask = ~isnan(window);
                if sum(window_mask(:)) > 0
                    weighted_values = window .* kernel;
                    weighted_values(~window_mask) = 0;
                    % Normalize by the valid weights
                    weight_sum = sum(kernel(window_mask));
                    if weight_sum > 0
                        smoothed_dff(i-floor(kernel_size/2), j-floor(kernel_size/2)) = sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end

    % Plot the smoothed heatmap
    imagesc(x_edges, y_edges, smoothed_dff, 'AlphaData', ~isnan(smoothed_dff));

    % Store the smoothed data for color scaling later
    masked_dff_before_smoothed = smoothed_dff;
else
    % No valid data for heatmap
    masked_dff_before_smoothed = NaN;
end

% Format heatmap before plot
axis equal tight
colormap(heatmap_before_axes, options.colormap);
c = colorbar;
c.Label.String = 'Z-scored \Delta F/F';
c.Label.FontSize = 12;
xlabel('X Position (cm)', 'FontSize', 12);
ylabel('Y Position (cm)', 'FontSize', 12);
title('Before Food Discovery: Z-scored DFF Heatmap', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box off;
axis off;

% Time series plot (before) - TOP MIDDLE
before_axes = subplot('Position', pos_time_before);

% Filter time series data
time_plot_before = time_before(valid_before);
dff_plot_before = dff_before(valid_before);
dist_plot_before = dist_before(valid_before);

if ~isempty(time_plot_before)
    % Plot data vs time (left y-axis)
    yyaxis left;
    plot(time_plot_before, dff_plot_before, 'LineWidth', 1.5, 'Color', [0, 0.4, 0.8]);
    ylabel('\Delta F/F', 'FontSize', 12, 'Color', [0, 0.4, 0.8]);
    % Set y-axis color
    ax = gca;
    ax.YAxis(1).Color = [0, 0.4, 0.8];

    % Plot data vs time (right y-axis)
    yyaxis right;
    plot(time_plot_before, dist_plot_before, 'LineWidth', 1.5, 'Color', [0.9, 0.3, 0.5]);
    ylabel('Distance to Food (cm)', 'FontSize', 12, 'Color', [0.9, 0.3, 0.5]);
    % Set y-axis color
    ax.YAxis(2).Color = [0.9, 0.3, 0.5];
end

% Format time series before plot
xlabel('Time (s)', 'FontSize', 12);
title('Before Food Discovery: df/f and Distance vs Time', 'FontSize', 14, 'FontWeight', 'bold');
grid off;
box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Correlation scatter plot (before) - TOP RIGHT
corr_before_axes = subplot('Position', pos_corr_before);
if ~isempty(dist_before_valid)
    scatter(dist_before_valid, dff_before_valid, 40, 'filled', 'MarkerFaceColor', [0.5, 0.7, 0.9], 'MarkerFaceAlpha', 0.7);
    hold on;

    % Add regression line
    if length(dist_before_valid) > 1
        p = polyfit(dist_before_valid, dff_before_valid, 1);
        x_range = linspace(min(dist_before_valid), max(dist_before_valid), 100);
        y_fit = polyval(p, x_range);
        plot(x_range, y_fit, 'LineWidth', 2, 'Color', [0.3, 0.5, 0.7]);
    end

    % Add correlation coefficient
    if ~isnan(rho_before)
        text(0.05, 0.95, sprintf('r = %.3f, p = %.3f', rho_before, pval_before), ...
            'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end
end

% Format before correlation plot
xlabel('Distance to Food (cm)', 'FontSize', 12);
ylabel('\Delta F/F', 'FontSize', 12);
title('Before Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

%% BOTTOM ROW: After discovery plots

% Heatmap of x,y activity (after) - BOTTOM LEFT
heatmap_after_axes = subplot('Position', pos_heatmap_after);

% Extract x,y position and dff data after discovery (FILTERED)
x_after_all = data(discovery:end_frame, COL_X);
y_after_all = data(discovery:end_frame, COL_Y);
dff_after_all = data(discovery:end_frame, COL_DFF);

% Apply the same filtering to heatmap data
x_after = x_after_all(valid_after);
y_after = y_after_all(valid_after);
dff_after_filtered = dff_after_all(valid_after);

% Z-score the filtered dF/F values for the after period
if ~isempty(dff_after_filtered)
    dff_after_zscore = (dff_after_filtered - mean(dff_after_filtered)) / std(dff_after_filtered);

    % Create a grid for visualization
    resolution = options.resolution; % Grid size
    x_min = min(x_after);
    x_max = max(x_after);
    y_min = min(y_after);
    y_max = max(y_after);
    x_edges = linspace(x_min, x_max, resolution);
    y_edges = linspace(y_min, y_max, resolution);

    % Initialize grid for averaged values and visit count
    dff_grid = zeros(resolution, resolution);
    visit_count = zeros(resolution, resolution);

    % Calculate which grid cell each data point belongs to
    x_idx = round((x_after - x_min) / (x_max - x_min) * (resolution-1)) + 1;
    y_idx = round((y_after - y_min) / (y_max - y_min) * (resolution-1)) + 1;

    % Enforce valid indices
    x_idx = max(1, min(resolution, x_idx));
    y_idx = max(1, min(resolution, y_idx));

    % Sum z-scored dff values and count visits for each grid cell
    for i = 1:length(x_after)
        xi = x_idx(i);
        yi = y_idx(i);
        dff_grid(yi, xi) = dff_grid(yi, xi) + dff_after_zscore(i);
        visit_count(yi, xi) = visit_count(yi, xi) + 1;
    end

    % Calculate mean z-scored dff for each cell and create a mask for visited locations
    visited_mask = visit_count > 0;
    mean_dff_grid = zeros(size(dff_grid));
    mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);

    % Create a masked array with NaNs for unvisited locations
    masked_dff = mean_dff_grid;
    masked_dff(~visited_mask) = NaN;

    % Apply Gaussian smoothing
    masked_dff_after = masked_dff;

    % Create a Gaussian kernel
    kernel_size = 7; % Must be odd
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*options.sigma^2));
    kernel = kernel / sum(kernel(:));

    % Apply the smoothing while preserving the mask
    smoothed_dff = NaN(size(masked_dff));
    valid_mask = ~isnan(masked_dff);

    % Pad the arrays to handle edges
    padded_dff = padarray(masked_dff, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);

    % Convolve with the kernel
    for i = 1+floor(kernel_size/2):size(padded_dff,1)-floor(kernel_size/2)
        for j = 1+floor(kernel_size/2):size(padded_dff,2)-floor(kernel_size/2)
            if padded_mask(i, j)
                window = padded_dff(i-floor(kernel_size/2):i+floor(kernel_size/2), j-floor(kernel_size/2):j+floor(kernel_size/2));
                window_mask = ~isnan(window);
                if sum(window_mask(:)) > 0
                    weighted_values = window .* kernel;
                    weighted_values(~window_mask) = 0;
                    % Normalize by the valid weights
                    weight_sum = sum(kernel(window_mask));
                    if weight_sum > 0
                        smoothed_dff(i-floor(kernel_size/2), j-floor(kernel_size/2)) = sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end

    % Plot the smoothed heatmap
    imagesc(x_edges, y_edges, smoothed_dff, 'AlphaData', ~isnan(smoothed_dff));

    % Store the smoothed data for color scaling
    masked_dff_after_smoothed = smoothed_dff;
else
    % No valid data for heatmap
    masked_dff_after_smoothed = NaN;
end

% Format heatmap after plot
axis equal tight
colormap(heatmap_after_axes, options.colormap);
c = colorbar;
c.Label.String = 'Z-scored \Delta F/F';
c.Label.FontSize = 12;
xlabel('X Position (cm)', 'FontSize', 12);
ylabel('Y Position (cm)', 'FontSize', 12);
title('After Food Discovery: Z-scored DFF Heatmap', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box off;
axis off;

% Time series plot (after) - BOTTOM MIDDLE
after_axes = subplot('Position', pos_time_after);

% Filter time series data
time_plot_after = time_after(valid_after);
dff_plot_after = dff_after(valid_after);
dist_plot_after = dist_after(valid_after);

if ~isempty(time_plot_after)
    % Plot data vs time (left y-axis)
    yyaxis left;
    plot(time_plot_after, dff_plot_after, 'LineWidth', 1.5, 'Color', [0, 0.4, 0.8]);
    ylabel('\Delta F/F', 'FontSize', 12, 'Color', [0, 0.4, 0.8]);
    % Set y-axis color
    ax = gca;
    ax.YAxis(1).Color = [0, 0.4, 0.8];

    % Plot data vs time (right y-axis)
    yyaxis right;
    plot(time_plot_after, dist_plot_after, 'LineWidth', 1.5, 'Color', [0.9, 0.3, 0.5]);
    ylabel('Distance to Food (cm)', 'FontSize', 12, 'Color', [0.9, 0.3, 0.5]);
    % Set y-axis color
    ax.YAxis(2).Color = [0.9, 0.3, 0.5];
end

% Format time series after plot
xlabel('Time (s)', 'FontSize', 12);
title('After Food Discovery: df/f and Distance vs Time', 'FontSize', 14, 'FontWeight', 'bold');
grid off;
box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Correlation scatter plot (after) - BOTTOM RIGHT
corr_after_axes = subplot('Position', pos_corr_after);
if ~isempty(dist_after_valid)
    scatter(dist_after_valid, dff_after_valid, 40, 'filled', 'MarkerFaceColor', [0.1, 0.3, 0.7], 'MarkerFaceAlpha', 0.7);
    hold on;

    % Add regression line
    if length(dist_after_valid) > 1
        p = polyfit(dist_after_valid, dff_after_valid, 1);
        x_range = linspace(min(dist_after_valid), max(dist_after_valid), 100);
        y_fit = polyval(p, x_range);
        plot(x_range, y_fit, 'LineWidth', 2, 'Color', [0, 0.2, 0.5]);
    end

    % Add correlation coefficient
    if ~isnan(rho_after)
        text(0.05, 0.95, sprintf('r = %.3f, p = %.3f', rho_after, pval_after), ...
            'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end
end

% Format after correlation plot
xlabel('Distance to Food (cm)', 'FontSize', 12);
ylabel('\Delta F/F', 'FontSize', 12);
title('After Food Discovery', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Create common axis limits for both heatmaps to ensure consistency
% First, we'll find the overall x/y limits
x_all = [x_before; x_after];
y_all = [y_before; y_after];
if ~isempty(x_all) && ~isempty(y_all)
    global_x_min = min(x_all);
    global_x_max = max(x_all);
    global_y_min = min(y_all);
    global_y_max = max(y_all);

    % Reset the before heatmap with global limits
    axes(heatmap_before_axes);
    xlim([global_x_min, global_x_max]);
    ylim([global_y_min, global_y_max]);

    % Reset the after heatmap with global limits
    axes(heatmap_after_axes);
    xlim([global_x_min, global_x_max]);
    ylim([global_y_min, global_y_max]);
end

% Set color limits for both DFF heatmaps
if isfield(options, 'colormap_limits') && length(options.colormap_limits) == 2
    % Use user-specified colormap limits
    global_clim = options.colormap_limits;
else
    % Use automatic limits based on data (existing behavior)
    all_values = [masked_dff_before_smoothed(:); masked_dff_after_smoothed(:)];
    all_values = all_values(~isnan(all_values)); % Remove NaNs
    if ~isempty(all_values)
        global_clim = [min(all_values), max(all_values)];
    else
        global_clim = [-1, 1]; % Default fallback
    end
end

% Apply the color limits to both heatmaps
axes(heatmap_before_axes);
caxis(global_clim);

axes(heatmap_after_axes);
caxis(global_clim);

% Make y-axis limits the same for both correlation plots
y_all = [dff_before_valid; dff_after_valid];
if ~isempty(y_all)
    y_min = min(y_all) - 0.1 * (max(y_all) - min(y_all));
    y_max = max(y_all) + 0.1 * (max(y_all) - min(y_all));
    axes(corr_before_axes); ylim([y_min, y_max]);
    axes(corr_after_axes); ylim([y_min, y_max]);
end

% Make x-axis limits the same for both correlation plots
x_all = [dist_before_valid; dist_after_valid];
if ~isempty(x_all)
    x_min = min(x_all) - 0.1 * (max(x_all) - min(x_all));
    x_max = max(x_all) + 0.1 * (max(x_all) - min(x_all));
    axes(corr_before_axes); xlim([x_min, x_max]);
    axes(corr_after_axes); xlim([x_min, x_max]);
end

% Set custom xlim for before time series plot if provided
if isfield(options, 'time_before_limits') && length(options.time_before_limits) == 2
    axes(before_axes);
    xlim(options.time_before_limits);
end

% Set custom xlim for after time series plot if provided
if isfield(options, 'time_after_limits') && length(options.time_after_limits) == 2
    axes(after_axes);
    xlim(options.time_after_limits);
end

% Add main title for the whole figure
filter_text = '';
if options.remove_grooming && options.remove_eating
    filter_text = ' (grooming & eating excluded)';
elseif options.remove_grooming
    filter_text = ' (grooming excluded)';
elseif options.remove_eating
    filter_text = ' (eating excluded)';
end

% ── Export time series data to Excel ───────────────────────────────────
excel_filename = sprintf('timeseries_%s_%s_%s_%s.xlsx', ...
    mouse_id(1:3), condition, stimulus, datestr(now, 'yyyymmdd'));

% Apply time window limits to before data (same as xlim applied to plot)
if isfield(options, 'time_before_limits') && length(options.time_before_limits) == 2
    t_mask_before = time_plot_before >= options.time_before_limits(1) & ...
        time_plot_before <= options.time_before_limits(2);
else
    t_mask_before = true(size(time_plot_before));
end

% Apply time window limits to after data
if isfield(options, 'time_after_limits') && length(options.time_after_limits) == 2
    t_mask_after = time_plot_after >= options.time_after_limits(1) & ...
        time_plot_after <= options.time_after_limits(2);
else
    t_mask_after = true(size(time_plot_after));
end

    before_table = table(...
        reshape(time_plot_before(t_mask_before),  [], 1), ...
        reshape(dff_plot_before(t_mask_before),   [], 1), ...
        reshape(dist_plot_before(t_mask_before),  [], 1), ...
        'VariableNames', {'Time_s', 'dFF_percent', 'Distance_to_food_cm'});
    
    after_table = table(...
        reshape(time_plot_after(t_mask_after),  [], 1), ...
        reshape(dff_plot_after(t_mask_after),   [], 1), ...
        reshape(dist_plot_after(t_mask_after),  [], 1), ...
        'VariableNames', {'Time_s', 'dFF_percent', 'Distance_to_food_cm'});

writetable(before_table, excel_filename, 'Sheet', 'before_discovery');
writetable(after_table,  excel_filename, 'Sheet', 'after_discovery');

fprintf('Time series data exported to: %s\n', excel_filename);

sgtitle({
    sprintf('Mouse %s: %s / %s%s', mouse_id(1:3), upper(condition), stimulus, filter_text), ...
    }, 'FontSize', 16, 'FontWeight', 'bold');
end