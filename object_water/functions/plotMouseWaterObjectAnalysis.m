function plotMouseWaterObjectAnalysis(mice_all, mouse_prefix, options)
% Plot comprehensive analysis for water and object sessions for a single mouse
% Creates four figures:
% 1. Water analysis: before, water, difference (water-before)
% 2. Object analysis: before, object, difference (object-before)
% 3. Detailed water analysis: before vs water (2 rows with heatmap, time series, correlation)
% 4. Detailed object analysis: before vs object (2 rows with heatmap, time series, correlation)
%
% Parameters:
%   mice_all: cell array with all mouse data
%   mouse_prefix: string prefix for mouse ID (e.g., 'F29', 'M28')
%   options: struct with options for filtering and plotting
%            - time_before_limits: [start_sec, end_sec] for before session time series
%            - time_water_limits: [start_sec, end_sec] for water session time series
%            - time_object_limits: [start_sec, end_sec] for object session time series

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
if ~isfield(options, 'sigma'), options.sigma = 4; end  % Increased from 2 to 4 for sparse data
if ~isfield(options, 'resolution'), options.resolution = 100; end
if ~isfield(options, 'time_before_limits'), options.time_before_limits = []; end
if ~isfield(options, 'time_water_limits'), options.time_water_limits = []; end
if ~isfield(options, 'time_object_limits'), options.time_object_limits = []; end

% Find sessions for this mouse
mouse_indices = find(contains(mice_all(:,1), mouse_prefix));

before_idx = []; water_idx = []; object_idx = [];

for i = mouse_indices'
    session_name = mice_all{i,1};
    if contains(session_name, 'before')
        before_idx = i;
    elseif contains(session_name, 'water')
        water_idx = i;
    elseif contains(session_name, 'object')
        object_idx = i;
    end
end

if isempty(before_idx) || isempty(water_idx) || isempty(object_idx)
    error('Could not find all required sessions (before, water, object) for mouse %s', mouse_prefix);
end

% Get session data
before_data = mice_all{before_idx, 4};
water_data = mice_all{water_idx, 4};
object_data = mice_all{object_idx, 4};
condition = mice_all{before_idx, 2};
stimulus = mice_all{before_idx, 3};

%% FIGURE 1: Water Analysis
figure('Position', [100, 100, 1600, 600]);
sgtitle(sprintf('Water Analysis: %s (%s / %s)', mouse_prefix, upper(condition), stimulus), ...
    'FontSize', 16, 'FontWeight', 'bold');

% Calculate common spatial bounds for water analysis
water_bounds = calculateCommonBounds(before_data, water_data, COL_X, COL_Y, COL_DIST_WATER, COL_GROOM, options);

% Before water heatmap
subplot(1, 3, 1);
before_grid_water = createSpatialGridWithBounds(before_data, COL_X, COL_Y, COL_DFF, COL_DIST_WATER, COL_GROOM, options, water_bounds);
plotHeatmap(before_grid_water, before_data, COL_X, COL_Y, options.colormap);
title('Before Water Discovery', 'FontSize', 14, 'FontWeight', 'bold');

% Water session heatmap
subplot(1, 3, 2);
water_grid = createSpatialGridWithBounds(water_data, COL_X, COL_Y, COL_DFF, COL_DIST_WATER, COL_GROOM, options, water_bounds);
plotHeatmap(water_grid, water_data, COL_X, COL_Y, options.colormap);
title('Water Session', 'FontSize', 14, 'FontWeight', 'bold');

% Difference (Water - Before)
subplot(1, 3, 3);
diff_grid_water = computeDifference(water_grid, before_grid_water);
plotDifference(diff_grid_water,  before_data, water_data,  COL_X, COL_Y, options.colormap, 'water',  mouse_prefix);
title('Difference (Water - Before)', 'FontSize', 14, 'FontWeight', 'bold');

% Sync color limits for before and water plots
syncColorLimits([1, 2], {before_grid_water, water_grid});

%% FIGURE 2: Object Analysis
figure('Position', [200, 200, 1600, 600]);
sgtitle(sprintf('Object Analysis: %s (%s / %s)', mouse_prefix, upper(condition), stimulus), ...
    'FontSize', 16, 'FontWeight', 'bold');

% Calculate common spatial bounds for object analysis
object_bounds = calculateCommonBounds(before_data, object_data, COL_X, COL_Y, COL_DIST_OBJECT, COL_GROOM, options);

% Before object heatmap
subplot(1, 3, 1);
before_grid_object = createSpatialGridWithBounds(before_data, COL_X, COL_Y, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, object_bounds);
plotHeatmap(before_grid_object, before_data, COL_X, COL_Y, options.colormap);
title('Before Object Discovery', 'FontSize', 14, 'FontWeight', 'bold');

% Object session heatmap
subplot(1, 3, 2);
object_grid = createSpatialGridWithBounds(object_data, COL_X, COL_Y, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, object_bounds);
plotHeatmap(object_grid, object_data, COL_X, COL_Y, options.colormap);
title('Object Session', 'FontSize', 14, 'FontWeight', 'bold');

% Difference (Object - Before)
subplot(1, 3, 3);
diff_grid_object = computeDifference(object_grid, before_grid_object);
plotDifference(diff_grid_object, before_data, object_data, COL_X, COL_Y, options.colormap, 'object', mouse_prefix);
title('Difference (Object - Before)', 'FontSize', 14, 'FontWeight', 'bold');

% Sync color limits for before and object plots
syncColorLimits([1, 2], {before_grid_object, object_grid});

%% FIGURE 3: Detailed Water Analysis
plotDetailedWaterAnalysis(before_data, water_data, mouse_prefix, condition, stimulus, options);

%% FIGURE 4: Detailed Object Analysis
plotDetailedObjectAnalysis(before_data, object_data, mouse_prefix, condition, stimulus, options);
end

function spatial_grid = createSpatialGrid(data, col_x, col_y, col_dff, col_dist, col_groom, options)
% Create spatial grid of z-scored dF/F values

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

% Create spatial bounds
if isempty(x_filt)
    spatial_grid = NaN(options.resolution, options.resolution);
    return;
end

x_min = min(x_filt); x_max = max(x_filt);
y_min = min(y_filt); y_max = max(y_filt);

% Initialize grids
resolution = options.resolution;
dff_grid = zeros(resolution, resolution);
visit_count = zeros(resolution, resolution);

% Map coordinates to grid indices
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

function plotCorrelation(data, col_dff, col_dist, col_groom, options, title_str, color, target_type)
% Plot correlation scatter

% Apply filtering
if options.remove_grooming && size(data, 2) >= col_groom
    valid_mask = data(:, col_dist) > options.dist_limit & data(:, col_groom) == 0;
else
    valid_mask = data(:, col_dist) > options.dist_limit;
end

dff_valid = data(valid_mask, col_dff);
dist_valid = data(valid_mask, col_dist);

if length(dff_valid) > 1
    [rho, pval] = corr(dff_valid, dist_valid, 'Type', 'Pearson');

    scatter(dist_valid, dff_valid, 40, 'filled', 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.7);
    hold on;

    % Regression line
    p = polyfit(dist_valid, dff_valid, 1);
    x_range = linspace(min(dist_valid), max(dist_valid), 100);
    y_fit = polyval(p, x_range);
    plot(x_range, y_fit, 'LineWidth', 2, 'Color', color*0.7);

    % Correlation info
    text(0.05, 0.95, sprintf('r = %.3f\np = %.3f', rho, pval), ...
        'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

if strcmp(target_type, 'Water')
    xlabel('Distance to Water (cm)', 'FontSize', 12);
else
    xlabel('Distance to Object (cm)', 'FontSize', 12);
end
ylabel('\Delta F/F', 'FontSize', 12);
title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
grid off; box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
end

function smoothed_grid = applySmoothing(input_grid, sigma)
% Apply Gaussian smoothing with kernel size scaling with sigma
kernel_size = 2*ceil(3*sigma) + 1; % Scale kernel with sigma
kernel_size = min(kernel_size, 21); % Cap at reasonable size

[X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), ...
    -floor(kernel_size/2):floor(kernel_size/2));
kernel = exp(-(X.^2 + Y.^2) / (2*sigma^2));
kernel = kernel / sum(kernel(:));

smoothed_grid = NaN(size(input_grid));
valid_mask = ~isnan(input_grid);

% Pad arrays
padded_grid = padarray(input_grid, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
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

function plotHeatmap(spatial_grid, data, col_x, col_y, colormap_choice)
% Plot spatial heatmap with crisp nearest neighbor interpolation
x_coords = data(:, col_x);
y_coords = data(:, col_y);
x_min = min(x_coords); x_max = max(x_coords);
y_min = min(y_coords); y_max = max(y_coords);

h = imagesc([x_min x_max], [y_min y_max], spatial_grid, 'AlphaData', ~isnan(spatial_grid));
h.Interpolation = 'nearest';  % Keep crisp appearance
axis equal tight
colormap(gca, colormap_choice);
c = colorbar;
c.Label.String = 'Z-scored \Delta F/F';
c.Label.FontSize = 12;
xlabel('X Position (cm)', 'FontSize', 12);
ylabel('Y Position (cm)', 'FontSize', 12);
set(gca, 'FontSize', 11);
box off; axis off;
end

function diff_grid = computeDifference(after_grid, before_grid)
% Compute difference between two spatial grids
valid_both = ~isnan(before_grid) & ~isnan(after_grid);
diff_grid = NaN(size(before_grid));
diff_grid(valid_both) = after_grid(valid_both) - before_grid(valid_both);
end

function plotDifference(diff_grid, before_data, after_data, col_x, col_y, colormap_choice, session_type, mouse_prefix)
% Plot difference grid with crisp nearest neighbor interpolation
x_before = before_data(:, col_x); y_before = before_data(:, col_y);
x_after = after_data(:, col_x); y_after = after_data(:, col_y);

x_min = min([x_before; x_after]); x_max = max([x_before; x_after]);
y_min = min([y_before; y_after]); y_max = max([y_before; y_after]);

h = imagesc([x_min x_max], [y_min y_max], diff_grid, 'AlphaData', ~isnan(diff_grid));
h.Interpolation = 'nearest';  % Keep crisp appearance
axis equal tight
colormap(gca, colormap_choice);

% Set symmetric color limits
diff_values = diff_grid(~isnan(diff_grid));
if ~isempty(diff_values)
    max_abs = max(abs(diff_values));
    caxis([-0.8, 0.8]);
end

c = colorbar;
c.Label.String = '\Delta(Z-scored \Delta F/F)';
c.Label.FontSize = 12;
xlabel('X Position (cm)', 'FontSize', 12);
ylabel('Y Position (cm)', 'FontSize', 12);
set(gca, 'FontSize', 11);
box off; axis off;

% ── Export difference grid data to Excel ──────────────────────────────
% Build x and y axis vectors
n_rows = size(diff_grid, 1);
n_cols = size(diff_grid, 2);
x_coords = linspace(x_min, x_max, n_cols);
y_coords = linspace(y_min, y_max, n_rows);

% Save diff_grid as matrix sheet + axis coordinates sheet
excel_filename = sprintf('difference_grid_%s_%s_%s.xlsx', mouse_prefix, session_type, datestr(now, 'yyyymmdd'));

% Sheet 1: the difference grid itself (rows = y, cols = x)
writematrix(diff_grid, excel_filename, 'Sheet', 'diff_grid');

% Sheet 2: x and y axis coordinates so grid can be reconstructed
axes_table = table(x_coords(:), y_coords(:), ...
    'VariableNames', {'x_coords_cm', 'y_coords_cm'});
writetable(axes_table, excel_filename, 'Sheet', 'axes_coords');

fprintf('Difference grid exported to: %s\n', excel_filename);
end

function syncColorLimits(subplot_nums, grids)
% Synchronize color limits across subplots
all_values = [];
for i = 1:length(grids)
    grid_vals = grids{i};
    all_values = [all_values; grid_vals(~isnan(grid_vals))];
end

if ~isempty(all_values)
    clim_range = [-0.8, 0.8];
    for i = 1:length(subplot_nums)
        subplot(1, 3, subplot_nums(i));
        caxis(clim_range);
    end
end
end

function plotDetailedWaterAnalysis(before_data, water_data, mouse_prefix, condition, stimulus, options)
% Create detailed 2-row water analysis plot

COL_TIME = 1; COL_X = 2; COL_Y = 3; COL_DFF = 6;
COL_DIST_WATER = 9; COL_GROOM = 13;

figure('Position', [300, 300, 1500, 800]);
excel_filename_water = sprintf('timeseries_%s_water_%s.xlsx', mouse_prefix, datestr(now, 'yyyymmdd'));

% Row positions for 2x3 layout
row_height = 0.35;
pos_before = [0.05, 0.55, 0.20, row_height; 0.30, 0.55, 0.45, row_height; 0.83, 0.55, 0.15, row_height];
pos_water = [0.05, 0.10, 0.20, row_height; 0.30, 0.10, 0.45, row_height; 0.83, 0.10, 0.15, row_height];

% Before session (Row 1)
subplot('Position', pos_before(1,:));
before_grid = createSpatialGrid(before_data, COL_X, COL_Y, COL_DFF, COL_DIST_WATER, COL_GROOM, options);
plotHeatmap(before_grid, before_data, COL_X, COL_Y, options.colormap);
title('Before: Z-scored dF/F', 'FontSize', 14, 'FontWeight', 'bold');

subplot('Position', pos_before(2,:));
plotTimeSeries(before_data, COL_TIME, COL_DFF, COL_DIST_WATER, COL_GROOM, options, ...
    'Before Session', options.time_before_limits, 'Water', mouse_prefix, excel_filename_water);


subplot('Position', pos_before(3,:));
plotCorrelation(before_data, COL_DFF, COL_DIST_WATER, COL_GROOM, options, ...
    'Before', [0.5, 0.7, 0.9], 'Water');

% Water session (Row 2)
subplot('Position', pos_water(1,:));
water_grid = createSpatialGrid(water_data, COL_X, COL_Y, COL_DFF, COL_DIST_WATER, COL_GROOM, options);
plotHeatmap(water_grid, water_data, COL_X, COL_Y, options.colormap);
title('Water: Z-scored dF/F', 'FontSize', 14, 'FontWeight', 'bold');

subplot('Position', pos_water(2,:));
plotTimeSeries(water_data, COL_TIME, COL_DFF, COL_DIST_WATER, COL_GROOM, options, ...
    'Water Session', options.time_water_limits, 'Water', mouse_prefix, excel_filename_water);

subplot('Position', pos_water(3,:));
plotCorrelation(water_data, COL_DFF, COL_DIST_WATER, COL_GROOM, options, ...
    'Water', [0.1, 0.3, 0.7], 'Water');

% Sync heatmap color limits
all_values = [before_grid(:); water_grid(:)];
all_values = all_values(~isnan(all_values));
if ~isempty(all_values)
    % clim_range = [min(all_values), max(all_values)];
    clim_range = [-0.8, 0.8];
    subplot('Position', pos_before(1,:)); caxis(clim_range);
    subplot('Position', pos_water(1,:)); caxis(clim_range);
end

sgtitle(sprintf('Detailed Water Analysis: %s (%s / %s)', mouse_prefix, upper(condition), stimulus), ...
    'FontSize', 16, 'FontWeight', 'bold');
end

function plotDetailedObjectAnalysis(before_data, object_data, mouse_prefix, condition, stimulus, options)
% Create detailed 2-row object analysis plot

COL_TIME = 1; COL_X = 2; COL_Y = 3; COL_DFF = 6;
COL_DIST_OBJECT = 10; COL_GROOM = 13;

figure('Position', [400, 400, 1500, 800]);
excel_filename_object = sprintf('timeseries_%s_object_%s.xlsx', mouse_prefix, datestr(now, 'yyyymmdd'));

% Row positions for 2x3 layout
row_height = 0.35;
pos_before = [0.05, 0.55, 0.20, row_height; 0.30, 0.55, 0.45, row_height; 0.83, 0.55, 0.15, row_height];
pos_object = [0.05, 0.10, 0.20, row_height; 0.30, 0.10, 0.45, row_height; 0.83, 0.10, 0.15, row_height];

% Before session (Row 1)
subplot('Position', pos_before(1,:));
before_grid = createSpatialGrid(before_data, COL_X, COL_Y, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options);
plotHeatmap(before_grid, before_data, COL_X, COL_Y, options.colormap);
title('Before: Z-scored dF/F', 'FontSize', 14, 'FontWeight', 'bold');

subplot('Position', pos_before(2,:));
plotTimeSeries(before_data, COL_TIME, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, ...
    'Before Session', options.time_before_limits, 'Object', mouse_prefix, excel_filename_object);


subplot('Position', pos_before(3,:));
plotCorrelation(before_data, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, ...
    'Before', [0.5, 0.7, 0.9], 'Object');

% Object session (Row 2)
subplot('Position', pos_object(1,:));
object_grid = createSpatialGrid(object_data, COL_X, COL_Y, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options);
plotHeatmap(object_grid, object_data, COL_X, COL_Y, options.colormap);
title('Object: Z-scored dF/F', 'FontSize', 14, 'FontWeight', 'bold');

subplot('Position', pos_object(2,:));
plotTimeSeries(object_data, COL_TIME, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, ...
    'Object Session', options.time_object_limits, 'Object', mouse_prefix, excel_filename_object);

subplot('Position', pos_object(3,:));
plotCorrelation(object_data, COL_DFF, COL_DIST_OBJECT, COL_GROOM, options, ...
    'Object', [0.8, 0.2, 0.4], 'Object');

% Sync heatmap color limits
all_values = [before_grid(:); object_grid(:)];
all_values = all_values(~isnan(all_values));
if ~isempty(all_values)
    % clim_range = [min(all_values), max(all_values)];
    clim_range = [-0.8, 0.8];
    subplot('Position', pos_before(1,:)); caxis(clim_range);
    subplot('Position', pos_object(1,:)); caxis(clim_range);
end

sgtitle(sprintf('Detailed Object Analysis: %s (%s / %s)', mouse_prefix, upper(condition), stimulus), ...
    'FontSize', 16, 'FontWeight', 'bold');
end

function plotTimeSeries(data, col_time, col_dff, col_dist, col_groom, options, title_str, time_limits, target_type, mouse_prefix, excel_filename)
% Plot time series with xlim for time limits

time_vec = data(:, col_time);
frame_rate = 1/mean(diff(time_vec));
time_plot = (1:length(time_vec)) / frame_rate;
dff_plot = data(:, col_dff);
dist_plot = data(:, col_dist);

% Apply filtering
if options.remove_grooming && size(data, 2) >= col_groom
    plot_mask = data(:, col_dist) > options.dist_limit & data(:, col_groom) == 0;
else
    plot_mask = data(:, col_dist) > options.dist_limit;
end

time_plot = time_plot(plot_mask);
dff_plot = dff_plot(plot_mask);
dist_plot = dist_plot(plot_mask);

if isempty(time_plot)
    text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Units', 'normalized', 'FontSize', 12);
    return;
end

% Plot data
yyaxis left;
plot(time_plot, dff_plot, 'LineWidth', 1.5, 'Color', [0, 0.4, 0.8]);
ylabel('\Delta F/F', 'FontSize', 12, 'Color', [0, 0.4, 0.8]);
ax = gca; ax.YAxis(1).Color = [0, 0.4, 0.8];

yyaxis right;
if strcmp(target_type, 'Water')
    plot(time_plot, dist_plot, 'LineWidth', 1.5, 'Color', [0.9, 0.3, 0.5]);
    ylabel('Distance to Water (cm)', 'FontSize', 12, 'Color', [0.9, 0.3, 0.5]);
    ax.YAxis(2).Color = [0.9, 0.3, 0.5];
else
    plot(time_plot, dist_plot, 'LineWidth', 1.5, 'Color', [0.8, 0.2, 0.4]);
    ylabel('Distance to Object (cm)', 'FontSize', 12, 'Color', [0.8, 0.2, 0.4]);
    ax.YAxis(2).Color = [0.8, 0.2, 0.4];
end

xlabel('Time (s)', 'FontSize', 12);

% Apply time limits using xlim
if ~isempty(time_limits) && length(time_limits) == 2
    xlim(time_limits);
    title(sprintf('%s (%.0f-%.0fs)', title_str, time_limits(1), time_limits(2)), ...
        'FontSize', 14, 'FontWeight', 'bold');
else
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
end

grid off; box off;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% ── Export time series data to Excel ───────────────────────────────────
if contains(title_str, 'Before')
    sheet_label = 'before';
elseif contains(title_str, 'Water')
    sheet_label = 'water';
elseif contains(title_str, 'Object')
    sheet_label = 'object';
else
    sheet_label = strrep(lower(title_str), ' ', '_');
end

% Apply time window to match what is visible in plot
if ~isempty(time_limits) && length(time_limits) == 2
    t_mask = time_plot >= time_limits(1) & time_plot <= time_limits(2);
else
    t_mask = true(size(time_plot));
end

if strcmp(target_type, 'Water')
    dist_label = 'Distance_to_water_cm';
else
    dist_label = 'Distance_to_object_cm';
end

export_table = table(...
    reshape(time_plot(t_mask), [], 1), ...
    reshape(dff_plot(t_mask),  [], 1), ...
    reshape(dist_plot(t_mask), [], 1), ...
    'VariableNames', {'Time_s', 'dFF', dist_label});

writetable(export_table, excel_filename, 'Sheet', sheet_label);
fprintf('Time series (%s) exported to: %s\n', sheet_label, excel_filename);
end

function bounds = calculateCommonBounds(data1, data2, col_x, col_y, col_dist, col_groom, options)
% Calculate common spatial bounds for two datasets

% Get coordinates from both datasets
x1 = data1(:, col_x); y1 = data1(:, col_y);
x2 = data2(:, col_x); y2 = data2(:, col_y);

% Apply filtering to both datasets
if options.remove_grooming && size(data1, 2) >= col_groom
    mask1 = data1(:, col_dist) > options.dist_limit & data1(:, col_groom) == 0;
    mask2 = data2(:, col_dist) > options.dist_limit & data2(:, col_groom) == 0;
else
    mask1 = data1(:, col_dist) > options.dist_limit;
    mask2 = data2(:, col_dist) > options.dist_limit;
end

% Get filtered coordinates
x1_filt = x1(mask1); y1_filt = y1(mask1);
x2_filt = x2(mask2); y2_filt = y2(mask2);

% Combine all coordinates to find common bounds
all_x = [x1_filt; x2_filt];
all_y = [y1_filt; y2_filt];

if isempty(all_x)
    bounds = struct('x_min', 0, 'x_max', 1, 'y_min', 0, 'y_max', 1);
else
    bounds = struct('x_min', min(all_x), 'x_max', max(all_x), ...
        'y_min', min(all_y), 'y_max', max(all_y));
end
end

function spatial_grid = createSpatialGridWithBounds(data, col_x, col_y, col_dff, col_dist, col_groom, options, bounds)
% Create spatial grid using specified bounds for coordinate alignment

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

% Use provided bounds instead of calculating from data
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