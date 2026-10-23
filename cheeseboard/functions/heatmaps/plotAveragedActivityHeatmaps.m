function plotAveragedActivityHeatmaps(mice_rotated, condition, options)
% PLOTAVERAGEDACTIVITYHEATMAPS - Plots averaged z-scored dF/F activity heatmaps across mice
%
% Creates heatmaps showing average neural activity patterns across all mice in a condition.
% Best used with rotated data from rotateMiceDataToAlignFood.
%
% Usage:
%   plotAveragedActivityHeatmaps(mice_rotated, 'saline')
%   plotAveragedActivityHeatmaps(mice_rotated, 'CNO', options)
%   
%   % With options:
%   options.time_limit = 10;            % Analyze only first 10 minutes
%   options.exclude_grooming = false;   % Exclude grooming periods
%   options.speed_threshold = 2;        % Minimum speed threshold
%   options.resolution = 100;           % Grid resolution
%   options.sigma = 2;                  % Gaussian smoothing sigma
%   options.colormap_name = 'jet';      % Colormap
%   options.clim_activity = [-2 2];     % Color limits for activity
%   options.clim_diff = [-1 1];         % Color limits for difference

%% Column definitions
COL_TIME = 1;
COL_X = 2;
COL_Y = 3;
COL_DFF = 6;
COL_SPEED = 7;
COL_GROOM = 10;

%% Input validation and defaults
if nargin < 2
    error('Both mice data and condition are required');
end

if nargin < 3
    options = struct();
end

% Set defaults
if ~isfield(options, 'exclude_grooming'), options.exclude_grooming = false; end
if ~isfield(options, 'resolution'), options.resolution = 100; end
if ~isfield(options, 'sigma'), options.sigma = 2; end
if ~isfield(options, 'colormap_name'), options.colormap_name = 'jet'; end
if ~isfield(options, 'figure_size'), options.figure_size = [100, 100, 1400, 450]; end

if isfield(options, 'time_limit') && ~isempty(options.time_limit)
    use_time_limit = true;
    time_limit_sec = options.time_limit * 60;
else
    use_time_limit = false;
    time_limit_sec = 0;
end

if isfield(options, 'speed_threshold') && ~isempty(options.speed_threshold)
    use_speed_threshold = true;
    speed_threshold = options.speed_threshold;
else
    use_speed_threshold = false;
    speed_threshold = 0;
end

%% Select mice based on condition
mouse_names = mice_rotated(:,1);
saline_idx = contains(mouse_names, 'saline');
cno_idx = contains(mouse_names, 'CNO');

has_memory_classification = size(mice_rotated, 2) >= 7 && ~all(cellfun(@isempty, mice_rotated(:,7)));
if has_memory_classification
    memory_classifications = mice_rotated(:,7);
    strong_memory_idx = strcmp(memory_classifications, 'strong memory');
    weak_memory_idx = strcmp(memory_classifications, 'weak memory');
end

switch lower(condition)
    case 'saline', selected_idx = saline_idx; condition_name = 'Saline';
    case 'cno', selected_idx = cno_idx; condition_name = 'CNO';
    case 'strong', selected_idx = strong_memory_idx; condition_name = 'Strong Memory';
    case 'weak', selected_idx = weak_memory_idx; condition_name = 'Weak Memory';
    otherwise, error('Invalid condition');
end

selected_mice = find(selected_idx);
n_mice = length(selected_mice);

fprintf('\n=== AVERAGED ACTIVITY HEATMAPS ===\n');
fprintf('Condition: %s (%d mice)\n', condition_name, n_mice);
fprintf('Resolution: %d x %d grid\n', options.resolution, options.resolution);

%% Determine common spatial bounds
all_coordinates = [];
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    pre_data = mice_rotated{mouse_idx, 3};
    test_data = mice_rotated{mouse_idx, 4};
    all_coordinates = [all_coordinates; pre_data(:, [COL_X, COL_Y]); test_data(:, [COL_X, COL_Y])];
end

x_min = min(all_coordinates(:,1)) - 10;
x_max = max(all_coordinates(:,1)) + 10;
y_min = min(all_coordinates(:,2)) - 10;
y_max = max(all_coordinates(:,2)) + 10;

bounds = struct('x_min', x_min, 'x_max', x_max, 'y_min', y_min, 'y_max', y_max);

fprintf('Coordinate ranges: X [%.1f, %.1f], Y [%.1f, %.1f]\n', x_min, x_max, y_min, y_max);

%% Accumulate activity maps across all mice
activity_pre_maps = [];
activity_test_maps = [];
food_positions = [];
center_positions = [];
n_mice_contributing = 0;

for i = 1:n_mice
    mouse_idx = selected_mice(i);
    
    pre_data = mice_rotated{mouse_idx, 3};
    test_data = mice_rotated{mouse_idx, 4};
    food_pos = mice_rotated{mouse_idx, 2};
    center_pos = mice_rotated{mouse_idx, 5};
    
    % Apply filters
    pre_filtered = apply_filters(pre_data, use_time_limit, time_limit_sec, ...
                                 use_speed_threshold, speed_threshold, options.exclude_grooming);
    test_filtered = apply_filters(test_data, use_time_limit, time_limit_sec, ...
                                  use_speed_threshold, speed_threshold, options.exclude_grooming);
    
    if size(pre_filtered,1) < 10 || size(test_filtered,1) < 10
        continue;
    end
    
    % Create activity maps for this mouse
    [activity_pre, activity_test] = createActivityMaps(pre_filtered, test_filtered, bounds, options);
    
    if isempty(activity_pre) || isempty(activity_test)
        continue;
    end
    
    % Store maps
    if isempty(activity_pre_maps)
        activity_pre_maps = zeros(size(activity_pre,1), size(activity_pre,2), n_mice);
        activity_test_maps = zeros(size(activity_test,1), size(activity_test,2), n_mice);
    end
    
    n_mice_contributing = n_mice_contributing + 1;
    activity_pre_maps(:,:,n_mice_contributing) = activity_pre;
    activity_test_maps(:,:,n_mice_contributing) = activity_test;
    
    food_positions = [food_positions; food_pos'];
    center_positions = [center_positions; center_pos'];
end

% Trim unused dimensions
activity_pre_maps = activity_pre_maps(:,:,1:n_mice_contributing);
activity_test_maps = activity_test_maps(:,:,1:n_mice_contributing);

% Average across mice (ignoring NaNs)
activity_pre_avg = nanmean(activity_pre_maps, 3);
activity_test_avg = nanmean(activity_test_maps, 3);

% Calculate difference
activity_diff = activity_test_avg - activity_pre_avg;

% Average food and center positions
avg_food = mean(food_positions, 1);
avg_center = mean(center_positions, 1);

fprintf('Contributing mice: %d\n', n_mice_contributing);
fprintf('Average food position: [%.1f, %.1f]\n', avg_food(1), avg_food(2));
fprintf('Average center position: [%.1f, %.1f]\n', avg_center(1), avg_center(2));

%% Create figure
figure('Name', sprintf('Averaged Activity: %s', condition_name), ...
       'Position', options.figure_size, 'Color', 'white');

% Create coordinate arrays
x_coords = linspace(bounds.x_min, bounds.x_max, options.resolution);
y_coords = linspace(bounds.y_min, bounds.y_max, options.resolution);

% Determine color limits
if isfield(options, 'clim_activity') && ~isempty(options.clim_activity)
    clim_act = options.clim_activity;
else
    all_vals = [activity_pre_avg(:); activity_test_avg(:)];
    all_vals_clean = all_vals(~isnan(all_vals));
    clim_act = [min(all_vals_clean), max(all_vals_clean)];
end

if isfield(options, 'clim_diff') && ~isempty(options.clim_diff)
    max_diff = max(abs(options.clim_diff));
else
    diff_clean = activity_diff(~isnan(activity_diff));
    max_diff = max(abs(diff_clean));
end

% Plot Before
subplot(1, 3, 1);
imagesc(x_coords, y_coords, activity_pre_avg, 'AlphaData', ~isnan(activity_pre_avg));
axis xy equal;
colormap(gca, options.colormap_name);
caxis(clim_act);
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Before Session', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Z-scored \DeltaF/F';
axis off;
hold off;

% Plot Test
subplot(1, 3, 2);
imagesc(x_coords, y_coords, activity_test_avg, 'AlphaData', ~isnan(activity_test_avg));
axis xy equal;
colormap(gca, options.colormap_name);
caxis(clim_act);
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Test Session', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Z-scored \DeltaF/F';
axis off;
hold off;

% Plot Difference
subplot(1, 3, 3);
imagesc(x_coords, y_coords, activity_diff, 'AlphaData', ~isnan(activity_diff));
axis xy equal;
colormap(gca, options.colormap_name);
caxis([-max_diff, max_diff]);
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Difference (Test - Before)', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = '\Delta(Z-scored \DeltaF/F)';
axis off;
hold off;

% Title
sgtitle(sprintf('Averaged Activity - %s (n=%d mice)', condition_name, n_mice_contributing), ...
        'FontSize', 16, 'FontWeight', 'bold');

% Print statistics
fprintf('\nActivity Statistics:\n');
fprintf('  Before - Mean: %.3f, Std: %.3f\n', nanmean(activity_pre_avg(:)), nanstd(activity_pre_avg(:)));
fprintf('  Test - Mean: %.3f, Std: %.3f\n', nanmean(activity_test_avg(:)), nanstd(activity_test_avg(:)));
fprintf('  Difference - Mean: %.3f, Max: %.3f, Min: %.3f\n', ...
        nanmean(activity_diff(:)), nanmax(activity_diff(:)), nanmin(activity_diff(:)));

end

%% Helper functions
function filtered_data = apply_filters(data, use_time_limit, time_limit_sec, ...
                                       use_speed_threshold, speed_threshold, exclude_grooming)
    COL_TIME = 1; COL_SPEED = 7; COL_GROOM = 10;
    filtered_data = data;
    
    if use_time_limit
        time = filtered_data(:, COL_TIME);
        time_mask = (time - min(time)) <= time_limit_sec;
        filtered_data = filtered_data(time_mask, :);
    end
    
    if use_speed_threshold
        speed_mask = filtered_data(:, COL_SPEED) >= speed_threshold;
        filtered_data = filtered_data(speed_mask, :);
    end
    
    if exclude_grooming
        non_grooming_mask = filtered_data(:, COL_GROOM) == 0;
        filtered_data = filtered_data(non_grooming_mask, :);
    end
end

function [activity_pre, activity_test] = createActivityMaps(pre_data, test_data, bounds, options)
    COL_X = 2; COL_Y = 3; COL_DFF = 6;
    
    % Extract data
    x_pre = pre_data(:, COL_X);
    y_pre = pre_data(:, COL_Y);
    dff_pre = pre_data(:, COL_DFF);
    
    x_test = test_data(:, COL_X);
    y_test = test_data(:, COL_Y);
    dff_test = test_data(:, COL_DFF);
    
    % Z-score separately
    if std(dff_pre) > 0
        dff_pre_zscore = (dff_pre - mean(dff_pre)) / std(dff_pre);
    else
        activity_pre = []; activity_test = []; return;
    end
    
    if std(dff_test) > 0
        dff_test_zscore = (dff_test - mean(dff_test)) / std(dff_test);
    else
        activity_pre = []; activity_test = []; return;
    end
    
    % Create spatial grids
    activity_pre = createSpatialGrid(x_pre, y_pre, dff_pre_zscore, bounds, options.resolution, options.sigma);
    activity_test = createSpatialGrid(x_test, y_test, dff_test_zscore, bounds, options.resolution, options.sigma);
end

function spatial_grid = createSpatialGrid(x_coords, y_coords, dff_values, bounds, resolution, sigma)
    dff_grid = zeros(resolution, resolution);
    visit_count = zeros(resolution, resolution);
    
    if isempty(x_coords) || bounds.x_max <= bounds.x_min || bounds.y_max <= bounds.y_min
        spatial_grid = NaN(resolution, resolution);
        return;
    end
    
    % Map to grid - x maps to columns, y maps to rows
    x_idx = round((x_coords - bounds.x_min) / (bounds.x_max - bounds.x_min) * (resolution-1)) + 1;
    y_idx = round((y_coords - bounds.y_min) / (bounds.y_max - bounds.y_min) * (resolution-1)) + 1;
    x_idx = max(1, min(resolution, x_idx));
    y_idx = max(1, min(resolution, y_idx));
    
    % Accumulate values using standard matrix indexing (row, col) = (y, x)
    for i = 1:length(x_coords)
        dff_grid(y_idx(i), x_idx(i)) = dff_grid(y_idx(i), x_idx(i)) + dff_values(i);
        visit_count(y_idx(i), x_idx(i)) = visit_count(y_idx(i), x_idx(i)) + 1;
    end
    
    % Calculate mean
    visited_mask = visit_count > 0;
    mean_dff_grid = zeros(size(dff_grid));
    mean_dff_grid(visited_mask) = dff_grid(visited_mask) ./ visit_count(visited_mask);
    masked_dff = mean_dff_grid;
    masked_dff(~visited_mask) = NaN;
    
    % Apply Gaussian smoothing
    kernel_size = 11;
    [X, Y] = meshgrid(-floor(kernel_size/2):floor(kernel_size/2), -floor(kernel_size/2):floor(kernel_size/2));
    kernel = exp(-(X.^2 + Y.^2) / (2*sigma^2));
    kernel = kernel / sum(kernel(:));
    
    smoothed_dff = NaN(size(masked_dff));
    valid_mask = ~isnan(masked_dff);
    padded_dff = padarray(masked_dff, [floor(kernel_size/2), floor(kernel_size/2)], NaN);
    padded_mask = padarray(valid_mask, [floor(kernel_size/2), floor(kernel_size/2)], 0);
    
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
                        smoothed_dff(i-floor(kernel_size/2), j-floor(kernel_size/2)) = sum(weighted_values(:)) / weight_sum;
                    end
                end
            end
        end
    end
    
    spatial_grid = smoothed_dff;
end