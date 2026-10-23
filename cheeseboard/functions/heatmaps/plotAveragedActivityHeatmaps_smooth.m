function plotAveragedActivityHeatmaps_smooth(mice_rotated, condition, options)
% PLOTAVERAGEDACTIVITYHEATMAPS_SMOOTH - Plots ultra-smooth averaged z-scored dF/F activity heatmaps
%
% Creates smooth heatmaps showing average neural activity patterns across all mice in a condition.
% Uses histcounts2 + imgaussfilt + optional interpolation for smooth, non-pixelated results.
% Best used with rotated data from rotateMiceDataToAlignFood.
%
% Usage:
%   plotAveragedActivityHeatmaps_smooth(mice_rotated, 'saline')
%   plotAveragedActivityHeatmaps_smooth(mice_rotated, 'CNO', options)
%   
%   % With options:
%   options.time_limit = 10;            % Analyze only first 10 minutes
%   options.exclude_grooming = false;   % Exclude grooming periods
%   options.speed_threshold = 2;        % Minimum speed threshold
%   options.bin_size = 4;               % Spatial bin size in cm (default: 4)
%   options.smooth_factor = 2.0;        % Gaussian smoothing factor (default: 2.0)
%   options.interpolate_display = true; % Interpolate to finer grid (default: true)
%   options.interp_factor = 4;          % Interpolation factor (default: 4)
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
if ~isfield(options, 'bin_size'), options.bin_size = 4; end
if ~isfield(options, 'smooth_factor'), options.smooth_factor = 2.0; end
if ~isfield(options, 'interpolate_display'), options.interpolate_display = true; end
if ~isfield(options, 'interp_factor'), options.interp_factor = 4; end
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

fprintf('\n=== AVERAGED ACTIVITY HEATMAPS (SMOOTH) ===\n');
fprintf('Condition: %s (%d mice)\n', condition_name, n_mice);
fprintf('Bin size: %.1f cm\n', options.bin_size);
fprintf('Smoothing factor: %.1f\n', options.smooth_factor);
if options.interpolate_display
    fprintf('Display interpolation: %dx finer grid\n', options.interp_factor);
end

%% Determine common spatial bounds
all_coordinates = [];
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    pre_data = mice_rotated{mouse_idx, 3};
    test_data = mice_rotated{mouse_idx, 4};
    all_coordinates = [all_coordinates; pre_data(:, [COL_X, COL_Y]); test_data(:, [COL_X, COL_Y])];
end

x_min = floor(min(all_coordinates(:,1)) / options.bin_size) * options.bin_size;
x_max = ceil(max(all_coordinates(:,1)) / options.bin_size) * options.bin_size;
y_min = floor(min(all_coordinates(:,2)) / options.bin_size) * options.bin_size;
y_max = ceil(max(all_coordinates(:,2)) / options.bin_size) * options.bin_size;

% Create spatial bins
x_edges = x_min:options.bin_size:x_max;
y_edges = y_min:options.bin_size:y_max;
x_centers = x_edges(1:end-1) + options.bin_size/2;
y_centers = y_edges(1:end-1) + options.bin_size/2;

fprintf('Spatial grid: %d x %d bins\n', length(x_centers), length(y_centers));
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
    
    % Create activity maps for this mouse using smooth method
    [activity_pre, activity_test] = createSmoothActivityMaps(pre_filtered, test_filtered, ...
                                                             x_edges, y_edges, x_centers, y_centers, options);
    
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
figure('Name', sprintf('Averaged Activity (Smooth): %s', condition_name), ...
       'Position', options.figure_size, 'Color', 'white');

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
    if ~isempty(diff_clean)
        max_diff = max(abs(diff_clean));
    else
        max_diff = 1;
    end
end

% Optional: Interpolate to finer grid for ultra-smooth display
if options.interpolate_display
    [X, Y] = meshgrid(x_centers, y_centers);
    x_fine = linspace(x_centers(1), x_centers(end), length(x_centers) * options.interp_factor);
    y_fine = linspace(y_centers(1), y_centers(end), length(y_centers) * options.interp_factor);
    [X_fine, Y_fine] = meshgrid(x_fine, y_fine);
    
    % Interpolate all maps
    activity_pre_fine = interp2(X, Y, activity_pre_avg', X_fine, Y_fine, 'cubic');
    activity_test_fine = interp2(X, Y, activity_test_avg', X_fine, Y_fine, 'cubic');
    activity_diff_fine = interp2(X, Y, activity_diff', X_fine, Y_fine, 'cubic');
    
    x_plot = x_fine;
    y_plot = y_fine;
else
    activity_pre_fine = activity_pre_avg';
    activity_test_fine = activity_test_avg';
    activity_diff_fine = activity_diff';
    
    x_plot = x_centers;
    y_plot = y_centers;
end

% Plot Before
subplot(1, 3, 1);
imagesc(x_plot, y_plot, activity_pre_fine, 'AlphaData', ~isnan(activity_pre_fine));
axis xy equal tight;
colormap(gca, options.colormap_name);
caxis(clim_act);
set(gca, 'Color', 'white');
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
imagesc(x_plot, y_plot, activity_test_fine, 'AlphaData', ~isnan(activity_test_fine));
axis xy equal tight;
colormap(gca, options.colormap_name);
caxis(clim_act);
set(gca, 'Color', 'white');
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
imagesc(x_plot, y_plot, activity_diff_fine, 'AlphaData', ~isnan(activity_diff_fine));
axis xy equal tight;
colormap(gca, options.colormap_name);
caxis([-max_diff, max_diff]);
set(gca, 'Color', 'white');
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
sgtitle(sprintf('Smooth Averaged Activity - %s (n=%d mice)', condition_name, n_mice_contributing), ...
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

function [activity_pre, activity_test] = createSmoothActivityMaps(pre_data, test_data, ...
                                                                  x_edges, y_edges, x_centers, y_centers, options)
    % Create smooth activity maps using histcounts2 + imgaussfilt
    COL_X = 2; COL_Y = 3; COL_DFF = 6;
    
    % Extract data
    x_pre = pre_data(:, COL_X);
    y_pre = pre_data(:, COL_Y);
    dff_pre = pre_data(:, COL_DFF);
    
    x_test = test_data(:, COL_X);
    y_test = test_data(:, COL_Y);
    dff_test = test_data(:, COL_DFF);
    
    % Z-score separately for each period
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
    
    % Create smooth spatial grids
    activity_pre = createSmoothHeatmap(x_pre, y_pre, dff_pre_zscore, ...
                                      x_edges, y_edges, x_centers, y_centers, options);
    activity_test = createSmoothHeatmap(x_test, y_test, dff_test_zscore, ...
                                       x_edges, y_edges, x_centers, y_centers, options);
end

function heatmap = createSmoothHeatmap(x_data, y_data, dff_data, x_edges, y_edges, x_centers, y_centers, options)
    % Create smooth heatmap using histcounts2 + imgaussfilt
    
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