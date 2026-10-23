function plotAveragedOccupancyHeatmaps(mice_rotated, condition, options)
% PLOTAVERAGEDOCCUPANCYHEATMAPS - Plots averaged occupancy heatmaps across mice
%
% Creates heatmaps showing average occupancy patterns across all mice in a condition.
% Best used with rotated data from rotateMiceDataToAlignFood.
% Uses masks to show only areas with data (transparent background elsewhere).
%
% Usage:
%   plotAveragedOccupancyHeatmaps(mice_rotated, 'saline')
%   plotAveragedOccupancyHeatmaps(mice_rotated, 'CNO', options)
%   
%   % With options:
%   options.time_limit = 10;            % Analyze only first 10 minutes
%   options.exclude_grooming = false;   % Exclude grooming periods
%   options.speed_threshold = 2;        % Minimum speed threshold
%   options.bin_size = 5;               % Spatial bin size in cm
%   options.smooth_factor = 1.5;        % Gaussian smoothing
%   options.colormap_name = 'hot';      % Colormap
%   options.clim_occupancy = [0 10];    % Color limits for occupancy
%   options.clim_diff = [-5 5];         % Color limits for difference

%% Column definitions
COL_TIME = 1;
COL_X = 2;
COL_Y = 3;
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
if ~isfield(options, 'bin_size'), options.bin_size = 5; end
if ~isfield(options, 'smooth_factor'), options.smooth_factor = 1.5; end
if ~isfield(options, 'colormap_name'), options.colormap_name = 'hot'; end
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

fprintf('\n=== AVERAGED OCCUPANCY HEATMAPS ===\n');
fprintf('Condition: %s (%d mice)\n', condition_name, n_mice);
fprintf('Bin size: %.1f cm\n', options.bin_size);

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

x_edges = x_min:options.bin_size:x_max;
y_edges = y_min:options.bin_size:y_max;
x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

fprintf('Spatial grid: %d x %d bins\n', length(x_edges)-1, length(y_edges)-1);

%% Accumulate occupancy across all mice
occupancy_pre_sum = zeros(length(x_edges)-1, length(y_edges)-1);
occupancy_test_sum = zeros(length(x_edges)-1, length(y_edges)-1);
mask_pre_sum = zeros(length(x_edges)-1, length(y_edges)-1);  % Track which bins have data
mask_test_sum = zeros(length(x_edges)-1, length(y_edges)-1);
n_mice_contributing = 0;

% Get average food and center positions
food_positions = [];
center_positions = [];

for i = 1:n_mice
    mouse_idx = selected_mice(i);
    
    pre_data = mice_rotated{mouse_idx, 3};
    test_data = mice_rotated{mouse_idx, 4};
    food_pos = mice_rotated{mouse_idx, 2};
    center_pos = mice_rotated{mouse_idx, 5};
    
    food_positions = [food_positions; food_pos'];
    center_positions = [center_positions; center_pos'];
    
    % Apply filters
    pre_filtered = apply_filters(pre_data, use_time_limit, time_limit_sec, ...
                                 use_speed_threshold, speed_threshold, options.exclude_grooming);
    test_filtered = apply_filters(test_data, use_time_limit, time_limit_sec, ...
                                  use_speed_threshold, speed_threshold, options.exclude_grooming);
    
    if size(pre_filtered,1) < 10 || size(test_filtered,1) < 10
        continue;
    end
    
    x_pre = pre_filtered(:, COL_X);
    y_pre = pre_filtered(:, COL_Y);
    x_test = test_filtered(:, COL_X);
    y_test = test_filtered(:, COL_Y);
    
    % Create histograms
    occ_pre = histcounts2(x_pre, y_pre, x_edges, y_edges);
    occ_test = histcounts2(x_test, y_test, x_edges, y_edges);
    
    % Create masks (bins with any data)
    mask_pre = occ_pre > 0;
    mask_test = occ_test > 0;
    
    % Normalize to percentage
    occ_pre_pct = (occ_pre / sum(occ_pre(:))) * 100;
    occ_test_pct = (occ_test / sum(occ_test(:))) * 100;
    
    occupancy_pre_sum = occupancy_pre_sum + occ_pre_pct;
    occupancy_test_sum = occupancy_test_sum + occ_test_pct;
    mask_pre_sum = mask_pre_sum + double(mask_pre);
    mask_test_sum = mask_test_sum + double(mask_test);
    n_mice_contributing = n_mice_contributing + 1;
end

% Average across mice
occupancy_pre_avg = occupancy_pre_sum / n_mice_contributing;
occupancy_test_avg = occupancy_test_sum / n_mice_contributing;

% Create combined masks (bins visited by any mouse)
mask_pre_combined = mask_pre_sum > 0;
mask_test_combined = mask_test_sum > 0;
mask_diff_combined = mask_pre_combined | mask_test_combined;  % Show where either session had data

% Set NaN for unvisited areas BEFORE smoothing
occupancy_pre_avg(~mask_pre_combined) = NaN;
occupancy_test_avg(~mask_test_combined) = NaN;

% Apply smoothing (NaN-aware)
if options.smooth_factor > 0
    % Smooth pre
    temp_pre = occupancy_pre_avg;
    temp_pre(isnan(temp_pre)) = 0;
    smoothed_pre = imgaussfilt(temp_pre, options.smooth_factor);
    mask_smooth_pre = imgaussfilt(double(mask_pre_combined), options.smooth_factor);
    smoothed_pre(mask_smooth_pre < 0.1) = NaN;
    occupancy_pre_avg = smoothed_pre;
    
    % Smooth test
    temp_test = occupancy_test_avg;
    temp_test(isnan(temp_test)) = 0;
    smoothed_test = imgaussfilt(temp_test, options.smooth_factor);
    mask_smooth_test = imgaussfilt(double(mask_test_combined), options.smooth_factor);
    smoothed_test(mask_smooth_test < 0.1) = NaN;
    occupancy_test_avg = smoothed_test;
end

% Calculate difference (only where both have data)
occupancy_diff = occupancy_test_avg - occupancy_pre_avg;

% Average food and center positions
avg_food = mean(food_positions, 1);
avg_center = mean(center_positions, 1);

fprintf('Contributing mice: %d\n', n_mice_contributing);
fprintf('Average food position: [%.1f, %.1f]\n', avg_food(1), avg_food(2));
fprintf('Average center position: [%.1f, %.1f]\n', avg_center(1), avg_center(2));

%% Create figure
figure('Name', sprintf('Averaged Occupancy: %s', condition_name), ...
       'Position', options.figure_size, 'Color', 'white');

% Determine color limits
if isfield(options, 'clim_occupancy') && ~isempty(options.clim_occupancy)
    clim_occ = options.clim_occupancy;
else
    all_occ_vals = [occupancy_pre_avg(:); occupancy_test_avg(:)];
    all_occ_vals = all_occ_vals(~isnan(all_occ_vals));
    clim_occ = [0, max(all_occ_vals)];
end

if isfield(options, 'clim_diff') && ~isempty(options.clim_diff)
    max_diff = max(abs(options.clim_diff));
else
    diff_vals = occupancy_diff(~isnan(occupancy_diff));
    if ~isempty(diff_vals)
        max_diff = max(abs(diff_vals));
    else
        max_diff = 1;
    end
end

% Plot Before
subplot(1, 3, 1);
imagesc(x_centers, y_centers, occupancy_pre_avg', 'AlphaData', ~isnan(occupancy_pre_avg'));
axis xy equal;
colormap(gca, options.colormap_name);
caxis(clim_occ);
set(gca, 'Color', 'white');  % Background color for areas with no data
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Before Session', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Occupancy (%)';
axis off;
hold off;

% Plot Test
subplot(1, 3, 2);
imagesc(x_centers, y_centers, occupancy_test_avg', 'AlphaData', ~isnan(occupancy_test_avg'));
axis xy equal;
colormap(gca, options.colormap_name);
caxis(clim_occ);
set(gca, 'Color', 'white');  % Background color for areas with no data
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Test Session', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Occupancy (%)';
axis off;
hold off;

% Plot Difference
subplot(1, 3, 3);
imagesc(x_centers, y_centers, occupancy_diff', 'AlphaData', ~isnan(occupancy_diff'));
axis xy equal;
colormap(gca, options.colormap_name);
caxis([-max_diff, max_diff]);
set(gca, 'Color', 'white');  % Background color for areas with no data
hold on;
plot(avg_food(1), avg_food(2), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'y', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 1);
% plot(avg_center(1), avg_center(2), 'x', 'MarkerSize', 12, 'Color', 'w', 'LineWidth', 3);
title('Difference (Test - Before)', 'FontSize', 14, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Occupancy Change (%)';
axis off;
hold off;

% Title
sgtitle(sprintf('Averaged Occupancy - %s (n=%d mice)', condition_name, n_mice_contributing), ...
        'FontSize', 16, 'FontWeight', 'bold');

end

%% Helper function
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