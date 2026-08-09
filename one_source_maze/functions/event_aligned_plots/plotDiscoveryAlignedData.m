function plotDiscoveryAlignedData(mice_all, options)
% PLOTDISCOVERYALIGNEDDATA Creates a discovery-aligned visualization with raster and average plots
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - time_window: two-element vector specifying time range around discovery [before after] in seconds (default: [-100 300])
%     - caxis_range: two-element vector specifying color axis range for DFF (default: [-2 2])
%     - sort_by: how to sort mice in raster ('discovery_time', 'condition', 'response_magnitude', default: 'condition')
%     - z_score: whether to z-score the DFF signals (default: true)
%     - smooth_window: window size for smoothing in frames (default: 30)
%     - plot_distance: plot distance to food instead of dFF (default: false)
%     - distance_ylim: y-axis limits for distance plot [min max] (default: [0 100])

% Constants
COL_TIME = 1;     % Time column
COL_DFF = 11;     % DFF data column
COL_DIST = 5;     % Distance to food

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'state')
    options.state = 'all';
end

if ~isfield(options, 'source')
    options.source = 'all';
end

if ~isfield(options, 'time_window')
    options.time_window = [-100, 300]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'sort_by')
    options.sort_by = 'condition';
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 30;
end

if ~isfield(options, 'plot_distance')
    options.plot_distance = false;
end

if ~isfield(options, 'distance_ylim')
    options.distance_ylim = [0, 150];
end

% Filter mice based on state and source
n = size(mice_all, 1);
mask = true(n, 1);

if ~strcmpi(options.state, 'all')
    mask = mask & strcmp(mice_all(:, 2), options.state);
end

if ~strcmpi(options.source, 'all')
    mask = mask & strcmp(mice_all(:, 3), options.source);
end

% Extract filtered mice
filtered_mice = mice_all(mask, :);
num_mice = size(filtered_mice, 1);

if num_mice == 0
    error('No mice found with the specified criteria');
end

% Initialize matrices to store aligned data and time
all_aligned_dff = cell(num_mice, 1);
all_aligned_time = cell(num_mice, 1);
all_aligned_dist = cell(num_mice, 1);

% Extract and align data for each mouse based on real time
for i = 1:num_mice
    data = filtered_mice{i, 4};
    discovery = filtered_mice{i, 6};
    
    % Get the discovery time in seconds
    discovery_time_sec = data(discovery, COL_TIME);
    
    % Extract all data and convert time to seconds from discovery
    time_sec = data(:, COL_TIME) - discovery_time_sec;
    dff = data(:, COL_DFF);
    dist = data(:, COL_DIST);
    
    % Filter data to be within the specified time window
    time_mask = time_sec >= options.time_window(1) & time_sec <= options.time_window(2);
    
    % Apply filter
    filtered_time = time_sec(time_mask);
    filtered_dff = dff(time_mask);
    filtered_dist = dist(time_mask);
    
    % Apply z-scoring if requested (only for DFF)
    if options.z_score && ~options.plot_distance
        filtered_dff = (filtered_dff - nanmean(filtered_dff)) / nanstd(filtered_dff);
    end
    
    % Apply smoothing if requested (for either DFF or distance)
    if options.smooth_window > 0
        if options.plot_distance
            filtered_dist = movmean(filtered_dist, options.smooth_window, 'omitnan');
        else
            filtered_dff = movmean(filtered_dff, options.smooth_window, 'omitnan');
        end
    end
    
    % Store aligned data
    all_aligned_dff{i} = filtered_dff;
    all_aligned_time{i} = filtered_time;
    all_aligned_dist{i} = filtered_dist;
end

% Determine the common time grid for plotting
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000); % 1000 time points for smooth visualization

% Interpolate data onto common time grid for each mouse
aligned_dff = nan(num_mice, length(time_grid));
aligned_dist = nan(num_mice, length(time_grid));

for i = 1:num_mice
    % Skip if no valid data
    if isempty(all_aligned_time{i}) || all(isnan(all_aligned_dff{i}))
        continue;
    end
    
    % Interpolate DFF onto common time grid
    aligned_dff(i, :) = interp1(all_aligned_time{i}, all_aligned_dff{i}, time_grid, 'linear', NaN);
    
    % Interpolate distance onto common time grid
    aligned_dist(i, :) = interp1(all_aligned_time{i}, all_aligned_dist{i}, time_grid, 'linear', NaN);
end

% Sort the mice based on sort_by option
switch lower(options.sort_by)
    case 'discovery_time'
        [~, sort_idx] = sort([filtered_mice{:, 6}]);
    
    case 'condition'
        % Sort by state first, then by source
        [~, sort_idx] = sortrows([filtered_mice(:, 2), filtered_mice(:, 3)]);
    
    case 'response_magnitude'
        % Find discovery index in time_grid
        [~, discovery_idx] = min(abs(time_grid));
        % Look at response 0-50 seconds after discovery
        response_window = discovery_idx:(discovery_idx + round(50 / (time_max-time_min) * length(time_grid)));
        response_window = response_window(response_window <= length(time_grid));
        
        % Sort based on the selected data type
        if options.plot_distance
            % Sort by maximum distance decrease
            [~, sort_idx] = sort(min(aligned_dist(:, response_window), [], 2), 'ascend');
        else
            % Sort by peak response magnitude for DFF
            [~, sort_idx] = sort(max(aligned_dff(:, response_window), [], 2), 'descend');
        end
        
    otherwise
        sort_idx = 1:num_mice;
end

% Apply sorting
aligned_dff = aligned_dff(sort_idx, :);
aligned_dist = aligned_dist(sort_idx, :);
filtered_mice = filtered_mice(sort_idx, :);

% Create figure with two subplots
figure('Position', [100, 100, 1000, 600]);

% 1. Average plot (left subplot)
subplot(1, 2, 1);
hold on;

% Select data type to plot based on options
if options.plot_distance
    plot_data = aligned_dist;
    data_label = 'Distance to food (cm)';
    plot_color = [0.2, 0.6, 0.2]; % Green for distance
else
    plot_data = aligned_dff;
    if options.z_score
        data_label = 'ΔF/F (z-scored)';
    else
        data_label = 'ΔF/F';
    end
    plot_color = [0.8, 0.1, 0.1]; % Red for neural activity
end

% Calculate average and SEM
avg_data = nanmean(plot_data, 1);
sem_data = nanstd(plot_data, 0, 1) ./ sqrt(sum(~isnan(plot_data), 1));

% Add discovery line
xline(0, '--k', 'LineWidth', 1.5);

% Plot SEM shading
ciplot(avg_data - sem_data, avg_data + sem_data, time_grid, plot_color, 0.3);

% Plot the mean line on top
plot(time_grid, avg_data, 'Color', plot_color, 'LineWidth', 2);

% Format plot
xlabel('Time from discovery (seconds)');
ylabel(data_label, 'FontSize', 14);

% Create title based on data type
if options.plot_distance
    title_prefix = 'Average distance to food';
else
    title_prefix = 'Average neural activity';
end

title({[title_prefix ' aligned to food discovery'], ...
       ['State: ' upper(options.state(1)) options.state(2:end), ...
        ', Source: ' upper(options.source(1)) options.source(2:end)]});
grid off;
box off;

% Set axis limits
xlim(options.time_window);
if options.plot_distance
    ylim(options.distance_ylim);
end

% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% 2. Raster plot (right subplot)
subplot(1, 2, 2);

% Plot data based on selected type
if options.plot_distance
    imagesc(time_grid, 1:num_mice, plot_data);
    % Use a colormap appropriate for distance (green-to-white)
    colormap(gca, bluewhitered(256)); % Flipped hot colormap works well for distance
    
    % Set color limits for distance
    if ~isfield(options, 'dist_caxis_range')
        caxis(options.distance_ylim);
    else
        caxis(options.dist_caxis_range);
    end
else
    % Standard DFF plot with blue-white-red colormap
    imagesc(time_grid, 1:num_mice, plot_data);
    colormap(gca, bluewhitered(256));
    caxis(options.caxis_range);
end

% Add discovery line
hold on;
xline(0, '--k', 'LineWidth', 1.5);

% Format plot
xlabel('Time from discovery (seconds)');
ylabel('Mouse #', 'FontSize', 14);

if options.plot_distance
    title_text = 'Distance to food across mice';
else
    title_text = 'Neural activity across mice';
end
title({title_text, ['n = ' num2str(num_mice) ' mice']});

% Add colorbar with appropriate label
h = colorbar;
ylabel(h, data_label, 'FontSize', 14);

% Add mouse IDs and conditions as y-tick labels if not too many mice
if num_mice <= 30
    mouse_labels = cell(num_mice, 1);
    for i = 1:num_mice
        if ischar(filtered_mice{i, 1})
            id_parts = strsplit(filtered_mice{i, 1}, '_');
            if length(id_parts) > 1
                mouse_id = id_parts{1};
            else
                mouse_id = filtered_mice{i, 1};
            end
        else
            mouse_id = num2str(i);
        end
        mouse_labels{i} = mouse_id;
    end
    set(gca, 'YTick', 1:num_mice, 'YTickLabel', mouse_labels);
else
    % Just show mouse numbers if there are too many
    set(gca, 'YTick', round(linspace(1, num_mice, 10)));
end

% Set x-axis limits
xlim(options.time_window);

% Add overall title
if options.plot_distance
    main_title = 'Distance to Food';
    subtitle = '';
else
    main_title = 'Neural Activity';
    if options.z_score
        subtitle = 'Z-scored ΔF/F values';
    else
        subtitle = 'Raw ΔF/F values';
    end
end

sgtitle({[main_title ' Aligned to Food Discovery: ', ...
         options.state, ' mice, ', options.source, ' stimulus'], ...
         subtitle}, 'FontWeight', 'bold');

%% Export data to Excel
% Determine data label for filename
if options.plot_distance
    data_type_label = 'distance';
else
    if options.z_score
        data_type_label = 'dFF_zscored';
    else
        data_type_label = 'dFF';
    end
end

excel_filename = sprintf('discovery_aligned_%s_%s_%s_%s.xlsx', ...
    data_type_label, options.state, options.source, datestr(now, 'yyyymmdd'));

% ── Sheet 1: Average trace (mean and SEM) ────────────────────────────────
avg_table = table(time_grid(:), avg_data(:), sem_data(:), ...
    'VariableNames', {'Time_s', 'Mean', 'SEM'});

writetable(avg_table, excel_filename, 'Sheet', 'average_trace');

% ── Sheet 2: Heatmap - one row per mouse ─────────────────────────────────
% Build mouse ID list (already sorted)
mouse_ids_sorted = cell(num_mice, 1);
for i = 1:num_mice
    mouse_ids_sorted{i} = filtered_mice{i, 1};
end

% Heatmap data: rows = mice, columns = time points
% Add mouse IDs as first column using a separate approach since
% writetable needs uniform types - write time axis first then data
time_header = table(time_grid(:)', 'VariableNames', {'Mouse_ID'});

% Build heatmap table: each row is a mouse, each column is a time point
% Column names: Mouse_ID, then t_1, t_2, ...
col_names = [{'Mouse_ID'}, arrayfun(@(t) sprintf('t_%.2f', t), time_grid, 'UniformOutput', false)];

% Write time axis as first row, then per-mouse data
% Use writecell for flexibility with mixed types
heatmap_cell = cell(num_mice + 1, length(time_grid) + 1);

% First row: header with time points
heatmap_cell{1, 1} = 'Mouse_ID';
for t = 1:length(time_grid)
    heatmap_cell{1, t+1} = time_grid(t);
end

% Fill mouse rows
for i = 1:num_mice
    heatmap_cell{i+1, 1} = mouse_ids_sorted{i};
    for t = 1:length(time_grid)
        if options.plot_distance
            heatmap_cell{i+1, t+1} = aligned_dist(i, t);
        else
            heatmap_cell{i+1, t+1} = aligned_dff(i, t);
        end
    end
end

writecell(heatmap_cell, excel_filename, 'Sheet', 'heatmap_per_mouse');

fprintf('Discovery-aligned data exported to: %s\n', excel_filename);

% Output some statistics
fprintf('Analysis complete: %d mice included\n', num_mice);
if options.plot_distance
    fprintf('Average pre-discovery distance: %.1f\n', nanmean(avg_data(time_grid < 0)));
    fprintf('Average post-discovery distance: %.1f\n', nanmean(avg_data(time_grid > 0)));
else
    fprintf('Average pre-discovery DFF: %.3f\n', nanmean(avg_data(time_grid < 0)));
    fprintf('Average post-discovery DFF: %.3f\n', nanmean(avg_data(time_grid > 0)));
end

% Optional: return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.dff = aligned_dff;
    aligned_data.time = time_grid;
    aligned_data.dist = aligned_dist;
    aligned_data.mouse_info = filtered_mice;
    varargout{1} = aligned_data;
end
end

function ciplot(lower, upper, x, color, alpha)
% ciplot creates a shaded area to visualize confidence intervals
% lower: lower bound
% upper: upper bound
% x: x values
% color: RGB color of the shaded area
% alpha: transparency (0-1)

% Create polygon vertices
x_polygon = [x, fliplr(x)];
y_polygon = [upper, fliplr(lower)];

% Remove any NaN values which can break the fill
nan_indices = isnan(x_polygon) | isnan(y_polygon);
x_polygon(nan_indices) = [];
y_polygon(nan_indices) = [];

% Create the polygon
h = fill(x_polygon, y_polygon, color);
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', alpha);
end

function cmap = bluewhitered(m)
% BLUEWHITERED creates a blue-to-white-to-red colormap
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a colormap
%   that transitions from dark blue through a small white region to true red

if nargin < 1
   m = 100;
end

% Define the color transitions: dark blue -> medium blue -> light blue -> white -> light red -> medium red -> true red
% These are the key colors in our gradient
colors = [
    0.0, 0.0, 0.5;  % Dark blue
    0.1, 0.3, 0.7;  % Softer blue
    0.4, 0.6, 0.9;  % Medium blue
    0.7, 0.85, 1.0; % Very light blue
    1.0, 1.0, 1.0;  % White (center point - reduced region)
    1.0, 0.85, 0.85; % Very light red
    0.9, 0.4, 0.4;  % Medium red
    0.8, 0.1, 0.1   % True red
];

n = size(colors, 1);
cmap = zeros(m, 3);

% Create non-uniform spacing to make white region smaller
% White will be at exactly the center but occupy fewer color indices
positions = [0, 0.15, 0.3, 0.45, 0.5, 0.55, 0.7, 1.0];
idx = round(positions * (m-1)) + 1;

% Interpolate colors between key points
for i = 1:n-1
    % Start and end indices for this segment
    i1 = idx(i);
    i2 = idx(i+1);
    
    % Number of colors in this segment
    ni = i2 - i1 + 1;
    
    % Linear interpolation for each RGB component
    r = linspace(colors(i,1), colors(i+1,1), ni);
    g = linspace(colors(i,2), colors(i+1,2), ni);
    b = linspace(colors(i,3), colors(i+1,3), ni);
    
    % Store in the colormap
    cmap(i1:i2, :) = [r', g', b'];
end
end