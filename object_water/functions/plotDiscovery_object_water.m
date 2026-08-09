function plotDiscovery_object_water(mice_all, options)
% PLOTDISCOVERY_OBJECT_WATER Creates discovery-aligned visualization for water and object
% discovery for all mice in the dataset, excluding grooming and climbing data
%
% Parameters:
%   mice_all - Cell array with mouse data
%   options - Struct with visualization parameters:
%     - time_window: two-element vector specifying time range around discovery [before after] in seconds (default: [-100 300])
%     - caxis_range: two-element vector specifying color axis range for DFF (default: [-2 2])
%     - z_score: whether to z-score the DFF signals (default: true)
%     - smooth_window: window size for smoothing in frames (default: 30)

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_dist_water = 9;  % Distance to water
col_dist_object = 10; % Distance to object
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'time_window')
    options.time_window = [-100, 300]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 30;
end

% Extract all unique mouse IDs from the dataset
mouse_ids = {};
for i = 1:size(mice_all, 1)
    current_id = mice_all{i, 1};
    if ischar(current_id)
        % Extract the mouse ID (e.g., 'F29' from 'F29_water_object_data_water')
        parts = strsplit(current_id, '_');
        mouse_id = parts{1};
        
        % Add to list if not already there
        if ~ismember(mouse_id, mouse_ids)
            mouse_ids{end+1} = mouse_id;
        end
    end
end

fprintf('Found %d unique mice in the dataset\n', length(mouse_ids));

% Process water sessions for all mice together
all_water_sessions = {};
all_water_mice_ids = {};
all_water_idx = [];

% Process object sessions for all mice together
all_object_sessions = {};
all_object_mice_ids = {};
all_object_idx = [];

% Collect sessions by type for all mice
for i = 1:size(mice_all, 1)
    session_id = mice_all{i, 1};
    if ischar(session_id)
        parts = strsplit(session_id, '_');
        mouse_id = parts{1}; % The mouse ID part
        
        if contains(session_id, 'water')
            all_water_sessions{end+1} = mice_all(i, :);
            all_water_mice_ids{end+1} = mouse_id;
            all_water_idx(end+1) = i;
        elseif contains(session_id, 'object')
            all_object_sessions{end+1} = mice_all(i, :);
            all_object_mice_ids{end+1} = mouse_id;
            all_object_idx(end+1) = i;
        end
        % Ignore 'before' sessions as they have no discovery events
    end
end

% Create water discovery plot for all mice
if ~isempty(all_water_sessions)
    createAllMiceDiscoveryPlot(mice_all(all_water_idx, :), all_water_mice_ids, 'water', options);
end

% Create object discovery plot for all mice
if ~isempty(all_object_sessions)
    createAllMiceDiscoveryPlot(mice_all(all_object_idx, :), all_object_mice_ids, 'object', options);
end

end

function createAllMiceDiscoveryPlot(sessions, mouse_ids, type, options)
% Create discovery aligned plots for either water or object for all mice
% sessions - Cell array with filtered sessions of one type
% mouse_ids - Cell array with corresponding mouse IDs
% type - 'water' or 'object'
% options - Plot options

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_dist_water = 9;  % Distance to water
col_dist_object = 10; % Distance to object
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)

num_sessions = size(sessions, 1);

% Initialize matrices to store aligned data and time
all_aligned_dff = cell(num_sessions, 1);
all_aligned_time = cell(num_sessions, 1);
all_mouse_ids = cell(num_sessions, 1);

% IMPORTANT: First collect and process all data BEFORE applying time window
% This ensures consistent data processing regardless of the window
for i = 1:num_sessions
    data = sessions{i, 4};
    discovery = sessions{i, 6};
    current_mouse_id = mouse_ids{i};
    
    % Get the appropriate discovery index based on type
    if strcmp(type, 'water')
        discovery_idx = discovery(1); % First element is water discovery
    else % object
        discovery_idx = discovery(2); % Second element is object discovery
    end
    
    % Skip if discovery index is NaN
    if isnan(discovery_idx)
        continue;
    end
    
    % Get the discovery time in seconds
    discovery_time_sec = data(discovery_idx, col_time);
    
    % Extract all data and convert time to seconds from discovery
    time_sec = data(:, col_time) - discovery_time_sec;
    dff = data(:, col_dff);
    
    % Create mask to exclude grooming and climbing data
    exclude_mask = false(size(data, 1), 1);
    
    % Check for grooming (column 13)
    if size(data, 2) >= col_grooming
        exclude_mask = exclude_mask | (data(:, col_grooming) > 0);
    end
    
    % Check for climbing (column 14) if it exists
    if size(data, 2) >= col_climbing
        exclude_mask = exclude_mask | (data(:, col_climbing) > 0);
    end
    
    % Apply mask to exclude grooming/climbing (but don't filter by time window yet)
    mask = ~exclude_mask;
    
    filtered_time = time_sec(mask);
    filtered_dff = dff(mask);
    
    % Store aligned data with all time points (not filtered by window)
    all_aligned_dff{i} = filtered_dff;
    all_aligned_time{i} = filtered_time;
    all_mouse_ids{i} = current_mouse_id;
end

% Apply z-scoring to each session's full dataset if requested
if options.z_score
    for i = 1:num_sessions
        if ~isempty(all_aligned_dff{i})
            % Z-score based on entire signal, not just window
            all_aligned_dff{i} = (all_aligned_dff{i} - nanmean(all_aligned_dff{i})) / nanstd(all_aligned_dff{i});
        end
    end
end

% Apply smoothing if requested (to entire signals)
if options.smooth_window > 0
    for i = 1:num_sessions
        if ~isempty(all_aligned_dff{i})
            all_aligned_dff{i} = movmean(all_aligned_dff{i}, options.smooth_window, 'omitnan');
        end
    end
end

% AFTER processing full signals, NOW determine the common time grid for plotting
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000); % 1000 time points for smooth visualization

% Interpolate data onto common time grid for each session
aligned_dff = nan(num_sessions, length(time_grid));

for i = 1:num_sessions
    % Skip if no valid data
    if isempty(all_aligned_time{i}) || all(isnan(all_aligned_dff{i}))
        continue;
    end
    
    % Interpolate DFF onto common time grid
    % Only interpolate within the valid range of the original data
    valid_times = time_grid >= min(all_aligned_time{i}) & time_grid <= max(all_aligned_time{i});
    
    if any(valid_times)
        aligned_dff(i, valid_times) = interp1(all_aligned_time{i}, all_aligned_dff{i}, time_grid(valid_times), 'linear', NaN);
    end
end

% Remove rows with all NaNs (no discovery or invalid data)
valid_rows = ~all(isnan(aligned_dff), 2);
aligned_dff = aligned_dff(valid_rows, :);
valid_mouse_ids = all_mouse_ids(valid_rows);
valid_sessions = sessions(valid_rows, :);
num_valid_sessions = sum(valid_rows);

if num_valid_sessions == 0
    warning(['No valid discovery events found for ' type]);
    return;
end

% Get unique mouse IDs for grouping
[unique_mouse_ids, ~, mouse_group_idx] = unique(valid_mouse_ids);
num_mice = length(unique_mouse_ids);

% Calculate per-mouse average
mouse_avg_dff = nan(num_mice, length(time_grid));
for i = 1:num_mice
    mouse_mask = (mouse_group_idx == i);
    mouse_avg_dff(i, :) = nanmean(aligned_dff(mouse_mask, :), 1);
end

% Create figure with two subplots
figure('Position', [100, 100, 1200, 600]);
sgtitle(['All Mice - ' upper(type(1)) type(2:end) ' Discovery Aligned Neural Activity (Excluding Grooming & Climbing)']);

% 1. Average DFF plot across all mice (left subplot)
subplot(1, 2, 1);
hold on;

% Calculate grand average and SEM across mice
grand_avg_dff = nanmean(mouse_avg_dff, 1);
sem_dff = nanstd(mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(mouse_avg_dff), 1));

% Add discovery line
xline(0, '--k', 'LineWidth', 1.5);

% Plot SEM shading
ciplot(grand_avg_dff - sem_dff, grand_avg_dff + sem_dff, time_grid,  [0.5, 0.5, 0.5], 0.3);

% Plot the mean line on top
plot(time_grid, grand_avg_dff, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);

% Format plot
xlabel('Time from discovery (seconds)');
if options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
ylim([-1.5 2]);

title({'Average neural activity across all mice', ...
       ['n = ' num2str(num_mice) ' mice, ' num2str(num_valid_sessions) ' total sessions']});
grid off;
box off;

% Set axis limits
xlim(options.time_window);

% 2. Raster plot by mouse (right subplot)
subplot(1, 2, 2);

% Create a grouped raster plot by mouse
raster_data = nan(num_mice, length(time_grid));
for i = 1:num_mice
    raster_data(i, :) = mouse_avg_dff(i, :);
end

% Plot mouse averages as a raster
imagesc(time_grid, 1:num_mice, raster_data);
colormap(gca, blueWhiteRed(256));
caxis(options.caxis_range);

% Add discovery line
hold on;
xline(0, '--k', 'LineWidth', 1.5);

% Format plot
xlabel('Time from discovery (seconds)');
ylabel('Mouse #', 'FontSize', 12);

title(['Individual mouse responses to ' type ' discovery']);

% Add colorbar with appropriate label
h = colorbar;
if options.z_score
    ylabel(h, 'ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel(h, 'ΔF/F', 'FontSize', 12);
end

% Set x-axis limits
xlim(options.time_window);

% Add mouse IDs as y-tick labels
set(gca, 'YTick', 1:num_mice, 'YTickLabel', unique_mouse_ids);


% Output some statistics
fprintf('%s discovery analysis complete: %d mice, %d total sessions\n', type, num_mice, num_valid_sessions);
fprintf('Average pre-discovery DFF: %.3f\n', nanmean(grand_avg_dff(time_grid < 0)));
fprintf('Average post-discovery DFF: %.3f\n', nanmean(grand_avg_dff(time_grid > 0)));
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
