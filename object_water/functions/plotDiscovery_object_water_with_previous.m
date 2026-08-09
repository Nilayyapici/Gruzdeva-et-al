function plotDiscovery_object_water_with_previous(mice_all, options)
% PLOTDISCOVERY_OBJECT_WATER_WITH_PREVIOUS Creates discovery-aligned visualization 
% including data from the previous session to have more baseline data
%
% For water discovery: includes "before" session data
% For object discovery: includes "water" session data
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

% Build a lookup structure for quick session finding
% Structure: session_map{mouse_id}{session_type} = index
session_map = struct();
for i = 1:size(mice_all, 1)
    session_id = mice_all{i, 1};
    if ischar(session_id)
        parts = strsplit(session_id, '_');
        mouse_id = parts{1};
        
        if contains(session_id, 'before')
            session_type = 'before';
        elseif contains(session_id, 'water')
            session_type = 'water';
        elseif contains(session_id, 'object')
            session_type = 'object';
        else
            continue;
        end
        
        % Store the index for this mouse/session combination
        if ~isfield(session_map, mouse_id)
            session_map.(mouse_id) = struct();
        end
        session_map.(mouse_id).(session_type) = i;
    end
end

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
        mouse_id = parts{1};
        
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

% Create water discovery plot for all mice (with before session prepended)
if ~isempty(all_water_sessions)
    createAllMiceDiscoveryPlot_with_previous(mice_all, all_water_idx, all_water_mice_ids, 'water', session_map, options);
end

% Create object discovery plot for all mice (with water session prepended)
if ~isempty(all_object_sessions)
    createAllMiceDiscoveryPlot_with_previous(mice_all, all_object_idx, all_object_mice_ids, 'object', session_map, options);
end

end

function createAllMiceDiscoveryPlot_with_previous(mice_all, session_indices, mouse_ids, type, session_map, options)
% Create discovery aligned plots for either water or object for all mice
% INCLUDING data from the previous session
%
% mice_all - Full cell array with all sessions
% session_indices - Indices of sessions to process
% mouse_ids - Cell array with corresponding mouse IDs
% type - 'water' or 'object'
% session_map - Structure for looking up session indices
% options - Plot options

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_dist_water = 9;  % Distance to water
col_dist_object = 10; % Distance to object
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)

num_sessions = length(session_indices);

% Initialize matrices to store aligned data and time
all_aligned_dff = cell(num_sessions, 1);
all_aligned_time = cell(num_sessions, 1);
all_mouse_ids = cell(num_sessions, 1);

% IMPORTANT: First collect and process all data BEFORE applying time window
% This ensures consistent data processing regardless of the window
for i = 1:num_sessions
    session_idx = session_indices(i);
    current_mouse_id = mouse_ids{i};
    
    % Get current session data
    current_data = mice_all{session_idx, 4};
    discovery = mice_all{session_idx, 6};
    
    % Get the appropriate discovery index based on type
    if strcmp(type, 'water')
        discovery_idx = discovery(1); % First element is water discovery
        prev_session_type = 'before';
    else % object
        discovery_idx = discovery(2); % Second element is object discovery
        prev_session_type = 'water';
    end
    
    % Skip if discovery index is NaN
    if isnan(discovery_idx)
        continue;
    end
    
    % Try to find and prepend previous session data
    prev_data = [];
    time_offset = 0;
    
    if isfield(session_map, current_mouse_id) && isfield(session_map.(current_mouse_id), prev_session_type)
        prev_session_idx = session_map.(current_mouse_id).(prev_session_type);
        prev_data = mice_all{prev_session_idx, 4};
        
        % Calculate time offset (end of previous session to start of current session)
        if ~isempty(prev_data)
            time_offset = prev_data(end, col_time) - current_data(1, col_time);
            fprintf('Mouse %s, %s session: prepending %d frames from %s session\n', ...
                current_mouse_id, type, size(prev_data, 1), prev_session_type);
        end
    end
    
    % Concatenate previous and current session data
    if ~isempty(prev_data)
        % Adjust current data time to be continuous with previous data
        adjusted_current_data = current_data;
        adjusted_current_data(:, col_time) = current_data(:, col_time) + time_offset;
        
        % Concatenate
        combined_data = [prev_data; adjusted_current_data];
        
        % Adjust discovery index to account for prepended data
        adjusted_discovery_idx = discovery_idx + size(prev_data, 1);
    else
        % No previous data available, use current session only
        combined_data = current_data;
        adjusted_discovery_idx = discovery_idx;
        fprintf('Mouse %s, %s session: no previous session found, using current session only\n', ...
            current_mouse_id, type);
    end
    
    % Get the discovery time in seconds
    discovery_time_sec = combined_data(adjusted_discovery_idx, col_time);
    
    % Extract all data and convert time to seconds from discovery
    time_sec = combined_data(:, col_time) - discovery_time_sec;
    dff = combined_data(:, col_dff);
    
    % Create mask to exclude grooming and climbing data
    exclude_mask = false(size(combined_data, 1), 1);
    
    % Check for grooming (column 13)
    if size(combined_data, 2) >= col_grooming
        exclude_mask = exclude_mask | (combined_data(:, col_grooming) > 0);
    end
    
    % Check for climbing (column 14) if it exists
    if size(combined_data, 2) >= col_climbing
        exclude_mask = exclude_mask | (combined_data(:, col_climbing) > 0);
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
    
    % Remove duplicate time points by averaging DFF values at same timepoints
    [unique_time, ~, ic] = unique(all_aligned_time{i});
    unique_dff = accumarray(ic, all_aligned_dff{i}, [], @nanmean);
    
    % Skip if we don't have enough unique points
    if length(unique_time) < 2
        continue;
    end
    
    % Interpolate DFF onto common time grid
    % Only interpolate within the valid range of the original data
    valid_times = time_grid >= min(unique_time) & time_grid <= max(unique_time);
    
    if any(valid_times)
        aligned_dff(i, valid_times) = interp1(unique_time, unique_dff, time_grid(valid_times), 'linear', NaN);
    end
end

% Remove rows with all NaNs (no discovery or invalid data)
valid_rows = ~all(isnan(aligned_dff), 2);
aligned_dff = aligned_dff(valid_rows, :);
valid_mouse_ids = all_mouse_ids(valid_rows);
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
sgtitle(['All Mice - ' upper(type(1)) type(2:end) ' Discovery Aligned Neural Activity (With Previous Session Data)']);

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

%% Export data to Excel
excel_filename = sprintf('discovery_aligned_%s_%s.xlsx', type, datestr(now, 'yyyymmdd'));

% ── Sheet 1: Average trace (mean and SEM across mice) ────────────────────
avg_table = table(time_grid(:), grand_avg_dff(:), sem_dff(:), ...
    'VariableNames', {'Time_s', 'Mean_dFF', 'SEM_dFF'});

writetable(avg_table, excel_filename, 'Sheet', 'average_trace');

% ── Sheet 2: Per-mouse average (one row per mouse) ───────────────────────
heatmap_cell = cell(num_mice + 1, length(time_grid) + 1);

% Header row: Mouse_ID then time points
heatmap_cell{1, 1} = 'Mouse_ID';
for t = 1:length(time_grid)
    heatmap_cell{1, t+1} = time_grid(t);
end

% One row per mouse
for i = 1:num_mice
    heatmap_cell{i+1, 1} = unique_mouse_ids{i};
    for t = 1:length(time_grid)
        heatmap_cell{i+1, t+1} = mouse_avg_dff(i, t);
    end
end

writecell(heatmap_cell, excel_filename, 'Sheet', 'heatmap_per_mouse');

fprintf('Discovery-aligned data exported to: %s\n', excel_filename);

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