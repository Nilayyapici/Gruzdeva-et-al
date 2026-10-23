function plotInteractionAligned(mice_all, options)
% PLOTINTERACTIONALIGNED Creates visualization of neural activity aligned to 
% interaction onsets for water and object interactions
%
% Parameters:
%   mice_all - Cell array with mouse data
%   options - Struct with visualization parameters:
%     - lim_dist: distance threshold for interaction in cm (default: 3)
%     - time_window: time range around interaction [before after] in seconds (default: [-5 10])
%     - caxis_range: color axis range for DFF (default: [-2 2])
%     - z_score: whether to z-score the DFF signals (default: true)
%     - smooth_window: window size for smoothing in frames (default: 10)
%     - min_duration: minimum interaction duration in frames (default: 5)

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

if ~isfield(options, 'lim_dist')
    options.lim_dist = 3; % Default distance threshold in cm
end

if ~isfield(options, 'time_window')
    options.time_window = [-5, 10]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 10;
end

if ~isfield(options, 'min_duration')
    options.min_duration = 5; % Minimum frames to consider a valid interaction
end

% Process water and object interactions separately
processInteractions(mice_all, 'water', options);
processInteractions(mice_all, 'object', options);

end

function processInteractions(mice_all, interaction_type, options)
% Process either water or object interactions
% interaction_type: 'water' or 'object'

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_dist_water = 9;  % Distance to water
col_dist_object = 10; % Distance to object
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)
lim_dist = options.lim_dist;

% Select appropriate distance column based on interaction type
if strcmp(interaction_type, 'water')
    dist_col = col_dist_water;
else
    dist_col = col_dist_object;
end

% Initialize data structures to collect all interactions
all_aligned_dff = {};
all_aligned_time = {};
all_mouse_ids = {};

% Process each session
for i = 1:size(mice_all, 1)
    % Get session data
    data = mice_all{i, 4};
    mouse_id = mice_all{i, 1};
    
    % Extract mouse ID (everything before first underscore)
    parts = strsplit(mouse_id, '_');
    mouse_id_short = parts{1};
    
    % Skip if not applicable session type
    if contains(mouse_id, 'before') || ...
       (strcmp(interaction_type, 'water') && contains(mouse_id, 'object')) || ...
       (strcmp(interaction_type, 'object') && contains(mouse_id, 'water'))
        continue;
    end
    
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
    
    % Get interaction vector
    interaction = zeros(size(data, 1), 1);
    interaction(data(:, dist_col) < lim_dist & ~exclude_mask) = 1;
    
    % Find all interaction onsets
    onsets = find(diff([0; interaction]) == 1);
    offsets = find(diff([interaction; 0]) == -1);
    
    % Calculate interaction durations
    durations = offsets - onsets;
    
    % Filter by minimum duration
    valid_idx = durations >= options.min_duration;
    onsets = onsets(valid_idx);
    
    if isempty(onsets)
        continue; % Skip if no valid interactions
    end
    
    fprintf('Found %d valid %s interactions for mouse %s\n', length(onsets), interaction_type, mouse_id_short);
    
    % Get data aligned to each interaction onset
    time = data(:, col_time);
    dff = data(:, col_dff);
    
    % Smooth DFF if requested
    if options.smooth_window > 0
        dff = movmean(dff, options.smooth_window, 'omitnan');
    end
    
    % Z-score full signal if requested
    if options.z_score
        dff = (dff - nanmean(dff(~exclude_mask))) / nanstd(dff(~exclude_mask));
    end
    
    % Process each interaction onset
    for j = 1:length(onsets)
        onset_time = time(onsets(j));
        
        % Align time to onset
        aligned_time = time - onset_time;
        
        % Filter data to be within the specified time window
        time_mask = aligned_time >= options.time_window(1) & ...
                    aligned_time <= options.time_window(2) & ...
                    ~exclude_mask;
        
        % Skip if not enough data points
        if sum(time_mask) < 10
            continue;
        end
        
        % Store aligned data
        all_aligned_dff{end+1} = dff(time_mask);
        all_aligned_time{end+1} = aligned_time(time_mask);
        all_mouse_ids{end+1} = mouse_id_short;  % Use the short mouse ID (everything before first underscore)
    end
end

% Check if we found any valid interactions
if isempty(all_aligned_dff)
    warning('No valid %s interactions found', interaction_type);
    return;
end

% Interpolate onto common time grid
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 500);

num_interactions = length(all_aligned_dff);
aligned_dff = nan(num_interactions, length(time_grid));

% Interpolate each interaction
for i = 1:num_interactions
    aligned_dff(i, :) = interp1(all_aligned_time{i}, all_aligned_dff{i}, time_grid, 'linear', NaN);
end

% Group by mouse ID
[unique_mouse_ids, ~, mouse_group_idx] = unique(all_mouse_ids);
num_mice = length(unique_mouse_ids);

% Calculate per-mouse average
mouse_avg_dff = nan(num_mice, length(time_grid));
for i = 1:num_mice
    mouse_mask = (mouse_group_idx == i);
    mouse_avg_dff(i, :) = nanmean(aligned_dff(mouse_mask, :), 1);
end

% Create figure
figure('Position', [100, 100, 1200, 600]);
sgtitle([upper(interaction_type(1)) interaction_type(2:end) ' Interaction Aligned Neural Activity']);

% 1. Average DFF plot (left subplot)
subplot(1, 2, 1);
hold on;

% Calculate grand average and SEM
grand_avg_dff = nanmean(mouse_avg_dff, 1);
sem_dff = nanstd(mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(mouse_avg_dff), 1));

% Add onset line
xline(0, '--k', 'LineWidth', 1.5);

% Plot SEM shading
ciplot(grand_avg_dff - sem_dff, grand_avg_dff + sem_dff, time_grid, [0.3, 0.7, 0.9], 0.3);

% Plot mean line
plot(time_grid, grand_avg_dff, 'Color', [0.3, 0.7, 0.9], 'LineWidth', 2);

% Format plot
title({'Average neural activity aligned to interaction onset', ...
       ['n = ' num2str(num_mice) ' mice, ' num2str(num_interactions) ' interactions']});
xlabel('Time from interaction onset (seconds)');
if options.z_score
    ylabel('ΔF/F (z-scored)');
else
    ylabel('ΔF/F');
end
grid off;
xlim(options.time_window);
ylim([-1.5 2]);

% 2. Raster plot (right subplot)
subplot(1, 2, 2);

% Plot mouse averages as a raster
imagesc(time_grid, 1:num_mice, mouse_avg_dff);
colormap(gca, blueWhiteRed(256));
caxis(options.caxis_range);

% Add onset line
hold on;
xline(0, '--k', 'LineWidth', 1.5);

% Count events per mouse for labels
events_per_mouse = zeros(num_mice, 1);
for i = 1:num_mice
    events_per_mouse(i) = sum(mouse_group_idx == i);
end

% Create labels with event counts
mouse_labels = cell(num_mice, 1);
for i = 1:num_mice
    mouse_labels{i} = sprintf('%s (%d)', unique_mouse_ids{i}, events_per_mouse(i));
end

% Format plot
title('Individual mouse responses');
xlabel('Time from interaction onset (seconds)');
ylabel('Mouse #');
set(gca, 'YTick', 1:num_mice, 'YTickLabel', mouse_labels);
xlim(options.time_window);
colorbar;

% Print statistics
fprintf('\n%s Interaction Analysis:\n', upper(interaction_type(1)) + interaction_type(2:end));
fprintf('Total mice: %d\n', num_mice);
fprintf('Total interactions: %d\n', num_interactions);
fprintf('Average pre-interaction DFF: %.3f\n', nanmean(grand_avg_dff(time_grid < 0)));
fprintf('Average during-interaction DFF: %.3f\n', nanmean(grand_avg_dff(time_grid >= 0)));

% =========================================================
% SAVE FIGURE DATA TO EXCEL
% =========================================================
excel_filename = sprintf('InteractionAligned_%s_dist%.1fcm.xlsx', ...
    interaction_type, options.lim_dist);

% --- Sheet 1: Average trace with SEM ---
trace_table = table(time_grid(:), grand_avg_dff(:), sem_dff(:), ...
    grand_avg_dff(:) - sem_dff(:), grand_avg_dff(:) + sem_dff(:), ...
    'VariableNames', {'Time_s', 'Mean_dFF', 'SEM_dFF', 'Lower_CI', 'Upper_CI'});

writetable(trace_table, excel_filename, 'Sheet', 'Average_Trace');

% --- Sheet 2: Per-mouse raster (one row per mouse) ---
time_header = arrayfun(@(t) sprintf('t=%.3f', t), time_grid, 'UniformOutput', false);
col_names   = [{'Mouse_ID', 'N_interactions'}, time_header];

raster_data = [unique_mouse_ids(:), num2cell(events_per_mouse(:)), num2cell(mouse_avg_dff)];
raster_table = cell2table(raster_data, 'VariableNames', col_names);

writetable(raster_table, excel_filename, 'Sheet', 'Raster_PerMouse');

fprintf('Figure data saved to: %s\n', excel_filename);
% =========================================================

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