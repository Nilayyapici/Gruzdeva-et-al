function plotGroomingAlignedData(mice_all, options)
% PLOTGROOMINGALIGNEDDATA Creates a grooming-aligned visualization with raster and average plots
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - time_window: two-element vector specifying time range around grooming onset [before after] in seconds (default: [-20 60])
%     - ylim_avg: y-axis limits for average plot [min max] (default: auto)
%     - caxis_range: two-element vector specifying color axis range for DFF (default: [-2 2])
%     - z_score: whether to z-score the DFF signals (default: true)
%     - smooth_window: window size for smoothing in frames (default: 15)
%     - min_grooming_duration: minimum duration of grooming events to include in seconds (default: 2)

% Constants
COL_TIME = 1;     % Time column
COL_DFF = 11;     % DFF data column
COL_GROOM = 10;   % Grooming column (0/1)

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
    options.time_window = [-20, 60]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_grooming_duration')
    options.min_grooming_duration = 2; % Default 2 seconds
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

% Initialize cell arrays to store aligned data
all_mouse_aligned_dff = cell(num_mice, 1);
all_mouse_event_counts = zeros(num_mice, 1);

% Process each mouse to find grooming events and align data
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    grooming = data(:, COL_GROOM);
    
    % Find onsets and offsets of grooming events
    groom_changes = diff([0; grooming; 0]); % Add padding and compute differences
    groom_onsets = find(groom_changes == 1);
    groom_offsets = find(groom_changes == -1) - 1;
    
    % Calculate durations in seconds
    groom_durations = time(groom_offsets) - time(groom_onsets);
    
    % Filter for minimum duration
    valid_events = groom_durations >= options.min_grooming_duration;
    groom_onsets = groom_onsets(valid_events);
    groom_offsets = groom_offsets(valid_events);
    groom_durations = groom_durations(valid_events);
    
    num_events = length(groom_onsets);
    all_mouse_event_counts(i) = num_events;
    
    if num_events == 0
        fprintf('No valid grooming events found for mouse %d\n', i);
        continue;
    end
    
    % Initialize cell array to store aligned data for each event
    all_events_aligned_dff = cell(num_events, 1);
    all_events_aligned_time = cell(num_events, 1);
    
    % Process each grooming event
    for j = 1:num_events
        onset = groom_onsets(j);
        
        % Calculate time range to extract (in indices)
        onset_time = time(onset);
        
        % Find indices within the time window
        time_mask = time >= (onset_time + options.time_window(1)) & ...
                    time <= (onset_time + options.time_window(2));
        
        % Extract aligned data
        aligned_time = time(time_mask) - onset_time;
        aligned_dff = dff(time_mask);
        
        % Apply z-scoring if requested
        if options.z_score
            aligned_dff = (aligned_dff - nanmean(aligned_dff)) / nanstd(aligned_dff);
        end
        
        % Apply smoothing if requested
        if options.smooth_window > 0
            aligned_dff = movmean(aligned_dff, options.smooth_window, 'omitnan');
        end
        
        % Store aligned data for this event
        all_events_aligned_dff{j} = aligned_dff;
        all_events_aligned_time{j} = aligned_time;
    end
    
    % Store all events for this mouse
    all_mouse_aligned_dff{i} = struct('events_dff', {all_events_aligned_dff}, ...
                                      'events_time', {all_events_aligned_time});
end

% Determine the common time grid for plotting
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000); % 1000 time points for smooth visualization

% Interpolate all events onto common time grid for each mouse
mouse_avg_dff = nan(num_mice, length(time_grid));

for i = 1:num_mice
    if all_mouse_event_counts(i) == 0
        continue;
    end
    
    % Initialize matrix to store interpolated data for all events
    events_interp = nan(all_mouse_event_counts(i), length(time_grid));
    
    % Interpolate each event onto common time grid
    for j = 1:all_mouse_event_counts(i)
        event_time = all_mouse_aligned_dff{i}.events_time{j};
        event_dff = all_mouse_aligned_dff{i}.events_dff{j};
        
        if isempty(event_time) || isempty(event_dff)
            continue;
        end
        
        % Interpolate this event
        events_interp(j, :) = interp1(event_time, event_dff, time_grid, 'linear', NaN);
    end
    
    % Calculate mouse average across all events
    mouse_avg_dff(i, :) = nanmean(events_interp, 1);
end

% Create figure with two subplots
figure('Position', [100, 100, 1000, 600]);

% 1. Average DFF plot (left subplot)
subplot(1, 2, 1);
hold on;

% Calculate overall average DFF across mice
avg_dff = nanmean(mouse_avg_dff, 1);
sem_dff = nanstd(mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(mouse_avg_dff), 1));

% Add grooming onset line
xline(0, '--k', 'LineWidth', 1.5);

% Add SEM shading
ciplot(avg_dff - sem_dff, avg_dff + sem_dff, time_grid, [0.5, 0, 0.5], 0.3); % Purple

% Plot the mean line on top
plot(time_grid, avg_dff, 'Color', [0.5, 0, 0.5], 'LineWidth', 2); % Purple line

% Format plot
xlabel('Time from grooming onset (seconds)', 'FontSize', 12);
if options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
title({'Average neural activity aligned to grooming onset', ...
       ['State: ' upper(options.state(1)) options.state(2:end), ...
        ', Source: ' upper(options.source(1)) options.source(2:end)]});
grid off;
box off;

% Set x-axis limits
xlim(options.time_window);

% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% 2. Raster plot (right subplot)
subplot(1, 2, 2);

% Create raster image
imagesc(time_grid, 1:num_mice, mouse_avg_dff);

% Add grooming onset line
hold on;
xline(0, '--k', 'LineWidth', 1.5);

% Format plot
xlabel('Time from grooming onset (seconds)', 'FontSize', 12);
ylabel('Mouse #', 'FontSize', 12);
title({'Neural activity across mice', ...
       ['n = ' num2str(sum(all_mouse_event_counts > 0)) ' mice, ' ...
        num2str(sum(all_mouse_event_counts)) ' grooming events']});
colormap(bluewhitered(256));
h = colorbar;
if options.z_score
    ylabel(h, 'ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel(h, 'ΔF/F', 'FontSize', 12);
end
caxis(options.caxis_range);

% Add mouse IDs and count of grooming events
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
        mouse_labels{i} = [mouse_id, ' (', num2str(all_mouse_event_counts(i)), ')'];
    end
    set(gca, 'YTick', 1:num_mice, 'YTickLabel', mouse_labels);
else
    % Just show mouse numbers if there are too many
    set(gca, 'YTick', round(linspace(1, num_mice, 10)));
end

% Set x-axis limits
xlim(options.time_window);

% Add overall title with z-score information
if options.z_score
    z_score_text = 'Z-scored';
else
    z_score_text = 'Raw';
end
sgtitle({['Neural Activity Aligned to Grooming Onset: ', ...
         options.state, ' mice, ', options.source, ' stimulus'], ...
         [z_score_text, ' ΔF/F values']}, 'FontWeight', 'bold');

% Output some statistics
fprintf('Analysis complete: %d mice, %d total grooming events\n', ...
        sum(all_mouse_event_counts > 0), sum(all_mouse_event_counts));
fprintf('Average pre-grooming DFF: %.3f\n', nanmean(avg_dff(time_grid < 0)));
fprintf('Average during-grooming DFF: %.3f\n', nanmean(avg_dff(time_grid >= 0)));

% Optional: return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.mouse_avg_dff = mouse_avg_dff;
    aligned_data.time_grid = time_grid;
    aligned_data.event_counts = all_mouse_event_counts;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.all_events = all_mouse_aligned_dff;
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