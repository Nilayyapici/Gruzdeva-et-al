function aligned_data = plotFoodEventsAligned(mice_all, options)
% PLOTFOODEEVENTSALIGNED Creates food visit/eating event-aligned visualization with raster and average plots
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - event_type: 'visit' or 'eating' (default: 'visit')
%     - time_window: two-element vector specifying time range around event onset [before after] in seconds (default: [-20 60])
%     - ylim_avg: y-axis limits for average plot [min max] (default: auto)
%     - caxis_range: two-element vector specifying color axis range for DFF (default: [-2 2])
%     - z_score: whether to z-score the DFF signals (default: false)
%     - smooth_window: window size for smoothing in frames (default: 15)
%     - min_event_duration: minimum duration of events to include in seconds (default: 1)
%
% NOTE on mouse filtering:
%   - The average trace (left plot + Excel 'average_trace' sheet) uses
%     ALL mice that pass the state/source filter, for maximum
%     statistical power.
%   - The raster plot (right plot) and the Excel 'heatmap_per_mouse'
%     sheet are both restricted to mice with MORE THAN 3 valid events,
%     and always show the exact same mouse set as each other.

% Constants
COL_TIME = 1;     % Time column
COL_VISIT = 8;    % Food visit column
COL_EATING = 9;   % Eating column
COL_DFF = 11;     % DFF data column

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

if ~isfield(options, 'session')
    options.session = 'all';
end

if ~isfield(options, 'event_type')
    options.event_type = 'visit';
end

if ~isfield(options, 'time_window')
    options.time_window = [-20, 60]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'z_score')
    options.z_score = false; % Default is now false (no z-scoring)
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_event_duration')
    options.min_event_duration = 1; % Default 1 second
end

% Minimum number of events a mouse needs to be included anywhere
% (plots AND Excel export)
MIN_EVENTS_TO_INCLUDE = 3;

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

% Determine which column to use based on event_type
if strcmpi(options.event_type, 'visit')
    event_col = COL_VISIT;
    event_name = 'Food Visit';
elseif strcmpi(options.event_type, 'eating')
    event_col = COL_EATING;
    event_name = 'Eating';
else
    error('Invalid event_type. Must be either "visit" or "eating"');
end

% Process each mouse to find events and align data
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    event_signal = data(:, event_col);

    % Get the food discovery frame
    discovery_frame = filtered_mice{i, 6};
    discovery_time_sec = data(discovery_frame, COL_TIME);

    % Find onsets and offsets of events
    event_changes = diff([0; event_signal; 0]); % Add padding and compute differences
    event_onsets = find(event_changes == 1);
    event_offsets = find(event_changes == -1) - 1;

    % Calculate durations in seconds
    event_durations = time(event_offsets) - time(event_onsets);

    % Filter for minimum duration and events after discovery
    valid_events = event_durations >= options.min_event_duration & ...
                   time(event_onsets) > discovery_time_sec; % Only include events after discovery

    event_onsets = event_onsets(valid_events);
    event_offsets = event_offsets(valid_events);
    event_durations = event_durations(valid_events);

    num_events = length(event_onsets);
    all_mouse_event_counts(i) = num_events;

    if num_events == 0
        fprintf('No valid %s events found for mouse %d after food discovery\n', lower(event_name), i);
        continue;
    end

    % Initialize cell array to store aligned data for each event
    all_events_aligned_dff = cell(num_events, 1);
    all_events_aligned_time = cell(num_events, 1);

    % Process each event
    for j = 1:num_events
        onset = event_onsets(j);

        % Calculate time range to extract (in indices)
        onset_time = time(onset);

        % Find indices within the time window
        time_mask = time >= (onset_time + options.time_window(1)) & ...
                    time <= (onset_time + options.time_window(2));

        % Extract aligned data
        aligned_time = time(time_mask) - onset_time;
        aligned_dff = dff(time_mask);

        % Apply z-scoring if requested
        if isfield(options, 'z_score') && options.z_score
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

% ── Apply the >3-events filter ONCE, up front, so every downstream
% output (average trace, raster, Excel) uses the exact same mouse set ──
mice_with_sufficient_events = all_mouse_event_counts > MIN_EVENTS_TO_INCLUDE;
filtered_mouse_avg_dff = mouse_avg_dff(mice_with_sufficient_events, :);
filtered_mice_info = filtered_mice(mice_with_sufficient_events, :);
filtered_event_counts = all_mouse_event_counts(mice_with_sufficient_events);
num_mice_raster = sum(mice_with_sufficient_events);

% Create figure with two subplots
figure('Position', [100, 100, 1000, 600]);

% 1. Average DFF plot (left subplot)
subplot(1, 2, 1);
hold on;

% Calculate overall average DFF across ALL mice that passed the
% state/source filter (not just the >3-events subset). This gives the
% average trace the full statistical power of all available mice, even
% though the raster plot and Excel export below are restricted to mice
% with >3 events.
avg_dff = nanmean(mouse_avg_dff, 1);
sem_dff = nanstd(mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(mouse_avg_dff), 1));

% Add event onset line
xline(0, '--k', 'LineWidth', 1.5);

% Determine colors based on event type
if strcmpi(options.event_type, 'visit')
    event_color = [0.0, 0.4, 0.8]; % Blue for visits
    sem_color = [0.6, 0.8, 1.0]; % Lighter blue for SEM
else % eating
    event_color = [0.8, 0.2, 0.2]; % Red for eating
    sem_color = [1.0, 0.7, 0.7]; % Lighter red for SEM
end

% Add SEM shading
ciplot(avg_dff - sem_dff, avg_dff + sem_dff, time_grid, sem_color, 0.5);

% Plot the mean line on top
plot(time_grid, avg_dff, 'Color', event_color, 'LineWidth', 2);

% Format plot
xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
if isfield(options, 'z_score') && options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
title({'ΔF/F aligned to onset', ...
       sprintf('%s', event_name)});
grid off;
box off;

% Set x-axis limits
xlim(options.time_window);

% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% 2. Raster plot (right subplot) - ONLY SHOWS MICE WITH >3 EVENTS
subplot(1, 2, 2);

if num_mice_raster == 0
    % If no mice have >3 events, show a message
    text(0.5, 0.5, {'No mice with >3 events', 'Cannot display raster plot'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'center', ...
         'FontSize', 14, 'Units', 'normalized');
    xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
    ylabel('Mouse #', 'FontSize', 12);
    title({'Activity across mice', ...
           'No mice with >3 events'});
else
    % Create raster image with filtered data
    imagesc(time_grid, 1:num_mice_raster, filtered_mouse_avg_dff);

    % Add event onset line
    hold on;
    xline(0, '--k', 'LineWidth', 1.5);

    % Format plot
    xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
    ylabel('Mouse #', 'FontSize', 12);
    title({'Activity across mice', ...
           ['n = ' num2str(num_mice_raster) ' mice, ' ...
            num2str(sum(filtered_event_counts)) ' events']});

    % Create colormap based on event type
    if strcmpi(options.event_type, 'visit')
        colormap(bluewhitered(256, true)); % Blue-white-red but emphasize blue
    else % eating
        colormap(bluewhitered(256, false)); % Blue-white-red but emphasize red
    end
    h = colorbar;
    if isfield(options, 'z_score') && options.z_score
        ylabel(h, 'ΔF/F (z-scored)', 'FontSize', 12);
    else
        ylabel(h, 'ΔF/F', 'FontSize', 12);
    end
    caxis(options.caxis_range);

    % Add mouse IDs and count of events for filtered mice
    if num_mice_raster <= 30
        mouse_labels = cell(num_mice_raster, 1);
        for i = 1:num_mice_raster
            if ischar(filtered_mice_info{i, 1})
                id_parts = strsplit(filtered_mice_info{i, 1}, '_');
                if length(id_parts) > 1
                    mouse_id = id_parts{1};
                else
                    mouse_id = filtered_mice_info{i, 1};
                end
            else
                mouse_id = num2str(i);
            end
            mouse_labels{i} = [mouse_id, ' (', num2str(filtered_event_counts(i)), ')'];
        end
        set(gca, 'YTick', 1:num_mice_raster, 'YTickLabel', mouse_labels);
    else
        % Just show mouse numbers if there are too many
        set(gca, 'YTick', round(linspace(1, num_mice_raster, 10)));
    end

    % Set x-axis limits
    xlim(options.time_window);
end

% Add overall title with z-score information
if isfield(options, 'z_score') && options.z_score
    z_score_text = 'Z-scored';
else
    z_score_text = 'Raw';
end
sgtitle({[z_score_text, ' ΔF/F values (post-discovery only) ', ...
         options.state, ' mice, ', options.source], ...
         }, 'FontWeight', 'bold');

%% Export data to Excel
% IMPORTANT: the two sheets below intentionally use DIFFERENT mouse sets,
% matching the two plots:
%   - Sheet 1 (average_trace) mirrors the LEFT panel: mean/SEM across ALL
%     mice that passed the state/source filter (num_mice), for maximum
%     statistical power.
%   - Sheet 2 (heatmap_per_mouse) mirrors the RIGHT panel (raster): only
%     the mice with >3 events (num_mice_raster), so it exactly matches
%     what's shown in the raster plot — no extra/excluded mice.
excel_filename = sprintf('event_aligned_%s_%s_%s_%s.xlsx', ...
    options.event_type, options.state, options.source, datestr(now, 'yyyymmdd'));

% ── Sheet 1: Average trace (mean and SEM across ALL filtered mice) ───────
avg_table = table(time_grid(:), avg_dff(:), sem_dff(:), ...
    'VariableNames', {'Time_s', 'Mean_dFF', 'SEM_dFF'});

writetable(avg_table, excel_filename, 'Sheet', 'average_trace');

% ── Sheet 2: Per-mouse average (one row per PLOTTED mouse only) ──────────
% Build heatmap cell: first row = time points, then one row per mouse
% that appears in the raster plot (num_mice_raster rows, not num_mice)
heatmap_cell = cell(num_mice_raster + 1, length(time_grid) + 1);

% Header row
heatmap_cell{1, 1} = 'Mouse_ID';
for t = 1:length(time_grid)
    heatmap_cell{1, t+1} = time_grid(t);
end

% One row per mouse that passed the >3-events filter
for i = 1:num_mice_raster
    % Get mouse ID
    id_parts = strsplit(filtered_mice_info{i, 1}, '_');
    if length(id_parts) > 1
        mouse_id = id_parts{1};
    else
        mouse_id = filtered_mice_info{i, 1};
    end
    heatmap_cell{i+1, 1} = [mouse_id, ' (', num2str(filtered_event_counts(i)), ')'];
    for t = 1:length(time_grid)
        heatmap_cell{i+1, t+1} = filtered_mouse_avg_dff(i, t);
    end
end

writecell(heatmap_cell, excel_filename, 'Sheet', 'heatmap_per_mouse');

fprintf('Event-aligned data exported to: %s\n', excel_filename);

% Output some statistics
fprintf('Analysis complete: %d total mice, %d mice with >%d events (plotted & exported), %d total %s events (after discovery, plotted mice only)\n', ...
        sum(all_mouse_event_counts > 0), num_mice_raster, MIN_EVENTS_TO_INCLUDE, sum(filtered_event_counts), lower(event_name));
fprintf('Average pre-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid < 0)));
fprintf('Average during-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid >= 0)));

% Return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.mouse_avg_dff = mouse_avg_dff; % all mice, unfiltered (kept for reference)
    aligned_data.filtered_mouse_avg_dff = filtered_mouse_avg_dff; % mice actually plotted/exported
    aligned_data.mice_with_sufficient_events = mice_with_sufficient_events; % filter mask
    aligned_data.time_grid = time_grid;
    aligned_data.event_counts = all_mouse_event_counts;
    aligned_data.filtered_event_counts = filtered_event_counts;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.filtered_mouse_info = filtered_mice_info;
    aligned_data.all_events = all_mouse_aligned_dff;
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

function cmap = bluewhitered(m, emphasize_blue)
% BLUEWHITERED creates a blue-to-white-to-red colormap
%   BLUEWHITERED(M, EMPHASIZE_BLUE) returns an M-by-3 matrix containing a colormap
%   that transitions from dark blue through a small white region to true red
%   If EMPHASIZE_BLUE is true, blue tones are emphasized; otherwise red tones are emphasized

if nargin < 1
   m = 100;
end

if nargin < 2
   emphasize_blue = true; % Default is to emphasize blue
end

% Define the color transitions
if emphasize_blue
    % Blue emphasized for visit events
    colors = [
        0.0, 0.2, 0.6;  % Dark blue
        0.1, 0.4, 0.8;  % Medium blue
        0.4, 0.7, 1.0;  % Light blue
        0.7, 0.9, 1.0;  % Very light blue
        1.0, 1.0, 1.0;  % White (center point)
        1.0, 0.8, 0.8;  % Light red
        0.9, 0.4, 0.4;  % Medium red
        0.7, 0.1, 0.1   % Dark red
    ];
    % More blue positions
    positions = [0, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 1.0];
else
    % Red emphasized for eating events
    colors = [
        0.0, 0.2, 0.6;  % Dark blue
        0.3, 0.5, 0.8;  % Medium blue
        0.6, 0.8, 1.0;  % Light blue
        1.0, 1.0, 1.0;  % White (center point)
        1.0, 0.8, 0.8;  % Light red
        1.0, 0.6, 0.6;  % Medium light red
        0.9, 0.3, 0.3;  % Medium red
        0.7, 0.0, 0.0   % Dark red
    ];
    % More red positions
    positions = [0, 0.15, 0.30, 0.45, 0.55, 0.7, 0.85, 1.0];
end

n = size(colors, 1);
cmap = zeros(m, 3);

% Interpolate colors between key points
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