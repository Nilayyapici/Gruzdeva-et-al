function plotFoodEventsAligned_sess(mice_all_reorganized, options)
% PLOTFOODEEVENTSALIGNED_SESSIONS Creates food visit/eating event-aligned visualization for specific sessions
% FIXED: Average trace now uses same mice as heatmap (mice with >2 events)

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
    options.session = 'sess1'; % Default to session 1
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
    options.z_score = false;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_event_duration')
    options.min_event_duration = 1; % Default 1 second
end

% Determine post_discovery_only default based on session
if ~isfield(options, 'post_discovery_only')
    if strcmpi(options.session, 'sess0')
        options.post_discovery_only = false; % Session 0 is before discovery
    else
        options.post_discovery_only = true; % All other sessions are after discovery
    end
end

% Filter mice based on state, source, and session
n = size(mice_all_reorganized, 1);
mask = true(n, 1);

% Filter by state
if ~strcmpi(options.state, 'all')
    mask = mask & strcmp(mice_all_reorganized(:, 2), options.state);
end

% Filter by source
if ~strcmpi(options.source, 'all')
    mask = mask & strcmp(mice_all_reorganized(:, 3), options.source);
end

% Filter by session
if ~strcmpi(options.session, 'all')
    session_mask = false(n, 1);
    for i = 1:n
        mouse_id = mice_all_reorganized{i, 1};
        if contains(mouse_id, ['_', options.session])
            session_mask(i) = true;
        end
    end
    mask = mask & session_mask;
end

% Extract filtered mice
filtered_mice = mice_all_reorganized(mask, :);
num_mice = size(filtered_mice, 1);

if num_mice == 0
    error('No mice found with the specified criteria (state: %s, source: %s, session: %s)', ...
          options.state, options.source, options.session);
end

fprintf('Found %d mice matching criteria (state: %s, source: %s, session: %s)\n', ...
        num_mice, options.state, options.source, options.session);

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

% Determine the common time grid for plotting (do this early for proper interpolation)
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000); % 1000 time points for smooth visualization

% Process each mouse to find events and align data
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    event_signal = data(:, event_col);
    
    % Get the food discovery frame (this is from the original data, not session-specific)
    discovery_frame = filtered_mice{i, 6};
    
    % Find onsets and offsets of events
    event_changes = diff([0; event_signal; 0]); % Add padding and compute differences
    event_onsets = find(event_changes == 1);
    event_offsets = find(event_changes == -1) - 1;
    
    % Calculate durations in seconds
    event_durations = time(event_offsets) - time(event_onsets);
    
    % Filter for minimum duration
    valid_events = event_durations >= options.min_event_duration;
    
    event_onsets = event_onsets(valid_events);
    event_offsets = event_offsets(valid_events);
    event_durations = event_durations(valid_events);
    
    num_events = length(event_onsets);
    all_mouse_event_counts(i) = num_events;
    
    if num_events == 0
        fprintf('No valid %s events found for mouse %d in session %s\n', ...
                lower(event_name), i, options.session);
        continue;
    end
    
    % Initialize matrix to store interpolated data for all events
    events_interp = nan(num_events, length(time_grid));
    
    % Process each event - FIXED: Interpolate FIRST, then smooth
    for j = 1:num_events
        onset_idx = event_onsets(j);
        onset_time = time(onset_idx);
        
        % Find indices within the time window, ensuring onset is included
        time_mask = time >= (onset_time + options.time_window(1)) & ...
                    time <= (onset_time + options.time_window(2));
        
        % Ensure the onset frame is included
        time_mask(onset_idx) = true;
        
        % Extract aligned data
        aligned_time = time(time_mask) - onset_time;
        aligned_dff = dff(time_mask);
        
        % Check if we have enough data points
        if length(aligned_time) < 3
            continue;
        end
        
        % Apply z-scoring if requested (before interpolation)
        if options.z_score
            aligned_dff = (aligned_dff - nanmean(aligned_dff)) / nanstd(aligned_dff);
        end
        
        % FIXED: Interpolate onto common time grid WITHOUT smoothing first
        try
            events_interp(j, :) = interp1(aligned_time, aligned_dff, time_grid, 'linear', NaN);
        catch
            % If interpolation fails, skip this event
            continue;
        end
    end
    
    % Store interpolated events for this mouse
    all_mouse_aligned_dff{i} = events_interp;
end

% Apply smoothing to the interpolated data if requested
if options.smooth_window > 0
    for i = 1:num_mice
        if all_mouse_event_counts(i) > 0
            for j = 1:size(all_mouse_aligned_dff{i}, 1)
                % Apply smoothing to each interpolated event
                all_mouse_aligned_dff{i}(j, :) = movmean(all_mouse_aligned_dff{i}(j, :), ...
                                                         options.smooth_window, 'omitnan');
            end
        end
    end
end

% Calculate mouse averages (average across events for each mouse)
mouse_avg_dff = nan(num_mice, length(time_grid));
for i = 1:num_mice
    if all_mouse_event_counts(i) > 0
        mouse_avg_dff(i, :) = nanmean(all_mouse_aligned_dff{i}, 1);
    end
end

% CRITICAL FIX: Filter for mice with sufficient events BEFORE averaging
% This ensures the average trace uses the SAME mice as the heatmap
mice_with_sufficient_events = all_mouse_event_counts > 0;
filtered_mouse_avg_dff = mouse_avg_dff(mice_with_sufficient_events, :);
filtered_mice_info = filtered_mice(mice_with_sufficient_events, :);
filtered_event_counts = all_mouse_event_counts(mice_with_sufficient_events);
num_mice_raster = sum(mice_with_sufficient_events);

% Create figure with two subplots
figure('Position', [100, 100, 1000, 600]);

% 1. Average DFF plot (left subplot)
subplot(1, 2, 1);
hold on;

% FIXED: Calculate average ONLY from mice with >2 events (same as heatmap)
avg_dff = nanmean(filtered_mouse_avg_dff, 1);
sem_dff = nanstd(filtered_mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(filtered_mouse_avg_dff), 1));

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
if options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
title({sprintf('ΔF/F aligned to %s onset', event_name), ...
       sprintf('Session: %s (n=%d mice)', options.session, num_mice_raster)});
grid off;
box off;

% Set x-axis limits
xlim(options.time_window);

% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% 2. Raster plot (right subplot) - Shows SAME mice as average
subplot(1, 2, 2);

if num_mice_raster == 0
    % If no mice have >2 events, show a message
    text(0.5, 0.5, {'No mice with >2 events', 'Cannot display raster plot'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'center', ...
         'FontSize', 14, 'Units', 'normalized');
    xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
    ylabel('Mouse #', 'FontSize', 12);
    title({'Activity across mice', ...
           sprintf('Session: %s - No mice with >2 events', options.session)});
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
           sprintf('Session: %s', options.session), ...
           sprintf('n = %d mice, %d events', num_mice_raster, sum(filtered_event_counts))});
    
    % Create colormap based on event type
    if strcmpi(options.event_type, 'visit')
        colormap(bluewhitered(256, true)); % Blue-white-red but emphasize blue
    else % eating
        colormap(bluewhitered(256, false)); % Blue-white-red but emphasize red
    end
    h = colorbar;
    if options.z_score
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
                % Get base mouse ID (without session suffix)
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

% Add overall title with session and z-score information
if options.z_score
    z_score_text = 'Z-scored';
else
    z_score_text = 'Raw';
end

% Create discovery text based on session
if strcmpi(options.session, 'sess0')
    discovery_text = '(before discovery)';
else
    discovery_text = '(post-discovery)';
end

sgtitle({[z_score_text, ' ΔF/F values ', discovery_text], ...
         [options.state, ' mice, ', options.source, ' - Session: ', options.session]}, ...
         'FontWeight', 'bold');

% Output some statistics
fprintf('\nAnalysis complete for session %s:\n', options.session);
fprintf('  %d total mice found\n', num_mice);
fprintf('  %d mice with >2 events (used for average and heatmap)\n', num_mice_raster);
fprintf('  %d total %s events\n', sum(all_mouse_event_counts), lower(event_name));
if num_mice_raster > 0
    fprintf('  Average pre-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid < 0)));
    fprintf('  Average during-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid >= 0)));
end

% Optional: return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.mouse_avg_dff = mouse_avg_dff;  % All mice
    aligned_data.filtered_mouse_avg_dff = filtered_mouse_avg_dff;  % Only mice with >2 events
    aligned_data.mice_with_sufficient_events = mice_with_sufficient_events;
    aligned_data.time_grid = time_grid;
    aligned_data.event_counts = all_mouse_event_counts;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.all_events = all_mouse_aligned_dff;
    aligned_data.session = options.session;
    varargout{1} = aligned_data;
end
end

function ciplot(lower, upper, x, color, alpha)
% ciplot creates a shaded area to visualize confidence intervals
x_polygon = [x, fliplr(x)];
y_polygon = [upper, fliplr(lower)];

nan_indices = isnan(x_polygon) | isnan(y_polygon);
x_polygon(nan_indices) = [];
y_polygon(nan_indices) = [];

h = fill(x_polygon, y_polygon, color);
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', alpha);
end

function cmap = bluewhitered(m, emphasize_blue)
if nargin < 1
   m = 100;
end

if nargin < 2
   emphasize_blue = true;
end

if emphasize_blue
    colors = [
        0.0, 0.2, 0.6;
        0.1, 0.4, 0.8;
        0.4, 0.7, 1.0;
        0.7, 0.9, 1.0;
        1.0, 1.0, 1.0;
        1.0, 0.8, 0.8;
        0.9, 0.4, 0.4;
        0.7, 0.1, 0.1
    ];
    positions = [0, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 1.0];
else
    colors = [
        0.0, 0.2, 0.6;
        0.3, 0.5, 0.8;
        0.6, 0.8, 1.0;
        1.0, 1.0, 1.0;
        1.0, 0.8, 0.8;
        1.0, 0.6, 0.6;
        0.9, 0.3, 0.3;
        0.7, 0.0, 0.0
    ];
    positions = [0, 0.15, 0.30, 0.45, 0.55, 0.7, 0.85, 1.0];
end

n = size(colors, 1);
cmap = zeros(m, 3);

idx = round(positions * (m-1)) + 1;

for i = 1:n-1
    i1 = idx(i);
    i2 = idx(i+1);
    ni = i2 - i1 + 1;
    
    r = linspace(colors(i,1), colors(i+1,1), ni);
    g = linspace(colors(i,2), colors(i+1,2), ni);
    b = linspace(colors(i,3), colors(i+1,3), ni);
    
    cmap(i1:i2, :) = [r', g', b'];
end
end