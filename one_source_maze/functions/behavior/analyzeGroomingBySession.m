function analyzeGroomingBySession(mice_all_reorganized, options)
% ANALYZEGROOOMINGBYSESSION Analyzes grooming events by session
%
% Parameters:
%   mice_all_reorganized - Cell array with reorganized mice data
%   options - Struct with visualization parameters (see plotGroomingAlignedData)
%     Additional options:
%     - sessions: which sessions to analyze, can be specified as:
%         - numeric array [0,1,2,3] 
%         - cell array of strings {'sess0', 'sess1', 'sess2', 'sess3'}
%         (default: all sessions)
%     - combined_figure: whether to create a combined figure comparison (default: true)

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

% Handle different formats for specifying sessions
if ~isfield(options, 'sessions')
    % Default: analyze all sessions
    options.sessions = [0, 1, 2, 3]; 
    options.session_suffixes = {'_sess0', '_sess1', '_sess2', '_sess3'};
else
    % Handle string format (e.g., {'sess0', 'sess1'})
    if iscell(options.sessions) && all(cellfun(@ischar, options.sessions))
        options.session_suffixes = cell(size(options.sessions));
        numeric_sessions = zeros(size(options.sessions));
        
        for i = 1:length(options.sessions)
            sess_str = options.sessions{i};
            if startsWith(sess_str, '_')
                options.session_suffixes{i} = sess_str;
            else
                options.session_suffixes{i} = ['_', sess_str];
            end
            
            % Extract numeric session ID
            if contains(sess_str, 'sess')
                numeric_sessions(i) = str2double(regexprep(sess_str, 'sess', ''));
            else
                numeric_sessions(i) = str2double(regexprep(sess_str, '_', ''));
            end
        end
        options.sessions = numeric_sessions;
    else
        % Handle numeric format (e.g., [0, 1])
        options.session_suffixes = cell(size(options.sessions));
        for i = 1:length(options.sessions)
            options.session_suffixes{i} = ['_sess', num2str(options.sessions(i))];
        end
    end
end

if ~isfield(options, 'combined_figure')
    options.combined_figure = true;
end

% Define session names for reference
session_names = {
    'Session 0: Before discovery, door closed',
    'Session 1: After discovery, door open',
    'Session 2: Second door closed period',
    'Session 3: Second door open period'
};

% Create subsets of data for each requested session
session_data = cell(1, length(options.sessions));
for i = 1:length(options.sessions)
    sess = options.sessions(i);
    suffix = options.session_suffixes{i};
    
    session_mask = false(size(mice_all_reorganized, 1), 1);
    for j = 1:size(mice_all_reorganized, 1)
        mouse_id = mice_all_reorganized{j, 1};
        if contains(mouse_id, suffix)
            session_mask(j) = true;
        end
    end
    session_data{i} = mice_all_reorganized(session_mask, :);
    fprintf('Session %d contains %d mice\n', sess, sum(session_mask));
end

% Analyze each requested session
results = cell(1, length(options.sessions));
for i = 1:length(options.sessions)
    sess = options.sessions(i);
    if sess < 0 || sess > 3
        warning('Invalid session number: %d. Skipping.', sess);
        continue;
    end
    
    sess_mice = session_data{i};
    if isempty(sess_mice)
        warning('No data found for session %d. Skipping.', sess);
        continue;
    end
    
    fprintf('\n==== Analyzing Session %d: %s ====\n', sess, session_names{sess+1});
    
    % Create a copy of options for this session
    sess_options = options;
    
    % Add session information to figure title
    if ~isfield(sess_options, 'figure_title')
        sess_options.figure_title = session_names{sess+1};
    end
    
    % Plot grooming-aligned data for this session
    try
        res = plotGroomingAlignedData(sess_mice, sess_options);
        results{i} = res;
        
        % Update figure title to include session information
        if ~isempty(findobj('type', 'figure'))
            fig = gcf;
            sgtitle({['Neural Activity Aligned to Grooming Onset: ', session_names{sess+1}], ...
                    ['State: ', sess_options.state, ', Source: ', sess_options.source]}, ...
                    'FontWeight', 'bold');
        end
    catch e
        warning('Error analyzing session %d: %s', sess, e.message);
    end
end

% Create a combined comparison figure if requested
if options.combined_figure && ~isempty(results) && any(~cellfun(@isempty, results))
    createCombinedFigure(results, options.sessions, session_names, options);
end

end

function createCombinedFigure(results, sessions, session_names, options)
% Create a combined figure comparing grooming responses across sessions

% Find valid results
valid_idx = find(~cellfun(@isempty, results));
if isempty(valid_idx)
    warning('No valid results to create combined figure');
    return;
end

% Create figure
figure('Position', [100, 100, 700, 800]);

% Extract data from each session
all_avg_dff = [];
all_sem_dff = [];
all_labels = {};
time_grid = [];

for i = 1:length(valid_idx)
    idx = valid_idx(i);
    res = results{idx};
    sess = sessions(idx);
    
    if isempty(res) || ~isfield(res, 'mouse_avg_dff')
        continue;
    end
    
    % Store time grid from first valid result
    if isempty(time_grid)
        time_grid = res.time_grid;
    end
    
    % Calculate average across mice
    avg_dff = nanmean(res.mouse_avg_dff, 1);
    sem_dff = nanstd(res.mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(res.mouse_avg_dff), 1));
    
    all_avg_dff = [all_avg_dff; avg_dff];
    all_sem_dff = [all_sem_dff; sem_dff];
    all_labels{end+1} = ['Session ', num2str(sess)];
end

% Plot averaged responses for each session (top subplot - full width)
subplot(2, 1, 1);
hold on;

% Define colors for sessions
session_colors = [
    0.2, 0.3, 0.8;  % Session 0: Blue
    0.8, 0.1, 0.1;  % Session 1: Dark red
    0.3, 0.7, 1.0;  % Session 2: Light blue
    1.0, 0.5, 0.3;  % Session 3: Light red
];

% Add grooming onset line
xline(0, '--k', 'LineWidth', 1.5);
% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% Plot each session
legend_handles = [];
for i = 1:size(all_avg_dff, 1)
    sess_idx = sessions(valid_idx(i)) + 1;
    color = session_colors(mod(sess_idx-1, size(session_colors, 1))+1, :);
    
    % Plot SEM as shaded area
    ciplot(all_avg_dff(i,:) - all_sem_dff(i,:), ...
           all_avg_dff(i,:) + all_sem_dff(i,:), ...
           time_grid, color, 0.2);
       
    % Plot mean line
    h = plot(time_grid, all_avg_dff(i,:), 'Color', color, 'LineWidth', 2);
    legend_handles = [legend_handles, h];
end

% Format plot
xlabel('Time from grooming onset (seconds)', 'FontSize', 12);
ylabel('ΔF/F (z-scored)', 'FontSize', 12);
title('Comparison of Neural Activity during Grooming across Sessions', 'FontSize', 14);
legend(legend_handles, all_labels, 'Location', 'best');
legend('boxoff');
grid on;
box off;

% Define constants for extraction
COL_TIME = 1;
COL_GROOM = 10;

% Get max grooming duration from options or use default
if ~isfield(options, 'max_grooming_duration')
    options.max_grooming_duration = 100000; % Default max 60 seconds
end

% Get min grooming duration from options or use default
if ~isfield(options, 'min_grooming_duration')
    options.min_grooming_duration = 0; % Default min 2 seconds
end

% Collect duration data and event counts for each session
all_session_durations = cell(1, length(valid_idx));
all_session_mice_ids = cell(1, length(valid_idx));
all_mouse_info = cell(1, length(valid_idx));
all_mouse_event_counts = cell(1, length(valid_idx));
all_unique_mice_with_events = cell(1, length(valid_idx));

for i = 1:length(valid_idx)
    idx = valid_idx(i);
    res = results{idx};
    
    if isempty(res) || ~isfield(res, 'mouse_info')
        continue;
    end
    
    % Store mouse info
    all_mouse_info{i} = res.mouse_info;
    
    % Initialize arrays for this session
    session_durations = [];
    session_mouse_ids = {};
    mouse_event_counts = zeros(size(res.mouse_info, 1), 1);
    mouse_ids_with_events = {}; % Track unique mice with grooming events
    
    % For each mouse, extract grooming durations from raw data
    for m = 1:size(res.mouse_info, 1)
        % Extract mouse data
        data = res.mouse_info{m, 4};
        mouse_id = res.mouse_info{m, 1};
        
        if isempty(data)
            continue;
        end
        
        % Extract time and grooming columns
        time = data(:, COL_TIME);
        grooming = data(:, COL_GROOM);
        
        % Find onsets and offsets of grooming events
        groom_changes = diff([0; grooming; 0]); % Add padding and compute differences
        groom_onsets = find(groom_changes == 1);
        groom_offsets = find(groom_changes == -1) - 1;
        
        % Count valid events for this mouse
        valid_event_count = 0;
        
        % Calculate grooming durations in seconds
        for g = 1:length(groom_onsets)
            if g <= length(groom_offsets)
                start_idx = groom_onsets(g);
                end_idx = groom_offsets(g);
                
                if start_idx <= length(time) && end_idx <= length(time)
                    duration = time(end_idx) - time(start_idx);
                    
                    % Only include durations between min and max range
                    if duration >= options.min_grooming_duration && duration <= options.max_grooming_duration
                        session_durations = [session_durations; duration];
                        session_mouse_ids = [session_mouse_ids; {mouse_id}];
                        valid_event_count = valid_event_count + 1;
                    end
                end
            end
        end
        
        % Store event count for this mouse
        mouse_event_counts(m) = valid_event_count;
        
        % Track mice with events
        if valid_event_count > 0
            mouse_ids_with_events{end+1} = mouse_id;
        end
    end
    
    % Store all durations and mouse-specific event counts
    all_session_durations{i} = session_durations;
    all_session_mice_ids{i} = session_mouse_ids;
    all_mouse_event_counts{i} = mouse_event_counts;
    all_unique_mice_with_events{i} = unique(mouse_ids_with_events);
end

% Calculate statistics for each session
mean_durations = zeros(1, length(valid_idx));
sem_durations = zeros(1, length(valid_idx));
event_counts = zeros(1, length(valid_idx));
mice_counts = zeros(1, length(valid_idx));
sess_labels = cell(1, length(valid_idx));

for i = 1:length(valid_idx)
    durations = all_session_durations{i};
    if ~isempty(durations)
        % Calculate statistics
        mean_durations(i) = mean(durations);
        sem_durations(i) = std(durations) / sqrt(length(durations));
        event_counts(i) = length(durations);
        
        % Count unique mice with grooming events
        mice_counts(i) = length(all_unique_mice_with_events{i});
        
        % Create session label
        sess_idx = sessions(valid_idx(i));
        sess_labels{i} = ['Session ', num2str(sess_idx)];
    end
end

% Second subplot: Grooming durations (bottom left)
subplot(2, 2, 3);
hold on;

% Plot individual data points with jitter
for i = 1:length(valid_idx)
    durations = all_session_durations{i};
    if ~isempty(durations)
        % Get session color
        sess_idx = sessions(valid_idx(i)) + 1;
        color = session_colors(mod(sess_idx-1, size(session_colors, 1))+1, :);
        
        % Create jitter
        jitter_width = 0.3;
        jitter = (rand(size(durations))-0.5) * jitter_width;
        
        % Plot points
        scatter(i + jitter, durations, 40, color, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

% Plot mean and error bars
errorbar(1:length(valid_idx), mean_durations, sem_durations, 'k.', 'LineWidth', 2, 'MarkerSize', 20);

% Plot mean durations as bar
bar_h = bar(1:length(valid_idx), mean_durations, 0.5, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);

% Add event counts as text
for i = 1:length(valid_idx)
    if event_counts(i) > 0
        text(i, max(all_session_durations{i}) + 2, ...
             sprintf('n=%d (%d)', mice_counts(i), event_counts(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
end

% Format plot
set(gca, 'XTick', 1:length(valid_idx));
set(gca, 'XTickLabel', sess_labels);
ylabel('Grooming Duration (seconds)', 'FontSize', 12);
title('Grooming Durations', 'FontSize', 14);
grid on;
box on;

% Set y-axis limits with padding
if ~isempty(all_session_durations) && any(~cellfun(@isempty, all_session_durations))
    % Filter out empty cells
    non_empty_cells = all_session_durations(~cellfun(@isempty, all_session_durations));
    
    % Calculate min and max, handling the case where all cells might be empty
    if ~isempty(non_empty_cells)
        y_min = min(cellfun(@min, non_empty_cells));
        y_max = max(cellfun(@max, non_empty_cells));
        
        % Set limits with padding
        ylim([max(0, y_min-5), y_max+5]);
    end
end

% Third subplot: Number of grooming events per mouse (bottom right)
subplot(2, 2, 4);
hold on;

% Find max event count for y-axis scaling
max_events = 0;
for i = 1:length(valid_idx)
    if ~isempty(all_mouse_event_counts{i})
        max_events = max(max_events, max(all_mouse_event_counts{i}));
    end
end

% Plot bar chart of mean events per mouse
mean_events_per_mouse = zeros(1, length(valid_idx));
sem_events_per_mouse = zeros(1, length(valid_idx));

for i = 1:length(valid_idx)
    if ~isempty(all_mouse_event_counts{i})
        counts = all_mouse_event_counts{i};
        counts = counts(counts > 0); % Only include mice with events
        if ~isempty(counts)
            mean_events_per_mouse(i) = mean(counts);
            sem_events_per_mouse(i) = std(counts) / sqrt(length(counts));
        end
    end
end

% Draw bars for mean events per mouse
bar_h2 = bar(1:length(valid_idx), mean_events_per_mouse, 0.5, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);

% Add error bars
errorbar(1:length(valid_idx), mean_events_per_mouse, sem_events_per_mouse, 'k.', 'LineWidth', 2, 'MarkerSize', 20);

% Plot individual data points with jitter
for i = 1:length(valid_idx)
    if ~isempty(all_mouse_event_counts{i})
        counts = all_mouse_event_counts{i};
        valid_counts = counts(counts > 0); % Only include mice with events
        
        if ~isempty(valid_counts)
            % Get session color
            sess_idx = sessions(valid_idx(i)) + 1;
            color = session_colors(mod(sess_idx-1, size(session_colors, 1))+1, :);
            
            % Create jitter
            jitter_width = 0.3;
            jitter = (rand(size(valid_counts))-0.5) * jitter_width;
            
            % Plot points
            scatter(i + jitter, valid_counts, 60, color, 'filled', 'MarkerFaceAlpha', 0.6);
        end
    end
end

% Format plot
set(gca, 'XTick', 1:length(valid_idx));
set(gca, 'XTickLabel', sess_labels);
ylabel('Grooming Events', 'FontSize', 12);
title('Grooming Events', 'FontSize', 14);
ylim([0, max_events*1.2]); % Add some padding to y-axis
grid on;
box on;

% Add number of mice as text
for i = 1:length(valid_idx)
    if ~isempty(all_mouse_event_counts{i}) && ~isempty(all_unique_mice_with_events{i})
        mice_with_events = length(all_unique_mice_with_events{i});
        if mice_with_events > 0
        end
    end
end

% Add overall title with state and source information
if isfield(options, 'state') && isfield(options, 'source')
    sgtitle({'Grooming-Related Activity across Sessions', ...
             ['State: ', upper(options.state(1)), options.state(2:end), ...
              ', Source: ', upper(options.source(1)), options.source(2:end)]}, ...
             'FontWeight', 'bold');
else
    sgtitle('Grooming-Related Activity across Sessions', 'FontWeight', 'bold');
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

% Wrapper function for plotGroomingAlignedData
function varargout = plotGroomingAlignedData(mice_all, options)
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
    warning('No mice found with the specified criteria');
    if nargout > 0
        varargout{1} = [];
    end
    return;
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
grid on;
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
    aligned_data.avg_dff = avg_dff;
    aligned_data.sem_dff = sem_dff;
    varargout{1} = aligned_data;
end
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