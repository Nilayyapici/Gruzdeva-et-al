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
%     - min_events: minimum number of valid events a mouse needs to be
%                   included (default: 1). A mouse is included if it has
%                   AT LEAST min_events events (>=). Applies to ALL panels
%                   (average, raster, bar plot) and all Excel sheets.
%     - bar_window: n seconds before/after event onset averaged for the
%                   pre vs post bar plot and paired t-test (default: 10)

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
    options.z_score = false;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_event_duration')
    options.min_event_duration = 1; % Default 1 second
end

if ~isfield(options, 'min_events')
    options.min_events = 1; % mouse needs >= this many events to be included
end

if ~isfield(options, 'bar_window')
    options.bar_window = 10; % seconds before/after onset for bar plot
end

if isfield(options, 'bar_mouse_set')
    warning(['options.bar_mouse_set is no longer used: all panels now use ' ...
             'the same mice, set by options.min_events.']);
end

MIN_EVENTS = options.min_events;
min_events_label = sprintf('\\geq%d events', MIN_EVENTS);   % for plot titles (TeX)
min_events_text  = sprintf('>=%d events', MIN_EVENTS);      % for console / Excel

% Make sure the bar window fits inside the time window
max_pre  = -options.time_window(1);
max_post =  options.time_window(2);
if options.bar_window > min(max_pre, max_post)
    warning('bar_window (%g s) exceeds time_window; clipping to %g s.', ...
        options.bar_window, min(max_pre, max_post));
    options.bar_window = min(max_pre, max_post);
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
    event_changes = diff([0; event_signal; 0]);
    event_onsets = find(event_changes == 1);
    event_offsets = find(event_changes == -1) - 1;

    % Calculate durations in seconds
    event_durations = time(event_offsets) - time(event_onsets);

    % Filter for minimum duration and events after discovery
    valid_events = event_durations >= options.min_event_duration & ...
                   time(event_onsets) > discovery_time_sec;

    event_onsets = event_onsets(valid_events);

    num_events = length(event_onsets);
    all_mouse_event_counts(i) = num_events;

    if num_events == 0
        fprintf('No valid %s events found for mouse %d after food discovery\n', lower(event_name), i);
        continue;
    end

    all_events_aligned_dff = cell(num_events, 1);
    all_events_aligned_time = cell(num_events, 1);

    for j = 1:num_events
        onset_time = time(event_onsets(j));

        time_mask = time >= (onset_time + options.time_window(1)) & ...
                    time <= (onset_time + options.time_window(2));

        aligned_time = time(time_mask) - onset_time;
        aligned_dff = dff(time_mask);

        if options.z_score
            aligned_dff = (aligned_dff - nanmean(aligned_dff)) / nanstd(aligned_dff);
        end

        if options.smooth_window > 0
            aligned_dff = movmean(aligned_dff, options.smooth_window, 'omitnan');
        end

        all_events_aligned_dff{j} = aligned_dff;
        all_events_aligned_time{j} = aligned_time;
    end

    all_mouse_aligned_dff{i} = struct('events_dff', {all_events_aligned_dff}, ...
                                      'events_time', {all_events_aligned_time});
end

% Determine the common time grid for plotting
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000);

% Interpolate all events onto common time grid for each mouse
mouse_avg_dff = nan(num_mice, length(time_grid));

for i = 1:num_mice
    if all_mouse_event_counts(i) == 0
        continue;
    end

    events_interp = nan(all_mouse_event_counts(i), length(time_grid));

    for j = 1:all_mouse_event_counts(i)
        event_time = all_mouse_aligned_dff{i}.events_time{j};
        event_dff = all_mouse_aligned_dff{i}.events_dff{j};

        if numel(event_time) < 2
            continue;
        end

        events_interp(j, :) = interp1(event_time, event_dff, time_grid, 'linear', NaN);
    end

    mouse_avg_dff(i, :) = nanmean(events_interp, 1);
end

% ── Apply the min_events filter ONCE: every panel and every Excel sheet
%    uses exactly this set of mice ──────────────────────────────────────
included = all_mouse_event_counts >= MIN_EVENTS;
inc_avg_dff     = mouse_avg_dff(included, :);
inc_mice_info   = filtered_mice(included, :);
inc_event_counts = all_mouse_event_counts(included);
num_inc         = sum(included);

n_excluded = sum(~included);
if n_excluded > 0
    fprintf('%d mouse/mice excluded (fewer than %d valid %s events).\n', ...
        n_excluded, MIN_EVENTS, lower(event_name));
end

if num_inc == 0
    error('No mice have %s %s events. Lower options.min_events.', ...
        min_events_text, lower(event_name));
end

% Short mouse labels (used for raster y-ticks and Excel)
inc_labels = cell(num_inc, 1);
for i = 1:num_inc
    inc_labels{i} = [shortId(inc_mice_info{i, 1}, i), ' (', num2str(inc_event_counts(i)), ')'];
end

% Create figure with three subplots
figure('Position', [100, 100, 1500, 550]);

%% 1. Average DFF plot (left subplot)
subplot(1, 3, 1);
hold on;

avg_dff = nanmean(inc_avg_dff, 1);
sem_dff = nanstd(inc_avg_dff, 0, 1) ./ sqrt(sum(~isnan(inc_avg_dff), 1));

xline(0, '--k', 'LineWidth', 1.5);
xline(-options.bar_window, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xline( options.bar_window, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

if strcmpi(options.event_type, 'visit')
    event_color = [0.0, 0.4, 0.8];
    sem_color = [0.6, 0.8, 1.0];
else
    event_color = [0.8, 0.2, 0.2];
    sem_color = [1.0, 0.7, 0.7];
end

ciplot(avg_dff - sem_dff, avg_dff + sem_dff, time_grid, sem_color, 0.5);
plot(time_grid, avg_dff, 'Color', event_color, 'LineWidth', 2);

xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
if options.z_score
    dff_label = 'ΔF/F (z-scored)';
else
    dff_label = 'ΔF/F';
end
ylabel(dff_label, 'FontSize', 12);
title({sprintf('ΔF/F aligned to %s onset', lower(event_name)), ...
       sprintf('n = %d mice (%s)', num_inc, min_events_label)});
grid off;
box off;

xlim(options.time_window);
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

%% 2. Raster plot (middle subplot)
subplot(1, 3, 2);

imagesc(time_grid, 1:num_inc, inc_avg_dff);
hold on;
xline(0, '--k', 'LineWidth', 1.5);

xlabel(sprintf('Time from %s onset (seconds)', event_name), 'FontSize', 12);
ylabel('Mouse #', 'FontSize', 12);
title({'Activity across mice', ...
       sprintf('n = %d mice, %d events (%s)', num_inc, sum(inc_event_counts), min_events_label)});

if strcmpi(options.event_type, 'visit')
    colormap(gca, bluewhitered(256, true));
else
    colormap(gca, bluewhitered(256, false));
end
h = colorbar;
ylabel(h, dff_label, 'FontSize', 12);
caxis(options.caxis_range);

if num_inc <= 30
    set(gca, 'YTick', 1:num_inc, 'YTickLabel', inc_labels, 'TickLabelInterpreter', 'none');
else
    set(gca, 'YTick', unique(round(linspace(1, num_inc, 10))));
end
xlim(options.time_window);

%% 3. Pre vs post onset bar plot with paired t-test (right subplot)
n_sec = options.bar_window;
pre_mask  = time_grid >= -n_sec & time_grid < 0;
post_mask = time_grid >  0      & time_grid <= n_sec;

% Per-mouse means for ALL state/source-filtered mice (for the Excel sheet)
pre_vals  = mean(mouse_avg_dff(:, pre_mask),  2, 'omitnan');
post_vals = mean(mouse_avg_dff(:, post_mask), 2, 'omitnan');

% Mice in the test = included AND have data in both windows
valid   = included & ~isnan(pre_vals) & ~isnan(post_vals);
pre_v   = pre_vals(valid);
post_v  = post_vals(valid);
n_valid = sum(valid);

if n_valid < num_inc
    warning('%d included mouse/mice dropped from the t-test (missing data in pre or post window).', ...
        num_inc - n_valid);
end

if n_valid >= 2
    [~, p_val, ~, t_stats] = ttest(post_v, pre_v);
else
    warning('Fewer than 2 mice with valid data; paired t-test not computed.');
    p_val = NaN;
    t_stats = struct('tstat', NaN, 'df', NaN, 'sd', NaN);
end

bar_means = [mean(pre_v), mean(post_v)];
bar_sems  = [std(pre_v), std(post_v)] ./ sqrt(n_valid);

subplot(1, 3, 3);
hold on;

b = bar(1:2, bar_means, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
b.CData = [0.65 0.65 0.65; event_color];
b.FaceAlpha = 0.7;

if n_valid > 0
    plot([1 2], [pre_v, post_v]', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75);
    scatter(ones(n_valid, 1),     pre_v,  30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2 * ones(n_valid, 1), post_v, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

errorbar(1:2, bar_means, bar_sems, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 12);

all_vals = [pre_v; post_v; bar_means(:) + bar_sems(:); bar_means(:) - bar_sems(:)];
y_top = max(all_vals);
y_bot = min([all_vals; 0]);
y_range = y_top - y_bot;
if isempty(y_range) || isnan(y_range) || y_range == 0
    y_range = 1;
end
y_brk = y_top + 0.08 * y_range;
plot([1 1 2 2], [y_brk - 0.03*y_range, y_brk, y_brk, y_brk - 0.03*y_range], ...
    'k', 'LineWidth', 1.2);
text(1.5, y_brk + 0.04 * y_range, sigStars(p_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

xlim([0.4 2.6]);
ylim([y_bot - 0.05 * y_range, y_brk + 0.18 * y_range]);
set(gca, 'XTick', 1:2, 'XTickLabel', ...
    {sprintf('Pre (-%g to 0 s)', n_sec), sprintf('Post (0 to %g s)', n_sec)});
ylabel(['Mean ' dff_label], 'FontSize', 12);
title({sprintf('Pre vs post %s onset (paired t-test)', lower(event_name)), ...
       sprintf('t(%d) = %.2f, p = %.4g, n = %d mice (%s)', ...
               t_stats.df, t_stats.tstat, p_val, n_valid, min_events_label)});
box off;

if options.z_score
    z_score_text = 'Z-scored';
else
    z_score_text = 'Raw';
end
sgtitle([z_score_text, ' ΔF/F values (post-discovery only) ', ...
         options.state, ' mice, ', options.source], 'FontWeight', 'bold');

%% Export data to Excel
% All sheets use the same included mice (>= min_events), except
% pre_post_per_mouse, which lists every mouse with an 'Included' column.
excel_filename = sprintf('event_aligned_%s_%s_%s_min%devents_%s.xlsx', ...
    options.event_type, options.state, options.source, MIN_EVENTS, datestr(now, 'yyyymmdd'));

% Sheet 1: Average trace
avg_table = table(time_grid(:), avg_dff(:), sem_dff(:), ...
    'VariableNames', {'Time_s', 'Mean_dFF', 'SEM_dFF'});
writetable(avg_table, excel_filename, 'Sheet', 'average_trace');

% Sheet 2: Per-mouse average (included mice, same as raster)
heatmap_cell = cell(num_inc + 1, length(time_grid) + 1);
heatmap_cell{1, 1} = 'Mouse_ID';
for t = 1:length(time_grid)
    heatmap_cell{1, t+1} = time_grid(t);
end
for i = 1:num_inc
    heatmap_cell{i+1, 1} = inc_labels{i};
    for t = 1:length(time_grid)
        heatmap_cell{i+1, t+1} = inc_avg_dff(i, t);
    end
end
writecell(heatmap_cell, excel_filename, 'Sheet', 'heatmap_per_mouse');

% Sheet 3: Pre/post means per mouse (all mice, with inclusion flags)
mouse_ids_all = cell(num_mice, 1);
for i = 1:num_mice
    mouse_ids_all{i} = shortId(filtered_mice{i, 1}, i);
end
pre_post_table = table(mouse_ids_all, all_mouse_event_counts, pre_vals, post_vals, ...
    post_vals - pre_vals, included, valid, ...
    'VariableNames', {'Mouse_ID', 'N_events', 'Pre_mean', 'Post_mean', ...
                      'Post_minus_Pre', 'Meets_min_events', 'Included_in_test'});
writetable(pre_post_table, excel_filename, 'Sheet', 'pre_post_per_mouse');

% Sheet 4: Paired t-test summary
stats_table = table(n_sec, MIN_EVENTS, n_valid, ...
    bar_means(1), bar_sems(1), bar_means(2), bar_sems(2), ...
    t_stats.tstat, t_stats.df, p_val, ...
    'VariableNames', {'Window_s', 'Min_events', 'N_mice', 'Pre_mean', 'Pre_SEM', ...
                      'Post_mean', 'Post_SEM', 't_stat', 'df', 'p_value'});
writetable(stats_table, excel_filename, 'Sheet', 'paired_ttest');

fprintf('Event-aligned data exported to: %s\n', excel_filename);

% Console summary
fprintf('Analysis complete: %d mice with events, %d mice with %s (used in all panels), %d %s events in those mice\n', ...
        sum(all_mouse_event_counts > 0), num_inc, min_events_text, sum(inc_event_counts), lower(event_name));
fprintf('Average pre-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid < 0)));
fprintf('Average during-%s DFF: %.3f\n', lower(event_name), nanmean(avg_dff(time_grid >= 0)));
fprintf('Paired t-test (%g s pre vs %g s post): pre = %.3f ± %.3f, post = %.3f ± %.3f (mean ± SEM)\n', ...
    n_sec, n_sec, bar_means(1), bar_sems(1), bar_means(2), bar_sems(2));
fprintf('  t(%d) = %.3f, p = %.4g, n = %d mice %s\n', ...
    t_stats.df, t_stats.tstat, p_val, n_valid, sigStars(p_val));

% Return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.mouse_avg_dff = mouse_avg_dff;         % all state/source-filtered mice
    aligned_data.included = included;                   % min_events mask
    aligned_data.included_mouse_avg_dff = inc_avg_dff;  % mice used in all panels
    aligned_data.time_grid = time_grid;
    aligned_data.event_counts = all_mouse_event_counts;
    aligned_data.included_event_counts = inc_event_counts;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.included_mouse_info = inc_mice_info;
    aligned_data.all_events = all_mouse_aligned_dff;
    aligned_data.min_events = MIN_EVENTS;
    aligned_data.pre_post.window_s = n_sec;
    aligned_data.pre_post.pre = pre_vals;
    aligned_data.pre_post.post = post_vals;
    aligned_data.pre_post.valid = valid;
    aligned_data.pre_post.p_value = p_val;
    aligned_data.pre_post.t_stats = t_stats;
end
end

function id = shortId(raw_id, fallback_idx)
% Returns the part of the mouse ID before the first underscore
if ischar(raw_id) || isstring(raw_id)
    parts = strsplit(char(raw_id), '_');
    id = parts{1};
else
    id = num2str(fallback_idx);
end
end

function s = sigStars(p)
if isnan(p)
    s = 'n/a';
elseif p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end

function ciplot(lower, upper, x, color, alpha)
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
% BLUEWHITERED creates a blue-to-white-to-red colormap
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
end
positions = [0, 0.15, 0.3, 0.45, 0.55, 0.7, 0.85, 1.0];

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