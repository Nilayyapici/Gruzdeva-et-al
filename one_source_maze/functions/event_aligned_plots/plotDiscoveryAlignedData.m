function varargout = plotDiscoveryAlignedData(mice_all, options)
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
%     - distance_ylim: y-axis limits for distance plot [min max] (default: [0 150])
%     - bar_window: n seconds before/after discovery averaged for the
%                   pre vs post bar plot and paired t-test (default: 30)

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

if ~isfield(options, 'bar_window')
    options.bar_window = 30; % seconds before/after discovery for bar plot
end

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

% Create figure with three subplots
figure('Position', [100, 100, 1500, 550]);

% 1. Average plot (left subplot)
subplot(1, 3, 1);
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

% Mark the pre/post windows used for the bar plot
xline(-options.bar_window, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xline( options.bar_window, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

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

% 2. Raster plot (middle subplot)
subplot(1, 3, 2);

% Plot data based on selected type
if options.plot_distance
    imagesc(time_grid, 1:num_mice, plot_data);
    colormap(gca, bluewhitered(256));
    
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

%% 3. Pre vs post discovery bar plot with paired t-test (right subplot)
n_sec = options.bar_window;
pre_mask  = time_grid >= -n_sec & time_grid < 0;
post_mask = time_grid >  0      & time_grid <= n_sec;

% Per-mouse mean in each window
pre_vals  = mean(plot_data(:, pre_mask),  2, 'omitnan');
post_vals = mean(plot_data(:, post_mask), 2, 'omitnan');

% Keep only mice with valid values in both windows (required for pairing)
valid   = ~isnan(pre_vals) & ~isnan(post_vals);
pre_v   = pre_vals(valid);
post_v  = post_vals(valid);
n_valid = sum(valid);

if n_valid < num_mice
    warning('%d mouse/mice excluded from the t-test (missing data in pre or post window).', ...
        num_mice - n_valid);
end

% Paired t-test (requires Statistics and Machine Learning Toolbox)
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

% Bars
b = bar(1:2, bar_means, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
b.CData = [0.65 0.65 0.65; plot_color];
b.FaceAlpha = 0.7;

% Individual mice: paired lines + dots
if n_valid > 0
    plot([1 2], [pre_v, post_v]', '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.75);
    scatter(ones(n_valid, 1),     pre_v,  30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2 * ones(n_valid, 1), post_v, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% SEM error bars
errorbar(1:2, bar_means, bar_sems, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 12);

% Significance bracket
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

% Format plot
xlim([0.4 2.6]);
ylim([y_bot - 0.05 * y_range, y_brk + 0.18 * y_range]);
set(gca, 'XTick', 1:2, 'XTickLabel', ...
    {sprintf('Pre (-%g to 0 s)', n_sec), sprintf('Post (0 to %g s)', n_sec)});
ylabel(['Mean ' data_label], 'FontSize', 14);
title({'Pre vs post discovery (paired t-test)', ...
       sprintf('t(%d) = %.2f, p = %.4g, n = %d', t_stats.df, t_stats.tstat, p_val, n_valid)});
box off;

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
mouse_ids_sorted = cell(num_mice, 1);
for i = 1:num_mice
    mouse_ids_sorted{i} = filtered_mice{i, 1};
end

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

% ── Sheet 3: Pre/post means per mouse ────────────────────────────────────
pre_post_table = table(mouse_ids_sorted, pre_vals, post_vals, post_vals - pre_vals, valid, ...
    'VariableNames', {'Mouse_ID', 'Pre_mean', 'Post_mean', 'Post_minus_Pre', 'Included_in_test'});
writetable(pre_post_table, excel_filename, 'Sheet', 'pre_post_per_mouse');

% ── Sheet 4: Paired t-test summary ───────────────────────────────────────
stats_table = table(n_sec, n_valid, bar_means(1), bar_sems(1), bar_means(2), bar_sems(2), ...
    t_stats.tstat, t_stats.df, p_val, ...
    'VariableNames', {'Window_s', 'N_mice', 'Pre_mean', 'Pre_SEM', ...
                      'Post_mean', 'Post_SEM', 't_stat', 'df', 'p_value'});
writetable(stats_table, excel_filename, 'Sheet', 'paired_ttest');

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
fprintf('Paired t-test (%g s pre vs %g s post): pre = %.3f ± %.3f, post = %.3f ± %.3f (mean ± SEM)\n', ...
    n_sec, n_sec, bar_means(1), bar_sems(1), bar_means(2), bar_sems(2));
fprintf('  t(%d) = %.3f, p = %.4g, n = %d mice %s\n', ...
    t_stats.df, t_stats.tstat, p_val, n_valid, sigStars(p_val));

% Optional: return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.dff = aligned_dff;
    aligned_data.time = time_grid;
    aligned_data.dist = aligned_dist;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.pre_post.window_s = n_sec;
    aligned_data.pre_post.pre = pre_vals;
    aligned_data.pre_post.post = post_vals;
    aligned_data.pre_post.valid = valid;
    aligned_data.pre_post.p_value = p_val;
    aligned_data.pre_post.t_stats = t_stats;
    varargout{1} = aligned_data;
end
end

function s = sigStars(p)
% SIGSTARS returns a significance label for a p-value
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
% ciplot creates a shaded area to visualize confidence intervals

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

if nargin < 1
   m = 100;
end

colors = [
    0.0, 0.0, 0.5;   % Dark blue
    0.1, 0.3, 0.7;   % Softer blue
    0.4, 0.6, 0.9;   % Medium blue
    0.7, 0.85, 1.0;  % Very light blue
    1.0, 1.0, 1.0;   % White
    1.0, 0.85, 0.85; % Very light red
    0.9, 0.4, 0.4;   % Medium red
    0.8, 0.1, 0.1    % True red
];

n = size(colors, 1);
cmap = zeros(m, 3);

positions = [0, 0.15, 0.3, 0.45, 0.5, 0.55, 0.7, 1.0];
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