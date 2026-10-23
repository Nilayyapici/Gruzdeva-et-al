function results = plot_dff_distance_correlation_clean(mice, condition, options)
% PLOT_DFF_DISTANCE_CORRELATION_CLEAN - dF/F vs distance-to-food correlations
%
% Key corrections vs previous version:
%   - Distance is now Euclidean distance to food well (from x,y + mice{i,2})
%     rather than column 9 (cumulative distance traveled)
%   - Speed threshold removed; grooming column used exclusively to filter
%     non-exploratory periods
%   - Paired t-test now correctly uses only mice valid in BOTH sessions
%   - Means/SEMs computed from the same paired subset used for statistics
%
% Usage:
%   results = plot_dff_distance_correlation_clean(mice, 'saline')
%   results = plot_dff_distance_correlation_clean(mice, 'CNO')
%   results = plot_dff_distance_correlation_clean(mice, 'strong')
%   results = plot_dff_distance_correlation_clean(mice, 'weak')
%
%   options.time_limit       = 10;    % First N minutes only (optional)
%   options.exclude_grooming = true;  % Exclude grooming frames (default: true)
%   options.min_points       = 50;    % Min frames for valid correlation (default: 50)
%   options.figure_size      = [100, 100, 600, 500];
%
% Data columns assumed:
%   1-time, 2-x, 3-y, 4-465nm, 5-405nm, 6-dff, 7-speed,
%   8-zones, 9-distance traveled, 10-grooming
% mice columns:
%   1-name, 2-[x,y] food position, 3-pre data, 4-test data,
%   5-start position, 6-latencies, 7-memory classification (optional)

%% Input validation and defaults
if nargin < 2
    error('Both mice data and condition are required');
end
if nargin < 3
    options = struct();
end

if ~isfield(options, 'exclude_grooming'), options.exclude_grooming = true;  end
if ~isfield(options, 'min_points'),       options.min_points       = 50;    end
if ~isfield(options, 'figure_size') || isempty(options.figure_size) || length(options.figure_size) ~= 4
    options.figure_size = [100, 100, 600, 500];
end

use_time_limit = isfield(options, 'time_limit') && ~isempty(options.time_limit);
if use_time_limit
    time_limit_min = options.time_limit;
    time_limit_sec = time_limit_min * 60;
    fprintf('Time limit: %.1f min\n', time_limit_min);
else
    fprintf('No time limit — full sessions\n');
end

% Speed threshold: use a low value (e.g. 2-5 cm/s) to exclude stops/sniffing.
% Grooming column handles grooming separately; these filters are complementary.
use_speed_threshold = isfield(options, 'speed_threshold') && ~isempty(options.speed_threshold);
if use_speed_threshold
    speed_threshold = options.speed_threshold;
    fprintf('Speed threshold: >= %.1f cm/s (excludes stops/sniffing)\n', speed_threshold);
else
    fprintf('No speed threshold\n');
end

%% Select mice by condition
mouse_names = mice(:, 1);

has_memory = size(mice, 2) >= 7 && ~all(cellfun(@isempty, mice(:, 7)));
if has_memory
    memory_class = mice(:, 7);
end

switch lower(condition)
    case 'saline'
        selected_idx = contains(mouse_names, 'saline');
        condition_name = 'Saline';
    case 'cno'
        selected_idx = contains(mouse_names, 'CNO');
        condition_name = 'CNO';
    case 'strong'
        if ~has_memory
            error('Memory classifications not found in mice column 7.');
        end
        selected_idx = strcmp(memory_class, 'strong memory');
        condition_name = 'Strong Memory';
    case 'weak'
        if ~has_memory
            error('Memory classifications not found in mice column 7.');
        end
        selected_idx = strcmp(memory_class, 'weak memory');
        condition_name = 'Weak Memory';
    otherwise
        error('Invalid condition. Use: saline, CNO, strong, or weak');
end

selected_mice = find(selected_idx);
n_mice        = length(selected_mice);

if n_mice == 0
    error('No mice found for condition: %s', condition);
end
fprintf('Found %d mice for condition: %s\n', n_mice, condition_name);

%% Calculate correlations
correlations_pre  = NaN(n_mice, 1);
correlations_test = NaN(n_mice, 1);
p_values_pre      = NaN(n_mice, 1);
p_values_test     = NaN(n_mice, 1);
valid_pre         = false(n_mice, 1);
valid_test        = false(n_mice, 1);

for i = 1:n_mice
    mouse_idx  = selected_mice(i);
    mouse_name = mice{mouse_idx, 1};

    % Food well position (pixels) — used to compute distance-to-food
    food_pos = mice{mouse_idx, 2};  % [x, y]

    pre_data  = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};

    fprintf('\nMouse %d/%d: %s\n', i, n_mice, mouse_name);

    % --- Apply time limit ---
    if use_time_limit
        pre_data  = apply_time_limit(pre_data,  time_limit_sec);
        test_data = apply_time_limit(test_data, time_limit_sec);
    end

    % --- Apply speed threshold (excludes stops/sniffing) ---
    if use_speed_threshold
        pre_data  = pre_data(pre_data(:,7)   >= speed_threshold, :);
        test_data = test_data(test_data(:,7) >= speed_threshold, :);
        fprintf('  After speed filter (>= %.1f cm/s): Pre=%d pts, Test=%d pts\n', ...
                speed_threshold, size(pre_data,1), size(test_data,1));
    end

    % --- Compute Euclidean distance to food well ---
    % Replaces column 9 (cumulative distance traveled) with true proximity
    pre_data  = add_distance_to_food(pre_data,  food_pos);
    test_data = add_distance_to_food(test_data, food_pos);

    % --- Compute correlations ---
    [correlations_pre(i),  p_values_pre(i),  valid_pre(i)]  = ...
        compute_correlation(pre_data,  options.exclude_grooming, options.min_points);
    [correlations_test(i), p_values_test(i), valid_test(i)] = ...
        compute_correlation(test_data, options.exclude_grooming, options.min_points);

    if valid_pre(i)
        fprintf('  Pre:  r = %.3f, p = %.3f\n', correlations_pre(i),  p_values_pre(i));
    else
        fprintf('  Pre:  insufficient data\n');
    end
    if valid_test(i)
        fprintf('  Test: r = %.3f, p = %.3f\n', correlations_test(i), p_values_test(i));
    else
        fprintf('  Test: insufficient data\n');
    end
end

%% Paired subset — only mice valid in BOTH sessions
valid_both = valid_pre & valid_test;
n_paired   = sum(valid_both);

if n_paired == 0
    warning('No mice have valid correlations in both sessions. Cannot plot or test.');
    results = struct();
    return;
end

pre_paired  = correlations_pre(valid_both);
test_paired = correlations_test(valid_both);

% Means and SEMs from the paired subset (consistent with the statistical test)
pre_mean  = mean(pre_paired);
test_mean = mean(test_paired);
pre_sem   = std(pre_paired)  / sqrt(n_paired);
test_sem  = std(test_paired) / sqrt(n_paired);

% Paired t-test
if n_paired > 1
    [~, p_val] = ttest(pre_paired, test_paired);
else
    p_val = NaN;
    warning('Only 1 paired mouse — cannot run t-test');
end

%% Figure
try
    fig = figure('Name', sprintf('dF/F-Distance Correlations: %s', condition_name), ...
                 'Position', options.figure_size, 'Color', 'white');
catch
    fig = figure('Name', sprintf('dF/F-Distance Correlations: %s', condition_name), ...
                 'Color', 'white');
end

means = [pre_mean, test_mean];
sems  = [pre_sem,  test_sem];

% Bars
bar(1, means(1), 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.7, ...
    'EdgeColor', 'none', 'BarWidth', 0.6);
hold on;
bar(2, means(2), 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.7, ...
    'EdgeColor', 'none', 'BarWidth', 0.6);

% Error bars
errorbar(1:2, means, sems, 'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');

% Individual paired data points + connecting lines
for i = 1:n_mice
    if valid_both(i)
        plot([1, 2], [correlations_pre(i), correlations_test(i)], ...
             'o-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
             'MarkerSize', 4, 'MarkerFaceColor', [0.3 0.3 0.3], ...
             'MarkerEdgeColor', 'none');
    end
end

% Significance annotation
all_vals = [pre_paired; test_paired];
y_range  = max(all_vals) - min(all_vals);
if y_range == 0, y_range = 0.1; end
max_y    = max([max(all_vals), means + sems]) + y_range * 0.15;

if ~isnan(p_val)
    if p_val < 0.001
        sig_text = '***';
    elseif p_val < 0.01
        sig_text = '**';
    elseif p_val < 0.05
        sig_text = '*';
    else
        sig_text = 'ns';
    end
    plot([1, 2], [max_y, max_y], 'k-', 'LineWidth', 1);
    plot([1, 1], [max_y - y_range*0.02, max_y], 'k-', 'LineWidth', 1);
    plot([2, 2], [max_y - y_range*0.02, max_y], 'k-', 'LineWidth', 1);
    text(1.5, max_y + y_range*0.04, sig_text, ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    text(1.5, max_y + y_range*0.12, sprintf('p = %.3f', p_val), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.3 0.3 0.3]);
end

% Axes formatting
ylim([min(all_vals) - y_range*0.15, max_y + y_range*0.25]);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Before', 'Test'}, ...
    'FontSize', 12, 'TickDir', 'out', 'LineWidth', 1.5);
ylabel('dF/F – Distance to Food Correlation (r)', 'FontSize', 13, 'FontWeight', 'bold');
box off; grid off;

% Title
title_str = sprintf('%s Group — dF/F vs Distance to Food', condition_name);
filter_parts = {};
if use_time_limit,       filter_parts{end+1} = sprintf('First %.1f min', time_limit_min);      end
if use_speed_threshold,  filter_parts{end+1} = sprintf('speed >= %.1f cm/s', speed_threshold); end
if ~isempty(filter_parts)
    title_str = [title_str, sprintf(' (%s, n=%d)', strjoin(filter_parts, ', '), n_paired)];
else
    title_str = [title_str, sprintf(' (Full session, n=%d)', n_paired)];
end
title(title_str, 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

% Sample size labels at bottom of plot
y_bot = min(ylim) + (max(ylim) - min(ylim)) * 0.03;
text(1, y_bot, sprintf('n=%d', n_paired), 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', [0.5 0.5 0.5]);
text(2, y_bot, sprintf('n=%d', n_paired), 'HorizontalAlignment', 'center', ...
     'FontSize', 10, 'Color', [0.5 0.5 0.5]);

hold off;

%% Results struct
results = struct();
results.condition          = condition_name;
results.correlations_pre   = correlations_pre;
results.correlations_test  = correlations_test;
results.p_values_pre       = p_values_pre;
results.p_values_test      = p_values_test;
results.valid_pre          = valid_pre;
results.valid_test         = valid_test;
results.valid_both         = valid_both;
results.mouse_names        = mouse_names(selected_mice);
results.n_mice             = n_mice;
results.n_paired           = n_paired;
results.pre_mean           = pre_mean;
results.test_mean          = test_mean;
results.pre_sem            = pre_sem;
results.test_sem           = test_sem;
results.ttest_p            = p_val;
results.significant        = ~isnan(p_val) && p_val < 0.05;
results.time_limit_used       = use_time_limit;
if use_time_limit,   results.time_limit_min  = time_limit_min;  end
results.speed_threshold_used  = use_speed_threshold;
if use_speed_threshold, results.speed_threshold = speed_threshold; end

%% Save bar plot data to xlsx
try
    % Build table with paired data
    mouse_names_paired = mouse_names(selected_mice(valid_both));
    
    % Per-mouse paired data
    T_mice = table(mouse_names_paired, pre_paired, test_paired, ...
        'VariableNames', {'Mouse', 'Correlation_Pre', 'Correlation_Test'});
    
    % Summary stats
    T_summary = table({'Before'; 'Test'; 'Paired_ttest_p'}, ...
        [pre_mean; test_mean; p_val], ...
        [pre_sem; test_sem; NaN], ...
        [n_paired; n_paired; NaN], ...
        'VariableNames', {'Condition', 'Mean_r', 'SEM_r', 'N'});
    
    % Filename
    xlsx_filename = sprintf('dff_distance_correlation_%s.xlsx', lower(condition_name));
    xlsx_filename = strrep(xlsx_filename, ' ', '_');
    
    writetable(T_mice,    xlsx_filename, 'Sheet', 'Per_Mouse_Data');
    writetable(T_summary, xlsx_filename, 'Sheet', 'Summary_Stats');
    
    fprintf('\nBar plot data saved to: %s\n', xlsx_filename);
catch ME
    warning('Could not save xlsx: %s', ME.message);
end

%% Console summary
fprintf('\n=== %s GROUP SUMMARY ===\n', upper(condition_name));
if use_time_limit,      fprintf('Time limit:       %.1f min\n',    time_limit_min);  end
if use_speed_threshold, fprintf('Speed threshold:  >= %.1f cm/s\n', speed_threshold); end
fprintf('Mice with valid pre correlation:  %d / %d\n', sum(valid_pre),  n_mice);
fprintf('Mice with valid test correlation: %d / %d\n', sum(valid_test), n_mice);
fprintf('Paired (valid in both):           %d / %d\n', n_paired, n_mice);
fprintf('Before mean ± SEM: %.3f ± %.3f\n', pre_mean,  pre_sem);
fprintf('Test   mean ± SEM: %.3f ± %.3f\n', test_mean, test_sem);
if ~isnan(p_val)
    fprintf('Paired t-test p = %.4f  %s\n', p_val, results.significant, '(significant)' : '(n.s.)');
end
fprintf('\nInterpretation:\n');
fprintf('  Negative r = higher dF/F when closer to food\n');
fprintf('  Positive r = higher dF/F when farther from food\n');

end

%% -------------------------------------------------------------------------
%  Helper: apply time limit to a data matrix
% --------------------------------------------------------------------------
function data_out = apply_time_limit(data, time_limit_sec)
    t      = data(:, 1);
    t_mask = (t - min(t)) <= time_limit_sec;
    data_out = data(t_mask, :);
end

%% -------------------------------------------------------------------------
%  Helper: compute Euclidean distance from x,y to food well
%  Adds/replaces column 9 with distance-to-food (in the same units as x,y)
% --------------------------------------------------------------------------
function data_out = add_distance_to_food(data, food_pos)
    x = data(:, 2);
    y = data(:, 3);
    data(:, 9) = sqrt((x - food_pos(1)).^2 + (y - food_pos(2)).^2);
    data_out = data;
end

%% -------------------------------------------------------------------------
%  Helper: Pearson correlation between dF/F and distance-to-food
% --------------------------------------------------------------------------
function [r, p, is_valid] = compute_correlation(data, exclude_grooming, min_points)
    dff      = data(:, 6);
    dist     = data(:, 9);   % now distance-to-food, set by add_distance_to_food
    grooming = data(:, 10);

    % Remove NaNs
    ok = ~isnan(dff) & ~isnan(dist);
    dff      = dff(ok);
    dist     = dist(ok);
    grooming = grooming(ok);

    % Exclude grooming frames
    if exclude_grooming
        not_grooming = grooming == 0;
        dff  = dff(not_grooming);
        dist = dist(not_grooming);
    end

    if length(dff) < min_points
        r = NaN; p = NaN; is_valid = false;
        return;
    end

    try
        [r, p]   = corr(dff, dist, 'type', 'Pearson');
        is_valid = true;
    catch
        r = NaN; p = NaN; is_valid = false;
    end
end