function results = plot_runs_dff_distance_correlation(run_data, condition, options)
% PLOT_RUNS_DFF_DISTANCE_CORRELATION - Plot dF/F vs distance correlations for runs only
%                                       Compares Before vs Test within Towards and Away runs
%
% Usage:
%   plot_runs_dff_distance_correlation(run_data, 'saline')
%   plot_runs_dff_distance_correlation(run_data, 'CNO')
%   plot_runs_dff_distance_correlation(run_data, 'strong')
%   plot_runs_dff_distance_correlation(run_data, 'weak')
%   
%   % With options:
%   options.min_runs = 1;               % Minimum runs per mouse per condition (default: 1)
%   options.figure_size = [100, 100, 600, 500]; % Figure position and size
%   options.session_filter = 'both';    % 'pre', 'test', or 'both' (default: 'both')
%                                       % Note: 'both' recommended for Before vs Test comparison
%   
%   results = plot_runs_dff_distance_correlation(run_data, condition, options);
%
% Input:
%   run_data - Output from analyze_cheeseboard_runs function
%   condition - 'saline', 'CNO', 'strong', 'weak', 'saline_strong', etc.
%   options - Optional parameters structure
%
% Output:
%   Creates a 2x2 grouped bar plot:
%   - X-axis: Towards Food (blue) vs Away from Food (red)
%   - Each group has Before (light) and Test (dark) bars
%   - Individual mouse trajectories connect Before to Test within each run type
%   - Statistical comparison: Before vs Test within each run type

%% Input validation and defaults
if nargin < 2
    error('Both run_data and condition are required');
end

if nargin < 3
    options = struct();
end

% Set defaults
if ~isfield(options, 'min_runs')
    options.min_runs = 1;  % Changed from 2 to 1 to include mice with single runs
end

if ~isfield(options, 'figure_size')
    options.figure_size = [100, 100, 600, 500];
end

if ~isfield(options, 'session_filter')
    options.session_filter = 'both';
end

if isempty(run_data)
    error('No run data provided');
end

fprintf('=== RUNS dF/F-DISTANCE CORRELATION ANALYSIS ===\n');
fprintf('Condition: %s\n', condition);
fprintf('Minimum runs per mouse per condition: %d\n', options.min_runs);
fprintf('Session filter: %s\n', options.session_filter);

%% Filter run_data by condition and session
filtered_data = filter_run_data_by_condition(run_data, condition, options.session_filter);

if isempty(filtered_data)
    error('No mice found for condition: %s', condition);
end

% Get unique mice
unique_mice = unique({filtered_data.mouse_id});
n_mice = length(unique_mice);

fprintf('Found %d mice for %s condition\n', n_mice, condition);

%% Calculate average correlations per mouse per run type and session
mouse_towards_before_corr = NaN(n_mice, 1);
mouse_towards_test_corr = NaN(n_mice, 1);
mouse_away_before_corr = NaN(n_mice, 1);
mouse_away_test_corr = NaN(n_mice, 1);

mouse_towards_before_valid = false(n_mice, 1);
mouse_towards_test_valid = false(n_mice, 1);
mouse_away_before_valid = false(n_mice, 1);
mouse_away_test_valid = false(n_mice, 1);

mouse_towards_before_nruns = zeros(n_mice, 1);
mouse_towards_test_nruns = zeros(n_mice, 1);
mouse_away_before_nruns = zeros(n_mice, 1);
mouse_away_test_nruns = zeros(n_mice, 1);

for i = 1:n_mice
    mouse_id = unique_mice{i};
    fprintf('\nProcessing mouse %d/%d: %s\n', i, n_mice, mouse_id);
    
    % Get all data for this mouse
    mouse_data = filtered_data(strcmp({filtered_data.mouse_id}, mouse_id));
    
    % Collect runs by session and type
    towards_before_runs = [];
    towards_test_runs = [];
    away_before_runs = [];
    away_test_runs = [];
    
    for j = 1:length(mouse_data)
        session = mouse_data(j).session;
        session_runs = mouse_data(j).runs;
        
        for k = 1:length(session_runs)
            if strcmp(session_runs(k).type, 'towards')
                if strcmp(session, 'pre')
                    towards_before_runs = [towards_before_runs; session_runs(k)];
                else
                    towards_test_runs = [towards_test_runs; session_runs(k)];
                end
            elseif strcmp(session_runs(k).type, 'away')
                if strcmp(session, 'pre')
                    away_before_runs = [away_before_runs; session_runs(k)];
                else
                    away_test_runs = [away_test_runs; session_runs(k)];
                end
            end
        end
    end
    
    mouse_towards_before_nruns(i) = length(towards_before_runs);
    mouse_towards_test_nruns(i) = length(towards_test_runs);
    mouse_away_before_nruns(i) = length(away_before_runs);
    mouse_away_test_nruns(i) = length(away_test_runs);
    
    fprintf('  Towards Before: %d, Towards Test: %d\n', mouse_towards_before_nruns(i), mouse_towards_test_nruns(i));
    fprintf('  Away Before: %d, Away Test: %d\n', mouse_away_before_nruns(i), mouse_away_test_nruns(i));
    
    % Calculate average correlations for each condition
    % Towards Before
    if length(towards_before_runs) >= options.min_runs
        correlations = [];
        for k = 1:length(towards_before_runs)
            [corr_val, is_valid] = calculate_run_correlation(towards_before_runs(k));
            if is_valid, correlations = [correlations; corr_val]; end
        end
        if ~isempty(correlations)
            mouse_towards_before_corr(i) = mean(correlations);
            mouse_towards_before_valid(i) = true;
            fprintf('    Towards Before: %.3f (avg of %d runs)\n', mouse_towards_before_corr(i), length(correlations));
        end
    end
    
    % Towards Test
    if length(towards_test_runs) >= options.min_runs
        correlations = [];
        for k = 1:length(towards_test_runs)
            [corr_val, is_valid] = calculate_run_correlation(towards_test_runs(k));
            if is_valid, correlations = [correlations; corr_val]; end
        end
        if ~isempty(correlations)
            mouse_towards_test_corr(i) = mean(correlations);
            mouse_towards_test_valid(i) = true;
            fprintf('    Towards Test: %.3f (avg of %d runs)\n', mouse_towards_test_corr(i), length(correlations));
        end
    end
    
    % Away Before
    if length(away_before_runs) >= options.min_runs
        correlations = [];
        for k = 1:length(away_before_runs)
            [corr_val, is_valid] = calculate_run_correlation(away_before_runs(k));
            if is_valid, correlations = [correlations; corr_val]; end
        end
        if ~isempty(correlations)
            mouse_away_before_corr(i) = mean(correlations);
            mouse_away_before_valid(i) = true;
            fprintf('    Away Before: %.3f (avg of %d runs)\n', mouse_away_before_corr(i), length(correlations));
        end
    end
    
    % Away Test
    if length(away_test_runs) >= options.min_runs
        correlations = [];
        for k = 1:length(away_test_runs)
            [corr_val, is_valid] = calculate_run_correlation(away_test_runs(k));
            if is_valid, correlations = [correlations; corr_val]; end
        end
        if ~isempty(correlations)
            mouse_away_test_corr(i) = mean(correlations);
            mouse_away_test_valid(i) = true;
            fprintf('    Away Test: %.3f (avg of %d runs)\n', mouse_away_test_corr(i), length(correlations));
        end
    end
end

%% Create clean plot with color-coded run types
condition_name = get_condition_display_name(condition);

fig = figure('Name', sprintf('Runs dF/F-Distance Correlations: %s', condition_name), ...
             'Position', options.figure_size, 'Color', 'white');

% Get valid data for each condition
towards_before_data = mouse_towards_before_corr(mouse_towards_before_valid);
towards_test_data = mouse_towards_test_corr(mouse_towards_test_valid);
away_before_data = mouse_away_before_corr(mouse_away_before_valid);
away_test_data = mouse_away_test_corr(mouse_away_test_valid);

if isempty(towards_before_data) && isempty(towards_test_data) && ...
   isempty(away_before_data) && isempty(away_test_data)
    text(0.5, 0.5, 'No valid correlations found', 'HorizontalAlignment', 'center', ...
         'Units', 'normalized', 'FontSize', 14, 'Color', [0.5 0.5 0.5]);
    return;
end

% Calculate means and SEMs for 2x2 design
group_means = [nanmean(towards_before_data), nanmean(towards_test_data); 
               nanmean(away_before_data), nanmean(away_test_data)];
group_errors = [nanstd(towards_before_data)/sqrt(length(towards_before_data)), ...
                nanstd(towards_test_data)/sqrt(length(towards_test_data));
                nanstd(away_before_data)/sqrt(length(away_before_data)), ...
                nanstd(away_test_data)/sqrt(length(away_test_data))];

% Handle NaN values
group_means(isnan(group_means)) = 0;
group_errors(isnan(group_errors)) = 0;

% Define color scheme: Blue for towards, Red for away; Light for before, Dark for test
towards_before_color = [0.7 0.8 1.0];    % Light blue
towards_test_color = [0.2 0.4 0.8];      % Dark blue
away_before_color = [1.0 0.7 0.7];       % Light red
away_test_color = [0.8 0.2 0.2];         % Dark red

% Create grouped bar plot with custom colors
b = bar(group_means, 'BarWidth', 0.9);  % or 1.0 for even wider bars
b(1).FaceColor = 'flat';  % Before colors
b(2).FaceColor = 'flat';  % Test colors

% Set individual bar colors
b(1).CData(1,:) = towards_before_color;  % Towards Before
b(1).CData(2,:) = away_before_color;     % Away Before
b(2).CData(1,:) = towards_test_color;    % Towards Test
b(2).CData(2,:) = away_test_color;       % Away Test

b(1).FaceAlpha = 0.8;
b(2).FaceAlpha = 0.8;

hold on;

% Add error bars
ngroups = size(group_means, 1);
nbars = size(group_means, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for j = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
    errorbar(x, group_means(:,j), group_errors(:,j), 'k.', 'LineWidth', 1.5);
end

% Add individual mouse trajectories for mice with both before and test data
% Get valid mice for each pair (mice with both before and test data)
towards_valid_both = mouse_towards_before_valid & mouse_towards_test_valid;
away_valid_both = mouse_away_before_valid & mouse_away_test_valid;

% Plot towards runs (before vs test) - Blue theme
if sum(towards_valid_both) > 0
    x_before = 1 - groupwidth/2 + groupwidth/(2*nbars);
    x_test = 1 - groupwidth/2 + 3*groupwidth/(2*nbars);
    
    % Individual point and line colors for towards runs
    towards_point_color = [0.1 0.3 0.6];  % Dark blue for points
    towards_line_color = [0.4 0.6 0.8];   % Medium blue for lines
    
    for i = 1:n_mice
        if towards_valid_both(i)
            % Add jitter for visibility
            jitter_before = x_before + (rand-0.5)*0;
            jitter_test = x_test + (rand-0.5)*0;
            
            scatter(jitter_before, mouse_towards_before_corr(i), 25, towards_point_color, 'filled', 'o', ...
                   'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
            scatter(jitter_test, mouse_towards_test_corr(i), 25, towards_point_color, 'filled', 'o', ...
                   'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
            plot([jitter_before, jitter_test], [mouse_towards_before_corr(i), mouse_towards_test_corr(i)], ...
                 '-', 'Color', towards_line_color, 'LineWidth', 1.5);
        end
    end
end

% Plot away runs (before vs test) - Red theme
if sum(away_valid_both) > 0
    x_before = 2 - groupwidth/2 + groupwidth/(2*nbars);
    x_test = 2 - groupwidth/2 + 3*groupwidth/(2*nbars);
    
    % Individual point and line colors for away runs
    away_point_color = [0.6 0.1 0.1];     % Dark red for points
    away_line_color = [0.8 0.4 0.4];      % Medium red for lines
    
    for i = 1:n_mice
        if away_valid_both(i)
            % Add jitter for visibility
            jitter_before = x_before + (rand-0.5)*0;
            jitter_test = x_test + (rand-0.5)*0;
            
            scatter(jitter_before, mouse_away_before_corr(i), 25, away_point_color, 'filled', 'o', ...
                   'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
            scatter(jitter_test, mouse_away_test_corr(i), 25, away_point_color, 'filled', 'o', ...
                   'MarkerEdgeColor', 'none', 'LineWidth', 0.5);
            plot([jitter_before, jitter_test], [mouse_away_before_corr(i), mouse_away_test_corr(i)], ...
                 '-', 'Color', away_line_color, 'LineWidth', 1.5);
        end
    end
end

% Perform paired statistical tests on mice with both conditions
towards_paired_before = mouse_towards_before_corr(towards_valid_both);
towards_paired_test = mouse_towards_test_corr(towards_valid_both);
away_paired_before = mouse_away_before_corr(away_valid_both);
away_paired_test = mouse_away_test_corr(away_valid_both);

% Statistical tests - only if we have paired data
p_towards = NaN;
p_away = NaN;

if length(towards_paired_before) > 1
    try
        [~, p_towards] = ttest(towards_paired_before, towards_paired_test);
    catch
        p_towards = NaN;
    end
end

if length(away_paired_before) > 1
    try
        [~, p_away] = ttest(away_paired_before, away_paired_test);
    catch
        p_away = NaN;
    end
end

% Get y-limits for significance indicators
all_vals = [towards_before_data; towards_test_data; away_before_data; away_test_data];
if ~isempty(all_vals)
    y_range = max(all_vals) - min(all_vals);
    max_y = max(towards_before_data) + y_range * 0.1;
    
    % Add significance indicators for within-run-type comparisons
    if ~isnan(p_towards) && p_towards < 0.05
        text(1, max_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, ...
             'Color', 'k', 'FontWeight', 'bold');
    end
    if ~isnan(p_away) && p_away < 0.05
        text(2, max_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, ...
             'Color', 'k', 'FontWeight', 'bold');
    end
    
    ylim([min(all_vals) - y_range * 0.1, max(all_vals) + y_range * 0.05]);
end

xlim([0.5, 2.5]);  % Adjust these values to get the desired spacing

% Clean up axes
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Towards Food', 'Away from Food'}, ...
    'FontSize', 12, 'FontWeight', 'normal');
ylabel('dF/F-Distance Correlation (r)', 'FontSize', 14, 'FontWeight', 'bold');

% Add custom legend with color-coded entries
legend_handles = [];
legend_labels = {};

% Create dummy handles for legend with appropriate colors
h1 = bar(NaN, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.8);
h2 = bar(NaN, 'FaceColor', [0.2 0.2 0.2], 'FaceAlpha', 0.8);
legend_handles = [h1, h2];
legend_labels = {'Before', 'Test'};

legend(legend_handles, legend_labels, 'Location', 'best');
legend('boxoff')

% Remove grid and clean up appearance
grid off;
box off;
set(gca, 'TickDir', 'out', 'LineWidth', 1.5);

% Add title with color information
title(sprintf('%s Group (Runs Only, n=%d mice)', condition_name, n_mice), ...
      'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

hold off;

%% Store results
results = struct();
results.condition = condition_name;
results.mouse_names = unique_mice;
results.mouse_towards_before_corr = mouse_towards_before_corr;
results.mouse_towards_test_corr = mouse_towards_test_corr;
results.mouse_away_before_corr = mouse_away_before_corr;
results.mouse_away_test_corr = mouse_away_test_corr;
results.mouse_towards_before_valid = mouse_towards_before_valid;
results.mouse_towards_test_valid = mouse_towards_test_valid;
results.mouse_away_before_valid = mouse_away_before_valid;
results.mouse_away_test_valid = mouse_away_test_valid;
results.mouse_towards_before_nruns = mouse_towards_before_nruns;
results.mouse_towards_test_nruns = mouse_towards_test_nruns;
results.mouse_away_before_nruns = mouse_away_before_nruns;
results.mouse_away_test_nruns = mouse_away_test_nruns;
results.n_mice = n_mice;

% Calculate means and SEMs
results.towards_before_mean = nanmean(towards_before_data);
results.towards_test_mean = nanmean(towards_test_data);
results.away_before_mean = nanmean(away_before_data);
results.away_test_mean = nanmean(away_test_data);

results.towards_before_sem = nanstd(towards_before_data)/sqrt(length(towards_before_data));
results.towards_test_sem = nanstd(towards_test_data)/sqrt(length(towards_test_data));
results.away_before_sem = nanstd(away_before_data)/sqrt(length(away_before_data));
results.away_test_sem = nanstd(away_test_data)/sqrt(length(away_test_data));

% Statistical tests
if length(towards_paired_before) > 1
    results.ttest_towards_p = p_towards;
    results.towards_significant = ~isnan(p_towards) && p_towards < 0.05;
    results.towards_paired_n = length(towards_paired_before);
else
    results.ttest_towards_p = NaN;
    results.towards_significant = false;
    results.towards_paired_n = 0;
end

if length(away_paired_before) > 1
    results.ttest_away_p = p_away;
    results.away_significant = ~isnan(p_away) && p_away < 0.05;
    results.away_paired_n = length(away_paired_before);
else
    results.ttest_away_p = NaN;
    results.away_significant = false;
    results.away_paired_n = 0;
end

%% Print summary
fprintf('\n=== %s GROUP RUNS CORRELATION SUMMARY ===\n', upper(condition_name));
fprintf('Total mice: %d\n', n_mice);
fprintf('\nVALID CORRELATIONS:\n');
fprintf('  Towards Before: %d (%.1f%%)\n', sum(mouse_towards_before_valid), sum(mouse_towards_before_valid)/n_mice*100);
fprintf('  Towards Test: %d (%.1f%%)\n', sum(mouse_towards_test_valid), sum(mouse_towards_test_valid)/n_mice*100);
fprintf('  Away Before: %d (%.1f%%)\n', sum(mouse_away_before_valid), sum(mouse_away_before_valid)/n_mice*100);
fprintf('  Away Test: %d (%.1f%%)\n', sum(mouse_away_test_valid), sum(mouse_away_test_valid)/n_mice*100);

fprintf('\nMEANS ± SEM:\n');
fprintf('  Towards Before: %.3f ± %.3f\n', results.towards_before_mean, results.towards_before_sem);
fprintf('  Towards Test: %.3f ± %.3f\n', results.towards_test_mean, results.towards_test_sem);
fprintf('  Away Before: %.3f ± %.3f\n', results.away_before_mean, results.away_before_sem);
fprintf('  Away Test: %.3f ± %.3f\n', results.away_test_mean, results.away_test_sem);

fprintf('\nSTATISTICAL TESTS (Paired Before vs Test):\n');
if ~isnan(results.ttest_towards_p)
    if results.towards_significant
        fprintf('  Towards: p = %.3f *** (significant, n=%d paired mice)\n', results.ttest_towards_p, results.towards_paired_n);
    else
        fprintf('  Towards: p = %.3f (not significant, n=%d paired mice)\n', results.ttest_towards_p, results.towards_paired_n);
    end
else
    fprintf('  Towards: Not enough paired data for comparison\n');
end

if ~isnan(results.ttest_away_p)
    if results.away_significant
        fprintf('  Away: p = %.3f *** (significant, n=%d paired mice)\n', results.ttest_away_p, results.away_paired_n);
    else
        fprintf('  Away: p = %.3f (not significant, n=%d paired mice)\n', results.ttest_away_p, results.away_paired_n);
    end
else
    fprintf('  Away: Not enough paired data for comparison\n');
end

% Print individual mouse data
fprintf('\nINDIVIDUAL MOUSE DATA:\n');
fprintf('Mouse\t\tT_Bef\tT_Test\tA_Bef\tA_Test\tT_B_n\tT_T_n\tA_B_n\tA_T_n\n');
fprintf('-----\t\t-----\t------\t-----\t------\t-----\t-----\t-----\t-----\n');
for i = 1:n_mice
    tb_str = 'N/A'; tt_str = 'N/A'; ab_str = 'N/A'; at_str = 'N/A';
    if mouse_towards_before_valid(i), tb_str = sprintf('%.3f', mouse_towards_before_corr(i)); end
    if mouse_towards_test_valid(i), tt_str = sprintf('%.3f', mouse_towards_test_corr(i)); end
    if mouse_away_before_valid(i), ab_str = sprintf('%.3f', mouse_away_before_corr(i)); end
    if mouse_away_test_valid(i), at_str = sprintf('%.3f', mouse_away_test_corr(i)); end
    
    fprintf('%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n', unique_mice{i}, ...
            tb_str, tt_str, ab_str, at_str, ...
            mouse_towards_before_nruns(i), mouse_towards_test_nruns(i), ...
            mouse_away_before_nruns(i), mouse_away_test_nruns(i));
end

fprintf('\nINTERPRETATION:\n');
fprintf('  Negative correlations = higher dF/F closer to food\n');
fprintf('  Positive correlations = higher dF/F farther from food\n');
fprintf('  Each point represents average correlation across all runs for that mouse/condition\n');
fprintf('  Blue = Towards runs, Red = Away runs\n');
fprintf('  Light colors = Before, Dark colors = Test\n');
fprintf('  * indicates significant Before vs Test difference within run type\n');

end

%% Helper functions

function filtered_data = filter_run_data_by_condition(run_data, condition, session_filter)
    % Filter run_data by condition and session
    % Note: For Before vs Test comparison, we need both sessions
    
    filtered_data = [];
    
    for i = 1:length(run_data)
        mouse_id = run_data(i).mouse_id;
        session = run_data(i).session;
        memory_strength = run_data(i).memory_strength;
        
        % Apply session filter - but warn if limiting to single session
        if ~strcmp(session_filter, 'both')
            if ~strcmp(session, session_filter)
                continue;
            end
            if i == 1  % Only warn once
                fprintf('Warning: Session filter "%s" selected - Before vs Test comparison will be limited\n', session_filter);
            end
        end
        
        % Apply condition filter
        include_mouse = false;
        
        switch lower(condition)
            case 'saline'
                if contains(mouse_id, 'saline')
                    include_mouse = true;
                end
            case 'cno'
                if contains(mouse_id, 'CNO')
                    include_mouse = true;
                end
            case 'strong'
                if strcmp(memory_strength, 'strong memory')
                    include_mouse = true;
                end
            case 'weak'
                if strcmp(memory_strength, 'weak memory')
                    include_mouse = true;
                end
            case 'saline_strong'
                if contains(mouse_id, 'saline') && strcmp(memory_strength, 'strong memory')
                    include_mouse = true;
                end
            case 'saline_weak'
                if contains(mouse_id, 'saline') && strcmp(memory_strength, 'weak memory')
                    include_mouse = true;
                end
            case 'cno_strong'
                if contains(mouse_id, 'CNO') && strcmp(memory_strength, 'strong memory')
                    include_mouse = true;
                end
            case 'cno_weak'
                if contains(mouse_id, 'CNO') && strcmp(memory_strength, 'weak memory')
                    include_mouse = true;
                end
            otherwise
                error('Invalid condition: %s', condition);
        end
        
        if include_mouse
            filtered_data = [filtered_data; run_data(i)];
        end
    end
end

function condition_name = get_condition_display_name(condition)
    % Convert condition string to display name
    
    switch lower(condition)
        case 'saline'
            condition_name = 'Saline';
        case 'cno'
            condition_name = 'CNO';
        case 'strong'
            condition_name = 'Strong Memory';
        case 'weak'
            condition_name = 'Weak Memory';
        case 'saline_strong'
            condition_name = 'Saline Strong Memory';
        case 'saline_weak'
            condition_name = 'Saline Weak Memory';
        case 'cno_strong'
            condition_name = 'CNO Strong Memory';
        case 'cno_weak'
            condition_name = 'CNO Weak Memory';
        otherwise
            condition_name = condition;
    end
end

function [correlation, is_valid] = calculate_run_correlation(run)
    % Calculate correlation for a single run
    
    % Minimum points for correlation - reduced to 3 for single runs
    min_points = 3;
    
    % Extract data
    distance = run.distance;
    dff = run.dff;
    
    % Remove NaN values
    valid_idx = ~isnan(distance) & ~isnan(dff) & isfinite(distance) & isfinite(dff);
    
    if sum(valid_idx) < min_points
        correlation = NaN;
        is_valid = false;
        return;
    end
    
    distance_clean = distance(valid_idx);
    dff_clean = dff(valid_idx);
    
    % Calculate correlation
    try
        correlation = corr(distance_clean, dff_clean, 'Type', 'Pearson');
        is_valid = ~isnan(correlation);
    catch
        correlation = NaN;
        is_valid = false;
    end
end