function mice = classify_memory_strength(mice, options)
% CLASSIFY_MEMORY_STRENGTH - Classify mice as having strong or weak memory
%
% Usage:
%   options.discovery_distance = 2;     % Distance threshold in cm (required)
%   options.time_limit = 10;            % Analyze only first X minutes (optional)
%   
%   mice = classify_memory_strength(mice, options);
%
% Input:
%   mice - cell array with structure:
%     mice{i,1} - mouse name
%     mice{i,2} - [x,y] food position
%     mice{i,3} - pre-test data (Nx10 double)
%     mice{i,4} - test data (Nx10 double)
%     mice{i,5} - [x0,y0] starting position
%     mice{i,6} - latencies [pre;test]
%
%   options - structure with fields:
%     .discovery_distance - distance threshold in cm (required)
%     .time_limit - (optional) analyze only first X minutes of each session
%
% Output:
%   mice - same cell array with added columns:
%     mice{i,7} - 'strong memory' or 'weak memory'
%     mice{i,8} - food zone change (test - pre) as percentage
%
% Memory Classification Criteria:
%   STRONG MEMORY: latency_pre > latency_test AND food_zone_pre < food_zone_test
%   WEAK MEMORY: all other cases
%
% Data columns:
%   1-time, 2,3-x,y position, 4-465nm, 5-405nm, 6-dff, 7-speed, 
%   8-zones (0-outside, 1-food, 2-area2, 3-area3), 9-distance, 10-grooming

%% Input validation
if nargin < 2
    error('Both mice data and options are required');
end

if ~isfield(options, 'discovery_distance') || isempty(options.discovery_distance)
    error('options.discovery_distance must be specified');
end

discovery_distance = options.discovery_distance;

% Check if time limit is specified
if isfield(options, 'time_limit') && ~isempty(options.time_limit)
    time_limit_min = options.time_limit;
    time_limit_sec = time_limit_min * 60; % Convert to seconds
    use_time_limit = true;
    fprintf('Using time limit: %.1f minutes (%.0f seconds)\n', time_limit_min, time_limit_sec);
else
    use_time_limit = false;
    fprintf('Analyzing full session duration\n');
end

fprintf('Food discovery threshold: %.1f cm\n', discovery_distance);
fprintf('Memory classification criteria:\n');
fprintf('  STRONG MEMORY: faster discovery (pre > test) AND more time in food zone (pre < test)\n');
fprintf('  WEAK MEMORY: all other cases\n\n');

%% Initialize storage
n_mice = size(mice, 1);
latencies_pre = zeros(n_mice, 1);
latencies_test = zeros(n_mice, 1);
food_zone_pre = zeros(n_mice, 1);
food_zone_test = zeros(n_mice, 1);
food_zone_change = zeros(n_mice, 1);  % New variable for food zone change
memory_classification = cell(n_mice, 1);

%% Analyze each mouse
fprintf('=== INDIVIDUAL MOUSE ANALYSIS ===\n');

for i = 1:n_mice
    fprintf('\nMouse %d/%d: %s\n', i, n_mice, mice{i,1});
    
    % Get pre and test data
    pre_data = mice{i, 3};
    test_data = mice{i, 4};
    
    % Apply time limit if specified
    if use_time_limit
        % Get time columns
        pre_time = pre_data(:, 1);
        test_time = test_data(:, 1);
        
        % Find indices within time limit (from start of session)
        pre_start_time = min(pre_time);
        test_start_time = min(test_time);
        
        pre_time_mask = (pre_time - pre_start_time) <= time_limit_sec;
        test_time_mask = (test_time - test_start_time) <= time_limit_sec;
        
        % Filter data
        pre_data_filtered = pre_data(pre_time_mask, :);
        test_data_filtered = test_data(test_time_mask, :);
    else
        % Use full data
        pre_data_filtered = pre_data;
        test_data_filtered = test_data;
    end
    
    %% Calculate Food Discovery Latencies
    % Extract distance and time data
    pre_distance = pre_data_filtered(:, 9);  % Use pre-calculated distance
    test_distance = test_data_filtered(:, 9);
    pre_time = pre_data_filtered(:, 1);
    test_time = test_data_filtered(:, 1);
    
    % Find first discovery
    pre_discovery_idx = find(pre_distance <= discovery_distance, 1, 'first');
    test_discovery_idx = find(test_distance <= discovery_distance, 1, 'first');
    
    % Calculate latencies
    if ~isempty(pre_discovery_idx)
        latencies_pre(i) = pre_time(pre_discovery_idx) - pre_time(1);
        fprintf('  Pre latency: %.1f seconds (discovered)\n', latencies_pre(i));
    else
        % Mouse never discovered food - use session duration
        if use_time_limit
            latencies_pre(i) = time_limit_sec;
        else
            latencies_pre(i) = max(pre_time) - min(pre_time);
        end
        fprintf('  Pre latency: %.1f seconds (NOT discovered - using session duration)\n', latencies_pre(i));
    end
    
    if ~isempty(test_discovery_idx)
        latencies_test(i) = test_time(test_discovery_idx) - test_time(1);
        fprintf('  Test latency: %.1f seconds (discovered)\n', latencies_test(i));
    else
        % Mouse never discovered food - use session duration
        if use_time_limit
            latencies_test(i) = time_limit_sec;
        else
            latencies_test(i) = max(test_time) - min(test_time);
        end
        fprintf('  Test latency: %.1f seconds (NOT discovered - using session duration)\n', latencies_test(i));
    end
    
    %% Calculate Food Zone Time
    % Extract zone information
    pre_zones = pre_data_filtered(:, 8);
    test_zones = test_data_filtered(:, 8);
    
    % Calculate percentage time in food zone (zone 1)
    food_zone_pre(i) = sum(pre_zones == 1) / length(pre_zones) * 100;
    food_zone_test(i) = sum(test_zones == 1) / length(test_zones) * 100;
    
    % Calculate food zone change (test - pre)
    food_zone_change(i) = food_zone_test(i) - food_zone_pre(i);
    
    fprintf('  Pre food zone time: %.1f%%\n', food_zone_pre(i));
    fprintf('  Test food zone time: %.1f%%\n', food_zone_test(i));
    fprintf('  Food zone change: %.1f%% (test - pre)\n', food_zone_change(i));
    
    %% Memory Classification
    % Criteria: STRONG if (latency improves) AND (food zone time improves)
    latency_improved = latencies_pre(i) > latencies_test(i);  % Faster discovery
    food_zone_improved = food_zone_pre(i) < food_zone_test(i);  % More time in food zone
    
    if latency_improved && food_zone_improved
        memory_classification{i} = 'strong memory';
        fprintf('  → STRONG MEMORY (faster discovery AND more food zone time)\n');
    else
        memory_classification{i} = 'weak memory';
        fprintf('  → WEAK MEMORY');
        if ~latency_improved
            fprintf(' (slower/same discovery)');
        end
        if ~food_zone_improved
            fprintf(' (less/same food zone time)');
        end
        fprintf('\n');
    end
    
    % Add detailed breakdown
    fprintf('    Latency change: %.1f → %.1f seconds (%.1f change)\n', ...
            latencies_pre(i), latencies_test(i), latencies_test(i) - latencies_pre(i));
    fprintf('    Food zone change: %.1f%% → %.1f%% (%.1f%% change)\n', ...
            food_zone_pre(i), food_zone_test(i), food_zone_change(i));
end

%% Add classifications and changes to mice array
fprintf('\n=== ADDING CLASSIFICATIONS AND CHANGES TO MICE ARRAY ===\n');

% Check if columns 7 and 8 already exist
if size(mice, 2) >= 7
    fprintf('Warning: Column 7 already exists. Overwriting with memory classifications.\n');
end
if size(mice, 2) >= 8
    fprintf('Warning: Column 8 already exists. Overwriting with food zone changes.\n');
end

% Add memory classifications to column 7 and food zone changes to column 8
for i = 1:n_mice
    mice{i, 7} = memory_classification{i};
    mice{i, 8} = food_zone_change(i);
end

fprintf('Column 7: Memory classification (''strong memory'' or ''weak memory'')\n');
fprintf('Column 8: Food zone change (test - pre) as percentage\n');

%% Summary Statistics
fprintf('\n=== SUMMARY ===\n');

strong_memory_count = sum(strcmp(memory_classification, 'strong memory'));
weak_memory_count = sum(strcmp(memory_classification, 'weak memory'));

fprintf('Total mice analyzed: %d\n', n_mice);
fprintf('Strong memory: %d mice (%.1f%%)\n', strong_memory_count, strong_memory_count/n_mice*100);
fprintf('Weak memory: %d mice (%.1f%%)\n', weak_memory_count, weak_memory_count/n_mice*100);

% Summary statistics for food zone changes
fprintf('\nFood Zone Change Statistics:\n');
fprintf('  Overall: %.1f ± %.1f%% (range: %.1f to %.1f%%)\n', ...
        mean(food_zone_change), std(food_zone_change), min(food_zone_change), max(food_zone_change));

if strong_memory_count > 0
    strong_indices = strcmp(memory_classification, 'strong memory');
    fprintf('  Strong memory mice: %.1f ± %.1f%%\n', ...
            mean(food_zone_change(strong_indices)), std(food_zone_change(strong_indices)));
end

if weak_memory_count > 0
    weak_indices = strcmp(memory_classification, 'weak memory');
    fprintf('  Weak memory mice: %.1f ± %.1f%%\n', ...
            mean(food_zone_change(weak_indices)), std(food_zone_change(weak_indices)));
end

% List mice by classification
fprintf('\nSTRONG MEMORY mice:\n');
strong_indices = find(strcmp(memory_classification, 'strong memory'));
for i = 1:length(strong_indices)
    idx = strong_indices(i);
    fprintf('  %s: %.1fs→%.1fs latency, %.1f%%→%.1f%% food zone (change: +%.1f%%)\n', ...
            mice{idx,1}, latencies_pre(idx), latencies_test(idx), ...
            food_zone_pre(idx), food_zone_test(idx), food_zone_change(idx));
end

fprintf('\nWEAK MEMORY mice:\n');
weak_indices = find(strcmp(memory_classification, 'weak memory'));
for i = 1:length(weak_indices)
    idx = weak_indices(i);
    fprintf('  %s: %.1fs→%.1fs latency, %.1f%%→%.1f%% food zone (change: %.1f%%)\n', ...
            mice{idx,1}, latencies_pre(idx), latencies_test(idx), ...
            food_zone_pre(idx), food_zone_test(idx), food_zone_change(idx));
end

%% Create visualization
figure('Name', 'Memory Classification Analysis', 'Position', [100, 100, 1200, 800]);

% Latency improvement plot
subplot(2, 3, 1);
strong_mask = strcmp(memory_classification, 'strong memory');
weak_mask = strcmp(memory_classification, 'weak memory');

% Plot strong memory mice
if sum(strong_mask) > 0
    plot([latencies_pre(strong_mask), latencies_test(strong_mask)]', 'g-o', 'LineWidth', 2);
    hold on;
end

% Plot weak memory mice
if sum(weak_mask) > 0
    plot([latencies_pre(weak_mask), latencies_test(weak_mask)]', 'r-o', 'LineWidth', 1);
end

% Add group means
plot([mean(latencies_pre), mean(latencies_test)], 'k-o', 'LineWidth', 3, 'MarkerSize', 8);

title('Food Discovery Latency');
ylabel('Latency (seconds)');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Pre', 'Test'});
legend('Strong Memory', 'Weak Memory', 'Group Mean', 'Location', 'best');
grid on;

% Food zone time plot
subplot(2, 3, 2);

% Plot strong memory mice
if sum(strong_mask) > 0
    plot([food_zone_pre(strong_mask), food_zone_test(strong_mask)]', 'g-o', 'LineWidth', 2);
    hold on;
end

% Plot weak memory mice
if sum(weak_mask) > 0
    plot([food_zone_pre(weak_mask), food_zone_test(weak_mask)]', 'r-o', 'LineWidth', 1);
end

% Add group means
plot([mean(food_zone_pre), mean(food_zone_test)], 'k-o', 'LineWidth', 3, 'MarkerSize', 8);

title('Food Zone Time');
ylabel('Time in Food Zone (%)');
set(gca, 'XTick', 1:2, 'XTickLabel', {'Pre', 'Test'});
legend('Strong Memory', 'Weak Memory', 'Group Mean', 'Location', 'best');
grid on;

% Classification pie chart
subplot(2, 3, 3);
pie([strong_memory_count, weak_memory_count], {'Strong Memory', 'Weak Memory'});
title(sprintf('Memory Classification\n(n=%d mice)', n_mice));

% Scatter plot: Latency change vs Food zone change
subplot(2, 3, 4);
latency_change = latencies_test - latencies_pre;  % Negative = improvement
food_zone_change_for_plot = food_zone_test - food_zone_pre;  % Positive = improvement

scatter(latency_change(strong_mask), food_zone_change_for_plot(strong_mask), 100, 'g', 'filled', 'o');
hold on;
scatter(latency_change(weak_mask), food_zone_change_for_plot(weak_mask), 100, 'r', 'filled', 's');

% Add quadrant lines
plot([0, 0], ylim, 'k--', 'LineWidth', 1);
plot(xlim, [0, 0], 'k--', 'LineWidth', 1);

% Add quadrant labels
text_props = {'FontSize', 10, 'FontWeight', 'bold'};
text(min(xlim)*0.8, max(ylim)*0.8, 'STRONG\n(faster + more food)', text_props{:}, 'Color', 'green');
text(max(xlim)*0.5, max(ylim)*0.8, 'WEAK\n(slower + more food)', text_props{:}, 'Color', 'red');
text(min(xlim)*0.8, min(ylim)*0.8, 'WEAK\n(faster + less food)', text_props{:}, 'Color', 'red');
text(max(xlim)*0.5, min(ylim)*0.8, 'WEAK\n(slower + less food)', text_props{:}, 'Color', 'red');

xlabel('Latency Change (Test - Pre, seconds)');
ylabel('Food Zone Change (Test - Pre, %)');
title('Memory Classification Criteria');
legend('Strong Memory', 'Weak Memory', 'Location', 'best');
grid on;

% Food zone change histogram
subplot(2, 3, 5);
histogram(food_zone_change(strong_mask), 'FaceColor', 'g', 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'BinWidth', 2);
hold on;
histogram(food_zone_change(weak_mask), 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'BinWidth', 2);

% Add vertical line at zero
plot([0, 0], ylim, 'k--', 'LineWidth', 2);
text(0, max(ylim)*0.9, 'No Change', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

xlabel('Food Zone Change (Test - Pre, %)');
ylabel('Number of Mice');
title('Distribution of Food Zone Changes');
legend('Strong Memory', 'Weak Memory', 'Location', 'best');
grid on;

% Individual mouse details
subplot(2, 3, 6);
axis off;

% Create text summary
text_str = {};
text_str{end+1} = sprintf('MEMORY CLASSIFICATION SUMMARY (Discovery threshold: %.1f cm)', discovery_distance);
text_str{end+1} = '';
text_str{end+1} = sprintf('Strong Memory: %d mice (%.1f%%) - Avg food zone change: +%.1f%%', ...
                         strong_memory_count, strong_memory_count/n_mice*100, ...
                         mean(food_zone_change(strong_mask)));
for i = 1:min(length(strong_indices), 4)  % Limit to first 4 to save space
    idx = strong_indices(i);
    text_str{end+1} = sprintf('  • %s (+%.1f%%)', mice{idx,1}, food_zone_change(idx));
end
if length(strong_indices) > 4
    text_str{end+1} = sprintf('  ... and %d more', length(strong_indices) - 4);
end

text_str{end+1} = '';
text_str{end+1} = sprintf('Weak Memory: %d mice (%.1f%%) - Avg food zone change: %.1f%%', ...
                         weak_memory_count, weak_memory_count/n_mice*100, ...
                         mean(food_zone_change(weak_mask)));
for i = 1:min(length(weak_indices), 4)  % Limit to first 4 to save space
    idx = weak_indices(i);
    text_str{end+1} = sprintf('  • %s (%.1f%%)', mice{idx,1}, food_zone_change(idx));
end
if length(weak_indices) > 4
    text_str{end+1} = sprintf('  ... and %d more', length(weak_indices) - 4);
end

text_str{end+1} = '';
text_str{end+1} = sprintf('Column 7: Memory classification');
text_str{end+1} = sprintf('Column 8: Food zone change (test-pre %%)');

text(0.05, 0.95, text_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 9);

fprintf('\nMemory classification completed! Columns 7 and 8 added to mice array.\n');
fprintf('Column 7: Memory classification (''strong memory'' or ''weak memory'')\n');
fprintf('Column 8: Food zone change (test - pre) as percentage\n');

end