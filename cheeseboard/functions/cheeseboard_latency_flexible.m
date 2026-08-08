function results = cheeseboard_latency_flexible(mice, options)
% CHEESEBOARD_LATENCY_FLEXIBLE - Paper-quality latency analysis for chosen group
%
% Usage:
%   options.group = 'saline';           % Choose 'saline' or 'CNO' (required)
%   options.time_limit = 10;            % Analyze only first 10 minutes (optional)
%   options.discovery_distance = 2;     % Distance threshold in cm (default: 2)
%   options.figure_size = [6, 4];       % Figure size in inches (optional)
%   options.font_size = 12;             % Font size for labels (optional)
%   options.save_figure = true;         % Save figure as PDF/PNG (optional)
%   options.output_dir = './figures';   % Output directory for saved figures (optional)
%   
%   results = cheeseboard_latency_flexible(mice, options);
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
%     .group - 'saline' or 'CNO' (required)
%     .time_limit - (optional) analyze only first X minutes of each session
%     .discovery_distance - (optional) distance threshold in cm (default: 2)
%     .figure_size - [width, height] in inches (default: [6, 4])
%     .font_size - font size for all text (default: 14)
%     .save_figure - whether to save figure (default: false)
%     .output_dir - directory to save figures (default: './figures')
%
% Data columns:
%   1-time, 2,3-x,y position, 4-465nm, 5-405nm, 6-dff, 7-speed, 
%   8-zones (0-outside, 1-food, 2-area2, 3-area3), 9-distance, 10-grooming

%% Input validation and default parameters
if nargin < 1
    error('mice data is required');
end

if nargin < 2
    options = struct();
end

% Check if group is specified
if ~isfield(options, 'group') || isempty(options.group)
    error('options.group must be specified as either "saline" or "CNO"');
end

% Validate group choice
if ~strcmpi(options.group, 'saline') && ~strcmpi(options.group, 'CNO')
    error('options.group must be either "saline" or "CNO"');
end

% Set default options
if ~isfield(options, 'time_limit'), options.time_limit = []; end
if ~isfield(options, 'discovery_distance'), options.discovery_distance = 2; end
if ~isfield(options, 'figure_size'), options.figure_size = [6, 4]; end
if ~isfield(options, 'font_size'), options.font_size = 14; end
if ~isfield(options, 'save_figure'), options.save_figure = false; end
if ~isfield(options, 'output_dir'), options.output_dir = './figures'; end

% Check if time limit is specified
if ~isempty(options.time_limit)
    time_limit_min = options.time_limit;
    time_limit_sec = time_limit_min * 60;
    use_time_limit = true;
    fprintf('Using time limit: %.1f minutes (%.0f seconds)\n', time_limit_min, time_limit_sec);
else
    use_time_limit = false;
    fprintf('Analyzing full session duration\n');
end

fprintf('Food discovery threshold: %.1f cm\n', options.discovery_distance);

%% Extract mice based on selected group
mouse_names = mice(:,1);

if strcmpi(options.group, 'saline')
    group_idx = contains(mouse_names, 'saline');
    group_name = 'Saline';
    fprintf('Analyzing SALINE group\n');
elseif strcmpi(options.group, 'CNO')
    group_idx = contains(mouse_names, 'CNO');
    group_name = 'CNO';
    fprintf('Analyzing CNO group\n');
end

selected_mice = find(group_idx);

if isempty(selected_mice)
    error('No %s mice found in the dataset', options.group);
end

n_mice = length(selected_mice);
fprintf('Found %d %s mice\n', n_mice, options.group);

%% Latency analysis
fprintf('\n=== FOOD DISCOVERY LATENCY ANALYSIS (%s GROUP) ===\n', upper(options.group));

% Initialize storage
latencies_pre = zeros(n_mice, 1);
latencies_test = zeros(n_mice, 1);
discovery_found_pre = false(n_mice, 1);
discovery_found_test = false(n_mice, 1);

% Calculate latencies for each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    fprintf('Processing mouse %d/%d: %s\n', i, n_mice, mice{mouse_idx,1});
    
    % Get pre and test data
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};
    food_pos = mice{mouse_idx, 2};
    
    % Apply time limit if specified
    if use_time_limit
        pre_time = pre_data(:, 1);
        test_time = test_data(:, 1);
        
        pre_start_time = min(pre_time);
        test_start_time = min(test_time);
        
        pre_time_mask = (pre_time - pre_start_time) <= time_limit_sec;
        test_time_mask = (test_time - test_start_time) <= time_limit_sec;
        
        pre_data_filtered = pre_data(pre_time_mask, :);
        test_data_filtered = test_data(test_time_mask, :);
        
        if size(pre_data_filtered, 1) < 10
            warning('Mouse %s: Less than 10 data points in pre session within time limit', mice{mouse_idx,1});
        end
        if size(test_data_filtered, 1) < 10
            warning('Mouse %s: Less than 10 data points in test session within time limit', mice{mouse_idx,1});
        end
    else
        pre_data_filtered = pre_data;
        test_data_filtered = test_data;
    end
    
    % Extract distance and time data
    pre_distance = pre_data_filtered(:, 9);  % Use pre-calculated distance
    test_distance = test_data_filtered(:, 9);
    pre_time = pre_data_filtered(:, 1);
    test_time = test_data_filtered(:, 1);
    
    % Find first discovery
    pre_discovery_idx = find(pre_distance <= options.discovery_distance, 1, 'first');
    test_discovery_idx = find(test_distance <= options.discovery_distance, 1, 'first');
    
    % Calculate latencies
    if ~isempty(pre_discovery_idx)
        latencies_pre(i) = pre_time(pre_discovery_idx) - pre_time(1);
        discovery_found_pre(i) = true;
        fprintf('  Pre: Food discovered at %.1f seconds\n', latencies_pre(i));
    else
        % Mouse never discovered food - use session duration as latency
        if use_time_limit
            latencies_pre(i) = time_limit_sec;
        else
            latencies_pre(i) = max(pre_time) - min(pre_time);
        end
        discovery_found_pre(i) = false;
        fprintf('  Pre: Food NOT discovered - latency set to %.1f seconds\n', latencies_pre(i));
    end
    
    if ~isempty(test_discovery_idx)
        latencies_test(i) = test_time(test_discovery_idx) - test_time(1);
        discovery_found_test(i) = true;
        fprintf('  Test: Food discovered at %.1f seconds\n', latencies_test(i));
    else
        % Mouse never discovered food - use session duration as latency
        if use_time_limit
            latencies_test(i) = time_limit_sec;
        else
            latencies_test(i) = max(test_time) - min(test_time);
        end
        discovery_found_test(i) = false;
        fprintf('  Test: Food NOT discovered - latency set to %.1f seconds\n', latencies_test(i));
    end
end

%% Create paper-quality figure
figure_width = options.figure_size(1);
figure_height = options.figure_size(2);

fig = figure('Position', [100, 100, figure_width*96, figure_height*96], ... % 96 DPI
             'Color', 'white', 'PaperUnits', 'inches', ...
             'PaperSize', [figure_width, figure_height], ...
             'PaperPosition', [0, 0, figure_width, figure_height]);

% Define colors - different for each group
colors = struct();
if strcmpi(options.group, 'saline')
    colors.test = [0.2, 0.4, 0.6];      % Dark blue for pre
    colors.pre = [0.4, 0.6, 0.8];     % Light blue for test
else % CNO
    colors.test = [0.6, 0.2, 0.2];      % Dark red for pre
    colors.pre = [0.8, 0.4, 0.4];     % Light red for test
end

% Calculate means and SEMs
latency_mean_pre = mean(latencies_pre);
latency_mean_test = mean(latencies_test);
latency_sem_pre = std(latencies_pre) / sqrt(n_mice);
latency_sem_test = std(latencies_test) / sqrt(n_mice);

% Create single plot for latencies
% Create grouped bar plot
bar_data = [latency_mean_pre, latency_mean_test];
bar_errors = [latency_sem_pre, latency_sem_test];

% Create bar plot
b = bar(1:2, bar_data, 'BarWidth', 0.7);
b.FaceColor = 'flat';
b.CData(1,:) = colors.pre;
b.CData(2,:) = colors.test;
b.EdgeColor = 'black';
b.LineWidth = 1.5;

hold on;

% Add error bars
errorbar(1:2, bar_data, bar_errors, 'k.', 'LineWidth', 2, 'CapSize', 10);

% Add individual data points with connecting lines
x_positions = [1, 2];
jitter_strength = 0;

for i = 1:n_mice
    % Add small jitter to x-positions for visibility
    x_jitter = x_positions + (rand(1,2)-0.5) * jitter_strength;
    
    % Different markers for discovered vs not discovered
    if discovery_found_pre(i) && discovery_found_test(i)
        % Both discovered - solid line, filled markers
        plot(x_jitter, [latencies_pre(i), latencies_test(i)], 'o-', ...
             'Color', [0.3, 0.3, 0.3], 'LineWidth', 1.5, ...
             'MarkerSize', 6, 'MarkerFaceColor', [0.6, 0.6, 0.6], ...
             'MarkerEdgeColor', 'black');
    elseif discovery_found_pre(i) && ~discovery_found_test(i)
        % Pre only - mixed markers
        plot(x_jitter(1), latencies_pre(i), 'o', 'Color', [0.3, 0.3, 0.3], ...
             'MarkerSize', 6, 'MarkerFaceColor', [0.6, 0.6, 0.6], 'MarkerEdgeColor', 'black');
        plot(x_jitter(2), latencies_test(i), 's', 'Color', [0.3, 0.3, 0.3], ...
             'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.8, 0.3, 0.3]);
        plot(x_jitter, [latencies_pre(i), latencies_test(i)], '--', 'Color', [0.6, 0.6, 0.6]);
    elseif ~discovery_found_pre(i) && discovery_found_test(i)
        % Test only - mixed markers
        plot(x_jitter(1), latencies_pre(i), 's', 'Color', [0.8, 0.3, 0.3], ...
             'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.8, 0.3, 0.3]);
        plot(x_jitter(2), latencies_test(i), 'o', 'Color', [0.3, 0.3, 0.3], ...
             'MarkerSize', 6, 'MarkerFaceColor', [0.6, 0.6, 0.6], 'MarkerEdgeColor', 'black');
        plot(x_jitter, [latencies_pre(i), latencies_test(i)], '-', 'Color', [0.6, 0.6, 0.6]);
    else
        % Neither discovered - open markers, dashed line
        plot(x_jitter, [latencies_pre(i), latencies_test(i)], 's--', ...
             'Color', [0.8, 0.3, 0.3], 'LineWidth', 1, ...
             'MarkerSize', 6, 'MarkerFaceColor', 'none', ...
             'MarkerEdgeColor', [0.8, 0.3, 0.3]);
    end
end

% Statistical test
[~, p_val] = ttest(latencies_pre, latencies_test);

% Calculate y-limits based on ALL data
all_data_points = [latencies_pre(:); latencies_test(:); (bar_data + bar_errors)'];
max_y = max(all_data_points) * 1.1;

% Add significance indicator
if p_val < 0.001
    sig_text = '***';
elseif p_val < 0.01
    sig_text = '**';
elseif p_val < 0.05
    sig_text = '*';
else
    sig_text = 'ns';
end

% Add significance bracket and text
bracket_y = max_y * 1.05;
plot([1, 2], [bracket_y, bracket_y], 'k-', 'LineWidth', 1.5);
plot([1, 1], [bracket_y-max_y*0.02, bracket_y], 'k-', 'LineWidth', 1.5);
plot([2, 2], [bracket_y-max_y*0.02, bracket_y], 'k-', 'LineWidth', 1.5);
text(1.5, bracket_y + max_y*0.05, sig_text, 'HorizontalAlignment', 'center', ...
     'FontSize', options.font_size+2, 'FontWeight', 'bold');

% Formatting
title(sprintf('Food Discovery Latency'), ...
      'FontSize', options.font_size+2, 'FontWeight', 'bold');
ylabel('Latency to Food Discovery (seconds)', 'FontSize', options.font_size);

% Set x-axis
set(gca, 'XTick', 1:2, 'XTickLabel', {'Before', 'Test'}, ...
         'FontSize', options.font_size);

% Set y-axis limits with padding
ylim([0, max_y * 1.2]);

% Clean up axes
box off;
set(gca, 'LineWidth', 1.5, 'TickDir', 'out');

discovery_rate_pre = sum(discovery_found_pre) / n_mice * 100;
discovery_rate_test = sum(discovery_found_test) / n_mice * 100;
% 
% % Add legend for marker types
% legend_x = 0.7;
% legend_y = 0.85;
% text(legend_x, legend_y, 'Legend:', 'Units', 'normalized', 'FontWeight', 'bold', ...
%      'FontSize', options.font_size-1);
% text(legend_x, legend_y-0.08, '● Food discovered', 'Units', 'normalized', ...
%      'FontSize', options.font_size-2);
% text(legend_x, legend_y-0.15, '□ Food missed', 'Units', 'normalized', ...
%      'FontSize', options.font_size-2);

% Print statistics to console
fprintf('\nLatency Results:\n');
fprintf('  Pre: %.1f±%.1f seconds (%.0f%% discovered)\n', ...
        latency_mean_pre, latency_sem_pre, discovery_rate_pre);
fprintf('  Test: %.1f±%.1f seconds (%.0f%% discovered)\n', ...
        latency_mean_test, latency_sem_test, discovery_rate_test);
fprintf('  p-value: %.4f\n', p_val);

%% Save figure if requested
if options.save_figure
    if ~exist(options.output_dir, 'dir')
        mkdir(options.output_dir);
    end
    
    % Create filename
    if use_time_limit
        filename = sprintf('%s_latency_%.0fmin', lower(options.group), time_limit_min);
    else
        filename = sprintf('%s_latency_full', lower(options.group));
    end
    
    % Save as both PDF and PNG
    pdf_file = fullfile(options.output_dir, [filename '.pdf']);
    png_file = fullfile(options.output_dir, [filename '.png']);
    
    % High-quality settings
    print(fig, pdf_file, '-dpdf', '-r300');
    print(fig, png_file, '-dpng', '-r300');
    
    fprintf('\nFigures saved:\n  PDF: %s\n  PNG: %s\n', pdf_file, png_file);
end

%% Store results
results = struct();
results.group = options.group;
results.group_name = group_name;
results.latencies_pre = latencies_pre;
results.latencies_test = latencies_test;
results.discovery_found_pre = discovery_found_pre;
results.discovery_found_test = discovery_found_test;
results.latency_mean_pre = latency_mean_pre;
results.latency_mean_test = latency_mean_test;
results.latency_sem_pre = latency_sem_pre;
results.latency_sem_test = latency_sem_test;
results.discovery_rate_pre = discovery_rate_pre;
results.discovery_rate_test = discovery_rate_test;
results.discovery_distance = options.discovery_distance;
results.selected_mice = selected_mice;
results.mouse_names = mouse_names(selected_mice);
results.n_mice = n_mice;
results.time_limit_used = use_time_limit;

if use_time_limit
    results.time_limit_min = time_limit_min;
    results.time_limit_sec = time_limit_sec;
end

% Add statistical results
results.statistics = struct(...
    'p_value', p_val, ...
    'significant', p_val < 0.05, ...
    'pre_mean', latency_mean_pre, ...
    'test_mean', latency_mean_test, ...
    'pre_sem', latency_sem_pre, ...
    'test_sem', latency_sem_test);

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS (%s GROUP) ===\n', upper(options.group));
if use_time_limit
    fprintf('Analysis period: First %.1f minutes of each session\n', time_limit_min);
else
    fprintf('Analysis period: Full session duration\n');
end
fprintf('Discovery threshold: %.1f cm from food\n', options.discovery_distance);
fprintf('Number of mice: %d\n\n', n_mice);

if results.statistics.significant
    sig_status = '(significant)';
else
    sig_status = '(n.s.)';
end

fprintf('Food Discovery Latency: Pre=%.1f±%.1f s, Test=%.1f±%.1f s, p=%.4f %s\n', ...
        results.statistics.pre_mean, results.statistics.pre_sem, ...
        results.statistics.test_mean, results.statistics.test_sem, ...
        results.statistics.p_value, sig_status);

fprintf('Discovery Success Rates: Pre=%.0f%% (%d/%d), Test=%.0f%% (%d/%d)\n', ...
        discovery_rate_pre, sum(discovery_found_pre), n_mice, ...
        discovery_rate_test, sum(discovery_found_test), n_mice);

fprintf('\nAnalysis completed!\n');

end