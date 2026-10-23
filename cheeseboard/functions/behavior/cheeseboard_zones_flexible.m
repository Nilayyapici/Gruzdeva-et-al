function results = cheeseboard_zones_flexible(mice, options)
% CHEESEBOARD_ZONES_FLEXIBLE - Paper-quality zone analysis for chosen group
%
% Usage:
%   options.group = 'saline';           % Choose 'saline' or 'CNO' (required)
%   options.time_limit = 10;            % Analyze only first 10 minutes (optional)
%   options.figure_size = [6, 4];       % Figure size in inches (optional)
%   options.font_size = 12;             % Font size for labels (optional)
%   options.save_figure = true;         % Save figure as PDF/PNG (optional)
%   options.output_dir = './figures';   % Output directory for saved figures (optional)
%
%   results = cheeseboard_zones_flexible(mice, options);
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
%     .figure_size - [width, height] in inches (default: [8, 6])
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
if ~isfield(options, 'figure_size'), options.figure_size = [8, 6]; end
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

%% Zone analysis
fprintf('\n=== ZONE ANALYSIS (%s GROUP) ===\n', upper(options.group));

% Initialize storage
zone_percent_pre = zeros(n_mice, 4);  % zones 0-3
zone_percent_test = zeros(n_mice, 4);
zone_names = {'Outside Arena', 'Food Zone', 'Area 2', 'Area 3'};
zone_labels = {'Outside', 'Food', 'Area 2', 'Area 3'}; % Shorter labels for plots

% Calculate zone percentages for each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    fprintf('Processing mouse %d/%d: %s\n', i, n_mice, mice{mouse_idx,1});

    % Get pre and test data
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};

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

    % Extract zone information
    pre_zones = pre_data_filtered(:, 8);
    test_zones = test_data_filtered(:, 8);

    % Calculate percentage time in each zone
    for zone = 0:3
        zone_percent_pre(i, zone+1) = sum(pre_zones == zone) / length(pre_zones) * 100;
        zone_percent_test(i, zone+1) = sum(test_zones == zone) / length(test_zones) * 100;
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
colors.gray = [0.7, 0.7, 0.7];     % Light gray for bars
colors.error = [0, 0, 0];          % Black for error bars

% Calculate means and SEMs
zone_means_pre = mean(zone_percent_pre, 1);
zone_means_test = mean(zone_percent_test, 1);
zone_sems_pre = std(zone_percent_pre, 0, 1) / sqrt(n_mice);
zone_sems_test = std(zone_percent_test, 0, 1) / sqrt(n_mice);

% Create subplot layout (1x4 for the four zones)
for zone = 1:4
    subplot(1, 4, zone);

    % Prepare data for this zone
    pre_vals = zone_percent_pre(:, zone);
    test_vals = zone_percent_test(:, zone);

    % Create grouped bar plot
    bar_data = [zone_means_pre(zone), zone_means_test(zone)];
    bar_errors = [zone_sems_pre(zone), zone_sems_test(zone)];

    % Create bar plot
    b = bar(1:2, bar_data, 'BarWidth', 0.7);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors.pre;
    b.CData(2,:) = colors.test;
    b.EdgeColor = 'black';
    b.LineWidth = 1;

    hold on;

    % Add error bars
    errorbar(1:2, bar_data, bar_errors, 'k.', 'LineWidth', 2, 'CapSize', 8);

    % Add individual data points with connecting lines
    x_positions = [1, 2];
    jitter_strength = 0;

    for i = 1:n_mice
        % Add small jitter to x-positions for visibility
        x_jitter = x_positions + (rand(1,2)-0.5) * jitter_strength;

        % Plot individual mouse data as connected points
        plot(x_jitter, [pre_vals(i), test_vals(i)], 'o-', ...
            'Color', [0.3, 0.3, 0.3], 'LineWidth', 1, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.6, 0.6, 0.6], ...
            'MarkerEdgeColor', 'black');
    end

    % Statistical test
    [~, p_val] = ttest(pre_vals, test_vals);

    % Calculate y-limits based on ALL data (bars, error bars, and individual points)
    all_data_points = [pre_vals(:); test_vals(:); (bar_data + bar_errors)'];
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
    plot([1, 2], [bracket_y, bracket_y], 'k-', 'LineWidth', 1);
    plot([1, 1], [bracket_y-max_y*0.02, bracket_y], 'k-', 'LineWidth', 1);
    plot([2, 2], [bracket_y-max_y*0.02, bracket_y], 'k-', 'LineWidth', 1);
    text(1.5, bracket_y + max_y*0.05, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', options.font_size, 'FontWeight', 'bold');

    % Formatting
    title(zone_labels{zone}, 'FontSize', options.font_size, 'FontWeight', 'bold');
    ylabel('Time (%)', 'FontSize', options.font_size);

    % Set x-axis
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Before', 'Test'}, ...
        'FontSize', options.font_size);

    % Set y-axis limits with padding to include all data points
    ylim([0, max_y * 1.2]);

    % Clean up axes
    box off;
    set(gca, 'LineWidth', 1, 'TickDir', 'out');

    % Print statistics to console
    fprintf('%s Zone: Pre=%.1f±%.1f%%, Test=%.1f±%.1f%%, p=%.4f\n', ...
        zone_names{zone}, zone_means_pre(zone), zone_sems_pre(zone), ...
        zone_means_test(zone), zone_sems_test(zone), p_val);
end

% Add overall figure title
if use_time_limit
    sgtitle(sprintf('%s Group Zone Analysis (First %.1f min, n=%d)', group_name, time_limit_min, n_mice), ...
        'FontSize', options.font_size+2, 'FontWeight', 'bold');
else
    sgtitle(sprintf('%s Group Zone Analysis (Full Session, n=%d)', group_name, n_mice), ...
        'FontSize', options.font_size+2, 'FontWeight', 'bold');
end

% Adjust spacing between subplots to make room for sgtitle and make them shorter
set(fig, 'Units', 'normalized');
subplots = findobj(fig, 'Type', 'axes');
for i = 1:length(subplots)
    pos = get(subplots(i), 'Position');
    % Make plots shorter by reducing height
    pos(4) = pos(4) * 0.8;  % Reduce height to 80% (make shorter)
    pos(2) = pos(2) + 0.05; % Move up slightly to center better
    set(subplots(i), 'Position', pos);
end

%% Store results
results = struct();
results.group = options.group;
results.group_name = group_name;
results.zone_percent_pre = zone_percent_pre;
results.zone_percent_test = zone_percent_test;
results.zone_means_pre = zone_means_pre;
results.zone_means_test = zone_means_test;
results.zone_sems_pre = zone_sems_pre;
results.zone_sems_test = zone_sems_test;
results.selected_mice = selected_mice;
results.mouse_names = mouse_names(selected_mice);
results.n_mice = n_mice;
results.time_limit_used = use_time_limit;

if use_time_limit
    results.time_limit_min = time_limit_min;
    results.time_limit_sec = time_limit_sec;
end

% Add statistical results
results.statistics = struct();
for zone = 1:4
    [~, p_val] = ttest(zone_percent_pre(:, zone), zone_percent_test(:, zone));
    results.statistics.(sprintf('zone_%d', zone-1)) = struct(...
        'p_value', p_val, ...
        'significant', p_val < 0.05, ...
        'pre_mean', zone_means_pre(zone), ...
        'test_mean', zone_means_test(zone), ...
        'pre_sem', zone_sems_pre(zone), ...
        'test_sem', zone_sems_test(zone));
end

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS (%s GROUP) ===\n', upper(options.group));
if use_time_limit
    fprintf('Analysis period: First %.1f minutes of each session\n', time_limit_min);
else
    fprintf('Analysis period: Full session duration\n');
end
fprintf('Number of mice: %d\n\n', n_mice);

for zone = 1:4
    stat = results.statistics.(sprintf('zone_%d', zone-1));
    if stat.significant
        sig_status = '(significant)';
    else
        sig_status = '(n.s.)';
    end

    fprintf('%s: Pre=%.1f±%.1f%%, Test=%.1f±%.1f%%, p=%.4f %s\n', ...
        zone_names{zone}, stat.pre_mean, stat.pre_sem, ...
        stat.test_mean, stat.test_sem, stat.p_value, sig_status);
end

%% Save bar plot data
save_data = struct();
save_data.zone_names = zone_names;
save_data.zone_labels = zone_labels;
save_data.group = options.group;
save_data.n_mice = n_mice;
save_data.mouse_names = mouse_names(selected_mice);

% Per-mouse data (rows = mice, cols = zones 0-3)
save_data.zone_percent_pre  = zone_percent_pre;
save_data.zone_percent_test = zone_percent_test;

% Summary stats
save_data.zone_means_pre  = zone_means_pre;
save_data.zone_means_test = zone_means_test;
save_data.zone_sems_pre   = zone_sems_pre;
save_data.zone_sems_test  = zone_sems_test;

% Save as .mat
if use_time_limit
    save_name = sprintf('zone_data_%s_%dmin', options.group, round(time_limit_min));
else
    save_name = sprintf('zone_data_%s_full', options.group);
end

if ~isfield(options, 'output_dir'), options.output_dir = './figures'; end
if ~exist(options.output_dir, 'dir'), mkdir(options.output_dir); end

save(fullfile(options.output_dir, [save_name '.mat']), 'save_data');
fprintf('Bar plot data saved to: %s\n', fullfile(options.output_dir, [save_name '.mat']));

% Save as Excel - all zones in one file, each zone on a separate sheet
xlsx_name = fullfile(options.output_dir, [save_name '.xlsx']);

for zone = 1:4
    T = table(mouse_names(selected_mice), ...
              zone_percent_pre(:, zone), ...
              zone_percent_test(:, zone), ...
              'VariableNames', {'Mouse', 'Pre_pct', 'Test_pct'});
    writetable(T, xlsx_name, 'Sheet', zone_labels{zone});
end

fprintf('Bar plot data saved to Excel: %s\n', xlsx_name);

fprintf('\nAnalysis completed!\n');

end