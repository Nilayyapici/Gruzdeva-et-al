function analyzeGroomingDurationByPrecedingBehavior(mice_all_reorganized, options)
% ANALYZEGROOMINGDURATIONBYPRECEDINGBEHAVIOR Analyzes how preceding behavior affects grooming duration
%
% Parameters:
%   mice_all_reorganized - Cell array with reorganized mice data
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - preceding_window: time window before grooming to check for behaviors [seconds] (default: 100)
%     - min_grooming_duration: minimum duration of grooming events (default: 2 seconds)
%     - max_grooming_duration: maximum duration of grooming events (default: 60 seconds)

% Constants
COL_TIME = 1;      % Time column
COL_GROOM = 10;    % Grooming column (0/1)
COL_FOOD = 8;      % Food interaction column
COL_EAT = 9;       % Eating column (0/1)

% Set default options
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'state')
    options.state = 'all';
end

if ~isfield(options, 'source')
    options.source = 'all';
end

if ~isfield(options, 'preceding_window')
    options.preceding_window = 100; % Look 100 seconds before grooming
end

if ~isfield(options, 'min_grooming_duration')
    options.min_grooming_duration = 2; % Default 2 seconds
end

if ~isfield(options, 'max_grooming_duration')
    options.max_grooming_duration = 60; % Default 60 seconds
end

% Filter mice based on state and source
n = size(mice_all_reorganized, 1);
mask = true(n, 1);

if ~strcmpi(options.state, 'all')
    mask = mask & strcmp(mice_all_reorganized(:, 2), options.state);
end

if ~strcmpi(options.source, 'all')
    mask = mask & strcmp(mice_all_reorganized(:, 3), options.source);
end

% Extract filtered mice
filtered_mice = mice_all_reorganized(mask, :);
num_mice = size(filtered_mice, 1);

if num_mice == 0
    warning('No mice found with the specified criteria');
    return;
end

% Categories for preceding behavior
categories = {'food interaction', 'eating', 'no interaction'};
num_categories = length(categories);

% Initialize data structures for each category
category_durations = cell(num_categories, 1);
category_mouse_ids = cell(num_categories, 1);
category_sessions = cell(num_categories, 1);

% Process each mouse
fprintf('Analyzing grooming durations with preceding behaviors...\n');
for m = 1:num_mice
    data = filtered_mice{m, 4};
    mouse_id = filtered_mice{m, 1};
    
    % Extract relevant columns
    time = data(:, COL_TIME);
    grooming = data(:, COL_GROOM);
    food_interaction = data(:, COL_FOOD);
    eating = data(:, COL_EAT);
    
    % Find session information
    if contains(mouse_id, '_sess')
        % Extract session from mouse_id
        session_parts = strsplit(mouse_id, '_sess');
        base_id = session_parts{1};
        session = str2double(session_parts{2});
    else
        base_id = mouse_id;
        session = NaN;
    end
    
    % Find onsets and offsets of grooming events
    groom_changes = diff([0; grooming; 0]); % Add padding and compute differences
    groom_onsets = find(groom_changes == 1);
    groom_offsets = find(groom_changes == -1) - 1;
    
    % Calculate durations in seconds
    groom_durations = zeros(length(groom_onsets), 1);
    for g = 1:length(groom_onsets)
        if g <= length(groom_offsets)
            onset_idx = groom_onsets(g);
            offset_idx = groom_offsets(g);
            
            if onset_idx <= length(time) && offset_idx <= length(time)
                groom_durations(g) = time(offset_idx) - time(onset_idx);
            end
        end
    end
    
    % Filter for valid durations
    valid_events = groom_durations >= options.min_grooming_duration & ...
                   groom_durations <= options.max_grooming_duration;
    groom_onsets = groom_onsets(valid_events);
    groom_offsets = groom_offsets(valid_events);
    groom_durations = groom_durations(valid_events);
    
    fprintf('Found %d valid grooming events for mouse %s\n', length(groom_onsets), mouse_id);
    
    % Process each grooming event
    for g = 1:length(groom_onsets)
        onset_idx = groom_onsets(g);
        onset_time = time(onset_idx);
        duration = groom_durations(g);
        
        % Check for preceding behaviors
        if onset_idx > 1
            % Find time window before grooming onset
            pre_time = onset_time - options.preceding_window;
            pre_indices = find(time >= pre_time & time < onset_time);
            
            % Check for different behaviors in the preceding window
            has_food_interaction = any(food_interaction(pre_indices) > 0);
            has_eating = any(eating(pre_indices) > 0);
            
            % Determine category
            if has_eating
                category_idx = 2; % Eating takes precedence
            elseif has_food_interaction
                category_idx = 1; % Food interaction but no eating
            else
                category_idx = 3; % No interaction
            end
            
            % Store data
            category_durations{category_idx}(end+1) = duration;
            category_mouse_ids{category_idx}{end+1} = base_id;
            category_sessions{category_idx}(end+1) = session;
        end
    end
end

% Check if we have data to analyze
if all(cellfun(@isempty, category_durations))
    warning('No valid grooming events found with preceding behaviors');
    return;
end

% Report the number of events and statistics for each category
fprintf('\nGrooming Duration Statistics by Preceding Behavior:\n');
fprintf('-------------------------------------------------\n');
for i = 1:num_categories
    durations = category_durations{i};
    if ~isempty(durations)
        unique_mice = unique(category_mouse_ids{i});
        fprintf('%s: %d events from %d mice\n', categories{i}, length(durations), length(unique_mice));
        fprintf('  Mean duration: %.2f ± %.2f seconds (SEM)\n', ...
                mean(durations), std(durations)/sqrt(length(durations)));
        fprintf('  Median duration: %.2f seconds\n', median(durations));
        fprintf('  Range: %.2f - %.2f seconds\n', min(durations), max(durations));
        fprintf('-------------------------------------------------\n');
    end
end

% Create figure for visualization
figure('Position', [100, 100, 1200, 800]);

% Define colors for categories
category_colors = [
    0.2, 0.6, 0.8;  % Food interaction: Blue
    0.8, 0.3, 0.3;  % Eating: Red
    0.3, 0.7, 0.3;  % No interaction: Green
];

% 1. Box plot of durations by category (top left)
subplot(2, 2, 1);
hold on;

% Prepare data for boxplot
boxplot_data = [];
boxplot_groups = [];
boxplot_labels = {};

for c = 1:num_categories
    durations = category_durations{c};
    if ~isempty(durations)
        boxplot_data = [boxplot_data; durations(:)];
        boxplot_groups = [boxplot_groups; c*ones(length(durations), 1)];
        boxplot_labels{end+1} = categories{c};
    end
end

% Create boxplot
boxplot(boxplot_data, boxplot_groups, 'Labels', boxplot_labels, 'Notch', 'on');

% Customize appearance
ylabel('Grooming Duration (seconds)', 'FontSize', 12);
title('Grooming Duration by Preceding Behavior', 'FontSize', 14);
grid on;

% 2. Bar plot with individual data points (top right)
subplot(2, 2, 2);
hold on;

% Calculate statistics
mean_durations = zeros(1, num_categories);
sem_durations = zeros(1, num_categories);
event_counts = zeros(1, num_categories);
mice_counts = zeros(1, num_categories);

for c = 1:num_categories
    durations = category_durations{c};
    if ~isempty(durations)
        mean_durations(c) = mean(durations);
        sem_durations(c) = std(durations) / sqrt(length(durations));
        event_counts(c) = length(durations);
        mice_counts(c) = length(unique(category_mouse_ids{c}));
    end
end

% Create bar plot
bar(1:num_categories, mean_durations, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Add error bars
errorbar(1:num_categories, mean_durations, sem_durations, 'k.', 'LineWidth', 1.5, 'MarkerSize', 15);

% Plot individual data points with jitter
for c = 1:num_categories
    durations = category_durations{c};
    if ~isempty(durations)
        jitter_width = 0.2;
        jitter = (rand(size(durations))-0.5) * jitter_width;
        scatter(c + jitter, durations, 40, category_colors(c,:), 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

% Format plot
set(gca, 'XTick', 1:num_categories);
set(gca, 'XTickLabel', categories);
ylabel('Grooming Duration (seconds)', 'FontSize', 12);
title('Mean Grooming Duration by Preceding Behavior', 'FontSize', 14);
grid on;
box on;

% Add event counts as text
for c = 1:num_categories
    if event_counts(c) > 0
        text(c, mean_durations(c) + sem_durations(c) + 2, ...
             sprintf('n=%d (%d)', mice_counts(c), event_counts(c)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

% 3. Violin plot (bottom left)
subplot(2, 2, 3);
hold on;

% Calculate kernel density for each category
x_points = linspace(0, options.max_grooming_duration, 100);
densities = cell(num_categories, 1);

for c = 1:num_categories
    durations = category_durations{c};
    if length(durations) > 5
        % Use kernel density estimation
        [density, xi] = ksdensity(durations, x_points);
        densities{c} = density;
    end
end

% Plot violin for each category
max_density = 0;
for c = 1:num_categories
    if ~isempty(densities{c})
        density = densities{c};
        density = density / max(density) * 0.3; % Scale for visualization
        max_density = max(max_density, max(density));
        
        % Plot violin
        patch([c-density, fliplr(c+density)], [x_points, fliplr(x_points)], ...
              category_colors(c,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        
        % Plot median line
        median_val = median(category_durations{c});
        plot([c-max_density, c+max_density], [median_val, median_val], 'k-', 'LineWidth', 1.5);
    end
end

% Format plot
set(gca, 'XTick', 1:num_categories);
set(gca, 'XTickLabel', categories);
xlabel('Preceding Behavior', 'FontSize', 12);
ylabel('Grooming Duration (seconds)', 'FontSize', 12);
title('Distribution of Grooming Durations', 'FontSize', 14);
ylim([0, options.max_grooming_duration]);
grid on;
box on;

% 4. Statistical comparison (bottom right)
subplot(2, 2, 4);

% Perform statistical tests
% ANOVA if more than 2 categories have data
valid_cats = find(event_counts > 0);
if length(valid_cats) > 2
    % Prepare data for ANOVA
    anova_data = [];
    anova_groups = [];
    
    for c = valid_cats
        durations = category_durations{c};
        anova_data = [anova_data; durations(:)];
        anova_groups = [anova_groups; c*ones(length(durations), 1)];
    end
    
    % Perform ANOVA
    [p, tbl, stats] = anova1(anova_data, anova_groups, 'off');
    
    % Run multiple comparisons
    [c, ~, ~, ~] = multcompare(stats, 'Display', 'off');
    
    % Create table of results
    anova_result = sprintf('ANOVA p = %.4f', p);
    if p < 0.05
        anova_result = [anova_result, ' *'];
    end
    
    % Format result text
    result_text = {anova_result, ''};
    for i = 1:size(c, 1)
        cat1 = categories{c(i,1)};
        cat2 = categories{c(i,2)};
        p_val = c(i,6);
        
        result_text{end+1} = sprintf('%s vs %s: p = %.4f', cat1, cat2, p_val);
        if p_val < 0.05
            result_text{end} = [result_text{end}, ' *'];
        end
    end
    
    % Display statistical results
    axis off;
    text(0.1, 0.9, 'Statistical Comparison:', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.1, 0.8, result_text, 'FontSize', 12);
    
elseif length(valid_cats) == 2
    % Perform t-test for two categories
    cat1 = valid_cats(1);
    cat2 = valid_cats(2);
    
    [h, p, ~, stats] = ttest2(category_durations{cat1}, category_durations{cat2});
    
    % Format result text
    if h == 1
        result = sprintf('t-test: %s vs %s\np = %.4f *\nt = %.2f, df = %d', ...
                        categories{cat1}, categories{cat2}, p, stats.tstat, stats.df);
    else
        result = sprintf('t-test: %s vs %s\np = %.4f (n.s.)\nt = %.2f, df = %d', ...
                        categories{cat1}, categories{cat2}, p, stats.tstat, stats.df);
    end
    
    % Display statistical results
    axis off;
    text(0.1, 0.9, 'Statistical Comparison:', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.1, 0.8, result, 'FontSize', 12);
    
else
    % Not enough categories for comparison
    axis off;
    text(0.1, 0.5, 'Not enough categories with data\nfor statistical comparison', ...
         'FontSize', 14, 'HorizontalAlignment', 'center');
end

% Add overall title with state and source information
if isfield(options, 'state') && isfield(options, 'source')
    sgtitle({['Grooming Duration by Preceding Behavior'], ...
             ['State: ', upper(options.state(1)), options.state(2:end), ...
              ', Source: ', upper(options.source(1)), options.source(2:end)]}, ...
             'FontWeight', 'bold');
else
    sgtitle('Grooming Duration by Preceding Behavior', 'FontWeight', 'bold');
end

% Return results if requested
if nargout > 0
    results = struct();
    results.categories = categories;
    results.category_durations = category_durations;
    results.category_mouse_ids = category_mouse_ids;
    results.category_sessions = category_sessions;
    results.mean_durations = mean_durations;
    results.sem_durations = sem_durations;
    results.event_counts = event_counts;
    results.mice_counts = mice_counts;
    varargout{1} = results;
end

end