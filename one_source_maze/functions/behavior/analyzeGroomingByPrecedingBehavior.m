function analyzeGroomingByPrecedingBehavior(mice_all_reorganized, options)
% ANALYZEGROOMBYPRECEDINGBEHAVIOR Analyzes neural response to grooming based on preceding behavior
%
% Parameters:
%   mice_all_reorganized - Cell array with reorganized mice data
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - preceding_window: time window before grooming to check for behaviors [seconds] (default: 100)
%     - time_window: time window around grooming onset [before after] in seconds (default: [-20 60])
%     - z_score: whether to z-score the DFF signals (default: true)
%     - smooth_window: window size for smoothing in frames (default: 15)
%     - min_grooming_duration: minimum duration of grooming events (default: 2 seconds)
%     - max_grooming_duration: maximum duration of grooming events (default: 60 seconds)

% Constants
COL_TIME = 1;      % Time column
COL_DFF = 11;      % DFF data column
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

if ~isfield(options, 'time_window')
    options.time_window = [-20, 60]; % Default in seconds
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_grooming_duration')
    options.min_grooming_duration = 1; % Default 2 seconds
end

if ~isfield(options, 'max_grooming_duration')
    options.max_grooming_duration = 100000; % Default 60 seconds
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

% Initialize data structures to store events by category
category_events = cell(num_categories, 1);
for i = 1:num_categories
    category_events{i} = struct('dff', {}, 'time', {}, 'mouse_id', {}, 'duration', {});
end

% Process each mouse
fprintf('Analyzing grooming events with preceding behaviors...\n');
for m = 1:num_mice
    data = filtered_mice{m, 4};
    mouse_id = filtered_mice{m, 1};
    
    % Extract relevant columns
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    grooming = data(:, COL_GROOM);
    food_interaction = data(:, COL_FOOD);
    eating = data(:, COL_EAT);
    
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
            
            % Extract time window around grooming onset
            time_mask = time >= (onset_time + options.time_window(1)) & ...
                        time <= (onset_time + options.time_window(2));
            
            if any(time_mask)
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
                
                % Store event data
                new_event = struct('dff', aligned_dff, 'time', aligned_time, ...
                                    'mouse_id', mouse_id, 'duration', duration);
                category_events{category_idx}(end+1) = new_event;
            end
        end
    end
end

% Report the number of events in each category
for i = 1:num_categories
    fprintf('%s: %d events\n', categories{i}, length(category_events{i}));
end

% Check if we have data to analyze
if all(cellfun(@isempty, category_events))
    warning('No valid grooming events found with preceding behaviors');
    return;
end

% Create common time grid for visualization
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000);

% Interpolate all events onto the common time grid
category_avg_dff = cell(num_categories, 1);
category_sem_dff = cell(num_categories, 1);
unique_mice_per_category = zeros(num_categories, 1);

for c = 1:num_categories
    events = category_events{c};
    num_events = length(events);
    
    if num_events == 0
        category_avg_dff{c} = nan(size(time_grid));
        category_sem_dff{c} = nan(size(time_grid));
        continue;
    end
    
    % Count unique mice
    mouse_ids = {events.mouse_id};
    unique_mice_per_category(c) = length(unique(mouse_ids));
    
    % Interpolate each event onto the time grid
    all_dff = nan(num_events, length(time_grid));
    
    for e = 1:num_events
        event_time = events(e).time;
        event_dff = events(e).dff;
        
        if ~isempty(event_time) && ~isempty(event_dff)
            all_dff(e, :) = interp1(event_time, event_dff, time_grid, 'linear', NaN);
        end
    end
    
    % Calculate average and SEM
    category_avg_dff{c} = nanmean(all_dff, 1);
    category_sem_dff{c} = nanstd(all_dff, 0, 1) ./ sqrt(sum(~isnan(all_dff), 1));
end

% Create figure for visualization
figure('Position', [100, 100, 500, 800]);

% 1. Average DFF plot (top subplot)
subplot(2, 1, 1);
hold on;

% Define colors for categories
category_colors = [
    0.2, 0.6, 0.8;  % Food interaction: Blue
    0.8, 0.3, 0.3;  % Eating: Red
    0.3, 0.7, 0.3;  % No interaction: Green
];

% Add grooming onset line
xline(0, '--k', 'LineWidth', 1.5);

% Plot each category
legend_handles = [];
legend_labels = {};

for c = 1:num_categories
    if ~isempty(category_events{c})
        color = category_colors(c, :);
        
        % Plot SEM as shaded area
        ciplot(category_avg_dff{c} - category_sem_dff{c}, ...
               category_avg_dff{c} + category_sem_dff{c}, ...
               time_grid, color, 0.2);
           
        % Plot mean line
        h = plot(time_grid, category_avg_dff{c}, 'Color', color, 'LineWidth', 2);
        legend_handles = [legend_handles, h];
        
        % Create label with event count
        label = sprintf('%s (n=%d events, %d mice)', ...
                categories{c}, length(category_events{c}), unique_mice_per_category(c));
        legend_labels{end+1} = label;
    end
end

% Format plot
xlabel('Time from grooming onset (seconds)', 'FontSize', 12);
if options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
title('Neural Activity during Grooming by Preceding Behavior', 'FontSize', 14);
legend(legend_handles, legend_labels, 'Location', 'best');
legend('boxoff');
grid on;
box off;

% Set x-axis limits
xlim(options.time_window);

% 2. Comparison of pre/during grooming activity (bottom subplot)
subplot(2, 1, 2);
hold on;

% Define time windows for comparison
pre_groom_mask = time_grid < 0;
during_groom_mask = time_grid >= 0 & time_grid <= 20; % First 20s of grooming

% Calculate average activity in each period
pre_activity = nan(num_categories, 1);
during_activity = nan(num_categories, 1);
increase = nan(num_categories, 1);
sem_increase = nan(num_categories, 1);
event_counts = zeros(num_categories, 1);

for c = 1:num_categories
    events = category_events{c};
    num_events = length(events);
    event_counts(c) = num_events;
    
    if num_events == 0
        continue;
    end
    
    % Calculate pre and during activity for each event
    all_pre = nan(num_events, 1);
    all_during = nan(num_events, 1);
    all_increase = nan(num_events, 1);
    
    for e = 1:num_events
        event_time = events(e).time;
        event_dff = events(e).dff;
        
        pre_mask = event_time < 0;
        during_mask = event_time >= 0 & event_time <= 20;
        
        if any(pre_mask) && any(during_mask)
            pre_avg = nanmean(event_dff(pre_mask));
            during_avg = nanmean(event_dff(during_mask));
            
            all_pre(e) = pre_avg;
            all_during(e) = during_avg;
            all_increase(e) = during_avg - pre_avg;
        end
    end
    
    % Store averages and SEM
    pre_activity(c) = nanmean(all_pre);
    during_activity(c) = nanmean(all_during);
    increase(c) = nanmean(all_increase);
    sem_increase(c) = nanstd(all_increase) / sqrt(sum(~isnan(all_increase)));
    
    % Plot individual data points with jitter
    jitter_width = 0.2;
    jitter = (rand(size(all_increase))-0.5) * jitter_width;
    scatter(c + jitter, all_increase, 40, category_colors(c,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

% Plot bar chart of increases
bar(1:num_categories, increase, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Add error bars
errorbar(1:num_categories, increase, sem_increase, 'k.', 'LineWidth', 1.5, 'MarkerSize', 15);

% Add zero reference line
yline(0, '--k');

% Format plot
set(gca, 'XTick', 1:num_categories);
set(gca, 'XTickLabel', categories);
ylabel('Increase in ΔF/F during grooming', 'FontSize', 12);
title('Change in Neural Activity: During vs. Before Grooming', 'FontSize', 14);
grid on;
box on;

% Add event counts as text
for c = 1:num_categories
    if event_counts(c) > 0
        text(c, increase(c) + sem_increase(c) + 0.1, ...
             sprintf('n=%d (%d)', unique_mice_per_category(c), event_counts(c)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

% Perform ANOVA if we have data in multiple categories
valid_categories = find(~isnan(increase));
if length(valid_categories) > 1
    % Prepare data for ANOVA
    anova_data = [];
    anova_groups = [];
    
    for c = valid_categories'
        events = category_events{c};
        
        for e = 1:length(events)
            event_time = events(e).time;
            event_dff = events(e).dff;
            
            pre_mask = event_time < 0;
            during_mask = event_time >= 0 & event_time <= 20;
            
            if any(pre_mask) && any(during_mask)
                pre_avg = nanmean(event_dff(pre_mask));
                during_avg = nanmean(event_dff(during_mask));
                
                anova_data = [anova_data; during_avg - pre_avg];
                anova_groups = [anova_groups; c];
            end
        end
    end
    
    % Perform ANOVA
    [p, ~, stats] = anova1(anova_data, anova_groups, 'off');
    
    % Add ANOVA p-value as text
    max_y = max(increase + sem_increase) + 0.3;
    text(mean(1:num_categories), max_y, sprintf('p = %.4f', p), ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Add significance stars if p < 0.05
    if p < 0.05
        [c, ~, ~, ~] = multcompare(stats, 'Display', 'off');
        
        % Add stars for significant pairwise comparisons
        line_height = max_y + 0.2;
        y_increment = 0.15;
        
        for i = 1:size(c, 1)
            group1 = c(i,1);
            group2 = c(i,2);
            p_val = c(i,6);
            
            if p_val < 0.05
                % Determine star level
                if p_val < 0.001
                    stars = '***';
                elseif p_val < 0.01
                    stars = '**';
                else
                    stars = '*';
                end
                
                % Draw comparison line with stars
                x1 = group1;
                x2 = group2;
                y = line_height;
                
                % Draw the line
                plot([x1, x2], [y, y], 'k-', 'LineWidth', 1);
                
                % Add stars in the middle
                text(mean([x1, x2]), y + 0.05, stars, ...
                     'HorizontalAlignment', 'center', 'FontSize', 14);
                
                % Increment height for next comparison
                line_height = line_height + y_increment;
            end
        end
        
        % Adjust y-axis to show stars
        ylim_current = ylim;
        ylim([ylim_current(1), line_height + 0.1]);
    end
end

% Return data if requested
if nargout > 0
    results = struct();
    results.category_events = category_events;
    results.categories = categories;
    results.time_grid = time_grid;
    results.category_avg_dff = category_avg_dff;
    results.category_sem_dff = category_sem_dff;
    results.pre_activity = pre_activity;
    results.during_activity = during_activity;
    results.increase = increase;
    results.sem_increase = sem_increase;
    results.unique_mice_per_category = unique_mice_per_category;
    varargout{1} = results;
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