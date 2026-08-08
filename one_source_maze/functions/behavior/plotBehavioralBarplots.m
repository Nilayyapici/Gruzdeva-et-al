function plotBehavioralBarplots(mice_all, options)
% Plot behavioral analysis barplots for eating, food visits, and grooming
%
% Inputs:
%   mice_all: cell array with mouse data (with session info in column 1)
%   options: struct with fields (all optional):
%     - behaviors: cell array of behaviors to analyze
%                  default: {'eating', 'food_visits', 'grooming'}
%     - metrics: cell array of metrics to plot
%                default: {'count', 'duration', 'total_time'}
%     - session: which session to analyze
%                ('sess0', 'sess1', 'sess2', 'sess3', or 'all')
%                default: 'all'
%     - same_ylim: if true, use same y-axis limits for all subplots
%                  default: false (auto-scale each subplot)
%
% Example usage:
%   options = struct();
%   options.session = 'sess1';  % Analyze only session 1
%   options.same_ylim = true;   % Use same y-axis for comparison
%   plotBehavioralBarplots(mice_all, options);

% Set defaults
if ~isfield(options, 'behaviors')
    options.behaviors = {'eating', 'food_visits', 'grooming'};
end
if ~isfield(options, 'metrics')
    options.metrics = {'count', 'duration', 'total_time'};
end
if ~isfield(options, 'session')
    options.session = 'all';
end
if ~isfield(options, 'same_ylim')
    options.same_ylim = false;  % Default: different y-axes for each subplot
end

% Column definitions
COL_TIME = 1;
COL_FOOD_INT = 8;  % Food interaction
COL_EATING = 9;    % Eating
COL_GROOMING = 10; % Grooming

% Filter data based on session
if ~strcmp(options.session, 'all')
    % Filter to only include the specified session
    session_mask = false(size(mice_all, 1), 1);
    for i = 1:size(mice_all, 1)
        mouse_id = mice_all{i, 1};
        if contains(mouse_id, ['_', options.session])
            session_mask(i) = true;
        end
    end
    mice_all = mice_all(session_mask, :);

    fprintf('Filtered to %d entries for session: %s\n', size(mice_all, 1), options.session);
end

% Initialize data structure
data_struct = struct();
conditions = {'fasted', 'fed'};
stimuli = {'food', 'gel'};

for cond = conditions
    for stim = stimuli
        field_name = [cond{1}, '_', stim{1}];
        data_struct.(field_name) = struct();
        for b_idx = 1:length(options.behaviors)
            behavior_name = options.behaviors{b_idx};
            data_struct.(field_name).(behavior_name) = struct(...
                'count', [], ...
                'duration', [], ...
                'total_time', []);
        end
    end
end

% Process each mouse
% First, group data by unique mouse ID (without session suffix)
mouse_data = struct();

for i = 1:size(mice_all, 1)
    mouse_id_full = mice_all{i, 1};
    condition = mice_all{i, 2};
    stimulus = mice_all{i, 3};
    data = mice_all{i, 4};

    % Extract base mouse ID (remove session suffix if present)
    mouse_id_base = mouse_id_full;
    if contains(mouse_id_full, '_sess')
        mouse_id_base = extractBefore(mouse_id_full, '_sess');
    end

    % Create unique key for this mouse
    field_name = [condition, '_', stimulus];
    mouse_key = [mouse_id_base, '_', field_name];

    % Initialize if this is a new mouse
    if ~isfield(mouse_data, mouse_key)
        mouse_data.(mouse_key) = struct('data', [], 'condition', condition, ...
            'stimulus', stimulus, 'mouse_id', mouse_id_base);
    end

    % Append data from this session
    mouse_data.(mouse_key).data = [mouse_data.(mouse_key).data; data];
end

% Now process aggregated data for each unique mouse
for b_idx = 1:length(options.behaviors)
    behavior_name = options.behaviors{b_idx};

    % Get the appropriate column
    switch behavior_name
        case 'eating'
            col = COL_EATING;
        case 'food_visits'
            col = COL_FOOD_INT;
        case 'grooming'
            col = COL_GROOMING;
        otherwise
            continue;
    end

    % Process each unique mouse
    mouse_keys = fieldnames(mouse_data);
    for k = 1:length(mouse_keys)
        mouse_key = mouse_keys{k};
        mouse_info = mouse_data.(mouse_key);

        data = mouse_info.data;
        condition = mouse_info.condition;
        stimulus = mouse_info.stimulus;
        field_name = [condition, '_', stimulus];

        % Find episodes
        episodes = find_episodes(data(:, col));

        if ~isempty(episodes)
            % Calculate metrics for this mouse
            num_episodes = size(episodes, 1);

            % Calculate durations (in seconds)
            durations = zeros(num_episodes, 1);
            for ep = 1:num_episodes
                start_idx = episodes(ep, 1);
                end_idx = episodes(ep, 2);
                durations(ep) = data(end_idx, COL_TIME) - data(start_idx, COL_TIME);
            end

            mean_duration = mean(durations);
            total_time = sum(durations);
        else
            num_episodes = 0;
            mean_duration = 0;
            total_time = 0;
        end

        % Store data (one point per unique mouse)
        data_struct.(field_name).(behavior_name).count(end+1) = num_episodes;
        data_struct.(field_name).(behavior_name).duration(end+1) = mean_duration;
        data_struct.(field_name).(behavior_name).total_time(end+1) = total_time;
    end
end

% ── Export individual mouse data to Excel ──────────────────────────────
% Determine session label for filename
if strcmp(options.session, 'all')
    sess_label = 'all_sessions';
else
    sess_label = options.session;
end

excel_filename = sprintf('behavioral_data_%s_%s.xlsx', sess_label, datestr(now, 'yyyymmdd'));

% Build one table per behavior x metric combination
for b_idx = 1:length(options.behaviors)
    behavior_name = options.behaviors{b_idx};

    for m_idx = 1:length(options.metrics)
        metric_name = options.metrics{m_idx};

        % Get data for each group
        fed_gel    = data_struct.fed_gel.(behavior_name).(metric_name)(:);
        fed_food   = data_struct.fed_food.(behavior_name).(metric_name)(:);
        fasted_gel  = data_struct.fasted_gel.(behavior_name).(metric_name)(:);
        fasted_food = data_struct.fasted_food.(behavior_name).(metric_name)(:);

        % Pad to equal length with NaN
        max_n = max([length(fed_gel), length(fed_food), ...
            length(fasted_gel), length(fasted_food)]);

        pad = @(x) [x(:); NaN(max_n - length(x), 1)];

        export_table = table(pad(fed_gel), pad(fed_food), pad(fasted_gel), pad(fasted_food), ...
            'VariableNames', {'Fed_Gel', 'Fed_Food', 'Fasted_Gel', 'Fasted_Food'});

        sheet_name = sprintf('%s_%s', behavior_name, metric_name);
        % Excel sheet names max 31 chars
        if length(sheet_name) > 31
            sheet_name = sheet_name(1:31);
        end

        writetable(export_table, excel_filename, 'Sheet', sheet_name);
    end
end

fprintf('Behavioral data exported to: %s\n', excel_filename);

% Create plots for each behavior and metric combination
behavior_labels = struct(...
    'eating', 'Eating', ...
    'food_visits', 'Food Visits', ...
    'grooming', 'Grooming');

metric_labels = struct(...
    'count', 'Number of Episodes', ...
    'duration', 'Mean Episode Duration (s)', ...
    'total_time', 'Total Time (s)');

% Plot each metric
for m = 1:length(options.metrics)
    metric_name = options.metrics{m};

    figure('Position', [100, 100, 1200, 400]);

    % Calculate common y-axis limits if requested
    common_ylim = [];
    if options.same_ylim
        % First pass: collect all y-values to determine common y-axis limits
        all_y_values = [];

        for b = 1:length(options.behaviors)
            behavior_name = options.behaviors{b};

            % Get data for each group
            fasted_gel = data_struct.fasted_gel.(behavior_name).(metric_name);
            fasted_food = data_struct.fasted_food.(behavior_name).(metric_name);
            fed_gel = data_struct.fed_gel.(behavior_name).(metric_name);
            fed_food = data_struct.fed_food.(behavior_name).(metric_name);

            % Collect all values
            all_y_values = [all_y_values; fasted_gel(:); fasted_food(:); fed_gel(:); fed_food(:)];
        end

        % Calculate common y-axis limits with some padding
        y_max = max(all_y_values);
        y_min = min(all_y_values);
        y_range = y_max - y_min;
        common_ylim = [min(0, y_min - 0.1*y_range), y_max + 0.15*y_range];
    end

    % Create subplots
    for b = 1:length(options.behaviors)
        behavior_name = options.behaviors{b};

        subplot(1, length(options.behaviors), b);

        % Extract data for plotting
        % Order: fed_gel, fed_food, fasted_gel, fasted_food
        means = zeros(2, 2); % [fed, fasted] x [gel, food]
        sems = zeros(2, 2);

        % Get data for each group
        fasted_gel = data_struct.fasted_gel.(behavior_name).(metric_name);
        fasted_food = data_struct.fasted_food.(behavior_name).(metric_name);
        fed_gel = data_struct.fed_gel.(behavior_name).(metric_name);
        fed_food = data_struct.fed_food.(behavior_name).(metric_name);

        % Calculate means and SEMs - Fed first, Fasted second
        if ~isempty(fed_gel)
            means(1, 1) = mean(fed_gel);
            sems(1, 1) = std(fed_gel) / sqrt(length(fed_gel));
        end
        if ~isempty(fed_food)
            means(1, 2) = mean(fed_food);
            sems(1, 2) = std(fed_food) / sqrt(length(fed_food));
        end
        if ~isempty(fasted_gel)
            means(2, 1) = mean(fasted_gel);
            sems(2, 1) = std(fasted_gel) / sqrt(length(fasted_gel));
        end
        if ~isempty(fasted_food)
            means(2, 2) = mean(fasted_food);
            sems(2, 2) = std(fasted_food) / sqrt(length(fasted_food));
        end

        % Create grouped bar plot
        bar_handle = bar(means);
        hold on;

        % Set colors: gel = blue, food = red
        bar_handle(1).FaceColor = [0.3, 0.6, 0.9]; % gel - blue
        bar_handle(2).FaceColor = [0.9, 0.3, 0.3]; % food - red

        % Add error bars
        % Get the x-coordinates for the bars
        x = 1:2; % fed, fasted
        for i = 1:2 % gel, food
            x_offset = bar_handle(i).XEndPoints;
            errorbar(x_offset, means(:, i), sems(:, i), 'k.', 'LineWidth', 1.5);
        end

        % Add individual data points
        jitter_amount = 0.1;
        all_data = {fed_gel, fed_food, fasted_gel, fasted_food};
        colors = {[0.3, 0.6, 0.9], [0.9, 0.3, 0.3], [0.3, 0.6, 0.9], [0.9, 0.3, 0.3]};
        positions = [bar_handle(1).XEndPoints(1), bar_handle(2).XEndPoints(1), ...
            bar_handle(1).XEndPoints(2), bar_handle(2).XEndPoints(2)];

        for i = 1:4
            if ~isempty(all_data{i})
                x_jitter = positions(i) + (rand(size(all_data{i})) - 0.5) * jitter_amount;
                scatter(x_jitter, all_data{i}, 20, 'o', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', colors{i}, ...
                    'MarkerFaceAlpha', 0.6);
            end
        end

        % Format plot
        set(gca, 'XTick', 1:2, 'XTickLabel', {'Fed', 'Fasted'});
        ylabel(metric_labels.(metric_name), 'FontSize', 11);
        title(behavior_labels.(behavior_name), 'FontSize', 12, 'FontWeight', 'bold');

        % Set common y-axis limits if requested
        if options.same_ylim && ~isempty(common_ylim)
            ylim(common_ylim);
        end

        % Perform statistical tests
        % Prepare data for analysis
        y_data = [];
        group_labels = {};
        group_condition = {};
        group_stimulus = {};

        % Fed-Gel (group 1)
        if ~isempty(fed_gel)
            y_data = [y_data; fed_gel(:)];
            group_labels = [group_labels; repmat({'Fed_Gel'}, length(fed_gel), 1)];
            group_condition = [group_condition; repmat({'Fed'}, length(fed_gel), 1)];
            group_stimulus = [group_stimulus; repmat({'Gel'}, length(fed_gel), 1)];
        end

        % Fed-Food (group 2)
        if ~isempty(fed_food)
            y_data = [y_data; fed_food(:)];
            group_labels = [group_labels; repmat({'Fed_Food'}, length(fed_food), 1)];
            group_condition = [group_condition; repmat({'Fed'}, length(fed_food), 1)];
            group_stimulus = [group_stimulus; repmat({'Food'}, length(fed_food), 1)];
        end

        % Fasted-Gel (group 3)
        if ~isempty(fasted_gel)
            y_data = [y_data; fasted_gel(:)];
            group_labels = [group_labels; repmat({'Fasted_Gel'}, length(fasted_gel), 1)];
            group_condition = [group_condition; repmat({'Fasted'}, length(fasted_gel), 1)];
            group_stimulus = [group_stimulus; repmat({'Gel'}, length(fasted_gel), 1)];
        end

        % Fasted-Food (group 4)
        if ~isempty(fasted_food)
            y_data = [y_data; fasted_food(:)];
            group_labels = [group_labels; repmat({'Fasted_Food'}, length(fasted_food), 1)];
            group_condition = [group_condition; repmat({'Fasted'}, length(fasted_food), 1)];
            group_stimulus = [group_stimulus; repmat({'Food'}, length(fasted_food), 1)];
        end

        % Run two-way ANOVA if we have enough data
        if length(y_data) >= 4 && length(unique(group_condition)) > 1 && length(unique(group_stimulus)) > 1
            [p_values, tbl, stats] = anovan(y_data, {group_condition, group_stimulus}, ...
                'model', 'interaction', 'display', 'off', 'varnames', {'Condition', 'Stimulus'});

            % Extract p-values
            p_condition = p_values(1);  % Main effect of Condition (Fed vs Fasted)
            p_stimulus = p_values(2);   % Main effect of Stimulus (Gel vs Food)
            p_interaction = p_values(3); % Interaction effect

            % Format p-values
            if p_condition < 0.001
                cond_str = 'p<0.001***';
            elseif p_condition < 0.01
                cond_str = 'p<0.01**';
            elseif p_condition < 0.05
                cond_str = 'p<0.05*';
            else
                cond_str = sprintf('p=%.3f', p_condition);
            end

            if p_stimulus < 0.001
                stim_str = 'p<0.001***';
            elseif p_stimulus < 0.01
                stim_str = 'p<0.01**';
            elseif p_stimulus < 0.05
                stim_str = 'p<0.05*';
            else
                stim_str = sprintf('p=%.3f', p_stimulus);
            end

            if p_interaction < 0.001
                int_str = 'p<0.001***';
            elseif p_interaction < 0.01
                int_str = 'p<0.01**';
            elseif p_interaction < 0.05
                int_str = 'p<0.05*';
            else
                int_str = sprintf('p=%.3f', p_interaction);
            end

            % Add ANOVA text annotation
            anova_text = sprintf('2-way ANOVA:\nState: %s\nSource: %s\nInt: %s', cond_str, stim_str, int_str);
            text(0.02, 0.98, anova_text, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontSize', 7, ...
                'BackgroundColor', 'white', 'EdgeColor', 'k');
        end

        % Run Kruskal-Wallis if we have enough data
        if length(y_data) >= 4 && length(unique(group_labels)) >= 2
            [p_kw, tbl, stats] = kruskalwallis(y_data, group_labels, 'off');

            % If significant, perform post-hoc pairwise comparisons
            if p_kw < 0.05
                % Perform pairwise Mann-Whitney U tests (uncorrected)
                group_names = unique(group_labels, 'stable');
                n_groups = length(group_names);

                % Get current y-axis limits for placing significance bars
                y_lim = ylim;
                y_range = y_lim(2) - y_lim(1);

                % Map group names to bar positions
                % Position mapping: Fed(1)-Gel, Fed(1)-Food, Fasted(2)-Gel, Fasted(2)-Food
                group_positions = struct();
                group_positions.Fed_Gel = bar_handle(1).XEndPoints(1);
                group_positions.Fed_Food = bar_handle(2).XEndPoints(1);
                group_positions.Fasted_Gel = bar_handle(1).XEndPoints(2);
                group_positions.Fasted_Food = bar_handle(2).XEndPoints(2);

                % Track bar heights to avoid overlapping
                bar_y_start = y_lim(2) + 0.02 * y_range;
                bar_y_increment = 0.08 * y_range;
                bar_count = 0;

                % Perform all pairwise comparisons
                for i = 1:n_groups-1
                    for j = i+1:n_groups
                        group1 = group_names{i};
                        group2 = group_names{j};

                        % Get data for these groups
                        data1 = y_data(strcmp(group_labels, group1));
                        data2 = y_data(strcmp(group_labels, group2));

                        % Mann-Whitney U test (ranksum in MATLAB)
                        [p_pair, ~] = ranksum(data1, data2);

                        % If significant (uncorrected p < 0.05), add bar
                        if p_pair < 0.05
                            x1 = group_positions.(group1);
                            x2 = group_positions.(group2);
                            y_bar = bar_y_start + bar_count * bar_y_increment;

                            % Draw the bar
                            plot([x1, x1, x2, x2], [y_bar-0.01*y_range, y_bar, y_bar, y_bar-0.01*y_range], ...
                                'k-', 'LineWidth', 1.5);

                            % Add star(s)
                            if p_pair < 0.001
                                star_text = '***';
                            elseif p_pair < 0.01
                                star_text = '**';
                            else
                                star_text = '*';
                            end

                            text(mean([x1, x2]), y_bar + 0.01*y_range, star_text, ...
                                'HorizontalAlignment', 'center', 'FontSize', 10);

                            bar_count = bar_count + 1;
                        end
                    end
                end

                % Adjust y-axis to show all bars
                if bar_count > 0
                    new_y_max = bar_y_start + bar_count * bar_y_increment + 0.05 * y_range;
                    ylim([y_lim(1), new_y_max]);
                end
            end
        end

        if b == length(options.behaviors)
            legend({'Gel', 'Food'}, 'Location', 'best');
        end

        grid off;
        box off;
        hold off;
    end

    % Add overall title with session info
    if strcmp(options.session, 'all')
        title_str = metric_labels.(metric_name);
    else
        title_str = sprintf('%s - %s', metric_labels.(metric_name), upper(options.session));
    end
    sgtitle(title_str, 'FontSize', 14, 'FontWeight', 'bold');
end
end

function episodes = find_episodes(binary_vector)
% Find continuous episodes where binary_vector is 1
% Returns: Nx2 matrix where each row is [start_index, end_index]

episodes = [];

if isempty(binary_vector)
    return;
end

% Find transitions
transitions = diff([0; binary_vector(:); 0]);
starts = find(transitions == 1);
ends = find(transitions == -1) - 1;

% Combine into episodes
if ~isempty(starts) && ~isempty(ends)
    episodes = [starts, ends];
end
end