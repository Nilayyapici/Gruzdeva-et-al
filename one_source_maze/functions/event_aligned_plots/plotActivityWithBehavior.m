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
    %
    % Example usage:
    %   options = struct();
    %   options.session = 'sess1';  % Analyze only session 1
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
    for i = 1:size(mice_all, 1)
        condition = mice_all{i, 2};
        stimulus = mice_all{i, 3};
        data = mice_all{i, 4};
        
        field_name = [condition, '_', stimulus];
        
        % Analyze each behavior
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
            
            % Store data
            data_struct.(field_name).(behavior_name).count(end+1) = num_episodes;
            data_struct.(field_name).(behavior_name).duration(end+1) = mean_duration;
            data_struct.(field_name).(behavior_name).total_time(end+1) = total_time;
        end
    end
    
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
        
        % Create subplot for each behavior
        for b = 1:length(options.behaviors)
            behavior_name = options.behaviors{b};
            
            subplot(1, length(options.behaviors), b);
            
            % Extract data for plotting
            % Order: fasted_gel, fasted_food, fed_gel, fed_food
            means = zeros(2, 2); % [fasted, fed] x [gel, food]
            sems = zeros(2, 2);
            
            % Get data for each group
            fasted_gel = data_struct.fasted_gel.(behavior_name).(metric_name);
            fasted_food = data_struct.fasted_food.(behavior_name).(metric_name);
            fed_gel = data_struct.fed_gel.(behavior_name).(metric_name);
            fed_food = data_struct.fed_food.(behavior_name).(metric_name);
            
            % Calculate means and SEMs
            if ~isempty(fasted_gel)
                means(1, 1) = mean(fasted_gel);
                sems(1, 1) = std(fasted_gel) / sqrt(length(fasted_gel));
            end
            if ~isempty(fasted_food)
                means(1, 2) = mean(fasted_food);
                sems(1, 2) = std(fasted_food) / sqrt(length(fasted_food));
            end
            if ~isempty(fed_gel)
                means(2, 1) = mean(fed_gel);
                sems(2, 1) = std(fed_gel) / sqrt(length(fed_gel));
            end
            if ~isempty(fed_food)
                means(2, 2) = mean(fed_food);
                sems(2, 2) = std(fed_food) / sqrt(length(fed_food));
            end
            
            % Create grouped bar plot
            bar_handle = bar(means);
            hold on;
            
            % Set colors: gel = blue, food = red
            bar_handle(1).FaceColor = [0.3, 0.6, 0.9]; % gel - blue
            bar_handle(2).FaceColor = [0.9, 0.3, 0.3]; % food - red
            
            % Add error bars
            % Get the x-coordinates for the bars
            x = 1:2; % fasted, fed
            for i = 1:2 % gel, food
                x_offset = bar_handle(i).XEndPoints;
                errorbar(x_offset, means(:, i), sems(:, i), 'k.', 'LineWidth', 1.5);
            end
            
            % Add individual data points
            jitter_amount = 0.1;
            all_data = {fasted_gel, fasted_food, fed_gel, fed_food};
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
            set(gca, 'XTick', 1:2, 'XTickLabel', {'Fasted', 'Fed'});
            ylabel(metric_labels.(metric_name), 'FontSize', 11);
            title(behavior_labels.(behavior_name), 'FontSize', 12, 'FontWeight', 'bold');
            
            if b == length(options.behaviors)
                legend({'Gel', 'Food'}, 'Location', 'best');
            end
            
            grid on;
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