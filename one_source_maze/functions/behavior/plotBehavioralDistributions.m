function plotBehavioralDistributions(mice_all, options)
    % Plot behavioral distributions with vertical shapes (violin-style)
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
    %   options.session = 'sess1';
    %   plotBehavioralDistributions(mice_all, options);
    
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
        options.same_ylim = false;
    end
    
    % Column definitions
    COL_TIME = 1;
    COL_FOOD_INT = 8;  % Food interaction
    COL_EATING = 9;    % Eating
    COL_GROOMING = 10; % Grooming
    
    % Filter data based on session
    if ~strcmp(options.session, 'all')
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
    
    % Process each mouse (group by unique mouse ID)
    mouse_data = struct();
    
    for i = 1:size(mice_all, 1)
        mouse_id_full = mice_all{i, 1};
        condition = mice_all{i, 2};
        stimulus = mice_all{i, 3};
        data = mice_all{i, 4};
        
        % Extract base mouse ID
        mouse_id_base = mouse_id_full;
        if contains(mouse_id_full, '_sess')
            mouse_id_base = extractBefore(mouse_id_full, '_sess');
        end
        
        field_name = [condition, '_', stimulus];
        mouse_key = [mouse_id_base, '_', field_name];
        
        if ~isfield(mouse_data, mouse_key)
            mouse_data.(mouse_key) = struct('data', [], 'condition', condition, ...
                                            'stimulus', stimulus, 'mouse_id', mouse_id_base);
        end
        
        mouse_data.(mouse_key).data = [mouse_data.(mouse_key).data; data];
    end
    
    % Process aggregated data
    for b_idx = 1:length(options.behaviors)
        behavior_name = options.behaviors{b_idx};
        
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
        
        mouse_keys = fieldnames(mouse_data);
        for k = 1:length(mouse_keys)
            mouse_key = mouse_keys{k};
            mouse_info = mouse_data.(mouse_key);
            
            data = mouse_info.data;
            condition = mouse_info.condition;
            stimulus = mouse_info.stimulus;
            field_name = [condition, '_', stimulus];
            
            episodes = find_episodes(data(:, col));
            
            if ~isempty(episodes)
                num_episodes = size(episodes, 1);
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
        
        % Calculate common y-axis limits if requested
        common_ylim = [];
        if options.same_ylim
            all_y_values = [];
            for b = 1:length(options.behaviors)
                behavior_name = options.behaviors{b};
                fasted_gel = data_struct.fasted_gel.(behavior_name).(metric_name);
                fasted_food = data_struct.fasted_food.(behavior_name).(metric_name);
                fed_gel = data_struct.fed_gel.(behavior_name).(metric_name);
                fed_food = data_struct.fed_food.(behavior_name).(metric_name);
                all_y_values = [all_y_values; fasted_gel(:); fasted_food(:); fed_gel(:); fed_food(:)];
            end
            y_max = max(all_y_values);
            y_min = max(0, min(all_y_values));  % Don't go below zero
            y_range = y_max - y_min;
            common_ylim = [0, y_max + 0.15*y_range];  % Start at 0
        end
        
        % Create subplots
        for b = 1:length(options.behaviors)
            behavior_name = options.behaviors{b};
            
            subplot(1, length(options.behaviors), b);
            
            % Get data for each group
            fed_gel = data_struct.fed_gel.(behavior_name).(metric_name);
            fed_food = data_struct.fed_food.(behavior_name).(metric_name);
            fasted_gel = data_struct.fasted_gel.(behavior_name).(metric_name);
            fasted_food = data_struct.fasted_food.(behavior_name).(metric_name);
            
            % Prepare data for plotting
            all_data = {fed_gel, fed_food, fasted_gel, fasted_food};
            group_labels = {'Fed-Gel', 'Fed-Food', 'Fasted-Gel', 'Fasted-Food'};
            colors = {[0.3, 0.6, 0.9], [0.9, 0.3, 0.3], [0.3, 0.6, 0.9], [0.9, 0.3, 0.3]};
            
            hold on;
            
            % Plot each group as a vertical distribution
            for g = 1:4
                if ~isempty(all_data{g})
                    plotDistribution(g, all_data{g}, colors{g});
                end
            end
            
            % Add horizontal lines for medians
            for g = 1:4
                if ~isempty(all_data{g})
                    med_val = median(all_data{g});
                    plot([g-0.35, g+0.35], [med_val, med_val], 'k-', 'LineWidth', 2);
                end
            end
            
            % Add baseline at y=0
            plot([0.5, 4.5], [0, 0], 'k-', 'LineWidth', 0.5);
            
            % Format plot
            xlim([0.5, 4.5]);
            set(gca, 'XTick', 1:4, 'XTickLabel', {'Fed\newlineGel', 'Fed\newlineFood', 'Fasted\newlineGel', 'Fasted\newlineFood'});
            ylabel(metric_labels.(metric_name), 'FontSize', 11);
            title(behavior_labels.(behavior_name), 'FontSize', 12, 'FontWeight', 'bold');
            
            if options.same_ylim && ~isempty(common_ylim)
                ylim(common_ylim);
            else
                % Auto-scale but ensure y-axis starts at 0
                current_ylim = ylim;
                ylim([0, current_ylim(2)]);
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

function plotDistribution(x_pos, data, color)
    % Plot a vertical distribution shape (simplified violin plot)
    
    if isempty(data) || all(isnan(data))
        return;
    end
    
    % Remove NaN values
    data = data(~isnan(data));
    
    if length(data) < 2
        % Just plot individual points if too few data points
        scatter(x_pos * ones(size(data)), data, 40, 'o', ...
               'MarkerEdgeColor', 'k', 'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.6);
        return;
    end
    
    % Calculate kernel density estimate
    [f, xi] = ksdensity(data, 'Function', 'pdf');
    
    % Truncate at zero (values can't be negative)
    valid_idx = xi >= 0;
    xi = xi(valid_idx);
    f = f(valid_idx);
    
    % If all values were below zero (shouldn't happen), return
    if isempty(xi)
        return;
    end
    
    % Add a point at zero if needed for clean boundary
    if xi(1) > 0
        xi = [0; xi(:)];
        f = [0; f(:)];
    end
    
    % Normalize the density for width (max width = 0.35)
    max_width = 0.35;
    f_norm = f / max(f) * max_width;
    
    % Create violin shape (mirrored distribution)
    x_violin = [x_pos - f_norm, fliplr(x_pos + f_norm)];
    y_violin = [xi, fliplr(xi)];
    
    % Plot the violin
    patch(x_violin, y_violin, color, 'FaceAlpha', 0.5, 'EdgeColor', color, 'LineWidth', 1);
    
    % Plot individual data points with jitter
    jitter_amount = 0.08;
    x_jitter = x_pos + (rand(size(data)) - 0.5) * jitter_amount;
    scatter(x_jitter, data, 20, 'o', 'MarkerEdgeColor', 'k', ...
           'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.6);
end

function episodes = find_episodes(binary_vector)
    % Find continuous episodes where binary_vector is 1
    episodes = [];
    
    if isempty(binary_vector)
        return;
    end
    
    transitions = diff([0; binary_vector(:); 0]);
    starts = find(transitions == 1);
    ends = find(transitions == -1) - 1;
    
    if ~isempty(starts) && ~isempty(ends)
        episodes = [starts, ends];
    end
end