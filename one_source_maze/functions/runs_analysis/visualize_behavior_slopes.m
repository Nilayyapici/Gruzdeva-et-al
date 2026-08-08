function visualize_behavior_slopes(slope_results, options)
    % visualize_behavior_slopes - Creates visualizations comparing eating vs visit slopes
    %
    % Inputs:
    %   slope_results - Output from analyze_behavior_slopes
    %   options - Visualization options:
    %       .plot_types - Cell array: {'barplot', 'scatter', 'heatmap'}
    %       .separate_directions - If true, separate plots for towards/away (default: true)
    %       .save_plots - Whether to save plots (default: false)
    %       .figure_path - Path to save figures (default: './')
    
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'plot_types')
        options.plot_types = {'barplot'};
    end
    if ~isfield(options, 'separate_directions'), options.separate_directions = true; end
    if ~isfield(options, 'save_plots'), options.save_plots = false; end
    if ~isfield(options, 'figure_path'), options.figure_path = './'; end
    
    % Extract data
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No slope analysis results to visualize.\n');
        return;
    end
    
    % Parse results
    [slope_data, sessions, directions, behaviors] = parse_behavior_slope_results(slope_results);
    
    % Create requested visualizations
    for i = 1:length(options.plot_types)
        plot_type = options.plot_types{i};
        
        switch plot_type
            case 'barplot'
                if options.separate_directions
                    create_behavior_barplot_separate(slope_data, sessions, directions, behaviors, slope_results, options);
                else
                    create_behavior_barplot_combined(slope_data, sessions, directions, behaviors, slope_results, options);
                end
            case 'scatter'
                create_behavior_scatter(slope_data, sessions, directions, behaviors, slope_results, options);
            case 'heatmap'
                create_behavior_heatmap(slope_data, sessions, directions, behaviors, options);
            otherwise
                warning('Unknown plot type: %s', plot_type);
        end
    end
end

function [slope_data, sessions, directions, behaviors] = parse_behavior_slope_results(slope_results)
    % Parse slope results into organized matrices
    
    detailed_fields = fieldnames(slope_results.detailed);
    
    sessions = [];
    directions = {};
    behaviors = {};
    
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        parts = split(field, '_');
        sess_num = str2double(parts{1}(5:end));
        direction = parts{2};
        behavior = parts{3};
        
        if ~ismember(sess_num, sessions)
            sessions = [sessions, sess_num];
        end
        if ~ismember(direction, directions)
            directions{end+1} = direction;
        end
        if ~ismember(behavior, behaviors)
            behaviors{end+1} = behavior;
        end
    end
    
    sessions = sort(sessions);
    
    % Initialize data structure
    slope_data = struct();
    for d = 1:length(directions)
        for b = 1:length(behaviors)
            key = sprintf('%s_%s', directions{d}, behaviors{b});
            slope_data.(key).slopes = nan(length(sessions), 1);
            slope_data.(key).p_values = nan(length(sessions), 1);
            slope_data.(key).r_squared = nan(length(sessions), 1);
            slope_data.(key).n_mice = nan(length(sessions), 1);
            slope_data.(key).is_significant = false(length(sessions), 1);
        end
    end
    
    % Fill data
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        
        parts = split(field, '_');
        sess_num = str2double(parts{1}(5:end));
        direction = parts{2};
        behavior = parts{3};
        
        s_idx = find(sessions == sess_num);
        key = sprintf('%s_%s', direction, behavior);
        
        if isfield(slope_data, key)
            slope_data.(key).slopes(s_idx) = stats.slope;
            slope_data.(key).p_values(s_idx) = stats.p_value;
            slope_data.(key).r_squared(s_idx) = stats.r_squared;
            slope_data.(key).n_mice(s_idx) = stats.n_mice;
            slope_data.(key).is_significant(s_idx) = stats.is_significant;
        end
    end
end

function create_behavior_barplot_separate(slope_data, sessions, directions, behaviors, slope_results, options)
    % Create separate bar plots for each direction, comparing eating vs visit
    
    direction_labels = {'Towards', 'Away'};
    behavior_labels = {'Eating', 'Visit'};
    
    % Colors for behaviors
    eating_color = [0.8, 0.2, 0.2];  % Dark red
    visit_color = [0.2, 0.4, 0.8];   % Dark blue
    behavior_colors = {eating_color, visit_color};
    
    % Calculate global y-limits first
    detailed_fields = fieldnames(slope_results.detailed);
    global_individual_data = [];
    global_mean_values = [];
    global_sem_values = [];
    
    for i = 1:length(detailed_fields)
        stats = slope_results.detailed.(detailed_fields{i});
        if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
            global_individual_data = [global_individual_data; stats.mouse_slopes];
            global_mean_values = [global_mean_values; mean(stats.mouse_slopes)];
            global_sem_values = [global_sem_values; std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes))];
        else
            global_mean_values = [global_mean_values; stats.slope];
            global_sem_values = [global_sem_values; stats.se_slope];
        end
    end
    
    % Calculate global y-range
    if ~isempty(global_individual_data)
        global_max = max([global_mean_values + global_sem_values; global_individual_data]);
        global_min = min([global_mean_values - global_sem_values; global_individual_data]);
    else
        global_max = max(global_mean_values + global_sem_values, [], 'omitnan');
        global_min = min(global_mean_values - global_sem_values, [], 'omitnan');
    end
    
    global_y_range = global_max - global_min;
    if global_y_range == 0 || isnan(global_y_range)
        global_y_range = 0.1;
    end
    
    global_bottom_margin = global_y_range * 0.15;
    global_top_margin = global_y_range * 0.35;
    global_ylim = [global_min - global_bottom_margin, global_max + global_top_margin];
    
    % Create one plot per direction
    for d = 1:length(directions)
        direction = directions{d};
        direction_label = direction_labels{d};
        
        figure('Name', sprintf('Behavior Slopes: %s Direction', direction_label), ...
               'Position', [100 + (d-1)*700, 100, 800, 600]);
        
        % Organize by session (x-axis) with behavior bars
        n_sessions = length(sessions);
        n_behaviors = length(behaviors);
        
        bar_width = 0.35;
        x_positions = 1:n_sessions;
        x_offset = [-0.2, 0.2];
        
        hold on;
        
        % Collect data for this direction
        individual_data = cell(n_behaviors, n_sessions);
        mouse_ids_data = cell(n_behaviors, n_sessions);
        mean_values = zeros(n_behaviors, n_sessions);
        sem_values = zeros(n_behaviors, n_sessions);
        
        for s = 1:n_sessions
            for b = 1:n_behaviors
                sess = sessions(s);
                behavior = behaviors{b};
                field_name = sprintf('sess%d_%s_%s', sess, direction, behavior);
                
                if isfield(slope_results.detailed, field_name)
                    stats = slope_results.detailed.(field_name);
                    
                    if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                        individual_data{b, s} = stats.mouse_slopes;
                        mean_values(b, s) = mean(stats.mouse_slopes);
                        sem_values(b, s) = std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes));
                        
                        % Try to get mouse IDs
                        if isfield(stats, 'mouse_ids')
                            mouse_ids_data{b, s} = stats.mouse_ids;
                        else
                            mouse_ids_data{b, s} = arrayfun(@(x) sprintf('M%d', x), ...
                                1:length(stats.mouse_slopes), 'UniformOutput', false);
                        end
                    else
                        mean_values(b, s) = stats.slope;
                        sem_values(b, s) = stats.se_slope;
                        individual_data{b, s} = [];
                        mouse_ids_data{b, s} = {};
                    end
                else
                    mean_values(b, s) = NaN;
                    sem_values(b, s) = NaN;
                    individual_data{b, s} = [];
                    mouse_ids_data{b, s} = {};
                end
            end
        end
        
        % Store mouse positions for connecting lines
        mouse_positions = cell(n_behaviors, n_sessions);
        
        % Plot bars
        for s = 1:n_sessions
            for b = 1:n_behaviors
                x_pos = x_positions(s) + x_offset(b);
                
                % Plot bar
                bar(x_pos, mean_values(b, s), bar_width, ...
                    'FaceColor', behavior_colors{b}, 'EdgeColor', behavior_colors{b} * 0.8, ...
                    'LineWidth', 1.5, 'FaceAlpha', 0.8);
                
                % Add error bars
                errorbar(x_pos, mean_values(b, s), sem_values(b, s), ...
                        'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');
                
                % Add individual points
                if ~isempty(individual_data{b, s})
                    mouse_slopes = individual_data{b, s};
                    mouse_names = mouse_ids_data{b, s};
                    n_mice = length(mouse_slopes);
                    
                    x_jitter = repmat(x_pos, n_mice, 1);
                    
                    mouse_positions{b, s} = struct('x_pos', x_jitter, 'y_pos', mouse_slopes, 'names', {mouse_names});
                    
                    scatter(x_jitter, mouse_slopes, 50, behavior_colors{b}, 'o', ...
                           'filled', 'MarkerFaceAlpha', 0.7, ...
                           'MarkerEdgeColor', behavior_colors{b} * 0.7, 'LineWidth', 1.5);
                end
            end
        end
        
        % Connect same mice across sessions for each behavior
        for b = 1:n_behaviors
            if n_sessions == 2 && ~isempty(mouse_positions{b, 1}) && ~isempty(mouse_positions{b, 2})
                sess0_data = mouse_positions{b, 1};
                sess1_data = mouse_positions{b, 2};
                
                for i = 1:length(sess0_data.names)
                    mouse_name = sess0_data.names{i};
                    sess1_idx = find(strcmp(sess1_data.names, mouse_name));
                    
                    if ~isempty(sess1_idx)
                        x_coords = [sess0_data.x_pos(i), sess1_data.x_pos(sess1_idx)];
                        y_coords = [sess0_data.y_pos(i), sess1_data.y_pos(sess1_idx)];
                        
                        plot(x_coords, y_coords, '-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1.0);
                    end
                end
            end
        end
        
        % Add significance indicators for individual slopes
        for s = 1:n_sessions
            for b = 1:n_behaviors
                sess = sessions(s);
                behavior = behaviors{b};
                key = sprintf('%s_%s', direction, behavior);
                
                if isfield(slope_data, key) && slope_data.(key).is_significant(s)
                    x_pos = x_positions(s) + x_offset(b);
                    y_pos = mean_values(b, s) + sem_values(b, s);
                    y_star = y_pos + global_y_range * 0.3;
                    
                    text(x_pos, y_star, '★', 'FontSize', 12, ...
                         'HorizontalAlignment', 'center', 'Color', 'k', ...
                         'FontName', 'Arial Unicode MS');
                end
            end
        end
        
        % Formatting
        xlabel('Sessions', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Slope (z-scored dF/F per distance unit)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('%s Direction: Eating vs Visit', direction_label), ...
              'FontSize', 18, 'FontWeight', 'bold');
        
        % Session labels
        session_labels = cell(1, n_sessions);
        for s = 1:n_sessions
            sess = sessions(s);
            if sess == 0
                session_labels{s} = 'Sess 0';
            elseif sess == 1
                session_labels{s} = 'Sess 1';
            else
                session_labels{s} = sprintf('Sess %d', sess);
            end
        end
        
        set(gca, 'XTick', x_positions, 'XTickLabel', session_labels, 'FontSize', 14);
        
        % Add zero line
        line([0.5, n_sessions + 0.5], [0, 0], 'Color', [0.3, 0.3, 0.3], ...
             'LineStyle', '--', 'LineWidth', 1.5);
        
        xlim([0.5, n_sessions + 0.5]);
        ylim(global_ylim);
        
        % Add sample sizes
        for s = 1:n_sessions
            for b = 1:n_behaviors
                if ~isempty(individual_data{b, s})
                    n_mice = length(individual_data{b, s});
                    x_pos = x_positions(s) + x_offset(b);
                    y_pos = global_min - global_y_range * 0.08;
                    
                    text(x_pos, y_pos, sprintf('n=%d', n_mice), ...
                         'FontSize', 10, 'HorizontalAlignment', 'center', ...
                         'Color', behavior_colors{b}, 'FontWeight', 'bold');
                end
            end
        end
        
        % Legend
        legend_handles = [];
        for b = 1:n_behaviors
            h = scatter(NaN, NaN, 100, behavior_colors{b}, 'o', 'filled', ...
                       'MarkerEdgeColor', behavior_colors{b} * 0.7, 'LineWidth', 1.5);
            legend_handles(b) = h;
        end
        legend(legend_handles, behavior_labels, 'Location', 'best', ...
               'FontSize', 14, 'Box', 'off');
        
        set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
        
        hold off;
        
        if options.save_plots
            filename = sprintf('behavior_slope_barplot_%s.png', lower(direction));
            saveas(gcf, fullfile(options.figure_path, filename));
        end
    end
end

function create_behavior_barplot_combined(slope_data, sessions, directions, behaviors, slope_results, options)
    % Create combined bar plot with both directions
    
    fprintf('Combined bar plot not yet implemented. Use separate_directions = true.\n');
end

function create_behavior_scatter(slope_data, sessions, directions, behaviors, slope_results, options)
    % Create scatter plots
    
    fprintf('Scatter plot not yet implemented.\n');
end

function create_behavior_heatmap(slope_data, sessions, directions, behaviors, options)
    % Create heatmap
    
    fprintf('Heatmap not yet implemented.\n');
end