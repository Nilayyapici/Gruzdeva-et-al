function visualize_combined_nonfood_slopes(slope_results, options)
    % visualize_combined_nonfood_slopes - Creates comprehensive visualizations of slope analysis with combined non-food arms
    %
    % Inputs:
    %   slope_results - Output from analyze_combined_nonfood_slopes function
    %   options - Visualization options (optional):
    %       .plot_types - Cell array of plot types to create (default: all)
    %                    {'heatmap', 'barplot', 'scatter', 'summary'}
    %       .group_by - How to group bars in barplot: 'session' (default) or 'arm'
    %       .save_plots - Whether to save plots (default: false)
    %       .figure_path - Path to save figures (default: current directory)
    %       .session_colors - Cell array of colors for sessions {Before, Learning, Test}
    %                        Default: {[0.6 0.6 0.6], [0.2 0.4 0.8], [0.1 0.7 0.9]}
    
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'plot_types')
        options.plot_types = {'heatmap','scatter'};
    end
    if ~isfield(options, 'save_plots'), options.save_plots = false; end
    if ~isfield(options, 'figure_path'), options.figure_path = './'; end
    if ~isfield(options, 'group_by'), options.group_by = 'session'; end % 'session' or 'arm'
    
    % Set default session colors if not provided
    if ~isfield(options, 'session_colors') || isempty(options.session_colors)
        % Default colors: Before (gray), Learning (blue), Test (light blue)
        options.session_colors = {[0.6 0.6 0.6], [0.2 0.4 0.8], [0.1 0.7 0.9]};
    end
    
    % Extract data from slope_results
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No slope analysis results to visualize.\n');
        return;
    end
    
    % Parse the results into organized data
    [slope_data, sessions, arms, directions] = parse_combined_slope_results(slope_results);
    
    % Create requested visualizations
    for i = 1:length(options.plot_types)
        plot_type = options.plot_types{i};
        
        switch plot_type
            case 'heatmap'
                create_combined_slope_heatmap(slope_data, sessions, arms, directions, options);
            case 'barplot'
                create_combined_slope_barplot(slope_data, sessions, arms, directions, slope_results, options);
            case 'scatter'
                create_combined_slope_scatter(slope_data, sessions, arms, directions, options);
            case 'summary'
                create_combined_summary_plot(slope_results, options);
            otherwise
                warning('Unknown plot type: %s', plot_type);
        end
    end
end

function [slope_data, sessions, arms, directions] = parse_combined_slope_results(slope_results)
    % Parse slope results into organized matrices for combined non-food arms
    
    detailed_fields = fieldnames(slope_results.detailed);
    
    % Extract unique sessions, arms, directions
    sessions = [];
    arms = {};
    directions = {};
    
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        parts = split(field, '_');
        sess_num = str2double(parts{1}(5:end));
        direction = parts{2};
        arm = parts{3};
        
        if ~ismember(sess_num, sessions)
            sessions = [sessions, sess_num];
        end
        if ~ismember(direction, directions)
            directions{end+1} = direction;
        end
        if ~ismember(arm, arms)
            arms{end+1} = arm;
        end
    end
    
    sessions = sort(sessions);
    directions = {'towards', 'away'}; % Force specific order
    arms = {'food', 'nonfood'}; % Should be only these two for combined analysis
    
    % Initialize data matrices
    slope_data = struct();
    slope_data.slopes = nan(length(sessions), length(directions), length(arms));
    slope_data.p_values = nan(length(sessions), length(directions), length(arms));
    slope_data.r_squared = nan(length(sessions), length(directions), length(arms));
    slope_data.n_mice = nan(length(sessions), length(directions), length(arms));
    slope_data.is_significant = false(length(sessions), length(directions), length(arms));
    slope_data.correct_direction = false(length(sessions), length(directions), length(arms));
    
    % Fill matrices
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        
        parts = split(field, '_');
        sess_num = str2double(parts{1}(5:end));
        direction = parts{2};
        arm = parts{3};
        
        s_idx = find(sessions == sess_num);
        d_idx = find(strcmp(directions, direction));
        a_idx = find(strcmp(arms, arm));
        
        if ~isempty(s_idx) && ~isempty(d_idx) && ~isempty(a_idx)
            slope_data.slopes(s_idx, d_idx, a_idx) = stats.slope;
            slope_data.p_values(s_idx, d_idx, a_idx) = stats.p_value;
            slope_data.r_squared(s_idx, d_idx, a_idx) = stats.r_squared;
            slope_data.n_mice(s_idx, d_idx, a_idx) = stats.n_mice;
            slope_data.is_significant(s_idx, d_idx, a_idx) = stats.is_significant;
            slope_data.correct_direction(s_idx, d_idx, a_idx) = stats.correct_direction;
        end
    end
end

function create_combined_slope_heatmap(slope_data, sessions, arms, directions, options)
    % Create heatmap showing slopes across conditions with combined non-food arms
    
    % Calculate figure dimensions for transposed heatmaps (sessions as columns, arms as rows)
    n_sessions = length(sessions);
    n_arms = length(arms); % Should be 2: food and nonfood
    
    % Calculate appropriate figure size for transposed layout
    cell_size = 100; % pixels per cell
    heatmap_height = n_arms * cell_size + 150; % Arms as rows
    heatmap_width = n_sessions * cell_size + 150; % Sessions as columns
    total_width = heatmap_width * 2 + 100; % two heatmaps plus spacing
    
    figure('Name', 'Slope Analysis Heatmap (Combined Non-Food)', 'Position', [100, 100, total_width, heatmap_height]);
    
    % Calculate global color limits for consistency
    all_slopes = slope_data.slopes(:);
    max_slope = max(abs(all_slopes), [], 'omitnan');
    if isnan(max_slope) || max_slope == 0
        color_limit = 0.01;
    else
        color_limit = max(max_slope, 0.005);
    end
    
    % Create session labels
    session_labels = cell(1, length(sessions));
    for i = 1:length(sessions)
        sess_num = sessions(i);
        switch sess_num
            case 0
                session_labels{i} = 'Before';
            case 1
                session_labels{i} = 'Learning';
            case 2
                session_labels{i} = 'Test';
            otherwise
                session_labels{i} = sprintf('Session %d', sess_num);
        end
    end
    
    % Create arm labels with food first
    arm_labels = cell(1, length(arms));
    arm_order = {}; % Reorder to put food first
    
    % Find food and nonfood indices
    food_idx = find(strcmp(arms, 'food'));
    nonfood_idx = find(strcmp(arms, 'nonfood'));
    
    if ~isempty(food_idx) && ~isempty(nonfood_idx)
        arm_order = {arms{food_idx}, arms{nonfood_idx}}; % Food first, then nonfood
        arm_labels = {'Food', 'Non-Food'};
    else
        % Fallback if arms don't match expected names
        arm_order = arms;
        for i = 1:length(arms)
            if strcmp(arms{i}, 'food')
                arm_labels{i} = 'Food';
            else
                arm_labels{i} = 'Non-Food';
            end
        end
    end
    
    for d = 1:length(directions)
        direction = directions{d};
        subplot(1, 2, d);
        
        % Get slope data for this direction and reorder arms (food first)
        slope_matrix = squeeze(slope_data.slopes(:, d, :)); % sessions x arms
        
        % Reorder arms: food first, then nonfood
        if ~isempty(food_idx) && ~isempty(nonfood_idx)
            slope_matrix_reordered = slope_matrix(:, [food_idx, nonfood_idx]);
        else
            slope_matrix_reordered = slope_matrix;
        end
        
        % Transpose: sessions become columns (x-axis), arms become rows (y-axis)
        slope_matrix_transposed = slope_matrix_reordered'; % Now: arms x sessions
        
        % Create heatmap with 3 decimal place formatting
        % X-axis: sessions, Y-axis: arms (food first)
        h = heatmap(session_labels, arm_labels, slope_matrix_transposed, ...
                   'Colormap', redblue(256), 'ColorLimits', [-color_limit, color_limit], ...
                   'CellLabelFormat', '%.3f');
        
        h.Title = sprintf('%s Runs - Slope Values', direction);
        h.XLabel = 'Sessions';
        h.YLabel = 'Arms';
        
        % Only show colorbar for the "away" direction
        if ~strcmp(direction, 'away')
            h.ColorbarVisible = 'off';
        else
            h.ColorbarVisible = 'on';
        end
    end
    
    sgtitle('Slope Analysis: Neural Activity vs Distance to Food (Combined Non-Food)', 'FontSize', 14, 'FontWeight', 'bold');
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'combined_slope_heatmap.png'));
    end
end

function create_combined_slope_scatter(slope_data, sessions, arms, directions, options)
    % Create scatter plot: Towards vs Away slopes for combined non-food analysis
    
    figure('Name', 'Slope Analysis: Towards vs Away Comparison (Combined Non-Food)', 'Position', [100, 100, 1200, 800]);
    
    % Define colors by SESSION and markers by ARM
    session_colors = options.session_colors;
    arm_markers = {'o', 's'}; % Circle (food), Square (non-food)
    arm_display_names = {'Food Arm', 'Non-Food Arms'};
    
    % Find direction indices
    towards_idx = find(strcmp(directions, 'towards'));
    away_idx = find(strcmp(directions, 'away'));
    
    if isempty(towards_idx) || isempty(away_idx)
        error('Both "towards" and "away" directions must be present in the data');
    end
    
    % Create main scatter plot
    subplot(1, 4, [1, 2, 3]); % Main plot takes 3/4 of width
    hold on;
    
    % Plot data points
    all_x_vals = [];
    all_y_vals = [];
    
    for a = 1:length(arms)
        for s = 1:length(sessions)
            x_val = slope_data.slopes(s, towards_idx, a); % Towards slope
            y_val = slope_data.slopes(s, away_idx, a);    % Away slope
            
            if ~isnan(x_val) && ~isnan(y_val)
                all_x_vals = [all_x_vals, x_val];
                all_y_vals = [all_y_vals, y_val];
                
                % Get color based on SESSION using options.session_colors
                sess_num = sessions(s);
                color_idx = sess_num + 1; % sess0->1, sess1->2, sess2->3
                
                if color_idx <= length(options.session_colors)
                    point_color = options.session_colors{color_idx};
                else
                    point_color = [0.5, 0.5, 0.5]; % Default gray for unknown sessions
                end
                
                % Check significance for both directions
                towards_sig = slope_data.is_significant(s, towards_idx, a);
                away_sig = slope_data.is_significant(s, away_idx, a);
                
                % Determine marker properties based on significance
                if towards_sig && away_sig
                    marker_size = 180;
                    edge_color = point_color;
                    line_width = 4;
                    face_alpha = 0.9;
                elseif towards_sig || away_sig
                    marker_size = 120;
                    edge_color = point_color;
                    line_width = 2;
                    face_alpha = 0.7;
                else
                    marker_size = 80;
                    edge_color = point_color;
                    line_width = 1;
                    face_alpha = 0.5;
                end
                
                % Plot the point
                h = scatter(x_val, y_val, marker_size, point_color, arm_markers{a}, ...
                           'filled', 'MarkerEdgeColor', edge_color, 'LineWidth', line_width, ...
                           'MarkerFaceAlpha', face_alpha);
                
                % Add session number as text label
                if towards_sig || away_sig
                    session_label = sprintf('S%d*', sess_num);
                else
                    session_label = sprintf('S%d', sess_num);
                end
                
                % Position label to avoid overlap
                text(x_val - max(abs(all_x_vals))*0.04, y_val + max(abs(all_y_vals))*0.04, ...
                     session_label, 'FontSize', 12, 'HorizontalAlignment', 'center', ...
                     'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.2]);
            end
        end
    end
    
    % Set axis limits with padding
    if ~isempty(all_x_vals) && ~isempty(all_y_vals)
        x_range = max(abs(all_x_vals));
        y_range = max(abs(all_y_vals));
        max_range = max(x_range, y_range) * 1.2;
        xlim([-max_range, max_range]);
        ylim([-max_range, max_range]);
    end
    
    % Add reference lines
    xl = xlim;
    yl = ylim;
    
    % Axes through origin
    line(xl, [0, 0], 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
    line([0, 0], yl, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
    
    % Diagonal line (y = x)
    diag_lim = min(max(xl), max(yl));
    line([-diag_lim, diag_lim], [-diag_lim, diag_lim], 'Color', [0.6, 0.6, 0.6], ...
         'LineStyle', ':', 'LineWidth', 1);
    
    % Expected pattern line (y = -x)  
    line([-diag_lim, diag_lim], [diag_lim, -diag_lim], 'Color', [0.8, 0.4, 0.4], ...
         'LineStyle', '-.', 'LineWidth', 2);
    
    % Add subtle quadrant shading for expected pattern
    fill([xl(1), 0, 0, xl(1)], [0, 0, yl(2), yl(2)], [0.9, 1, 0.9], ...
         'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % Formatting
    xlabel('Towards Food Slope (z-score/cm)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Away from Food Slope (z-score/cm)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Neural Activity Slopes: Towards vs Away Food (Combined Non-Food)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Create dual legend: sessions and arms
    dummy_x = xl(1) - 1000; % Plot outside visible area
    dummy_y = yl(1) - 1000;
    
    % Session legend
    session_legend_handles = [];
    session_legend_labels = {};
    session_names = {'Before', 'Learning', 'Test'};
    
    for s = 1:length(sessions)
        sess_num = sessions(s);
        color_idx = sess_num + 1; % sess0->1, sess1->2, sess2->3
        
        if color_idx <= length(options.session_colors)
            h_sess = scatter(dummy_x, dummy_y, 120, options.session_colors{color_idx}, 'o', 'filled');
            session_legend_handles(end+1) = h_sess;
            if color_idx <= length(session_names)
                session_legend_labels{end+1} = session_names{color_idx};
            else
                session_legend_labels{end+1} = sprintf('Session %d', sess_num);
            end
        end
    end
    
    % Arm legend  
    arm_legend_handles = [];
    arm_legend_labels = {};
    for a = 1:length(arms)
        % Create empty (unfilled) markers with black edges
        h_arm = scatter(dummy_x, dummy_y, 120, [1, 1, 1], arm_markers{a}, ...
                       'MarkerEdgeColor', [0, 0, 0], 'LineWidth', 2);
        arm_legend_handles(end+1) = h_arm;
        arm_legend_labels{end+1} = arm_display_names{a};
    end
    
    % Create combined legend
    combined_handles = [session_legend_handles, arm_legend_handles];
    combined_labels = [session_legend_labels, arm_legend_labels];
    
    legend(combined_handles, combined_labels, 'Location', 'northwest', 'FontSize', 14, ...
           'Box', 'off', 'NumColumns', 2);
    
    % Add quadrant interpretation with two-line labels
    text(xl(2)*0.95, yl(2)*0.95, {'Towards↑', 'Away↑'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'Margin', 2);
    text(xl(1)*0.95, yl(2)*0.95, {'Towards↓', 'Away↑'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [0.9, 1, 0.9], 'Margin', 2);
    text(xl(1)*0.95, yl(1)*0.95, {'Towards↓', 'Away↓'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
         'BackgroundColor', 'white', 'Margin', 2);
    text(xl(2)*0.95, yl(1)*0.95, {'Towards↑', 'Away↓'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
         'BackgroundColor', 'white','Margin', 2);
    
    % Create interpretation panel (right side)
    subplot(1, 4, 4);
    axis off;
    
    % Title
    text(0.05, 0.95, 'INTERPRETATION GUIDE', 'FontSize', 14, 'FontWeight', 'bold', ...
         'Color', [0.2, 0.2, 0.2]);
    
    % Quadrant explanations
    text(0.05, 0.88, 'QUADRANTS:', 'FontSize', 12, 'FontWeight', 'bold');
    
    text(0.05, 0.82, 'I (↑↑):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.8]);
    text(0.05, 0.79, 'Activity increases both', 'FontSize', 10);
    text(0.05, 0.76, 'approaching & leaving food', 'FontSize', 10);
    
    text(0.05, 0.70, 'II (↓↑):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0.7, 0]);
    text(0.05, 0.67, 'Activity decreases approaching,', 'FontSize', 10);
    text(0.05, 0.64, 'increases leaving food', 'FontSize', 10);
    text(0.05, 0.61, '(EXPECTED PATTERN)', 'FontSize', 9, 'FontWeight', 'bold', ...
         'Color', [0, 0.7, 0]);
    
    text(0.05, 0.55, 'III (↓↓):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.6, 0, 0.6]);
    text(0.05, 0.52, 'Activity decreases both', 'FontSize', 10);
    text(0.05, 0.49, 'approaching & leaving food', 'FontSize', 10);
    
    text(0.05, 0.43, 'IV (↑↓):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.6, 0, 0.6]);
    text(0.05, 0.40, 'Activity increases approaching,', 'FontSize', 10);
    text(0.05, 0.37, 'decreases leaving food', 'FontSize', 10);
    
    % Visual elements guide
    text(0.05, 0.30, 'VISUAL ELEMENTS:', 'FontSize', 12, 'FontWeight', 'bold');
    
    text(0.05, 0.26, 'Colors (Sessions):', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Display session colors dynamically based on options.session_colors
    session_names = {'Before', 'Learning', 'Test'};
    y_pos = [0.23, 0.20, 0.17];
    
    for s = 1:min(length(sessions), 3) % Show up to 3 sessions
        sess_num = sessions(s);
        color_idx = sess_num + 1;
        
        if color_idx <= length(options.session_colors) && s <= length(y_pos)
            if color_idx <= length(session_names)
                label_text = sprintf('• %s = %s', session_names{color_idx}, session_names{color_idx});
            else
                label_text = sprintf('• Session %d', sess_num);
            end
            
            text(0.05, y_pos(s), label_text, 'FontSize', 10, 'Color', options.session_colors{color_idx});
        end
    end
    
    text(0.05, 0.13, 'Shapes (Arms):', 'FontSize', 11, 'FontWeight', 'bold');
    text(0.05, 0.10, '• Circle = Food arm', 'FontSize', 10);
    text(0.05, 0.07, '• Square = Non-food arms', 'FontSize', 10);
    
    text(0.05, 0.03, 'Size = Significance, * = Significant', 'FontSize', 10, 'FontWeight', 'bold');
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'combined_slope_scatter.png'));
    end
end

function create_combined_slope_barplot(slope_data, sessions, arms, directions, slope_results, options)
    % Create separate bar plots for towards and away directions
    % Can group by session (default) or by arm (alternative view)
    
    % Check grouping option
    if ~isfield(options, 'group_by')
        options.group_by = 'session'; % default
    end
    
    if strcmp(options.group_by, 'arm')
        % NEW: Group by arm (Food, Non-Food), bars are sessions (Before, Learning, Test)
        create_combined_slope_barplot_by_arm(slope_data, sessions, arms, directions, slope_results, options);
    else
        % ORIGINAL: Group by session (Before, Learning, Test), bars are arms (Food, Non-Food)
        create_combined_slope_barplot_by_session(slope_data, sessions, arms, directions, slope_results, options);
    end
end

function create_combined_summary_plot(slope_results, options)
    % Create summary visualization
    fprintf('Summary plot visualization not implemented for combined non-food analysis.\n');
end

% ========== HELPER FUNCTIONS FOR BARPLOT ==========

function create_combined_slope_barplot_by_session(slope_data, sessions, arms, directions, slope_results, options)
    % Original function: Groups by session, bars are arms
    
    direction_labels = {'Towards', 'Away'};
    towards_colors = {[0.8, 0.5, 0.2], [0.4, 0.7, 0.3]}; % Orange (food), Green (non-food)
    away_colors = {[0.8, 0.5, 0.2], [0.4, 0.7, 0.3]};     % Same colors for consistency
    
    direction_color_sets = {towards_colors, away_colors};
    
    % CALCULATE GLOBAL Y-LIMITS FIRST (for consistent scaling across both plots)
    detailed_fields = fieldnames(slope_results.detailed);
    global_individual_data = [];
    global_mean_values = [];
    global_sem_values = [];
    
    % Collect all data across all directions to determine global y-range
    for d = 1:length(directions)
        for s = 1:length(sessions)
            for a = 1:length(arms)
                sess_num = sessions(s);
                direction = directions{d};
                arm = arms{a};
                field_name = sprintf('sess%d_%s_%s', sess_num, direction, arm);
                
                if any(strcmp(detailed_fields, field_name))
                    stats = slope_results.detailed.(field_name);
                    
                    if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                        global_individual_data = [global_individual_data; stats.mouse_slopes];
                        global_mean_values = [global_mean_values; mean(stats.mouse_slopes)];
                        global_sem_values = [global_sem_values; std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes))];
                    else
                        global_mean_values = [global_mean_values; stats.slope];
                        global_sem_values = [global_sem_values; stats.se_slope];
                    end
                end
            end
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
    
    % Apply margins for global limits
    global_bottom_margin = global_y_range * 0.15;
    global_top_margin = global_y_range * 0.35;
    global_ylim = [global_min - global_bottom_margin, global_max + global_top_margin];
    
    for d = 1:length(directions)
        direction = directions{d};
        direction_label = direction_labels{d};
        colors = direction_color_sets{d};
        
        figure('Name', sprintf('Combined Non-Food Slope Analysis - %s Direction (Grouped by Session)', direction_label), ...
               'Position', [100 + (d-1)*800, 100, 700, 700]);
        
        % Create session labels
        session_labels = cell(1, length(sessions));
        for i = 1:length(sessions)
            sess_num = sessions(i);
            switch sess_num
                case 0
                    session_labels{i} = 'Before';
                case 1
                    session_labels{i} = 'Learning';
                case 2
                    session_labels{i} = 'Test';
                otherwise
                    session_labels{i} = sprintf('Session %d', sess_num);
            end
        end
        
        % Arm labels
        arm_labels = {'Food', 'Non-Food'};
        
        % Get individual mouse data for this direction
        individual_data = cell(length(sessions), length(arms));
        mouse_ids = cell(length(sessions), length(arms));
        mean_values = zeros(length(sessions), length(arms));
        sem_values = zeros(length(sessions), length(arms));
        
        for s = 1:length(sessions)
            for a = 1:length(arms)
                sess_num = sessions(s);
                arm = arms{a};
                field_name = sprintf('sess%d_%s_%s', sess_num, direction, arm);
                
                if any(strcmp(detailed_fields, field_name))
                    stats = slope_results.detailed.(field_name);
                    
                    if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                        individual_data{s, a} = stats.mouse_slopes;
                        mean_values(s, a) = mean(stats.mouse_slopes);
                        sem_values(s, a) = std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes));
                        
                        if isfield(stats, 'mouse_ids') && ~isempty(stats.mouse_ids)
                            mouse_ids{s, a} = stats.mouse_ids;
                        else
                            mouse_ids{s, a} = arrayfun(@(x) sprintf('Mouse_%d', x), 1:length(stats.mouse_slopes), 'UniformOutput', false);
                        end
                    else
                        mean_values(s, a) = stats.slope;
                        sem_values(s, a) = stats.se_slope;
                        individual_data{s, a} = [];
                        mouse_ids{s, a} = {};
                    end
                else
                    mean_values(s, a) = NaN;
                    sem_values(s, a) = NaN;
                    individual_data{s, a} = [];
                    mouse_ids{s, a} = {};
                end
            end
        end
        
        % Use global y-range for all positioning calculations
        y_range = global_y_range;
        max_val = global_max;
        min_val = global_min;
        
        % Create bar plot
        bar_width = 0.3;
        group_spacing = 1.2; % Adjust this to increase/decrease space between groups (default: 1)
        x_positions = (1:length(sessions)) * group_spacing;
        x_offset = [-0.2, 0.2]; % Food left, non-food right
        
        hold on;
        
        % Store mouse positions for connecting lines
        mouse_positions = cell(length(sessions), length(arms));
        
        for s = 1:length(sessions)
            for a = 1:length(arms)
                x_pos = x_positions(s) + x_offset(a);
                
                % Plot bar
                h_bar = bar(x_pos, mean_values(s, a), bar_width, ...
                           'FaceColor', colors{a}, 'EdgeColor', colors{a} * 0.8, ...
                           'LineWidth', 1.5, 'FaceAlpha', 0.8);
                
                % Add error bars
                errorbar(x_pos, mean_values(s, a), sem_values(s, a), ...
                        'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');
                
                % Add individual mouse data points
                if ~isempty(individual_data{s, a})
                    mouse_slopes = individual_data{s, a};
                    mouse_names = mouse_ids{s, a};
                    n_mice = length(mouse_slopes);
                    
                    % No jitter for x-coordinates (for connecting lines)
                    x_jitter = repmat(x_pos, n_mice, 1);
                    
                    % Store positions for connecting lines
                    mouse_positions{s, a} = struct('x_pos', x_jitter, 'y_pos', mouse_slopes, 'names', {mouse_names});
                    
                    % Plot individual points
                    scatter(x_jitter, mouse_slopes, 50, colors{a}, 'o', ...
                           'filled', 'MarkerFaceAlpha', 0.7, ...
                           'MarkerEdgeColor', colors{a} * 0.7, 'LineWidth', 1.5);
                end
            end
        end
        
        % Connect same mice across sessions for each arm separately
        for a = 1:length(arms)
            arm_name = arms{a};
            
            % Connect across all three sessions if data available
            for s = 1:(length(sessions)-1)
                if ~isempty(mouse_positions{s, a}) && ~isempty(mouse_positions{s+1, a})
                    sess_curr = mouse_positions{s, a};
                    sess_next = mouse_positions{s+1, a};
                    
                    % Find matching mice between consecutive sessions
                    for i = 1:length(sess_curr.names)
                        mouse_name = sess_curr.names{i};
                        next_idx = find(strcmp(sess_next.names, mouse_name));
                        
                        if ~isempty(next_idx)
                            % Connect this mouse across sessions
                            x_coords = [sess_curr.x_pos(i), sess_next.x_pos(next_idx)];
                            y_coords = [sess_curr.y_pos(i), sess_next.y_pos(next_idx)];
                            
                            plot(x_coords, y_coords, '-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1.0);
                        end
                    end
                end
            end
            
            % Perform paired t-test between Before and Test (sess 0 vs sess 2)
            if length(sessions) >= 3 && ~isempty(mouse_positions{1, a}) && ~isempty(mouse_positions{3, a})
                sess0_data = mouse_positions{1, a}; % Before
                sess2_data = mouse_positions{3, a}; % Test
                
                paired_sess0 = [];
                paired_sess2 = [];
                
                for i = 1:length(sess0_data.names)
                    mouse_name = sess0_data.names{i};
                    sess2_idx = find(strcmp(sess2_data.names, mouse_name));
                    
                    if ~isempty(sess2_idx)
                        paired_sess0 = [paired_sess0; sess0_data.y_pos(i)];
                        paired_sess2 = [paired_sess2; sess2_data.y_pos(sess2_idx)];
                    end
                end
                
                % Perform paired t-test
                if length(paired_sess0) >= 3
                    [~, p_value, ~, stats] = ttest(paired_sess2, paired_sess0);
                    
                    % Add significance indicator if significant
                    if p_value < 0.05
                        % Position bracket between Before and Test sessions
                        bracket_x1 = x_positions(1) + x_offset(a);
                        bracket_x2 = x_positions(3) + x_offset(a);
                        
                        arm_max_vals = [mean_values(1, a) + sem_values(1, a); 
                                       mean_values(3, a) + sem_values(3, a)];
                        if ~isempty(individual_data{1, a})
                            arm_max_vals = [arm_max_vals; max(individual_data{1, a})];
                        end
                        if ~isempty(individual_data{3, a})
                            arm_max_vals = [arm_max_vals; max(individual_data{3, a})];
                        end
                        max_val_arm = max(arm_max_vals);
                        
                        bracket_y = max_val_arm + y_range * 0.06;
                        star_y = max_val_arm + y_range * 0.08;
                        star_x = mean([bracket_x1, bracket_x2]);
                        
                        % Draw bracket
                        plot([bracket_x1, bracket_x1, bracket_x2, bracket_x2], ...
                             [bracket_y - y_range*0.01, bracket_y, bracket_y, bracket_y - y_range*0.01], ...
                             'k-', 'LineWidth', 1.5);
                        
                        % Add significance star
                        if p_value < 0.001
                            star_text = '***';
                            p_text = 'p<0.001';
                        elseif p_value < 0.01
                            star_text = '**';
                            p_text = 'p<0.01';
                        else
                            star_text = '*';
                            p_text = sprintf('p=%.3f', p_value);
                        end
                        
                        text(star_x, star_y, star_text, 'FontSize', 14, 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center', 'Color', 'k');
                        text(star_x, star_y - y_range * 0.04, p_text, 'FontSize', 9, ...
                             'HorizontalAlignment', 'center', 'Color', 'k');
                        
                        % Print results
                        fprintf('\n===== %s DIRECTION - %s ARM: BEFORE vs TEST =====\n', ...
                                upper(direction), upper(arm_name));
                        fprintf('Paired t-test: t(%d) = %.3f, p = %.4f\n', ...
                                length(paired_sess0)-1, stats.tstat, p_value);
                        fprintf('Mean difference: %.4f\n', mean(paired_sess2 - paired_sess0));
                        fprintf('==========================================\n');
                    end
                end
            end
        end
        
        % COMPARE FOOD vs NON-FOOD within each session (paired t-test)
        for s = 1:length(sessions)
            if length(arms) >= 2 && ~isempty(mouse_positions{s, 1}) && ~isempty(mouse_positions{s, 2})
                food_data = mouse_positions{s, 1}; % Food arm (a=1)
                nonfood_data = mouse_positions{s, 2}; % Non-food arm (a=2)
                
                paired_food = [];
                paired_nonfood = [];
                
                % Find matching mice between food and non-food arms
                for i = 1:length(food_data.names)
                    mouse_name = food_data.names{i};
                    nonfood_idx = find(strcmp(nonfood_data.names, mouse_name));
                    
                    if ~isempty(nonfood_idx)
                        paired_food = [paired_food; food_data.y_pos(i)];
                        paired_nonfood = [paired_nonfood; nonfood_data.y_pos(nonfood_idx)];
                    end
                end
                
                % Perform paired t-test if we have enough pairs
                if length(paired_food) >= 3
                    [~, p_value, ~, stats] = ttest(paired_food, paired_nonfood);
                    
                    % Add significance indicator if significant
                    if p_value < 0.05
                        % Position bracket between food and non-food bars within this session
                        bracket_x1 = x_positions(s) + x_offset(1); % Food
                        bracket_x2 = x_positions(s) + x_offset(2); % Non-food
                        
                        % Find highest point in this session group
                        session_max_vals = [mean_values(s, 1) + sem_values(s, 1); 
                                          mean_values(s, 2) + sem_values(s, 2)];
                        if ~isempty(individual_data{s, 1})
                            session_max_vals = [session_max_vals; max(individual_data{s, 1})];
                        end
                        if ~isempty(individual_data{s, 2})
                            session_max_vals = [session_max_vals; max(individual_data{s, 2})];
                        end
                        max_val_session = max(session_max_vals);
                        
                        % Position bracket lower than cross-session comparisons
                        bracket_y = max_val_session + y_range * 0.02;
                        star_y = max_val_session + y_range * 0.04;
                        star_x = x_positions(s); % Center of session
                        
                        % Draw bracket (shorter, within session)
                        plot([bracket_x1, bracket_x1, bracket_x2, bracket_x2], ...
                             [bracket_y - y_range*0.005, bracket_y, bracket_y, bracket_y - y_range*0.005], ...
                             'Color', [0.3, 0.3, 0.7], 'LineWidth', 1.2);
                        
                        % Add significance star
                        if p_value < 0.001
                            star_text = '***';
                            p_text = 'p<0.001';
                        elseif p_value < 0.01
                            star_text = '**';
                            p_text = 'p<0.01';
                        else
                            star_text = '*';
                            p_text = sprintf('p=%.3f', p_value);
                        end
                        
                        text(star_x, star_y, star_text, 'FontSize', 11, 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.7]);
                        text(star_x, star_y - y_range * 0.02, p_text, 'FontSize', 8, ...
                             'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.7]);
                        
                        % Print results
                        fprintf('\n===== %s DIRECTION - %s: FOOD vs NON-FOOD =====\n', ...
                                upper(direction), session_labels{s});
                        fprintf('Paired t-test: t(%d) = %.3f, p = %.4f\n', ...
                                length(paired_food)-1, stats.tstat, p_value);
                        fprintf('Mean difference (Food - Non-Food): %.4f\n', mean(paired_food - paired_nonfood));
                        fprintf('==========================================\n');
                    end
                end
            end
        end
        
        % Add individual slope significance indicators
        for s = 1:length(sessions)
            for a = 1:length(arms)
                if slope_data.is_significant(s, d, a)
                    x_pos = x_positions(s) + x_offset(a);
                    y_pos = mean_values(s, a) + sem_values(s, a);
                    y_star = y_pos + y_range * 0.3;
                    
                    text(x_pos, y_star, '★', 'FontSize', 12, ...
                         'HorizontalAlignment', 'center', 'Color', 'k', ...
                         'FontName', 'Arial Unicode MS');
                end
            end
        end
        
        % Formatting
        xlabel('Sessions', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Slope (z-scored dF/F per distance unit)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('Neural Activity Slopes: %s Direction (Combined Non-Food)', direction_label), ...
              'FontSize', 18, 'FontWeight', 'bold');
        
        % Set x-axis
        set(gca, 'XTick', x_positions, 'XTickLabel', session_labels, 'FontSize', 14);
        
        % Add horizontal line at zero
        line([x_positions(1) - 0.5, x_positions(end) + 0.5], [0, 0], 'Color', [0.3, 0.3, 0.3], ...
             'LineStyle', '--', 'LineWidth', 1.5);
        
        % Set limits (use global limits for consistency)
        xlim([x_positions(1) - 0.5, x_positions(end) + 0.5]);
        ylim(global_ylim);
        
        % Add legend
        legend_handles = [];
        for a = 1:length(arms)
            h = scatter(NaN, NaN, 100, colors{a}, 'o', 'filled', ...
                       'MarkerEdgeColor', colors{a} * 0.7, 'LineWidth', 1.5);
            legend_handles(a) = h;
        end
        legend(legend_handles, arm_labels, 'Location', 'best', 'FontSize', 14, ...
               'Box', 'off');
        
        % Add sample size information
        for s = 1:length(sessions)
            for a = 1:length(arms)
                if ~isempty(individual_data{s, a})
                    n_mice = length(individual_data{s, a});
                    x_pos = x_positions(s) + x_offset(a);
                    y_pos = global_min - global_y_range * 0.08;
                    
                    text(x_pos, y_pos, sprintf('n=%d', n_mice), ...
                         'FontSize', 10, 'HorizontalAlignment', 'center', ...
                         'Color', colors{a}, 'FontWeight', 'bold');
                end
            end
        end
        
        % Formatting
        set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
        
        hold off;
        
        % ===== PERFORM TWO-WAY REPEATED MEASURES ANOVA =====
        perform_repeated_measures_anova(individual_data, mouse_ids, sessions, arms, direction, session_labels, arm_labels);
        
        if options.save_plots
            filename = sprintf('combined_slope_barplot_%s_by_session.png', lower(direction));
            saveas(gcf, fullfile(options.figure_path, filename));
        end
    end
end

function create_combined_slope_barplot_by_arm(slope_data, sessions, arms, directions, slope_results, options)
    % NEW function: Groups by arm, bars are sessions
    
    direction_labels = {'Towards', 'Away'};
    
    % Session colors: Before (gray), Learning (blue), Test (light blue)
    session_colors = {[0.6, 0.6, 0.6], [0.2, 0.4, 0.8], [0.1, 0.7, 0.9]};
    
    % CALCULATE GLOBAL Y-LIMITS FIRST
    detailed_fields = fieldnames(slope_results.detailed);
    global_individual_data = [];
    global_mean_values = [];
    global_sem_values = [];
    
    for d = 1:length(directions)
        for s = 1:length(sessions)
            for a = 1:length(arms)
                sess_num = sessions(s);
                direction = directions{d};
                arm = arms{a};
                field_name = sprintf('sess%d_%s_%s', sess_num, direction, arm);
                
                if any(strcmp(detailed_fields, field_name))
                    stats = slope_results.detailed.(field_name);
                    
                    if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                        global_individual_data = [global_individual_data; stats.mouse_slopes];
                        global_mean_values = [global_mean_values; mean(stats.mouse_slopes)];
                        global_sem_values = [global_sem_values; std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes))];
                    else
                        global_mean_values = [global_mean_values; stats.slope];
                        global_sem_values = [global_sem_values; stats.se_slope];
                    end
                end
            end
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
    
    % Apply margins
    global_bottom_margin = global_y_range * 0.15;
    global_top_margin = global_y_range * 0.35;
    global_ylim = [global_min - global_bottom_margin, global_max + global_top_margin];
    
    for d = 1:length(directions)
        direction = directions{d};
        direction_label = direction_labels{d};
        
        figure('Name', sprintf('Combined Non-Food Slope Analysis - %s Direction (Grouped by Arm)', direction_label), ...
               'Position', [100 + (d-1)*800, 100, 600, 700]);
        
        % Create session labels
        session_labels = cell(1, length(sessions));
        for i = 1:length(sessions)
            sess_num = sessions(i);
            switch sess_num
                case 0
                    session_labels{i} = 'Before';
                case 1
                    session_labels{i} = 'Learning';
                case 2
                    session_labels{i} = 'Test';
                otherwise
                    session_labels{i} = sprintf('Session %d', sess_num);
            end
        end
        
        % Arm labels for x-axis groups
        arm_labels = {'Food', 'Non-Food'};
        
        % Get individual mouse data for this direction
        individual_data = cell(length(arms), length(sessions)); % Swapped dimensions
        mouse_ids = cell(length(arms), length(sessions));
        mean_values = zeros(length(arms), length(sessions));
        sem_values = zeros(length(arms), length(sessions));
        
        for a = 1:length(arms)
            for s = 1:length(sessions)
                sess_num = sessions(s);
                arm = arms{a};
                field_name = sprintf('sess%d_%s_%s', sess_num, direction, arm);
                
                if any(strcmp(detailed_fields, field_name))
                    stats = slope_results.detailed.(field_name);
                    
                    if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                        individual_data{a, s} = stats.mouse_slopes;
                        mean_values(a, s) = mean(stats.mouse_slopes);
                        sem_values(a, s) = std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes));
                        
                        if isfield(stats, 'mouse_ids') && ~isempty(stats.mouse_ids)
                            mouse_ids{a, s} = stats.mouse_ids;
                        else
                            mouse_ids{a, s} = arrayfun(@(x) sprintf('Mouse_%d', x), 1:length(stats.mouse_slopes), 'UniformOutput', false);
                        end
                    else
                        mean_values(a, s) = stats.slope;
                        sem_values(a, s) = stats.se_slope;
                        individual_data{a, s} = [];
                        mouse_ids{a, s} = {};
                    end
                else
                    mean_values(a, s) = NaN;
                    sem_values(a, s) = NaN;
                    individual_data{a, s} = [];
                    mouse_ids{a, s} = {};
                end
            end
        end
        
        y_range = global_y_range;
        max_val = global_max;
        min_val = global_min;
        
        % Create bar plot - NOW GROUPED BY ARM
        bar_width = 0.25;
        group_spacing = 1.2; % Adjust this to increase/decrease space between groups (default: 1)
        x_positions = (1:length(arms)) * group_spacing; % Two groups: Food, Non-Food
        n_sessions = length(sessions);
        
        % Calculate offsets for 3 bars (Before, Learning, Test)
        if n_sessions == 3
            x_offset = [-0.3, 0, 0.3]; % Before left, Learning center, Test right
        elseif n_sessions == 2
            x_offset = [-0.15, 0.15];
        else
            x_offset = 0;
        end
        
        hold on;
        
        % Store mouse positions for connecting lines
        mouse_positions = cell(length(arms), length(sessions)); % [arm, session]
        
        % Plot bars: for each arm group, plot session bars
        for a = 1:length(arms)
            for s = 1:length(sessions)
                x_pos = x_positions(a) + x_offset(s);
                
                % Plot bar with session-specific color
                h_bar = bar(x_pos, mean_values(a, s), bar_width, ...
                           'FaceColor', session_colors{s}, 'EdgeColor', session_colors{s} * 0.8, ...
                           'LineWidth', 1.5, 'FaceAlpha', 0.8);
                
                % Add error bars
                errorbar(x_pos, mean_values(a, s), sem_values(a, s), ...
                        'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');
                
                % Add individual mouse data points
                if ~isempty(individual_data{a, s})
                    mouse_slopes = individual_data{a, s};
                    mouse_names = mouse_ids{a, s};
                    n_mice = length(mouse_slopes);
                    
                    x_jitter = repmat(x_pos, n_mice, 1);
                    
                    % Store positions
                    mouse_positions{a, s} = struct('x_pos', x_jitter, 'y_pos', mouse_slopes, 'names', {mouse_names});
                    
                    % Plot individual points
                    scatter(x_jitter, mouse_slopes, 50, session_colors{s}, 'o', ...
                           'filled', 'MarkerFaceAlpha', 0.7, ...
                           'MarkerEdgeColor', session_colors{s} * 0.7, 'LineWidth', 1.5);
                end
            end
        end
        
        % Connect same mice across SESSIONS within each arm
        for a = 1:length(arms)
            for s = 1:(length(sessions)-1)
                if ~isempty(mouse_positions{a, s}) && ~isempty(mouse_positions{a, s+1})
                    sess_curr = mouse_positions{a, s};
                    sess_next = mouse_positions{a, s+1};
                    
                    for i = 1:length(sess_curr.names)
                        mouse_name = sess_curr.names{i};
                        next_idx = find(strcmp(sess_next.names, mouse_name));
                        
                        if ~isempty(next_idx)
                            x_coords = [sess_curr.x_pos(i), sess_next.x_pos(next_idx)];
                            y_coords = [sess_curr.y_pos(i), sess_next.y_pos(next_idx)];
                            
                            plot(x_coords, y_coords, '-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1.0);
                        end
                    end
                end
            end
            
            % Within-arm comparison: Before vs Test (paired t-test)
            if length(sessions) >= 3 && ~isempty(mouse_positions{a, 1}) && ~isempty(mouse_positions{a, 3})
                sess0_data = mouse_positions{a, 1}; % Before
                sess2_data = mouse_positions{a, 3}; % Test
                
                paired_sess0 = [];
                paired_sess2 = [];
                
                for i = 1:length(sess0_data.names)
                    mouse_name = sess0_data.names{i};
                    sess2_idx = find(strcmp(sess2_data.names, mouse_name));
                    
                    if ~isempty(sess2_idx)
                        paired_sess0 = [paired_sess0; sess0_data.y_pos(i)];
                        paired_sess2 = [paired_sess2; sess2_data.y_pos(sess2_idx)];
                    end
                end
                
                if length(paired_sess0) >= 3
                    [~, p_value, ~, stats] = ttest(paired_sess2, paired_sess0);
                    
                    if p_value < 0.05
                        bracket_x1 = x_positions(a) + x_offset(1); % Before
                        bracket_x2 = x_positions(a) + x_offset(3); % Test
                        
                        arm_max_vals = [mean_values(a, 1) + sem_values(a, 1);
                                       mean_values(a, 3) + sem_values(a, 3)];
                        if ~isempty(individual_data{a, 1})
                            arm_max_vals = [arm_max_vals; max(individual_data{a, 1})];
                        end
                        if ~isempty(individual_data{a, 3})
                            arm_max_vals = [arm_max_vals; max(individual_data{a, 3})];
                        end
                        max_val_arm = max(arm_max_vals);
                        
                        bracket_y = max_val_arm + y_range * 0.06;
                        star_y = max_val_arm + y_range * 0.08;
                        star_x = mean([bracket_x1, bracket_x2]);
                        
                        % Draw bracket
                        plot([bracket_x1, bracket_x1, bracket_x2, bracket_x2], ...
                             [bracket_y - y_range*0.01, bracket_y, bracket_y, bracket_y - y_range*0.01], ...
                             'k-', 'LineWidth', 1.5);
                        
                        % Add significance star
                        if p_value < 0.001
                            star_text = '***';
                            p_text = 'p<0.001';
                        elseif p_value < 0.01
                            star_text = '**';
                            p_text = 'p<0.01';
                        else
                            star_text = '*';
                            p_text = sprintf('p=%.3f', p_value);
                        end
                        
                        text(star_x, star_y, star_text, 'FontSize', 14, 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center', 'Color', 'k');
                        text(star_x, star_y - y_range * 0.04, p_text, 'FontSize', 9, ...
                             'HorizontalAlignment', 'center', 'Color', 'k');
                        
                        fprintf('\n===== %s DIRECTION - %s ARM: BEFORE vs TEST =====\n', ...
                                upper(direction), upper(arms{a}));
                        fprintf('Paired t-test: t(%d) = %.3f, p = %.4f\n', ...
                                length(paired_sess0)-1, stats.tstat, p_value);
                        fprintf('Mean difference: %.4f\n', mean(paired_sess2 - paired_sess0));
                        fprintf('==========================================\n');
                    end
                end
            end
        end
        
        % Between-arm comparison: Food vs Non-Food for each session
        for s = 1:length(sessions)
            if length(arms) >= 2 && ~isempty(mouse_positions{1, s}) && ~isempty(mouse_positions{2, s})
                food_data = mouse_positions{1, s}; % Food
                nonfood_data = mouse_positions{2, s}; % Non-food
                
                paired_food = [];
                paired_nonfood = [];
                
                for i = 1:length(food_data.names)
                    mouse_name = food_data.names{i};
                    nonfood_idx = find(strcmp(nonfood_data.names, mouse_name));
                    
                    if ~isempty(nonfood_idx)
                        paired_food = [paired_food; food_data.y_pos(i)];
                        paired_nonfood = [paired_nonfood; nonfood_data.y_pos(nonfood_idx)];
                    end
                end
                
                if length(paired_food) >= 3
                    [~, p_value, ~, stats] = ttest(paired_food, paired_nonfood);
                    
                    if p_value < 0.05
                        bracket_x1 = x_positions(1) + x_offset(s); % Food
                        bracket_x2 = x_positions(2) + x_offset(s); % Non-food
                        
                        session_max_vals = [mean_values(1, s) + sem_values(1, s);
                                          mean_values(2, s) + sem_values(2, s)];
                        if ~isempty(individual_data{1, s})
                            session_max_vals = [session_max_vals; max(individual_data{1, s})];
                        end
                        if ~isempty(individual_data{2, s})
                            session_max_vals = [session_max_vals; max(individual_data{2, s})];
                        end
                        max_val_session = max(session_max_vals);
                        
                        bracket_y = max_val_session + y_range * 0.02;
                        star_y = max_val_session + y_range * 0.04;
                        star_x = mean([bracket_x1, bracket_x2]);
                        
                        % Draw bracket
                        plot([bracket_x1, bracket_x1, bracket_x2, bracket_x2], ...
                             [bracket_y - y_range*0.005, bracket_y, bracket_y, bracket_y - y_range*0.005], ...
                             'Color', [0.3, 0.3, 0.7], 'LineWidth', 1.2);
                        
                        if p_value < 0.001
                            star_text = '***';
                            p_text = 'p<0.001';
                        elseif p_value < 0.01
                            star_text = '**';
                            p_text = 'p<0.01';
                        else
                            star_text = '*';
                            p_text = sprintf('p=%.3f', p_value);
                        end
                        
                        text(star_x, star_y, star_text, 'FontSize', 11, 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.7]);
                        text(star_x, star_y - y_range * 0.02, p_text, 'FontSize', 8, ...
                             'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.7]);
                        
                        fprintf('\n===== %s DIRECTION - %s: FOOD vs NON-FOOD =====\n', ...
                                upper(direction), session_labels{s});
                        fprintf('Paired t-test: t(%d) = %.3f, p = %.4f\n', ...
                                length(paired_food)-1, stats.tstat, p_value);
                        fprintf('Mean difference (Food - Non-Food): %.4f\n', mean(paired_food - paired_nonfood));
                        fprintf('==========================================\n');
                    end
                end
            end
        end
        
        % Add individual slope significance indicators
        for a = 1:length(arms)
            for s = 1:length(sessions)
                if slope_data.is_significant(s, d, a)
                    x_pos = x_positions(a) + x_offset(s);
                    y_pos = mean_values(a, s) + sem_values(a, s);
                    y_star = y_pos + y_range * 0.3;
                    
                    text(x_pos, y_star, '★', 'FontSize', 12, ...
                         'HorizontalAlignment', 'center', 'Color', 'k', ...
                         'FontName', 'Arial Unicode MS');
                end
            end
        end
        
        % Formatting
        xlabel('Arms', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Slope (z-scored dF/F per distance unit)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('Neural Activity Slopes: %s Direction (Combined Non-Food)', direction_label), ...
              'FontSize', 18, 'FontWeight', 'bold');
        
        % Set x-axis
        set(gca, 'XTick', x_positions, 'XTickLabel', arm_labels, 'FontSize', 14);
        
        % Add horizontal line at zero
        line([x_positions(1) - 0.5, x_positions(end) + 0.5], [0, 0], 'Color', [0.3, 0.3, 0.3], ...
             'LineStyle', '--', 'LineWidth', 1.5);
        
        % Set limits
        xlim([x_positions(1) - 0.5, x_positions(end) + 0.5]);
        ylim(global_ylim);
        
        % Add legend for sessions
        legend_handles = [];
        for s = 1:length(sessions)
            h = scatter(NaN, NaN, 100, session_colors{s}, 'o', 'filled', ...
                       'MarkerEdgeColor', session_colors{s} * 0.7, 'LineWidth', 1.5);
            legend_handles(s) = h;
        end
        legend(legend_handles, session_labels, 'Location', 'best', 'FontSize', 14, ...
               'Box', 'off');
        
        % Add sample size information
        for a = 1:length(arms)
            for s = 1:length(sessions)
                if ~isempty(individual_data{a, s})
                    n_mice = length(individual_data{a, s});
                    x_pos = x_positions(a) + x_offset(s);
                    y_pos = global_min - global_y_range * 0.08;
                    
                    text(x_pos, y_pos, sprintf('n=%d', n_mice), ...
                         'FontSize', 10, 'HorizontalAlignment', 'center', ...
                         'Color', session_colors{s}, 'FontWeight', 'bold');
                end
            end
        end
        
        set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
        
        hold off;
        
        if options.save_plots
            filename = sprintf('combined_slope_barplot_%s_by_arm.png', lower(direction));
            saveas(gcf, fullfile(options.figure_path, filename));
        end
    end
end

function perform_repeated_measures_anova(individual_data, mouse_ids, sessions, arms, direction, session_labels, arm_labels)
    % Perform two-way repeated measures ANOVA: Session × Arm
    % individual_data: cell array (varies by grouping - could be sessions×arms or arms×sessions)
    % mouse_ids: cell array of mouse identifiers
    
    fprintf('\n\n========================================\n');
    fprintf('TWO-WAY REPEATED MEASURES ANOVA: %s DIRECTION\n', upper(direction));
    fprintf('========================================\n');
    
    % Determine data organization (by session or by arm)
    [dim1, dim2] = size(individual_data);
    if dim1 == length(sessions) && dim2 == length(arms)
        % Data organized as: individual_data{session, arm}
        n_sessions = dim1;
        n_arms = dim2;
    elseif dim1 == length(arms) && dim2 == length(sessions)
        % Data organized as: individual_data{arm, session}
        n_arms = dim1;
        n_sessions = dim2;
        % Transpose to standardize
        individual_data_transposed = cell(n_sessions, n_arms);
        mouse_ids_transposed = cell(n_sessions, n_arms);
        for s = 1:n_sessions
            for a = 1:n_arms
                individual_data_transposed{s, a} = individual_data{a, s};
                mouse_ids_transposed{s, a} = mouse_ids{a, s};
            end
        end
        individual_data = individual_data_transposed;
        mouse_ids = mouse_ids_transposed;
    else
        fprintf('ERROR: Cannot determine data organization for ANOVA\n');
        return;
    end
    
    % Find common mice across all conditions
    all_mice = {};
    for s = 1:length(sessions)
        for a = 1:length(arms)
            if ~isempty(mouse_ids{s, a})
                all_mice = [all_mice; mouse_ids{s, a}(:)];
            end
        end
    end
    unique_mice = unique(all_mice);
    
    % Build data matrix: each row = one mouse, each column = one condition
    % Columns order: sess1_arm1, sess1_arm2, sess2_arm1, sess2_arm2, ...
    data_matrix = [];
    valid_mice = {};
    
    for m = 1:length(unique_mice)
        mouse_name = unique_mice{m};
        mouse_data = nan(1, n_sessions * n_arms);
        
        % Extract data for this mouse across all conditions
        has_all_data = true;
        col_idx = 1;
        for s = 1:n_sessions
            for a = 1:n_arms
                if ~isempty(mouse_ids{s, a})
                    mouse_idx = find(strcmp(mouse_ids{s, a}, mouse_name));
                    if ~isempty(mouse_idx)
                        mouse_data(col_idx) = individual_data{s, a}(mouse_idx);
                    else
                        has_all_data = false;
                    end
                else
                    has_all_data = false;
                end
                col_idx = col_idx + 1;
            end
        end
        
        % Only include mice with complete data
        if has_all_data && ~any(isnan(mouse_data))
            data_matrix = [data_matrix; mouse_data];
            valid_mice{end+1} = mouse_name;
        end
    end
    
    if size(data_matrix, 1) < 3
        fprintf('ERROR: Not enough mice with complete data (n=%d, need ≥3 for ANOVA)\n', size(data_matrix, 1));
        return;
    end
    
    fprintf('Complete data available for %d mice\n', size(data_matrix, 1));
    
    % Create table for repeated measures ANOVA
    % Column names: Sess0_Food, Sess0_NonFood, Sess1_Food, Sess1_NonFood, ...
    col_names = {};
    for s = 1:n_sessions
        for a = 1:n_arms
            % Remove hyphens and spaces from labels for valid variable names
            sess_label_clean = strrep(session_labels{s}, ' ', '');
            arm_label_clean = strrep(strrep(arm_labels{a}, '-', ''), ' ', '');
            col_names{end+1} = sprintf('%s_%s', sess_label_clean, arm_label_clean);
        end
    end
    
    data_table = array2table(data_matrix, 'VariableNames', col_names);
    
    % Define within-subjects factors
    % Create a table defining the design
    Session = [];
    Arm = [];
    for s = 1:n_sessions
        for a = 1:n_arms
            Session = [Session; s];
            Arm = [Arm; a];
        end
    end
    within_design = table(Session, Arm, 'VariableNames', {'Session', 'Arm'});
    
    % Fit repeated measures model
    try
        rm = fitrm(data_table, sprintf('%s-%s~1', col_names{1}, col_names{end}), ...
                   'WithinDesign', within_design);
        
        % Run repeated measures ANOVA
        ranova_results = ranova(rm, 'WithinModel', 'Session*Arm');
        
        % Display results
        fprintf('\n--- ANOVA Results ---\n');
        disp(ranova_results);
        
        % Extract p-values
        p_session = ranova_results.pValue(1);
        p_arm = ranova_results.pValue(2);
        p_interaction = ranova_results.pValue(3);
        
        % Report main effects and interaction
        fprintf('\n--- Interpretation ---\n');
        fprintf('Main Effect of Session (Before, Learning, Test): ');
        if p_session < 0.001
            fprintf('p < 0.001 ***\n');
        elseif p_session < 0.01
            fprintf('p = %.4f **\n', p_session);
        elseif p_session < 0.05
            fprintf('p = %.4f *\n', p_session);
        else
            fprintf('p = %.4f (n.s.)\n', p_session);
        end
        
        fprintf('Main Effect of Arm (Food vs Non-Food): ');
        if p_arm < 0.001
            fprintf('p < 0.001 ***\n');
        elseif p_arm < 0.01
            fprintf('p = %.4f **\n', p_arm);
        elseif p_arm < 0.05
            fprintf('p = %.4f *\n', p_arm);
        else
            fprintf('p = %.4f (n.s.)\n', p_arm);
        end
        
        fprintf('Session × Arm Interaction: ');
        if p_interaction < 0.001
            fprintf('p < 0.001 ***\n');
        elseif p_interaction < 0.01
            fprintf('p = %.4f **\n', p_interaction);
        elseif p_interaction < 0.05
            fprintf('p = %.4f *\n', p_interaction);
        else
            fprintf('p = %.4f (n.s.)\n', p_interaction);
        end
        
        % Post-hoc tests if main effects are significant
        if p_session < 0.05
            fprintf('\n--- Post-hoc: Pairwise Session Comparisons (Bonferroni corrected) ---\n');
            % Average across arms for each session
            session_means = zeros(n_sessions, 1);
            for s = 1:n_sessions
                cols = (s-1)*n_arms + (1:n_arms);
                session_means(s) = mean(data_matrix(:, cols), 'all');
            end
            
            % Pairwise comparisons with Bonferroni correction
            n_comparisons = nchoosek(n_sessions, 2);
            alpha_corrected = 0.05 / n_comparisons;
            
            comparison_pairs = nchoosek(1:n_sessions, 2);
            for i = 1:size(comparison_pairs, 1)
                s1 = comparison_pairs(i, 1);
                s2 = comparison_pairs(i, 2);
                
                % Get data for these two sessions (averaged across arms)
                data_s1 = mean(data_matrix(:, (s1-1)*n_arms + (1:n_arms)), 2);
                data_s2 = mean(data_matrix(:, (s2-1)*n_arms + (1:n_arms)), 2);
                
                [~, p, ~, stats] = ttest(data_s1, data_s2);
                
                fprintf('%s vs %s: t(%d) = %.3f, p = %.4f', ...
                    session_labels{s1}, session_labels{s2}, ...
                    length(data_s1)-1, stats.tstat, p);
                
                if p < alpha_corrected
                    fprintf(' ** (Bonferroni corrected p < %.4f)\n', alpha_corrected);
                else
                    fprintf(' (n.s.)\n');
                end
            end
        end
        
        if p_arm < 0.05
            fprintf('\n--- Post-hoc: Food vs Non-Food (averaged across sessions) ---\n');
            % Average across sessions for each arm
            data_food = data_matrix(:, 1:n_arms:end); % Every n_arms columns starting from 1
            data_nonfood = data_matrix(:, 2:n_arms:end); % Every n_arms columns starting from 2
            
            [~, p, ~, stats] = ttest(mean(data_food, 2), mean(data_nonfood, 2));
            
            fprintf('Food vs Non-Food: t(%d) = %.3f, p = %.4f', ...
                length(data_food)-1, stats.tstat, p);
            if p < 0.05
                fprintf(' *\n');
            else
                fprintf(' (n.s.)\n');
            end
        end
        
        if p_interaction < 0.05
            fprintf('\n--- Post-hoc: Simple Effects Tests (Interaction significant) ---\n');
            fprintf('Testing Food vs Non-Food separately for each session:\n');
            
            alpha_corrected = 0.05 / n_sessions; % Bonferroni for multiple sessions
            for s = 1:n_sessions
                cols = (s-1)*n_arms + (1:n_arms);
                data_food_sess = data_matrix(:, cols(1));
                data_nonfood_sess = data_matrix(:, cols(2));
                
                [~, p, ~, stats] = ttest(data_food_sess, data_nonfood_sess);
                
                fprintf('  %s: t(%d) = %.3f, p = %.4f', ...
                    session_labels{s}, length(data_food_sess)-1, stats.tstat, p);
                
                if p < alpha_corrected
                    fprintf(' ** (Bonferroni corrected p < %.4f)\n', alpha_corrected);
                else
                    fprintf(' (n.s.)\n');
                end
            end
        end
        
    catch ME
        fprintf('ERROR running repeated measures ANOVA: %s\n', ME.message);
    end
    
    fprintf('========================================\n\n');
end