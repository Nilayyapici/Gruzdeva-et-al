function visualize_cheeseboard_slopes(slope_results, options)
    % visualize_cheeseboard_slopes - Creates comprehensive visualizations of slope analysis for cheeseboard maze
    %
    % Inputs:
    %   slope_results - Output from analyze_cheeseboard_slopes function
    %   options - Visualization options (optional):
    %       .plot_types - Cell array of plot types to create (default: {'heatmap', 'scatter', 'barplot'})
    %                    {'heatmap', 'barplot', 'scatter', 'summary'}
    %       .save_plots - Whether to save plots (default: false)
    %       .figure_path - Path to save figures (default: current directory)
    %       .session_colors - Cell array of colors for sessions {Before, Test}
    %                        Default: {[0.6, 0.6, 0.6], [0.3, 0.3, 0.3]}
    
    if nargin < 2
        options = struct();
    end
    
    if ~isfield(options, 'plot_types')
        options.plot_types = {'heatmap', 'scatter', 'barplot'};
    end
    if ~isfield(options, 'save_plots'), options.save_plots = false; end
    if ~isfield(options, 'figure_path'), options.figure_path = './'; end
    
    % Set default session colors if not provided
    if ~isfield(options, 'session_colors') || isempty(options.session_colors)
        % Default colors: Before (light gray), Test (dark gray)
        options.session_colors = {[0.6, 0.6, 0.6], [0.3, 0.3, 0.3]};
    end
    
    % Extract data from slope_results
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No slope analysis results to visualize.\n');
        return;
    end
    
    % Parse the results into organized data
    [slope_data, sessions, directions] = parse_cheeseboard_slope_results(slope_results);
    
    % Create requested visualizations
    for i = 1:length(options.plot_types)
        plot_type = options.plot_types{i};
        
        switch plot_type
            case 'heatmap'
                create_cheeseboard_slope_heatmap(slope_data, sessions, directions, options);
            case 'barplot'
                create_cheeseboard_slope_barplot(slope_data, sessions, directions, slope_results, options);
            case 'scatter'
                create_cheeseboard_slope_scatter(slope_data, sessions, directions, options);
            case 'summary'
                create_cheeseboard_summary_plot(slope_results, options);
            otherwise
                warning('Unknown plot type: %s', plot_type);
        end
    end
end

function [slope_data, sessions, directions] = parse_cheeseboard_slope_results(slope_results)
    % Parse slope results into organized matrices for cheeseboard analysis
    
    detailed_fields = fieldnames(slope_results.detailed);
    
    % Extract unique sessions and directions
    sessions = {};
    directions = {};
    
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        parts = split(field, '_');
        session_name = parts{1};
        direction = parts{2};
        
        if ~ismember(session_name, sessions)
            sessions{end+1} = session_name;
        end
        if ~ismember(direction, directions)
            directions{end+1} = direction;
        end
    end
    
    % Ensure specific order for consistency
    sessions = {'pre', 'test'};
    directions = {'towards', 'away'};
    
    % Initialize data matrices
    slope_data = struct();
    slope_data.slopes = nan(length(sessions), length(directions));
    slope_data.p_values = nan(length(sessions), length(directions));
    slope_data.r_squared = nan(length(sessions), length(directions));
    slope_data.n_mice = nan(length(sessions), length(directions));
    slope_data.is_significant = false(length(sessions), length(directions));
    slope_data.correct_direction = false(length(sessions), length(directions));
    
    % Fill matrices
    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        
        parts = split(field, '_');
        session_name = parts{1};
        direction = parts{2};
        
        s_idx = find(strcmp(sessions, session_name));
        d_idx = find(strcmp(directions, direction));
        
        if ~isempty(s_idx) && ~isempty(d_idx)
            slope_data.slopes(s_idx, d_idx) = stats.slope;
            slope_data.p_values(s_idx, d_idx) = stats.p_value;
            slope_data.r_squared(s_idx, d_idx) = stats.r_squared;
            slope_data.n_mice(s_idx, d_idx) = stats.n_mice;
            slope_data.is_significant(s_idx, d_idx) = stats.is_significant;
            slope_data.correct_direction(s_idx, d_idx) = stats.correct_direction;
        end
    end
end

function create_cheeseboard_slope_heatmap(slope_data, sessions, directions, options)
    % Create heatmap showing slopes across sessions and directions
    
    % Calculate figure dimensions
    n_sessions = length(sessions);
    n_directions = length(directions);
    
    % Calculate appropriate figure size
    cell_size = 150; % pixels per cell
    heatmap_height = n_directions * cell_size + 150;
    heatmap_width = n_sessions * cell_size + 150;
    
    figure('Name', 'Cheeseboard Slope Analysis Heatmap', 'Position', [100, 100, heatmap_width, heatmap_height]);
    
    % Calculate global color limits for consistency
    all_slopes = slope_data.slopes(:);
    max_slope = max(abs(all_slopes), [], 'omitnan');
    if isnan(max_slope) || max_slope == 0
        color_limit = 0.01;
    else
        color_limit = max(max_slope, 0.005);
    end
    
    % Create session labels
    session_labels = {'Before', 'Test'};
    
    % Create direction labels
    direction_labels = cell(1, length(directions));
    for i = 1:length(directions)
        direction_labels{i} = [upper(directions{i}(1)), directions{i}(2:end)];
    end
    
    % Transpose matrix: sessions as columns (x-axis), directions as rows (y-axis)
    slope_matrix_transposed = slope_data.slopes';
    
    % Create heatmap with 3 decimal place formatting
    h = heatmap(session_labels, direction_labels, slope_matrix_transposed, ...
               'Colormap', redblue(256), 'ColorLimits', [-color_limit, color_limit], ...
               'CellLabelFormat', '%.3f');
    
    h.Title = 'Neural Activity Slopes: dF/F vs Distance to Food (Cheeseboard)';
    h.XLabel = 'Sessions';
    h.YLabel = 'Movement Direction';
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'cheeseboard_slope_heatmap.png'));
    end
end

function create_cheeseboard_slope_scatter(slope_data, sessions, directions, options)
    % Create scatter plot: Towards vs Away slopes for cheeseboard analysis
    
    figure('Name', 'Cheeseboard Slope Analysis: Towards vs Away Comparison', 'Position', [100, 100, 1200, 800]);
    
    % Find direction indices
    towards_idx = find(strcmp(directions, 'towards'));
    away_idx = find(strcmp(directions, 'away'));
    
    if isempty(towards_idx) || isempty(away_idx)
        error('Both "towards" and "away" directions must be present in the data');
    end
    
    % Create main scatter plot
    subplot(1, 4, [1, 2, 3]);
    hold on;
    
    % Plot data points
    all_x_vals = [];
    all_y_vals = [];
    
    for s = 1:length(sessions)
        x_val = slope_data.slopes(s, towards_idx); % Towards slope
        y_val = slope_data.slopes(s, away_idx);    % Away slope
        
        if ~isnan(x_val) && ~isnan(y_val)
            all_x_vals = [all_x_vals, x_val];
            all_y_vals = [all_y_vals, y_val];
            
            % Get color based on SESSION
            point_color = options.session_colors{s};
            
            % Check significance for both directions
            towards_sig = slope_data.is_significant(s, towards_idx);
            away_sig = slope_data.is_significant(s, away_idx);
            
            % Determine marker properties based on significance
            if towards_sig && away_sig
                marker_size = 200;
                edge_color = [0, 0, 0];
                line_width = 4;
                face_alpha = 0.9;
            elseif towards_sig || away_sig
                marker_size = 150;
                edge_color = point_color * 0.7;
                line_width = 3;
                face_alpha = 0.7;
            else
                marker_size = 100;
                edge_color = point_color * 0.5;
                line_width = 2;
                face_alpha = 0.5;
            end
            
            % Plot the point
            h = scatter(x_val, y_val, marker_size, point_color, 'o', ...
                       'filled', 'MarkerEdgeColor', edge_color, 'LineWidth', line_width, ...
                       'MarkerFaceAlpha', face_alpha);
            
            % Add session label
            session_names = {'Before', 'Test'};
            if towards_sig || away_sig
                session_label = sprintf('%s*', session_names{s});
            else
                session_label = session_names{s};
            end
            
            % Position label to avoid overlap
            if ~isempty(all_x_vals) && ~isempty(all_y_vals)
                x_offset = max(abs(all_x_vals)) * 0.05;
                y_offset = max(abs(all_y_vals)) * 0.05;
            else
                x_offset = 0.01;
                y_offset = 0.01;
            end
            
            text(x_val + x_offset, y_val + y_offset, session_label, ...
                 'FontSize', 12, 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.2]);
        end
    end
    
    % Set axis limits with padding
    if ~isempty(all_x_vals) && ~isempty(all_y_vals)
        x_range = max(abs(all_x_vals));
        y_range = max(abs(all_y_vals));
        max_range = max(x_range, y_range) * 1.3;
        xlim([-max_range, max_range]);
        ylim([-max_range, max_range]);
    else
        xlim([-0.1, 0.1]);
        ylim([-0.1, 0.1]);
    end
    
    % Add reference lines
    xl = xlim;
    yl = ylim;
    
    % Axes through origin
    line(xl, [0, 0], 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
    line([0, 0], yl, 'Color', [0.3, 0.3, 0.3], 'LineStyle', '-', 'LineWidth', 1.5);
    
    % Diagonal line (y = x)
    diag_lim = min(max(abs(xl)), max(abs(yl)));
    line([-diag_lim, diag_lim], [-diag_lim, diag_lim], 'Color', [0.6, 0.6, 0.6], ...
         'LineStyle', ':', 'LineWidth', 1);
    
    % Expected pattern line (y = x for cheeseboard - both should increase with food proximity)
    line([-diag_lim, diag_lim], [-diag_lim, diag_lim], 'Color', [0.4, 0.8, 0.4], ...
         'LineStyle', '-.', 'LineWidth', 2);
    
    % Add subtle quadrant shading for expected pattern
    fill([0, xl(2), xl(2), 0], [0, yl(2), yl(2), 0], [0.9, 1, 0.9], ...
         'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % Formatting
    xlabel('Towards Food Slope (z-score/distance)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Away from Food Slope (z-score/distance)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Neural Activity Slopes: Towards vs Away Food (Cheeseboard)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Create legend for sessions
    session_legend_handles = [];
    session_legend_labels = {};
    session_names = {'Before', 'Test'};
    
    % Plot dummy points outside visible area for legend
    dummy_x = xl(1) - 1000;
    dummy_y = yl(1) - 1000;
    
    for s = 1:length(sessions)
        h_sess = scatter(dummy_x, dummy_y, 120, options.session_colors{s}, 'o', 'filled');
        session_legend_handles(end+1) = h_sess;
        session_legend_labels{end+1} = session_names{s};
    end
    
    legend(session_legend_handles, session_legend_labels, 'Location', 'northwest', ...
           'FontSize', 14, 'Box', 'off');
    
    % Add quadrant interpretation
    text(xl(2)*0.95, yl(2)*0.95, {'Towards↑', 'Away↑'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [0.9, 1, 0.9], 'Margin', 2);
    text(xl(1)*0.95, yl(2)*0.95, {'Towards↓', 'Away↑'}, 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'Margin', 2);
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
    
    text(0.05, 0.82, 'I (↑↑):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0.7, 0]);
    text(0.05, 0.79, 'Activity increases both', 'FontSize', 10);
    text(0.05, 0.76, 'approaching & leaving food', 'FontSize', 10);
    text(0.05, 0.73, '(EXPECTED FOR CHEESEBOARD)', 'FontSize', 9, 'FontWeight', 'bold', ...
         'Color', [0, 0.7, 0]);
    
    text(0.05, 0.67, 'II (↓↑):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0, 0, 0.8]);
    text(0.05, 0.64, 'Activity decreases approaching,', 'FontSize', 10);
    text(0.05, 0.61, 'increases leaving food', 'FontSize', 10);
    
    text(0.05, 0.55, 'III (↓↓):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.6, 0, 0.6]);
    text(0.05, 0.52, 'Activity decreases both', 'FontSize', 10);
    text(0.05, 0.49, 'approaching & leaving food', 'FontSize', 10);
    
    text(0.05, 0.43, 'IV (↑↓):', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.6, 0, 0.6]);
    text(0.05, 0.40, 'Activity increases approaching,', 'FontSize', 10);
    text(0.05, 0.37, 'decreases leaving food', 'FontSize', 10);
    
    % Visual elements guide
    text(0.05, 0.30, 'VISUAL ELEMENTS:', 'FontSize', 12, 'FontWeight', 'bold');
    
    text(0.05, 0.26, 'Colors (Sessions):', 'FontSize', 11, 'FontWeight', 'bold');
    
    % Display session colors
    session_display_names = {'Before', 'Test'};
    y_pos = [0.23, 0.20];
    
    for s = 1:min(length(sessions), 2)
        if s <= length(y_pos)
            label_text = sprintf('• %s', session_display_names{s});
            text(0.05, y_pos(s), label_text, 'FontSize', 10, 'Color', options.session_colors{s});
        end
    end
    
    text(0.05, 0.15, 'Size & Edge = Significance', 'FontSize', 10, 'FontWeight', 'bold');
    text(0.05, 0.12, '* = Significant slope', 'FontSize', 10, 'FontWeight', 'bold');
    
    text(0.05, 0.06, 'Expected: Upper-right quadrant', 'FontSize', 10, 'FontWeight', 'bold', ...
         'Color', [0, 0.7, 0]);
    text(0.05, 0.03, '(Activity ↑ both directions)', 'FontSize', 9, 'Color', [0, 0.7, 0]);
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'cheeseboard_slope_scatter.png'));
    end
end

function create_cheeseboard_slope_barplot(slope_data, sessions, directions, slope_results, options)
    % Create enhanced bar plot with individual mouse data points and SEM
    
    figure('Name', 'Cheeseboard Slope Analysis: Enhanced Bar Plot', 'Position', [100, 100, 600, 700]);
    
    % Create session labels
    session_labels = {'Before', 'Test'};
    
    % Direction labels for x-axis
    direction_labels = {'Towards', 'Away'};
    
    % Enhanced colors based on direction and session
    % Towards: bluish (before lighter, test darker)
    % Away: reddish (before lighter, test darker)
    towards_colors = {[0.7, 0.8, 1.0], [0.3, 0.5, 0.9]}; % Light blue, Dark blue
    away_colors = {[1.0, 0.7, 0.7], [0.9, 0.3, 0.3]};     % Light red, Dark red
    
    % Create direction-specific color scheme
    direction_session_colors = cell(length(directions), length(sessions));
    direction_session_colors_light = cell(length(directions), length(sessions));
    
    for d = 1:length(directions)
        if strcmp(directions{d}, 'towards')
            direction_session_colors(d, :) = towards_colors;
            % Lighter versions for bar faces
            direction_session_colors_light{d, 1} = towards_colors{1} * 0.8 + [0.2, 0.2, 0.2];
            direction_session_colors_light{d, 2} = towards_colors{2} * 0.8 + [0.2, 0.2, 0.2];
        else % away
            direction_session_colors(d, :) = away_colors;
            % Lighter versions for bar faces  
            direction_session_colors_light{d, 1} = away_colors{1} * 0.8 + [0.2, 0.2, 0.2];
            direction_session_colors_light{d, 2} = away_colors{2} * 0.8 + [0.2, 0.2, 0.2];
        end
    end
    
    % Get individual mouse data and calculate SEM
    individual_data = cell(length(sessions), length(directions));
    mouse_ids = cell(length(sessions), length(directions));
    mean_values = zeros(length(sessions), length(directions));
    sem_values = zeros(length(sessions), length(directions));
    
    % Extract individual mouse slopes from slope_results
    detailed_fields = fieldnames(slope_results.detailed);
    
    for s = 1:length(sessions)
        for d = 1:length(directions)
            session_name = sessions{s};
            direction = directions{d};
            field_name = sprintf('%s_%s', session_name, direction);
            
            if any(strcmp(detailed_fields, field_name))
                stats = slope_results.detailed.(field_name);
                
                % Get individual mouse slopes if available
                if isfield(stats, 'mouse_slopes') && ~isempty(stats.mouse_slopes)
                    individual_data{s, d} = stats.mouse_slopes;
                    mean_values(s, d) = mean(stats.mouse_slopes);
                    sem_values(s, d) = std(stats.mouse_slopes) / sqrt(length(stats.mouse_slopes));
                    
                    % Try to get mouse IDs if available
                    if isfield(stats, 'mouse_ids') && ~isempty(stats.mouse_ids)
                        mouse_ids{s, d} = stats.mouse_ids;
                    else
                        mouse_ids{s, d} = arrayfun(@(x) sprintf('Mouse_%d', x), 1:length(stats.mouse_slopes), 'UniformOutput', false);
                    end
                else
                    % Fallback to overall slope if individual data not available
                    mean_values(s, d) = stats.slope;
                    sem_values(s, d) = stats.se_slope;
                    individual_data{s, d} = []; 
                    mouse_ids{s, d} = {};
                end
            else
                mean_values(s, d) = NaN;
                sem_values(s, d) = NaN;
                individual_data{s, d} = [];
                mouse_ids{s, d} = {};
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % XLSX EXPORT
    % Saves barplot data to 'cheeseboard_slope_barplot_data.xlsx' with two sheets:
    %   1. "Summary" - mean slopes, SEM, significance, n per group
    %   2. "Individual Mouse Data" - per-mouse slopes in wide format (one column per group)
    % -------------------------------------------------------------------------
    save_barplot_data_to_xlsx(slope_data, sessions, directions, session_labels, ...
        direction_labels, mean_values, sem_values, individual_data, mouse_ids, options);
    % -------------------------------------------------------------------------
    
    % Create beautiful grouped bar plot with increased spacing
    bar_width = 0.3;
    x_positions = 1:length(directions);
    x_offset = [-0.2, 0.2];
    
    hold on;
    
    % Calculate proper y-range for positioning elements
    all_individual_values = [];
    for s = 1:length(sessions)
        for d = 1:length(directions)
            if ~isempty(individual_data{s, d})
                all_individual_values = [all_individual_values; individual_data{s, d}];
            end
        end
    end
    
    % Calculate comprehensive min/max
    group_max = max(mean_values(:) + sem_values(:), [], 'omitnan');
    group_min = min(mean_values(:) - sem_values(:), [], 'omitnan');
    
    if ~isempty(all_individual_values)
        individual_max = max(all_individual_values);
        individual_min = min(all_individual_values);
        max_val = max(group_max, individual_max);
        min_val = min(group_min, individual_min);
    else
        max_val = group_max;
        min_val = group_min;
    end
    
    y_range = max_val - min_val;
    if y_range == 0 || isnan(y_range)
        y_range = 0.1;
    end
    
    % Plot bars and collect mouse positions for connections
    bar_handles = [];
    mouse_positions = cell(length(directions), length(sessions));
    
    for s = 1:length(sessions)
        for d = 1:length(directions)
            x_pos = x_positions(d) + x_offset(s);
            
            session_direction_mean = mean_values(s, d);
            session_direction_sem = sem_values(s, d);
            
            % Create bars with direction-specific colors
            h_bar = bar(x_pos, session_direction_mean, bar_width, ...
                       'FaceColor', direction_session_colors_light{d, s}, ...
                       'EdgeColor', direction_session_colors{d, s}, ...
                       'LineWidth', 1.5, ...
                       'FaceAlpha', 0.8);
            
            if d == 1
                bar_handles(s) = h_bar;
            end
            
            % Add error bars (SEM)
            errorbar(x_pos, session_direction_mean, session_direction_sem, ...
                    'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');
            
            % Add individual mouse data points and store positions
            if ~isempty(individual_data{s, d})
                mouse_slopes = individual_data{s, d};
                mouse_names = mouse_ids{s, d};
                n_mice = length(mouse_slopes);
                
                jitter_amount = bar_width * 0;
                x_jitter = x_pos + (rand(n_mice, 1) - 0.5) * jitter_amount;
                
                mouse_positions{d, s} = struct('x_pos', x_jitter, 'y_pos', mouse_slopes, 'names', {mouse_names});
                
                scatter(x_jitter, mouse_slopes, 50, direction_session_colors{d, s}, 'o', ...
                       'filled', 'MarkerFaceAlpha', 0.7, ...
                       'MarkerEdgeColor', direction_session_colors{d, s} * 0.7, ...
                       'LineWidth', 1.5);
            end
        end
    end
    
    % Connect same mice across sessions within each direction group
    session_comparison_results = struct();
    
    for d = 1:length(directions)
        if length(sessions) == 2 && ~isempty(mouse_positions{d, 1}) && ~isempty(mouse_positions{d, 2})
            before_data = mouse_positions{d, 1};
            test_data = mouse_positions{d, 2};
            
            paired_before = [];
            paired_test = [];
            
            % Find matching mice between sessions
            for i = 1:length(before_data.names)
                mouse_name = before_data.names{i};
                test_idx = find(strcmp(test_data.names, mouse_name));
                
                if ~isempty(test_idx)
                    paired_before = [paired_before; before_data.y_pos(i)];
                    paired_test = [paired_test; test_data.y_pos(test_idx)];
                    
                    x_coords = [before_data.x_pos(i), test_data.x_pos(test_idx)];
                    y_coords = [before_data.y_pos(i), test_data.y_pos(test_idx)];
                    
                    plot(x_coords, y_coords, '-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1.0);
                end
            end
            
            % Perform paired t-test
            if length(paired_before) >= 3
                [~, p_value, ~, stats] = ttest(paired_test, paired_before);
                
                direction_name = directions{d};
                session_comparison_results.(direction_name) = struct();
                session_comparison_results.(direction_name).p_value = p_value;
                session_comparison_results.(direction_name).t_stat = stats.tstat;
                session_comparison_results.(direction_name).n_pairs = length(paired_before);
                session_comparison_results.(direction_name).is_significant = p_value < 0.05;
                session_comparison_results.(direction_name).mean_diff = mean(paired_test - paired_before);
                
                % Add significance indicator
                if p_value < 0.05
                    direction_max_vals = mean_values(:, d) + sem_values(:, d);
                    
                    if ~isempty(individual_data{1, d})
                        direction_max_vals = [direction_max_vals; max(individual_data{1, d})];
                    end
                    if ~isempty(individual_data{2, d})
                        direction_max_vals = [direction_max_vals; max(individual_data{2, d})];
                    end
                    
                    max_val_direction = max(direction_max_vals);
                    
                    star_x = x_positions(d);
                    star_y = max_val_direction + y_range * 0.08;
                    
                    bracket_y = max_val_direction + y_range * 0.06;
                    bracket_x1 = x_positions(d) + x_offset(1);
                    bracket_x2 = x_positions(d) + x_offset(2);
                    
                    if abs(bracket_x2 - bracket_x1) > 0.1
                        plot([bracket_x1, bracket_x1, bracket_x2, bracket_x2], ...
                             [bracket_y - y_range*0.01, bracket_y, bracket_y, bracket_y - y_range*0.01], ...
                             'k-', 'LineWidth', 1.5);
                    end
                    
                    if p_value < 0.001
                        star_text = '***';
                    elseif p_value < 0.01
                        star_text = '**';
                    else
                        star_text = '*';
                    end
                    
                    text(star_x, star_y, star_text, 'FontSize', 14, 'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'center', 'Color', 'k');
                    
                    if p_value < 0.001
                        p_text = 'p<0.001';
                    elseif p_value < 0.01
                        p_text = 'p<0.01';
                    else
                        p_text = sprintf('p=%.3f', p_value);
                    end
                    
                    text(star_x, star_y - y_range * 0.04, p_text, 'FontSize', 9, ...
                         'HorizontalAlignment', 'center', 'Color', 'k');
                end
            end
        end
    end
    
    % Add individual slope significance indicators
    for d = 1:length(directions)
        for s = 1:length(sessions)
            if slope_data.is_significant(s, d)
                x_pos = x_positions(d) + x_offset(s);
                y_pos = mean_values(s, d) + sem_values(s, d);
                y_star = y_pos + y_range * 0.3;
                
                text(x_pos, y_star, '★', 'FontSize', 12, ...
                     'HorizontalAlignment', 'center', 'Color', 'k', ...
                     'FontName', 'Arial Unicode MS');
            end
        end
    end
    
    % Enhanced formatting
    xlabel('Movement Direction', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Slope (z-scored dF/F per cm)', 'FontSize', 16, 'FontWeight', 'bold');
    
    set(gca, 'XTick', x_positions, 'XTickLabel', direction_labels, 'FontSize', 14);
    
    line([0.5, length(directions) + 0.5], [0, 0], 'Color', [0.3, 0.3, 0.3], 'LineStyle', '--', 'LineWidth', 1.5);
    
    % Enhanced legend
    legend_handles = [];
    legend_labels = {};
    
    for s = 1:length(sessions)
        legend_color = options.session_colors{s};
        h_legend = scatter(NaN, NaN, 100, legend_color, 'o', 'filled', ...
                          'MarkerEdgeColor', legend_color * 0.7, 'LineWidth', 1.5);
        legend_handles(s) = h_legend;
        legend_labels{s} = session_labels{s};
    end
    
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 14, ...
           'Box', 'off', 'EdgeColor', 'none');
    
    xlim([0.5, length(directions) + 0.5]);
    
    set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
    
    bottom_margin = y_range * 0.15;
    top_margin = y_range * 0.30;
    
    ylim([min_val - bottom_margin, max_val + top_margin]);
    
    % Add sample size information
    for d = 1:length(directions)
        for s = 1:length(sessions)
            if ~isempty(individual_data{s, d})
                n_mice = length(individual_data{s, d});
                x_pos = x_positions(d) + x_offset(s);
                y_pos = min_val - y_range * 0.08;
                
                text(x_pos, y_pos, sprintf('n=%d', n_mice), ...
                     'FontSize', 10, 'HorizontalAlignment', 'center', ...
                     'Color', direction_session_colors{d, s}, 'FontWeight', 'bold');
            end
        end
    end
    
    hold off;
    
    % Print session comparison results
    if ~isempty(fieldnames(session_comparison_results))
        fprintf('\n===== SESSION COMPARISON RESULTS =====\n');
        direction_names = fieldnames(session_comparison_results);
        for i = 1:length(direction_names)
            dir_name = direction_names{i};
            results = session_comparison_results.(dir_name);
            
            significance_text = '';
            if results.is_significant
                if results.p_value < 0.001
                    significance_text = ' ***';
                elseif results.p_value < 0.01
                    significance_text = ' **';
                else
                    significance_text = ' *';
                end
            end
            
            fprintf('%s: t(%.0f) = %.3f, p = %.4f, mean difference = %.4f%s\n', ...
                upper(dir_name), results.n_pairs-1, results.t_stat, results.p_value, ...
                results.mean_diff, significance_text);
        end
        fprintf('==========================================\n');
    end
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'cheeseboard_slope_barplot_enhanced.png'));
    end
end

function save_barplot_data_to_xlsx(slope_data, sessions, directions, session_labels, ...
        direction_labels, mean_values, sem_values, individual_data, mouse_ids, options)
    % save_barplot_data_to_xlsx - Exports barplot data to a two-sheet Excel file.
    %
    % Sheet 1 "Summary": group-level statistics (mean, SEM, n, p-value, significance).
    % Sheet 2 "Individual Mouse Data": per-mouse slopes in wide format, one column per group.
    
    xlsx_path = fullfile(options.figure_path, 'cheeseboard_slope_barplot_data.xlsx');
    
    % ---- Sheet 1: Summary ------------------------------------------------
    % Columns: Direction | Session | Mean Slope | SEM | N | P-value | Significant
    summary_header = {'Direction', 'Session', 'Mean Slope', 'SEM', 'N', 'P-value', 'Significant'};
    summary_rows = {};
    
    for d = 1:length(directions)
        for s = 1:length(sessions)
            n_mice = sum(~isnan(individual_data{s, d}));
            if isempty(individual_data{s, d})
                n_mice = slope_data.n_mice(s, d);
            end
            
            sig_str = 'No';
            if slope_data.is_significant(s, d)
                sig_str = 'Yes';
            end
            
            summary_rows{end+1, 1} = direction_labels{d};
            summary_rows{end,   2} = session_labels{s};
            summary_rows{end,   3} = mean_values(s, d);
            summary_rows{end,   4} = sem_values(s, d);
            summary_rows{end,   5} = n_mice;
            summary_rows{end,   6} = slope_data.p_values(s, d);
            summary_rows{end,   7} = sig_str;
        end
    end
    
    summary_table = [summary_header; summary_rows];
    writecell(summary_table, xlsx_path, 'Sheet', 'Summary');
    
    % ---- Sheet 2: Individual Mouse Data ----------------------------------
    % Build one column per (direction x session) group.
    % Header row: "Towards - Before", "Towards - Test", "Away - Before", "Away - Test"
    % Subsequent rows: individual mouse slopes (mouse ID | slope value pairs)
    
    % Determine column headers and gather data per group
    n_groups = length(directions) * length(sessions);
    group_headers_id   = cell(1, n_groups);  % e.g. "Towards - Before (Mouse ID)"
    group_headers_slope = cell(1, n_groups); % e.g. "Towards - Before (Slope)"
    group_mouse_ids    = cell(1, n_groups);
    group_slopes       = cell(1, n_groups);
    
    col = 0;
    for d = 1:length(directions)
        for s = 1:length(sessions)
            col = col + 1;
            label = sprintf('%s - %s', direction_labels{d}, session_labels{s});
            group_headers_id{col}    = sprintf('%s (Mouse ID)', label);
            group_headers_slope{col} = sprintf('%s (Slope)',    label);
            
            if ~isempty(individual_data{s, d})
                group_slopes{col}   = num2cell(individual_data{s, d}(:));
                group_mouse_ids{col} = mouse_ids{s, d}(:);
            else
                group_slopes{col}   = {};
                group_mouse_ids{col} = {};
            end
        end
    end
    
    % Interleave ID and Slope columns: ID_1, Slope_1, ID_2, Slope_2, ...
    n_header_cols = n_groups * 2;
    interleaved_header = cell(1, n_header_cols);
    for col = 1:n_groups
        interleaved_header{2*col - 1} = group_headers_id{col};
        interleaved_header{2*col}     = group_headers_slope{col};
    end
    
    % Find max number of mice across all groups
    max_mice = 0;
    for col = 1:n_groups
        max_mice = max(max_mice, length(group_slopes{col}));
    end
    
    % Build data rows
    indiv_rows = cell(max_mice, n_header_cols);
    for row = 1:max_mice
        for col = 1:n_groups
            id_col    = 2*col - 1;
            slope_col = 2*col;
            
            if row <= length(group_mouse_ids{col})
                indiv_rows{row, id_col}    = group_mouse_ids{col}{row};
                indiv_rows{row, slope_col} = group_slopes{col}{row};
            else
                indiv_rows{row, id_col}    = '';
                indiv_rows{row, slope_col} = NaN;
            end
        end
    end
    
    indiv_table = [interleaved_header; indiv_rows];
    writecell(indiv_table, xlsx_path, 'Sheet', 'Individual Mouse Data');
    
    fprintf('Barplot data saved to: %s\n', xlsx_path);
end

function create_cheeseboard_summary_plot(slope_results, options)
    % Create summary visualization showing key statistics
    
    figure('Name', 'Cheeseboard Slope Analysis: Summary', 'Position', [100, 100, 400, 600]);
    
    summary = slope_results.summary;
    
    subplot(2, 3, [1, 2]);
    axis off;
    
    % Title
    text(0.5, 0.95, 'CHEESEBOARD SLOPE SUMMARY', 'FontSize', 16, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
    
    % Key statistics
    y_pos = 0.85;
    line_spacing = 0.08;
    
    text(0.05, y_pos, sprintf('Total Analyses: %d', summary.total_analyses), ...
         'FontSize', 12, 'FontWeight', 'bold');
    y_pos = y_pos - line_spacing;
    
    text(0.05, y_pos, sprintf('Significant Slopes: %d (%.1f%%)', ...
         summary.significant_slopes, 100*summary.significant_slopes/summary.total_analyses), ...
         'FontSize', 12);
    y_pos = y_pos - line_spacing;
    
    text(0.05, y_pos, sprintf('Expected Direction: %d (%.1f%%)', ...
         summary.correct_direction, 100*summary.correct_direction/summary.total_analyses), ...
         'FontSize', 12);
    y_pos = y_pos - line_spacing;
    
    text(0.05, y_pos, sprintf('Mean Slope: %.4f ± %.4f', ...
         summary.mean_slope, summary.std_slope), 'FontSize', 12);
    y_pos = y_pos - line_spacing;
    
    text(0.05, y_pos, sprintf('Mean R²: %.3f', summary.mean_r_squared), ...
         'FontSize', 12);
    y_pos = y_pos - line_spacing;
    
    text(0.05, y_pos, sprintf('Median p-value: %.4f', summary.median_p_value), ...
         'FontSize', 12);
    
    % Individual results
    subplot(2, 3, [4, 5, 6]);
    
    detailed_fields = fieldnames(slope_results.detailed);
    n_results = length(detailed_fields);
    
    if n_results > 0
        y_positions = 1:n_results;
        colors = [];
        slopes = [];
        labels = {};
        
        for i = 1:n_results
            field = detailed_fields{i};
            stats = slope_results.detailed.(field);
            
            slopes(i) = stats.slope;
            labels{i} = strrep(field, '_', ' ');
            
            % Color by significance and direction
            if stats.is_significant && stats.correct_direction
                colors(i, :) = [0, 0.7, 0]; % Green
            elseif stats.is_significant
                colors(i, :) = [0.7, 0.7, 0]; % Yellow
            else
                colors(i, :) = [0.7, 0, 0]; % Red
            end
        end
        
        barh(y_positions, slopes, 'FaceColor', 'flat', 'CData', colors);
        set(gca, 'YTick', y_positions, 'YTickLabel', labels);
        xlabel('Slope Value');
        title('Individual Slope Results');
        grid on;
        
        % Add vertical line at zero
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    end
    
    % Legend for colors
    subplot(2, 3, 3);
    axis off;
    
    text(0.1, 0.8, 'LEGEND:', 'FontSize', 12, 'FontWeight', 'bold');
    
    rectangle('Position', [0.1, 0.65, 0.05, 0.05], 'FaceColor', [0, 0.7, 0]);
    text(0.2, 0.675, 'Significant & Expected', 'FontSize', 10);
    
    rectangle('Position', [0.1, 0.55, 0.05, 0.05], 'FaceColor', [0.7, 0.7, 0]);
    text(0.2, 0.575, 'Significant & Unexpected', 'FontSize', 10);
    
    rectangle('Position', [0.1, 0.45, 0.05, 0.05], 'FaceColor', [0.7, 0, 0]);
    text(0.2, 0.475, 'Not Significant', 'FontSize', 10);
    
    if options.save_plots
        saveas(gcf, fullfile(options.figure_path, 'cheeseboard_slope_summary.png'));
    end
end