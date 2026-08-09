function plot_combined_nonfood_compare_sessions(run_data, options)
    % plot_combined_nonfood_compare_sessions - Creates a figure comparing sessions for food vs combined non-food arms
    % with flexible layout based on number of sessions
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_mouse_runs
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to include (e.g., {'sess0', 'sess1', 'sess2'})
    %       .ylim - Y-axis limits [min max], default [-2 2]
    %       .smoothing - Window size for smoothing, default 5
    %       .figure_width - Width of the figure in pixels (default: 1200)
    %       .figure_height - Height of the figure in pixels (default: 800)
    %       .title - Main figure title (optional, auto-generated if not provided)
    %       .show_title - Whether to show the overall figure title (default: true)
    %       .plot_sem - Whether to plot SEM shading, default true
    %       .max_distance - Maximum distance to plot, default 200
    %       .run_stats - Whether to perform statistical tests, default true
    %       .stat_alpha - Alpha level for significance, default 0.05
    %       .min_mice_per_bin - Minimum mice per bin for stats, default 3
    %       .colors - Colors for sessions {Before, Learning, Test}, default gray/blue scheme
    %       .axis_label_font_size - Font size for axis labels (default: 12)
    %       .tick_font_size - Font size for axis tick numbers (default: 10)
    %       .subplot_title_font_size - Font size for subplot titles (default: 14)
    %       .main_title_font_size - Font size for main figure title (default: 16)
    %       .x_tick_spacing - Spacing between x-axis ticks in cm (default: 20)
    %       .y_tick_spacing - Spacing between y-axis ticks (default: automatic)
    
    % Default options
    if nargin < 2
        options = struct();
    end
    
    % Set default options if not provided
    if ~isfield(options, 'sessions')
        % Find all unique session numbers in run_data
        all_sessions = [];
        for m = 1:length(run_data)
            if ~ismember(run_data(m).session, all_sessions)
                all_sessions = [all_sessions, run_data(m).session];
            end
        end
        sessions_to_plot_nums = sort(all_sessions);
        
        % Convert numbers to strings
        options.sessions = cell(1, length(sessions_to_plot_nums));
        for i = 1:length(sessions_to_plot_nums)
            options.sessions{i} = ['sess' num2str(sessions_to_plot_nums(i))];
        end
    end
    
    if ~isfield(options, 'ylim'), options.ylim = [-2 2]; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'figure_width'), options.figure_width = 1200; end
    if ~isfield(options, 'figure_height'), options.figure_height = 800; end
    if ~isfield(options, 'show_title'), options.show_title = true; end
    if ~isfield(options, 'plot_sem'), options.plot_sem = true; end
    if ~isfield(options, 'max_distance') || isempty(options.max_distance)
        options.max_distance = 200; % Default max distance
    end
    if ~isfield(options, 'run_stats') || isempty(options.run_stats)
        options.run_stats = true; % Default: run statistics
    end
    if ~isfield(options, 'stat_alpha') || isempty(options.stat_alpha)
        options.stat_alpha = 0.05; % Default alpha level
    end
    if ~isfield(options, 'min_mice_per_bin') || isempty(options.min_mice_per_bin)
        options.min_mice_per_bin = 3; % Minimum mice per bin for statistics
    end
    
    % Font size options
    if ~isfield(options, 'axis_label_font_size') || isempty(options.axis_label_font_size)
        options.axis_label_font_size = 12; % Default axis labels font size
    end
    if ~isfield(options, 'tick_font_size') || isempty(options.tick_font_size)
        options.tick_font_size = 10; % Default tick labels font size
    end
    if ~isfield(options, 'subplot_title_font_size') || isempty(options.subplot_title_font_size)
        options.subplot_title_font_size = 14; % Default subplot title font size
    end
    if ~isfield(options, 'main_title_font_size') || isempty(options.main_title_font_size)
        options.main_title_font_size = 16; % Default main title font size
    end
    
    % Tick spacing options
    if ~isfield(options, 'x_tick_spacing') || isempty(options.x_tick_spacing)
        options.x_tick_spacing = 20; % Default x-axis tick spacing in cm
    end
    if ~isfield(options, 'y_tick_spacing') || isempty(options.y_tick_spacing)
        options.y_tick_spacing = []; % Default: auto y-axis tick spacing
    end
    
    % Create a mapping from session string to index and vice versa
    session_map = containers.Map();
    session_numbers = [];
    for i = 1:length(options.sessions)
        session_str = options.sessions{i};
        session_num = str2double(session_str(5:end));
        session_map(session_str) = i; % Map session string to sequential index
        session_numbers = [session_numbers, session_num];
    end
    
    % Sort session numbers for consistent ordering
    [session_numbers, sort_idx] = sort(session_numbers);
    options.sessions = options.sessions(sort_idx);
    
    % Set default colors for sessions if not provided
    if ~isfield(options, 'colors')
        % Default colors: Before (gray), Learning (blue), Test (light blue)
        default_colors = {[0.6 0.6 0.6], [0.2 0.4 0.8], [0.1 0.7 0.9]};
        options.colors = default_colors(1:length(options.sessions));
    else
        % Make sure we have enough colors for all sessions
        if length(options.colors) < length(options.sessions)
            warning('Not enough colors provided for all sessions. Using defaults for missing colors.');
            default_colors = {[0.6 0.6 0.6], [0.2 0.4 0.8], [0.1 0.7 0.9], [0.8 0.4 0.4], [0.4 0.8 0.4]};
            for i = length(options.colors)+1:length(options.sessions)
                options.colors{i} = default_colors{mod(i-1, length(default_colors))+1};
            end
        end
    end
    
    % Set default title
    if ~isfield(options, 'title') || isempty(options.title)
        options.title = 'Z-scored dF/F vs Distance - Session Comparison (Food vs Non-Food)';
    end
    
    % Display which sessions are being plotted
    fprintf('Comparing dF/F across sessions with combined non-food arms: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf(' (max distance: %.0f)\n', options.max_distance);
    
    % Create session name mapping for plot labels
    session_labels = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        sess_num = session_numbers(i);
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

    % Create the figure
    fig = figure('Name', 'dF/F vs Distance - Session Comparison (Combined Non-Food)', ...
        'Position', [100, 100, options.figure_width, options.figure_height]);

    % Subplot titles for 2x2 layout (Food/Non-Food x Towards/Away)
    arm_titles = {'Food Arm', 'Non-Food Arms'};
    direction_titles = {'Towards', 'Away'};
    
    % Define arms and run types (simplified to food vs combined non-food)
    arm_types = {'food', 'nonfood'};
    type_names = {'towards', 'away'};
    
    % Use the specified max distance
    max_distance = options.max_distance;
    fprintf('Using maximum distance: %.0f\n', max_distance);
    
    % Create distance bins
    bin_width = 1; % 1 unit bins for better resolution
    dist_bins = 0:bin_width:max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialize data structure to collect per-mouse averages for each bin
    all_data = struct();
    mouse_bin_data = struct(); % For statistical testing - store per-mouse averages
    
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        for t = 1:length(type_names)
            for a = 1:length(arm_types)
                type = type_names{t};
                arm = arm_types{a};
                key = sprintf('sess%d_%s_%s', sess, type, arm);
                all_data.(key) = cell(length(dist_bins)-1, 1);
                mouse_bin_data.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    all_data.(key){i} = [];
                    mouse_bin_data.(key){i} = containers.Map(); % Mouse ID -> mean signal values
                end
            end
        end
    end
    
    % Process each mouse
    for m = 1:length(run_data)
        mouse_id = run_data(m).mouse_id;
        mouse_sess = run_data(m).session;
        
        % Skip if session not in sessions to plot
        sess_str = ['sess' num2str(mouse_sess)];
        if ~ismember(sess_str, options.sessions)
            continue;
        end
        
        runs = run_data(m).runs;
        
        % Initialize temporary storage for this mouse's binned data
        mouse_temp_bins = struct();
        for t = 1:length(type_names)
            for a = 1:length(arm_types)
                type = type_names{t};
                arm = arm_types{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, type, arm);
                mouse_temp_bins.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_temp_bins.(key){i} = [];
                end
            end
        end
        
        % Collect all dF/F values for this mouse to calculate z-score parameters
        all_dff_values = [];
        for r = 1:length(runs)
            all_dff_values = [all_dff_values; runs(r).dff];
        end
        
        % Calculate mean and std for z-scoring
        dff_mean = mean(all_dff_values);
        dff_std = std(all_dff_values);
        
        % Skip mouse if std is zero or very small (would cause division by zero)
        if dff_std < 1e-10
            warning('Mouse %s has near-zero dF/F standard deviation, skipping.', mouse_id);
            continue;
        end
        
        % Process runs for this mouse
        for r = 1:length(runs)
            run = runs(r);
            
            % Determine if this is food or non-food arm
            if strcmp(run.arm, 'food')
                combined_arm = 'food';
            elseif ismember(run.arm, {'nonfood1', 'nonfood2'})
                combined_arm = 'nonfood'; % Combine both non-food arms
            else
                continue; % Skip unrecognized arms
            end
            
            % Skip if not a recognized type
            if ~ismember(run.type, type_names)
                continue;
            end
            
            % Apply distance limit - only include data points within max_distance
            valid_distance_idx = run.distance <= max_distance;
            if ~any(valid_distance_idx)
                continue; % Skip this run if no points within distance limit
            end
            
            % Filter the run data to only include points within distance limit
            filtered_distance = run.distance(valid_distance_idx);
            filtered_dff = run.dff(valid_distance_idx);
            
            % Z-score this run's dF/F
            z_scored_dff = (filtered_dff - dff_mean) / dff_std;
            
            % Bin the z-scored data
            for i = 1:length(dist_bins)-1
                indices = filtered_distance >= dist_bins(i) & filtered_distance < dist_bins(i+1);
                if any(indices)
                    key = sprintf('sess%d_%s_%s', mouse_sess, run.type, combined_arm);
                    all_data.(key){i} = [all_data.(key){i}; z_scored_dff(indices)];
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; z_scored_dff(indices)];
                end
            end
        end
        
        % Calculate per-mouse averages for statistical testing
        for t = 1:length(type_names)
            for a = 1:length(arm_types)
                type = type_names{t};
                arm = arm_types{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, type, arm);
                
                for i = 1:length(dist_bins)-1
                    if ~isempty(mouse_temp_bins.(key){i})
                        mouse_bin_avg = mean(mouse_temp_bins.(key){i});
                        
                        % If this mouse already has data in this bin, average with existing
                        if mouse_bin_data.(key){i}.isKey(mouse_id)
                            existing_avg = mouse_bin_data.(key){i}(mouse_id);
                            mouse_bin_data.(key){i}(mouse_id) = (existing_avg + mouse_bin_avg) / 2;
                        else
                            mouse_bin_data.(key){i}(mouse_id) = mouse_bin_avg;
                        end
                    end
                end
            end
        end
    end
    
    % Create 2×2 subplot layout (2 arms, 2 directions) and plot data
    for a = 1:length(arm_types)
        arm = arm_types{a};
        
        for t = 1:length(type_names)
            type = type_names{t};
            
            % Calculate subplot index
            sp_idx = (a-1)*2 + t;
            subplot(2, 2, sp_idx);
            hold on;
            
            % Plot all selected sessions
            for s = 1:length(session_numbers)
                sess = session_numbers(s);
                
                % Create the key for this combination
                key = sprintf('sess%d_%s_%s', sess, type, arm);
                
                % Skip if this key doesn't exist in all_data
                if ~isfield(all_data, key)
                    continue;
                end
                
                % Calculate mean and SEM for each distance bin
                means = nan(length(bin_centers), 1);
                sems = nan(length(bin_centers), 1);
                
                for i = 1:length(dist_bins)-1
                    bin_data = all_data.(key){i};
                    if ~isempty(bin_data)
                        means(i) = mean(bin_data);
                        sems(i) = std(bin_data) / sqrt(length(bin_data));
                    end
                end
                
                % Apply smoothing to the means and SEMs
                valid = ~isnan(means);
                if sum(valid) > options.smoothing
                    % Get the valid data only
                    x_valid = bin_centers(valid);
                    y_valid = means(valid);
                    sem_valid = sems(valid);
                    
                    % For 'towards' runs, make x-coordinates negative (distance to food)
                    if strcmp(type, 'towards')
                        x_valid = -x_valid; % Convert to negative distances
                    end
                    
                    % Apply moving average smoothing
                    y_smoothed = movmean(y_valid, options.smoothing);
                    sem_smoothed = movmean(sem_valid, options.smoothing);
                    
                    % Get session label and color
                    sess_label = session_labels{s};
                    sess_color = options.colors{s};
                    
                    % Plot SEM shading if requested
                    if options.plot_sem
                        for i = 1:length(x_valid)-1
                            x_patch = [x_valid(i), x_valid(i+1), x_valid(i+1), x_valid(i)];
                            y_patch = [y_smoothed(i)-sem_smoothed(i), y_smoothed(i+1)-sem_smoothed(i+1), ...
                                      y_smoothed(i+1)+sem_smoothed(i+1), y_smoothed(i)+sem_smoothed(i)];
                            patch(x_patch, y_patch, sess_color, 'EdgeColor', 'none', ...
                                  'FaceAlpha', 0.2, 'HandleVisibility', 'off');
                        end
                    end
                    
                    % Plot the smoothed line with session label
                    plot(x_valid, y_smoothed, 'Color', sess_color, 'LineWidth', 3, 'DisplayName', sess_label);
                end
            end
            
            % Perform point-by-point statistical tests and add significance bars if requested
            if options.run_stats && length(session_numbers) >= 2
                addCombinedSessionSignificanceBars(session_numbers, type, arm, mouse_bin_data, bin_centers, options);
            end
            
            % Set subplot title and labels
            t = title(sprintf('%s - %s', arm_titles{a}, direction_titles{t}), 'FontWeight', 'bold');
            t.FontSize = options.subplot_title_font_size;
            
            if strcmp(type, 'towards')
                if a == length(arm_types) % Only add xlabel for bottom row
                    xlabel('Distance to Food', 'FontSize', options.axis_label_font_size);
                end
            else
                if a == length(arm_types) % Only add xlabel for bottom row
                    xlabel('Distance from Food', 'FontSize', options.axis_label_font_size);
                end
            end
            ylabel('Z-scored dF/F', 'FontSize', options.axis_label_font_size);
            
            % Set y-axis limits
            ylim(options.ylim);
            
            % Set y-axis ticks if specified
            if ~isempty(options.y_tick_spacing)
                y_min = options.ylim(1);
                y_max = options.ylim(2);
                y_ticks = y_min:options.y_tick_spacing:y_max;
                yticks(y_ticks);
            end
            
            % Set x-axis limits and ticks based on direction
            if strcmp(type, 'towards')
                xlim([-max_distance, 0]);
                % Set x-axis ticks
                x_ticks = -max_distance:options.x_tick_spacing:0;
                xticks(x_ticks);
            else
                xlim([0, max_distance]);
                % Set x-axis ticks
                x_ticks = 0:options.x_tick_spacing:max_distance;
                xticks(x_ticks);
            end
            
            % Set axis tick font size and ensure axis labels keep their font size
            set(gca, 'FontSize', options.tick_font_size);
            ax = gca;
            ax.XLabel.FontSize = options.axis_label_font_size;
            ax.YLabel.FontSize = options.axis_label_font_size;
            
            % Add horizontal line at y=0
            line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            % Add legend
            legend('show', 'Location', 'best');
            legend('boxoff');
            
            % Make legend lines shorter
            leg = legend;
            leg.ItemTokenSize = [15, 8];
            grid off;
        end
    end
    
    % Add overall super title
    if options.show_title
        sgtitle(options.title, 'FontSize', options.main_title_font_size, 'FontWeight', 'bold');
    end
end

function addCombinedSessionSignificanceBars(session_numbers, type, arm, mouse_bin_data, bin_centers, options)
    % Helper function to add point-by-point significance bars comparing sessions
    % Uses ANOVA + post-hoc tests for proper multi-session comparisons
    
    % Get current y-limits to position significance bars
    yl = ylim;
    y_range = yl(2) - yl(1);
    
    % Position bars at the top of the plot
    bar_y_base = yl(2) - 0.05 * y_range;
    bar_height = 0.01 * y_range;
    
    % Initialize arrays to store significance results for each bin
    num_bins = length(bin_centers);
    
    % Generate all pairwise comparisons
    comparisons = {};
    comparison_names = {};
    comparison_colors = {};
    color_options = {[0.8, 0, 0.8], [0, 0.8, 0.8], [0.8, 0.8, 0], [0.8, 0.4, 0], [0.4, 0.8, 0], [0.4, 0, 0.8]};
    
    % Build all pairwise comparisons with specific colors
    for i = 1:length(session_numbers)-1
        for j = i+1:length(session_numbers)
            sess1 = session_numbers(i);
            sess2 = session_numbers(j);
            comparisons{end+1} = {sess1, sess2};
            
            % Assign specific colors based on session comparison
            if (sess1 == 0 && sess2 == 2) || (sess1 == 2 && sess2 == 0)
                % S0 vs S2 - cian
                comparison_colors{end+1} = [0, 0, 1];
            elseif (sess1 == 0 && sess2 == 1) || (sess1 == 1 && sess2 == 0)
                % S0 vs S1 - Magenta
                comparison_colors{end+1} = [0.8, 0, 0.8];
            elseif (sess1 == 1 && sess2 == 2) || (sess1 == 2 && sess2 == 1)
                % S1 vs S2 - Gray
                comparison_colors{end+1} = [0.5, 0.5, 0.5];
            else
                % Default color for any other combinations
                comparison_colors{end+1} = [0.3, 0.3, 0.3];
            end
            
            % No comparison names needed since we won't display text
            comparison_names{end+1} = '';
        end
    end
    
    if isempty(comparisons)
        return;
    end
    
    % Store significant bins for each comparison
    all_sig_bins = cell(length(comparisons), 1);
    all_sig_x_coords = cell(length(comparisons), 1);
    
    % Test each bin independently using ANOVA + post-hoc
    for bin_idx = 1:num_bins
        % Collect data for all sessions at this bin
        session_data = cell(1, length(session_numbers));
        session_counts = zeros(1, length(session_numbers));
        
        for s = 1:length(session_numbers)
            sess = session_numbers(s);
            key = sprintf('sess%d_%s_%s', sess, type, arm);
            
            if isfield(mouse_bin_data, key) && bin_idx <= length(mouse_bin_data.(key))
                data_map = mouse_bin_data.(key){bin_idx};
                if data_map.Count >= options.min_mice_per_bin
                    session_data{s} = cell2mat(values(data_map));
                    session_counts(s) = length(session_data{s});
                else
                    session_data{s} = [];
                end
            else
                session_data{s} = [];
            end
        end
        
        % Check if we have enough data for ANOVA
        sessions_with_data = find(session_counts >= options.min_mice_per_bin);
        
        if length(sessions_with_data) >= 2
            % Prepare data for ANOVA
            all_values = [];
            group_labels = [];
            
            for s = sessions_with_data
                all_values = [all_values; session_data{s}(:)];
                group_labels = [group_labels; s * ones(length(session_data{s}), 1)];
            end
            
            % Perform one-way ANOVA for this bin
            try
                [p_anova, ~, stats] = anova1(all_values, group_labels, 'off');
                
                % If ANOVA is significant, perform post-hoc tests
                if p_anova < options.stat_alpha && length(sessions_with_data) >= 2
                    [comparisons_results, ~, ~, ~] = multcompare(stats, 'Alpha', options.stat_alpha, 'Display', 'off');
                    
                    % Check each specified comparison
                    for comp_idx = 1:length(comparisons)
                        sess1 = comparisons{comp_idx}{1};
                        sess2 = comparisons{comp_idx}{2};
                        
                        % Find indices of these sessions in our data
                        sess1_idx = find(session_numbers == sess1);
                        sess2_idx = find(session_numbers == sess2);
                        
                        % Check if both sessions have data and are in the ANOVA
                        sess1_anova_idx = find(sessions_with_data == sess1_idx);
                        sess2_anova_idx = find(sessions_with_data == sess2_idx);
                        
                        if ~isempty(sess1_anova_idx) && ~isempty(sess2_anova_idx)
                            % Find this comparison in the multcompare results
                            comp_row = find((comparisons_results(:,1) == sess1_anova_idx & comparisons_results(:,2) == sess2_anova_idx) | ...
                                          (comparisons_results(:,1) == sess2_anova_idx & comparisons_results(:,2) == sess1_anova_idx));
                            
                            if ~isempty(comp_row)
                                p_posthoc = comparisons_results(comp_row(1), 6); % p-value column
                                
                                if p_posthoc < options.stat_alpha
                                    % This comparison is significant at this bin
                                    all_sig_bins{comp_idx} = [all_sig_bins{comp_idx}, bin_idx];
                                    
                                    % Calculate x coordinate for this bin
                                    x_coord = bin_centers(bin_idx);
                                    if strcmp(type, 'towards')
                                        x_coord = -x_coord;
                                    end
                                    all_sig_x_coords{comp_idx} = [all_sig_x_coords{comp_idx}, x_coord];
                                end
                            end
                        end
                    end
                end
            catch ME
                % Skip if statistical test fails
                continue;
            end
        end
    end
    
    % Draw significance bars for each comparison
    max_bar_level = 0;
    
    for comp_idx = 1:length(comparisons)
        sig_bins = all_sig_bins{comp_idx};
        sig_x_coords = all_sig_x_coords{comp_idx};
        
        if ~isempty(sig_bins)
            % Find consecutive runs of significant bins
            if length(sig_bins) == 1
                run_starts = 1;
                run_ends = 1;
            else
                diff_bins = diff(sig_bins);
                break_points = find(diff_bins > 1);
                
                if isempty(break_points)
                    run_starts = 1;
                    run_ends = length(sig_bins);
                else
                    run_starts = [1, break_points + 1];
                    run_ends = [break_points, length(sig_bins)];
                end
            end
            
            % Draw a bar for each consecutive run
            for run_idx = 1:length(run_starts)
                start_idx = run_starts(run_idx);
                end_idx = run_ends(run_idx);
                
                if start_idx <= length(sig_x_coords) && end_idx <= length(sig_x_coords)
                    x_start = sig_x_coords(start_idx);
                    x_end = sig_x_coords(end_idx);
                    
                    % Extend the bar slightly beyond the actual significant points
                    bin_width = abs(bin_centers(2) - bin_centers(1));
                    if strcmp(type, 'towards')
                        bin_width = -bin_width;
                    end
                    x_start = x_start - bin_width/2;
                    x_end = x_end + bin_width/2;
                    
                    % Position this comparison's bar
                    bar_y = bar_y_base - (comp_idx - 1) * (bar_height + 0.005 * y_range);
                    max_bar_level = max(max_bar_level, comp_idx);
                    
                    % Draw the significance bar
                    plot([x_start, x_end], [bar_y, bar_y], 'Color', comparison_colors{comp_idx}, ...
                         'LineWidth', 3, 'HandleVisibility', 'off');
                    
                    % No text labels - just the colored bars
                end
            end
            
            % Print results to console
            sess1 = comparisons{comp_idx}{1};
            sess2 = comparisons{comp_idx}{2};
            fprintf('%s, %s: Session %d vs Session %d significant in %d bins (ANOVA + post-hoc)\n', ...
                arm, type, sess1, sess2, length(sig_bins));
        end
    end
    
    % Adjust y-limits to accommodate all significance bars
    if max_bar_level > 0
        new_ylim = [yl(1), yl(2) + 0.1 * y_range + max_bar_level * (bar_height + 0.005 * y_range)];
        ylim(new_ylim);
    end

end


