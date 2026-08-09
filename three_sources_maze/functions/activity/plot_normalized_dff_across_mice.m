function plot_normalized_dff_across_mice(run_data, options)
    % plot_normalized_dff_across_mice - Creates plots showing signals vs distance across all mice
    % with point-by-point statistical comparisons along trace length
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_mouse_runs
    %   options - Structure with the following fields (all optional):
    %       .sessions - Cell array of sessions to analyze (e.g., {'sess1', 'sess2'})
    %       .smoothing - Window size for smoothing
    %       .y_limits - Limits for y-axis [min max]
    %       .title - Custom title for the figure
    %       .subplot_titles - Custom subplot titles (cell array)
    %       .bin_width - Width of distance bins (default: 1)
    %       .signal_type - Signal to analyze ('dff', '465', '405'), default 'dff'
    %       .max_distance - Maximum distance to plot (default: 200)
    %       .apply_zscore - Whether to apply z-scoring (default: true)
    %       .run_stats - Whether to perform statistical tests (default: true)
    %       .stat_alpha - Alpha level for significance (default: 0.05)
    %       .min_mice_per_bin - Minimum mice per bin for stats (default: 3)
    %       .stat_window - Window size for statistical smoothing (default: 5)
    
    % Set default options if not provided
    if nargin < 2
        options = struct();
    end
    
    % Extract and validate options, set defaults if not specified
    if ~isfield(options, 'smoothing') || isempty(options.smoothing)
        options.smoothing = 5; % Default smoothing window
    end
    
    if ~isfield(options, 'y_limits') || isempty(options.y_limits)
        options.y_limits = [-1, 1]; % Default y-limits
    end
    
    if ~isfield(options, 'bin_width') || isempty(options.bin_width)
        options.bin_width = 1; % Default bin width
    end

    if ~isfield(options, 'signal_type') || isempty(options.signal_type)
        options.signal_type = 'dff'; % Default signal type
    end

    if ~isfield(options, 'max_distance') || isempty(options.max_distance)
        options.max_distance = 200; % Default max distance
    end
    
    if ~isfield(options, 'apply_zscore') || isempty(options.apply_zscore)
        options.apply_zscore = true; % Default: apply z-scoring
    end
    
    % Statistical options
    if ~isfield(options, 'run_stats') || isempty(options.run_stats)
        options.run_stats = true; % Default: run statistics
    end
    
    if ~isfield(options, 'stat_alpha') || isempty(options.stat_alpha)
        options.stat_alpha = 0.05; % Default alpha level
    end
    
    if ~isfield(options, 'min_mice_per_bin') || isempty(options.min_mice_per_bin)
        options.min_mice_per_bin = 3; % Minimum mice per bin for statistics
    end
    
    if ~isfield(options, 'stat_window') || isempty(options.stat_window)
        options.stat_window = 5; % Statistical smoothing window
    end
    
    % Process signal type option
    signal_type = lower(options.signal_type);
    switch signal_type
        case 'dff'
            signal_field = 'dff';
            signal_label = 'dF/F';
        case '465'
            signal_field = 'signal465';
            signal_label = '465nm';
        case '405'
            signal_field = 'signal405';
            signal_label = '405nm';
        otherwise
            warning('Unrecognized signal type "%s". Using dF/F instead.', signal_type);
            signal_field = 'dff';
            signal_label = 'dF/F';
    end
    
    % Set default title based on signal type and z-scoring
    if ~isfield(options, 'title') || isempty(options.title)
        if options.apply_zscore
            options.title = sprintf('Z-scored %s vs Distance - Population Average', signal_label);
        else
            options.title = sprintf('%s vs Distance - Population Average', signal_label);
        end
    end
    
    % Check if sessions option is provided, otherwise detect from data
    if ~isfield(options, 'sessions') || isempty(options.sessions)
        % Find all unique session numbers in run_data
        all_sessions = [];
        for m = 1:length(run_data)
            if ~ismember(run_data(m).session, all_sessions)
                all_sessions = [all_sessions, run_data(m).session];
            end
        end
        sessions_to_plot_nums = sort(all_sessions);
        
        % Convert numbers to session names
        options.sessions = cell(1, length(sessions_to_plot_nums));
        for i = 1:length(sessions_to_plot_nums)
            options.sessions{i} = ['sess' num2str(sessions_to_plot_nums(i))];
        end
    end
    
    % Create a mapping from session string to number and vice versa
    session_map = containers.Map();
    session_numbers = [];
    for i = 1:length(options.sessions)
        session_str = options.sessions{i};
        session_num = str2double(session_str(5:end));
        session_map(session_str) = session_num;
        session_numbers = [session_numbers, session_num];
    end
    
    % Sort session numbers for consistent ordering
    session_numbers = sort(session_numbers);
    
    % Display which sessions are being plotted with z-scoring information
    if options.apply_zscore
        fprintf('Plotting z-scored %s across mice for sessions: ', signal_label);
    else
        fprintf('Plotting raw %s across mice for sessions: ', signal_label);
    end
    
    for i = 1:length(session_numbers)
        fprintf('sess%d ', session_numbers(i));
    end
    fprintf('\n');
    
    % Create session name mapping for plot titles
    session_names = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        sess_num = session_numbers(i);
        session_names{i} = sprintf('Session %d', sess_num);
    end
    
    % Calculate number of rows and columns for subplots
    num_sessions = length(session_numbers);
    num_types = 2; % 'towards' and 'away'
    num_subplots = num_sessions * num_types;
    
    % Create subplot layout
    rows = num_sessions;
    cols = num_types;
    
    % Create subplot titles if not provided
    if ~isfield(options, 'subplot_titles') || isempty(options.subplot_titles)
        options.subplot_titles = cell(num_subplots, 1);
        for s = 1:num_sessions
            for t = 1:2
                subplot_idx = (s-1)*2 + t;
                type_label = {'Towards', 'Away'};
                options.subplot_titles{subplot_idx} = sprintf('%s %s', session_names{s}, type_label{t});
            end
        end
    end
    
    % Create the figure
    figure_width = min(1200, 600 * cols);
    figure_height = min(1000, 400 * rows);
    figure('Name', options.title, 'Position', [100, 100, figure_width, figure_height]);
    
    % Define color gradients for each arm and type
    arm_colors = struct();
    
    % Blue gradient for 'towards' runs
    arm_colors.towards_food = [0, 0, 0.8];         % Dark blue
    arm_colors.towards_nonfood1 = [0.1, 0.7, 1];   % Medium blue
    arm_colors.towards_nonfood2 = [0.6, 0.8, 1];   % Light blue
    
    % Red gradient for 'away' runs
    arm_colors.away_food = [0.8, 0, 0];           % Dark red
    arm_colors.away_nonfood1 = [1, 0.5, 0.5];     % Medium red
    arm_colors.away_nonfood2 = [1, 0.8, 0.8];     % Light red
    
    % Define arms in plotting order (nonfood first, food last for better visibility)
    arm_names = {'food', 'nonfood1', 'nonfood2'};
    arm_plot_order = {'nonfood2', 'nonfood1', 'food'}; % Plot order (food last to be on top)
    type_names = {'towards', 'away'};
    
    % Use user-specified max distance
    max_distance = options.max_distance;
    fprintf('Using maximum distance: %.2f\n', max_distance);
    
    % Create distance bins
    bin_width = options.bin_width;
    dist_bins = 0:bin_width:max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialize data structure to collect per-mouse averages for each bin
    mouse_bin_data = struct();
    
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        for t = 1:length(type_names)
            for a = 1:length(arm_names)
                type = type_names{t};
                arm = arm_names{a};
                key = sprintf('sess%d_%s_%s', sess, type, arm);
                % Each bin will contain a containers.Map: mouse_id -> mean_value_for_that_bin
                mouse_bin_data.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_bin_data.(key){i} = containers.Map();
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
        
        % Collect all signal values for this mouse for potential z-scoring
        all_signal_values = [];
        for r = 1:length(runs)
            % Skip if run doesn't have the requested signal field
            if ~isfield(runs(r), signal_field)
                warning('Run %d for mouse %s does not have signal field "%s". Skipping.', ...
                    r, mouse_id, signal_field);
                continue;
            end
            
            all_signal_values = [all_signal_values; runs(r).(signal_field)];
        end
        
        % Skip mouse if no valid runs found
        if isempty(all_signal_values)
            warning('Mouse %s has no valid runs with signal field "%s". Skipping.', ...
                mouse_id, signal_field);
            continue;
        end
        
        % Calculate mean and std for z-scoring if needed
        if options.apply_zscore
            signal_mean = mean(all_signal_values);
            signal_std = std(all_signal_values);
            
            % Skip mouse if std is zero or very small (would cause division by zero)
            if signal_std < 1e-10
                warning('Mouse %s has near-zero %s standard deviation, skipping.', mouse_id, signal_label);
                continue;
            end
        end
        
        % Initialize temporary storage for this mouse's binned data
        mouse_temp_bins = struct();
        for t = 1:length(type_names)
            for a = 1:length(arm_names)
                type = type_names{t};
                arm = arm_names{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, type, arm);
                mouse_temp_bins.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    mouse_temp_bins.(key){i} = [];
                end
            end
        end
        
        % Process runs for this mouse
        for r = 1:length(runs)
            run = runs(r);
            
            % Skip if not a recognized arm or type
            if ~ismember(run.arm, arm_names) || ~ismember(run.type, type_names)
                continue;
            end
            
            % Skip if run doesn't have the requested signal field
            if ~isfield(run, signal_field)
                continue;
            end
            
            % Apply z-scoring if requested, otherwise use raw signal
            if options.apply_zscore
                % Z-score this run's signal
                processed_signal = (run.(signal_field) - signal_mean) / signal_std;
            else
                % Use raw signal values
                processed_signal = run.(signal_field);
            end
            
            % Bin the processed data for this mouse
            for i = 1:length(dist_bins)-1
                indices = run.distance >= dist_bins(i) & run.distance < dist_bins(i+1);
                if any(indices)
                    key = sprintf('sess%d_%s_%s', mouse_sess, run.type, run.arm);
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; processed_signal(indices)];
                end
            end
        end
        
        % Calculate per-mouse, per-bin averages and store them
        for t = 1:length(type_names)
            for a = 1:length(arm_names)
                type = type_names{t};
                arm = arm_names{a};
                key = sprintf('sess%d_%s_%s', mouse_sess, type, arm);
                
                for i = 1:length(dist_bins)-1
                    if ~isempty(mouse_temp_bins.(key){i})
                        mouse_bin_avg = mean(mouse_temp_bins.(key){i});
                        mouse_bin_data.(key){i}(mouse_id) = mouse_bin_avg;
                    end
                end
            end
        end
    end
    
    % Calculate means and SEMs for each condition and plot
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        
        for t = 1:length(type_names)
            type = type_names{t};
            
            % Calculate subplot index
            sp_idx = (s-1)*2 + t;
            subplot(rows, cols, sp_idx);
            hold on;
            
            % Store plot data and mouse data for statistics
            plot_data = struct();
            x_plot_data = struct();
            
            % Plot each arm in the specified order (nonfood first, food last)
            for a = 1:length(arm_plot_order)
                arm = arm_plot_order{a};
                key = sprintf('sess%d_%s_%s', sess, type, arm);
                
                % Get the appropriate color based on type and arm
                color_key = [type '_' arm];
                line_color = arm_colors.(color_key);
                
                % Calculate mean and SEM for each distance bin
                means = nan(length(bin_centers), 1);
                sems = nan(length(bin_centers), 1);
                
                for i = 1:length(dist_bins)-1
                    if isfield(mouse_bin_data, key) && ~isempty(mouse_bin_data.(key){i})
                        bin_values = values(mouse_bin_data.(key){i});
                        if ~isempty(bin_values)
                            bin_data = cell2mat(bin_values);
                            means(i) = mean(bin_data);
                            sems(i) = std(bin_data) / sqrt(length(bin_data));
                        end
                    end
                end
                
                % Apply smoothing to the means and SEMs
                valid = ~isnan(means);
                if sum(valid) > options.smoothing
                    % Get the valid data only
                    x_valid = bin_centers(valid);
                    y_valid = means(valid);
                    sem_valid = sems(valid);
                    
                    % Apply moving average smoothing
                    y_smoothed = movmean(y_valid, options.smoothing);
                    sem_smoothed = movmean(sem_valid, options.smoothing);
                    
                    % Create proper arm labels for the legend
                    arm_label = arm;
                    if strcmp(arm, 'food')
                        arm_label = 'Food';
                    elseif strcmp(arm, 'nonfood1')
                        arm_label = 'Non-food 1';
                    elseif strcmp(arm, 'nonfood2')
                        arm_label = 'Non-food 2';
                    end
                    
                    % For 'towards' runs, negate x values to represent distance to target
                    if strcmp(type, 'towards')
                        x_valid = -x_valid; % Negate x values for plotting
                    end
                    
                    % Store plot data for statistics
                    plot_data.(arm) = y_smoothed;
                    x_plot_data.(arm) = x_valid;
                    
                    % Plot SEM shading using patches - set HandleVisibility to 'off'
                    % to exclude them from the legend
                    for i = 1:length(x_valid)-1
                        x_patch = [x_valid(i), x_valid(i+1), x_valid(i+1), x_valid(i)];
                        y_patch = [y_smoothed(i)-sem_smoothed(i), y_smoothed(i+1)-sem_smoothed(i+1), ...
                                  y_smoothed(i+1)+sem_smoothed(i+1), y_smoothed(i)+sem_smoothed(i)];
                        patch(x_patch, y_patch, line_color, 'EdgeColor', 'none', ...
                              'FaceAlpha', 0.2, 'HandleVisibility', 'off');
                    end
                    
                    % Plot the smoothed line on top with a DisplayName for the legend
                    plot(x_valid, y_smoothed, 'Color', line_color, 'LineWidth', 2, 'DisplayName', arm_label);
                end
            end
            
            % Perform point-by-point statistical tests and add significance bars if requested
            if options.run_stats
                addPointwiseSignificanceBars(sess, type, mouse_bin_data, bin_centers, options, arm_names);
            end

            % Set subplot properties
            title(options.subplot_titles{sp_idx});
            
            % Adjust y-axis label based on whether z-scoring was applied
            if options.apply_zscore
                ylabel(sprintf('Z-scored %s', signal_label));
            else
                ylabel(signal_label);
            end

            % Apply user-specified y-limits
            ylim(options.y_limits);

            % Set x-axis limits based on run type
            if strcmp(type, 'towards')
                xlabel('Distance to Food');
                xlim([-max_distance, 0]);  % Negative values for towards runs
            else
                xlabel('Distance from Food');
                xlim([0, max_distance]);   % Positive values for away runs
            end

            % Add a horizontal dashed line at y=0
            line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

            % Create a clean legend with only the line plots
            legend('show', 'Location', 'best');
            legend('boxoff');
            grid on;
        end
    end
    
    % Add an overall title
    sgtitle(options.title, 'FontWeight', 'bold');
end

function addPointwiseSignificanceBars(sess, type, mouse_bin_data, bin_centers, options, arm_names)
    % Helper function to add point-by-point significance bars to the plot
    % Uses ANOVA + post-hoc tests for proper 3-way comparisons
    
    % Get current y-limits to position significance bars
    yl = ylim;
    y_range = yl(2) - yl(1);
    
    % Position bars at the top of the plot
    bar_y_base = yl(2) - 0.05 * y_range;
    bar_height = 0.01 * y_range;
    
    % Get x-axis limits for bar positioning
    xl = xlim;
    
    % Initialize arrays to store significance results for each bin
    num_bins = length(bin_centers);
    
    % All possible pairwise comparisons
    comparisons = {{'food', 'nonfood1'}, {'food', 'nonfood2'}, {'nonfood1', 'nonfood2'}};
    comparison_colors = {[0.8, 0, 0.8], [0, 0.8, 0.8], [0.8, 0.8, 0]}; % Different colors for different comparisons
    comparison_names = {'F≠N1', 'F≠N2', 'N1≠N2'};
    
    % Store significant bins for each comparison
    all_sig_bins = cell(length(comparisons), 1);
    all_sig_x_coords = cell(length(comparisons), 1);
    
    % Test each bin independently using ANOVA + post-hoc
    for bin_idx = 1:num_bins
        % Collect data for all three arms at this bin
        arm_data = cell(1, 3);
        arm_counts = zeros(1, 3);
        
        for a = 1:length(arm_names)
            arm = arm_names{a};
            key = sprintf('sess%d_%s_%s', sess, type, arm);
            
            if isfield(mouse_bin_data, key) && bin_idx <= length(mouse_bin_data.(key))
                data_map = mouse_bin_data.(key){bin_idx};
                if data_map.Count >= options.min_mice_per_bin
                    arm_data{a} = cell2mat(values(data_map));
                    arm_counts(a) = length(arm_data{a});
                else
                    arm_data{a} = [];
                end
            else
                arm_data{a} = [];
            end
        end
        
        % Check if we have enough data for ANOVA (at least 2 arms with sufficient mice)
        arms_with_data = find(arm_counts >= options.min_mice_per_bin);
        
        if length(arms_with_data) >= 2
            % Prepare data for ANOVA
            all_values = [];
            group_labels = [];
            
            for a = arms_with_data
                all_values = [all_values; arm_data{a}(:)];
                group_labels = [group_labels; a * ones(length(arm_data{a}), 1)];
            end
            
            % Perform one-way ANOVA for this bin
            try
                [p_anova, ~, stats] = anova1(all_values, group_labels, 'off');
                
                % If ANOVA is significant, perform post-hoc tests
                if p_anova < options.stat_alpha && length(arms_with_data) >= 2
                    [comparisons_results, ~, ~, ~] = multcompare(stats, 'Alpha', options.stat_alpha, 'Display', 'off');
                    
                    % Check each pairwise comparison
                    for comp_idx = 1:length(comparisons)
                        arm1 = comparisons{comp_idx}{1};
                        arm2 = comparisons{comp_idx}{2};
                        
                        % Find indices of these arms in our data
                        arm1_idx = find(strcmp(arm_names, arm1));
                        arm2_idx = find(strcmp(arm_names, arm2));
                        
                        % Check if both arms have data and are in the ANOVA
                        arm1_anova_idx = find(arms_with_data == arm1_idx);
                        arm2_anova_idx = find(arms_with_data == arm2_idx);
                        
                        if ~isempty(arm1_anova_idx) && ~isempty(arm2_anova_idx)
                            % Find this comparison in the multcompare results
                            comp_row = find((comparisons_results(:,1) == arm1_anova_idx & comparisons_results(:,2) == arm2_anova_idx) | ...
                                          (comparisons_results(:,1) == arm2_anova_idx & comparisons_results(:,2) == arm1_anova_idx));
                            
                            if ~isempty(comp_row)
                                p_posthoc = comparisons_results(comp_row(1), 6); % p-value column
                                
                                if p_posthoc < options.stat_alpha
                                    % This comparison is significant at this bin
                                    all_sig_bins{comp_idx} = [all_sig_bins{comp_idx}, bin_idx];
                                    
                                    % Calculate x coordinate for this bin
                                    x_coord = bin_centers(bin_idx);
                                    if strcmp(type, 'towards')
                                        x_coord = -x_coord; % Negate for towards runs
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
                % Single significant bin
                run_starts = 1;
                run_ends = 1;
            else
                % Multiple bins - find consecutive runs
                diff_bins = diff(sig_bins);
                break_points = find(diff_bins > 1);
                
                if isempty(break_points)
                    % All bins are consecutive
                    run_starts = 1;
                    run_ends = length(sig_bins);
                else
                    % Multiple runs
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
                    
                    % Add text label at the center of the bar
                    text_x = (x_start + x_end) / 2;
                    text_y = bar_y + bar_height;
                    
                    % Only add text label for the first bar of each comparison to avoid clutter
                    if run_idx == 1
                        text(text_x, text_y, comparison_names{comp_idx}, ...
                             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                             'FontSize', 8, 'FontWeight', 'bold', 'Color', comparison_colors{comp_idx});
                    end
                end
            end
            
            % Print results to console
            arm1_display = getArmDisplayName(comparisons{comp_idx}{1});
            arm2_display = getArmDisplayName(comparisons{comp_idx}{2});
            fprintf('Session %d, %s: %s vs %s significant in %d bins (ANOVA + post-hoc)\n', ...
                sess, type, arm1_display, arm2_display, length(sig_bins));
        end
    end
    
    % Adjust y-limits to accommodate all significance bars
    if max_bar_level > 0
        new_ylim = [yl(1), yl(2) + 0.1 * y_range + max_bar_level * (bar_height + 0.005 * y_range)];
        ylim(new_ylim);
    end
end

function display_name = getArmDisplayName(arm_name)
    % Helper function to convert arm names to display names
    switch arm_name
        case 'food'
            display_name = 'Food';
        case 'nonfood1'
            display_name = 'Non-food 1';
        case 'nonfood2'
            display_name = 'Non-food 2';
        otherwise
            display_name = arm_name;
    end
end