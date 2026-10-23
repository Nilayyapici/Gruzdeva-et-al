function plot_combined_nonfood_dff_across_mice(run_data, options)
    % plot_combined_nonfood_dff_across_mice - Creates plots showing signals vs distance across all mice
    % with non-food arms combined together for clearer food vs non-food comparison
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_mouse_runs
    %   options - Structure with the following fields (all optional):
    %       .sessions - Cell array of sessions to analyze (e.g., {'sess1', 'sess2'})
    %       .smoothing - Window size for smoothing
    %       .y_limits - Limits for y-axis [min max]
    %       .title - Custom title for the figure (optional, auto-generated if not provided)
    %       .subplot_titles - Custom subplot titles (cell array)
    %       .bin_width - Width of distance bins (default: 1)
    %       .signal_type - Signal to analyze ('dff', '465', '405'), default 'dff'
    %       .max_distance - Maximum distance to plot (default: 200)
    %       .apply_zscore - Whether to apply z-scoring (default: true)
    %       .run_stats - Whether to perform statistical tests (default: true)
    %       .stat_alpha - Alpha level for significance (default: 0.05)
    %       .min_mice_per_bin - Minimum mice per bin for stats (default: 3)
    %       .stat_window - Window size for statistical smoothing (default: 5)
    %       .figure_width - Width of the figure in pixels (default: 1200)
    %       .figure_height - Height of the figure in pixels (default: 800)
    %       .show_title - Whether to show the overall figure title (default: true)
    %       .axis_label_font_size - Font size for axis labels (xlabel, ylabel) (default: 12)
    %       .tick_font_size - Font size for axis tick numbers (default: 10)
    %       .subplot_title_font_size - Font size for subplot titles (default: 14)
    %       .main_title_font_size - Font size for main figure title (default: 16)
    %       .x_tick_spacing - Spacing between x-axis ticks in cm (default: 20)
    %       .y_tick_spacing - Spacing between y-axis ticks (default: automatic)
    
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
    
    % Figure size options
    if ~isfield(options, 'figure_width') || isempty(options.figure_width)
        options.figure_width = 1200; % Default figure width
    end
    
    if ~isfield(options, 'figure_height') || isempty(options.figure_height)
        options.figure_height = 800; % Default figure height
    end
    
    % Title display option
    if ~isfield(options, 'show_title') || isempty(options.show_title)
        options.show_title = true; % Default: show title
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
    
    % X-axis and Y-axis tick options
    if ~isfield(options, 'x_tick_spacing') || isempty(options.x_tick_spacing)
        options.x_tick_spacing = 20; % Default x-axis tick spacing in cm
    end
    
    if ~isfield(options, 'y_tick_spacing') || isempty(options.y_tick_spacing)
        options.y_tick_spacing = []; % Default: auto y-axis tick spacing (empty means automatic)
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
            options.title = sprintf('Z-scored %s vs Distance - Food vs Non-Food Combined', signal_label);
        else
            options.title = sprintf('%s vs Distance - Food vs Non-Food Combined', signal_label);
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
        fprintf('Plotting z-scored %s with combined non-food arms for sessions: ', signal_label);
    else
        fprintf('Plotting raw %s with combined non-food arms for sessions: ', signal_label);
    end
    
    for i = 1:length(session_numbers)
        fprintf('sess%d ', session_numbers(i));
    end
    fprintf('\n');
    
    % Create session name mapping for plot titles
    session_names = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        sess_num = session_numbers(i);
        switch sess_num
            case 0
                session_names{i} = 'Before';
            case 1
                session_names{i} = 'Learning';
            case 2
                session_names{i} = 'Test';
            otherwise
                session_names{i} = sprintf('Session %d', sess_num);
        end
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
    figure('Name', options.title, 'Position', [100, 100, options.figure_width, options.figure_height]);
    
    % Define colors for food vs non-food
    arm_colors = struct();
    
    % Blue for 'towards' runs
    arm_colors.towards_food = [0, 0, 0.8];         % Dark blue for food
    arm_colors.towards_nonfood = [0.4, 0.7, 1];    % Light blue for non-food
    
    % Red for 'away' runs
    arm_colors.away_food = [0.8, 0, 0];           % Dark red for food
    arm_colors.away_nonfood = [1, 0.6, 0.6];      % Light red for non-food
    
    % Define arm types (now simplified to food vs non-food)
    arm_types = {'food', 'nonfood'}; % Combined non-food
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
            for a = 1:length(arm_types)
                type = type_names{t};
                arm = arm_types{a};
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
        
        % Process runs for this mouse
        for r = 1:length(runs)
            run = runs(r);
            
            % Skip if run doesn't have the requested signal field
            if ~isfield(run, signal_field)
                continue;
            end
            
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
                    key = sprintf('sess%d_%s_%s', mouse_sess, run.type, combined_arm);
                    mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; processed_signal(indices)];
                end
            end
        end
        
        % Calculate per-mouse, per-bin averages and store them
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
    
    % Calculate means and SEMs for each condition and plot
    for s = 1:length(session_numbers)
        sess = session_numbers(s);
        
        for t = 1:length(type_names)
            type = type_names{t};
            
            % Calculate subplot index
            sp_idx = (s-1)*2 + t;
            subplot(rows, cols, sp_idx);
            hold on;
            
            % Plot each arm type (food vs non-food)
            for a = 1:length(arm_types)
                arm = arm_types{a};
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
                    if strcmp(arm, 'food')
                        arm_label = 'Food';
                    else
                        arm_label = 'Non-Food';
                    end
                    
                    % For 'towards' runs, negate x values to represent distance to target
                    if strcmp(type, 'towards')
                        x_valid = -x_valid; % Negate x values for plotting
                    end
                    
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
                    plot(x_valid, y_smoothed, 'Color', line_color, 'LineWidth', 3, 'DisplayName', arm_label);
                end
            end
            
            % Perform statistical tests and add significance bars if requested
            if options.run_stats
                addFoodVsNonfoodSignificanceBars(sess, type, mouse_bin_data, bin_centers, options);
            end

            % Set subplot properties
            title(options.subplot_titles{sp_idx}, 'FontSize', options.subplot_title_font_size, 'FontWeight', 'bold');
            
            % Adjust y-axis label based on whether z-scoring was applied
            if options.apply_zscore
                ylabel(sprintf('Z-scored %s', signal_label), 'FontSize', options.axis_label_font_size);
            else
                ylabel(signal_label, 'FontSize', options.axis_label_font_size);
            end

            % Apply user-specified y-limits
            ylim(options.y_limits);
            
            % Set y-axis ticks if specified
            if ~isempty(options.y_tick_spacing)
                y_min = options.y_limits(1);
                y_max = options.y_limits(2);
                y_ticks = y_min:options.y_tick_spacing:y_max;
                yticks(y_ticks);
            end

            % Set x-axis limits based on run type
            if strcmp(type, 'towards')
                % Only add xlabel for the last row
                if s == length(session_numbers)
                    xlabel('Distance to Food', 'FontSize', options.axis_label_font_size);
                end
                xlim([-max_distance, 0]);  % Negative values for towards runs
                % Set x-axis ticks every x_tick_spacing cm
                x_ticks = -max_distance:options.x_tick_spacing:0;
                xticks(x_ticks);
            else
                % Only add xlabel for the last row
                if s == length(session_numbers)
                    xlabel('Distance from Food', 'FontSize', options.axis_label_font_size);
                end
                xlim([0, max_distance]);   % Positive values for away runs
                % Set x-axis ticks every x_tick_spacing cm
                x_ticks = 0:options.x_tick_spacing:max_distance;
                xticks(x_ticks);
            end

            % Set axis tick font size (separate from axis labels)
            set(gca, 'FontSize', options.tick_font_size);
            
            % Ensure axis labels keep their specified font size (this overrides the gca FontSize for labels)
            ax = gca;
            ax.XLabel.FontSize = options.axis_label_font_size;
            ax.YLabel.FontSize = options.axis_label_font_size;

            % Add a horizontal dashed line at y=0
            line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

            % Create a clean legend with only the line plots
            legend('show', 'Location', 'best');
            legend('boxoff');
            
            % Make legend lines shorter
            leg = legend;
            leg.ItemTokenSize = [15, 8]; % [width, height] - default is usually [30, 18]
            grid off;
        end
    end
    
    % Add an overall title if requested
    if options.show_title
        sgtitle(options.title, 'FontWeight', 'bold', 'FontSize', options.main_title_font_size);
    end
end

function addFoodVsNonfoodSignificanceBars(sess, type, mouse_bin_data, bin_centers, options)
    % Helper function to add significance bars comparing food vs non-food arms
    
    % Get current y-limits to position significance bars
    yl = ylim;
    y_range = yl(2) - yl(1);
    
    % Position bars at the top of the plot
    bar_y = yl(2) - 0.05 * y_range;
    bar_height = 0.01 * y_range;
    
    % Initialize arrays to store significance results for each bin
    num_bins = length(bin_centers);
    sig_bins = [];
    sig_x_coords = [];
    
    % Test each bin independently using t-test
    for bin_idx = 1:num_bins
        % Get data for food and non-food arms at this bin
        food_key = sprintf('sess%d_%s_food', sess, type);
        nonfood_key = sprintf('sess%d_%s_nonfood', sess, type);
        
        food_data = [];
        nonfood_data = [];
        
        % Extract food arm data
        if isfield(mouse_bin_data, food_key) && bin_idx <= length(mouse_bin_data.(food_key))
            food_map = mouse_bin_data.(food_key){bin_idx};
            if food_map.Count >= options.min_mice_per_bin
                food_data = cell2mat(values(food_map));
            end
        end
        
        % Extract non-food arm data
        if isfield(mouse_bin_data, nonfood_key) && bin_idx <= length(mouse_bin_data.(nonfood_key))
            nonfood_map = mouse_bin_data.(nonfood_key){bin_idx};
            if nonfood_map.Count >= options.min_mice_per_bin
                nonfood_data = cell2mat(values(nonfood_map));
            end
        end
        
        % Perform t-test if we have enough data for both groups
        if ~isempty(food_data) && ~isempty(nonfood_data) && ...
           length(food_data) >= options.min_mice_per_bin && length(nonfood_data) >= options.min_mice_per_bin
            
            try
                [~, p_ttest] = ttest2(food_data, nonfood_data);
                
                if p_ttest < options.stat_alpha
                    % This bin shows significant difference
                    sig_bins = [sig_bins, bin_idx];
                    
                    % Calculate x coordinate for this bin
                    x_coord = bin_centers(bin_idx);
                    if strcmp(type, 'towards')
                        x_coord = -x_coord; % Negate for towards runs
                    end
                    sig_x_coords = [sig_x_coords, x_coord];
                end
            catch ME
                % Skip if statistical test fails
                continue;
            end
        end
    end
    
    % Draw significance bars for consecutive significant bins
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
                
                % Draw the significance bar
                plot([x_start, x_end], [bar_y, bar_y], 'Color', [0.2, 0.2, 0.2], ...
                     'LineWidth', 4, 'HandleVisibility', 'off');
                
                % Add asterisks to indicate significance level
                % text_x = (x_start + x_end) / 2;
                % text_y = bar_y + bar_height;
                % text(text_x, text_y, '*', 'HorizontalAlignment', 'center', ...
                %      'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
            end
        end
        
        % Adjust y-limits to accommodate significance bars
        new_ylim = [yl(1), yl(2) + 0.08 * y_range];
        ylim(new_ylim);
        
        % Print results to console
        fprintf('Session %d, %s: Food vs Non-Food significant in %d bins (t-test)\n', ...
            sess, type, length(sig_bins));
    end
end