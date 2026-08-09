function plot_correlation_comparison(run_data, options)
    % PLOT_CORRELATION_COMPARISON Creates a plot comparing towards vs away correlations
    %
    % This function plots the correlation between distance and dF/F for a single 
    % condition, comparing towards and away runs for selected sessions with
    % statistical comparisons between sessions.
    %
    % Inputs:
    %   run_data - Run data structure from analyze_single_arm_runs
    %   options - Structure with configuration parameters:
    %       .sessions - Cell array of session names (e.g. {'sess0', 'sess1'})
    %       .run_category - Which runs to include: 'food', 'not_food', or 'all' (default: 'food')
    %       .title - Optional custom title for the plot
    %       .separate_plots - Boolean, if true creates separate plots for towards/away (default: false)
    %       .use_fisher_z - Boolean, if true uses Fisher z-transformation for stats (default: true)
    %
    % Example Usage:
    %   options = struct();
    %   options.sessions = {'sess0', 'sess1'};
    %   options.run_category = 'food';  % or 'not_food' or 'all'
    %   options.title = 'Fasted Food: Towards vs Away Correlations';
    %   options.separate_plots = true;  % Create separate plots
    %   options.use_fisher_z = false;   % Use raw correlations for stats
    %   plot_correlation_comparison(runs_fasted_food, options);
    
    % Set default options if not provided
    if nargin < 2 || isempty(options)
        options = struct();
    end
    
    % Set default run_category
    if ~isfield(options, 'run_category')
        options.run_category = 'food';  % Default to food runs only
    end
    
    % Validate run_category
    valid_categories = {'food', 'not_food', 'all'};
    if ~ismember(options.run_category, valid_categories)
        error('run_category must be one of: ''food'', ''not_food'', or ''all''');
    end
    
    % Set default sessions if not provided
    if ~isfield(options, 'sessions') || isempty(options.sessions)
        options.sessions = {'sess0', 'sess1'};
        warning('No sessions specified. Using default: sess0 and sess1');
    end
    
    % Set default for separate_plots
    if ~isfield(options, 'separate_plots')
        options.separate_plots = false;
    end
    
    % Set default for Fisher z-transformation
    if ~isfield(options, 'use_fisher_z')
        options.use_fisher_z = true;
    end
    
    % Convert session names to numbers for processing
    session_nums = [];
    for i = 1:length(options.sessions)
        session_str = char(options.sessions{i});
        if contains(session_str, 'sess')
            session_num = str2double(extractAfter(session_str, 'sess'));
            if ~isnan(session_num)
                session_nums = [session_nums, session_num];
            end
        end
    end
    
    if isempty(session_nums)
        error('No valid session numbers found in options.sessions');
    end
    
    % Display which sessions and run category are being analyzed
    fprintf('Analyzing sessions: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf('\nRun category: %s\n\n', options.run_category);
    
    % 1. Calculate mouse-level averages
    avg_corr = calculate_mouse_averages(run_data, session_nums, options.run_category);
    
    % 2. Plot the comparison based on separate_plots option
    if options.separate_plots
        plot_separate_directions(avg_corr, options, session_nums);
    else
        plot_towards_vs_away(avg_corr, options, session_nums);
    end
end

function avg_corr = calculate_mouse_averages(runs, session_nums, run_category)
    % Calculate average correlations for each mouse by session and run type
    
    % Initialize output structure
    avg_corr = struct('mouse_id', {}, 'session', {}, 'run_type', {}, ...
                     'avg_r', {}, 'n_runs', {}, 'z_avg', {});
    
    % Exit if runs is empty
    if isempty(runs)
        return;
    end
    
    % Determine which run types to include based on run_category
    switch run_category
        case 'food'
            valid_types = {'towards', 'away'};
            direction_labels = {'towards', 'away'};
        case 'not_food'
            valid_types = {'not_food_towards', 'not_food_away'};
            direction_labels = {'towards', 'away'};  % Use simplified labels for plotting
        case 'all'
            valid_types = {'towards', 'away', 'not_food_towards', 'not_food_away'};
            direction_labels = {'towards', 'away'};  % Combine into two groups
    end
    
    % Process each mouse
    for i = 1:length(runs)
        mouse_id = runs(i).mouse_id;
        session = runs(i).session;
        
        % Skip if not in selected sessions
        if ~ismember(session, session_nums)
            continue;
        end
        
        % Get runs for this mouse
        mouse_runs = runs(i).runs;
        
        % Skip if no runs
        if isempty(mouse_runs)
            continue;
        end
        
        % Process by direction (towards/away) - combining food and not_food if 'all'
        for dir_idx = 1:length(direction_labels)
            direction = direction_labels{dir_idx};
            
            % Determine which specific run types to include for this direction
            if strcmp(run_category, 'food')
                types_to_include = {direction};
            elseif strcmp(run_category, 'not_food')
                types_to_include = {['not_food_' direction]};
            else  % 'all'
                types_to_include = {direction, ['not_food_' direction]};
            end
            
            % Find runs matching these types
            type_indices = [];
            if isstruct(mouse_runs) && isfield(mouse_runs, 'type')
                for type_idx = 1:length(types_to_include)
                    type_indices = [type_indices, find(strcmp({mouse_runs.type}, types_to_include{type_idx}))];
                end
            elseif iscell(mouse_runs)
                % Handle case where runs might be stored differently
                for k = 1:length(mouse_runs)
                    if isfield(mouse_runs{k}, 'type')
                        for type_idx = 1:length(types_to_include)
                            if strcmp(mouse_runs{k}.type, types_to_include{type_idx})
                                type_indices = [type_indices, k];
                                break;
                            end
                        end
                    end
                end
            else
                % If structure is different, try to find type field
                fprintf('Warning: Unexpected runs structure for mouse %s, session %d\n', mouse_id, session);
                continue;
            end
            
            % Skip if no runs of this type
            if isempty(type_indices)
                continue;
            end
            
            % Calculate correlation for each run
            r_values = [];
            z_values = [];
            
            for k = 1:length(type_indices)
                run_idx = type_indices(k);
                if iscell(mouse_runs)
                    run = mouse_runs{run_idx};
                else
                    run = mouse_runs(run_idx);
                end
                
                % Check if run has required fields
                if ~isfield(run, 'distance') || ~isfield(run, 'dff')
                    fprintf('Warning: Run missing distance or dff field for mouse %s\n', mouse_id);
                    continue;
                end
                
                % Calculate correlation between distance and dF/F
                [r, p] = corr(run.distance, run.dff, 'type', 'Pearson', 'rows', 'complete');
                
                % Skip if correlation couldn't be calculated
                if isnan(r) || isnan(p)
                    continue;
                end
                
                % Store correlation value
                r_values = [r_values, r];
                
                % Calculate Fisher's z-transform
                if r < 1 && r > -1  % Avoid issues with perfect correlations
                    z = 0.5 * log((1 + r) / (1 - r));
                    z_values = [z_values, z];
                end
            end
            
            % Skip if no valid correlations
            if isempty(r_values)
                continue;
            end
            
            % Calculate average correlation
            avg_r = mean(r_values);
            
            % Calculate average z-transform (for statistical comparisons)
            avg_z = mean(z_values);
            
            % Add to output structure (use direction label for consistency)
            entry = struct('mouse_id', mouse_id, ...
                          'session', session, ...
                          'run_type', direction, ...
                          'avg_r', avg_r, ...
                          'n_runs', length(r_values), ...
                          'z_avg', avg_z);
            
            avg_corr = [avg_corr; entry];
        end
    end
end

function plot_separate_directions(avg_corr, options, session_nums)
    % Plot towards and away runs in separate figures
    
    % Exit if no data
    if isempty(avg_corr)
        fprintf('No data available for the selected sessions and run category.\n');
        return;
    end
    
    run_types = {'towards', 'away'};
    run_type_labels = {'Towards', 'Away'};
    
    % Add run category label
    switch options.run_category
        case 'food'
            category_label = '(Food Runs)';
        case 'not_food'
            category_label = '(Non-Food Runs)';
        case 'all'
            category_label = '(All Runs)';
    end
    
    % Define colors for each run type
    % Towards: Blue colors (light blue for sess0, dark blue for sess1, etc.)
    towards_colors = [0.6, 0.8, 1.0; 0.2, 0.4, 0.8; 0.4, 0.6, 0.9; 0.1, 0.3, 0.7];
    % Away: Red colors (light red for sess0, dark red for sess1, etc.)
    away_colors = [1.0, 0.6, 0.6; 0.8, 0.2, 0.2; 0.9, 0.3, 0.3; 0.7, 0.1, 0.1];
    
    % Create separate figure for each run type
    for type_idx = 1:length(run_types)
        run_type = run_types{type_idx};
        run_type_label = run_type_labels{type_idx};
        
        % Select colors for this run type
        if strcmp(run_type, 'towards')
            colors = towards_colors;
        else
            colors = away_colors;
        end
        
        % Create figure
        figure('Position', [100 + (type_idx-1)*600, 100, 500, 600]);
        
        % Extract data for this run type across sessions
        type_data = [];
        
        for i = 1:length(session_nums)
            sess = session_nums(i);
            
            % Get data for this session and run type
            type_idx_data = [avg_corr.session] == sess & strcmp({avg_corr.run_type}, run_type);
            
            if sum(type_idx_data) > 0
                session_mean = mean([avg_corr(type_idx_data).avg_r]);
                session_sem = std([avg_corr(type_idx_data).avg_r]) / sqrt(sum(type_idx_data));
                session_count = sum(type_idx_data);
                session_values = [avg_corr(type_idx_data).avg_r];
                session_mice = {avg_corr(type_idx_data).mouse_id};
            else
                session_mean = NaN;
                session_sem = NaN;
                session_count = 0;
                session_values = [];
                session_mice = {};
            end
            
            % Store data for this session
            type_data(i).mean = session_mean;
            type_data(i).sem = session_sem;
            type_data(i).values = session_values;
            type_data(i).mice = session_mice;
            type_data(i).count = session_count;
            type_data(i).session = sess;
        end
        
        % Plot bars for each session
        bar_width = 0.7;
        x_positions = 1:length(session_nums);
        
        hold on;
        
        % Plot bars and data points
        for i = 1:length(session_nums)
            if type_data(i).count > 0
                % Plot bar
                bar(x_positions(i), type_data(i).mean, bar_width, 'FaceColor', colors(i,:));
                
                % Add error bar
                errorbar(x_positions(i), type_data(i).mean, type_data(i).sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
                
                % Add individual data points
                if ~isempty(type_data(i).values)
                    x_jitter = repmat(x_positions(i), 1, length(type_data(i).values));
                    scatter(x_jitter, type_data(i).values, 25, 'k', 'filled');
                end
            end
        end
        
        % Add connecting lines for individual mice if we have 2 sessions
        if length(session_nums) == 2 && type_data(1).count > 0 && type_data(2).count > 0
            common_mice = intersect(type_data(1).mice, type_data(2).mice);
            
            for i = 1:length(common_mice)
                mouse_id = common_mice{i};
                
                % Find values for this mouse in both sessions
                idx1 = find(strcmp(type_data(1).mice, mouse_id));
                idx2 = find(strcmp(type_data(2).mice, mouse_id));
                
                if ~isempty(idx1) && ~isempty(idx2)
                    x1 = x_positions(1);
                    x2 = x_positions(2);
                    y1 = type_data(1).values(idx1);
                    y2 = type_data(2).values(idx2);
                    
                    plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.7]);
                end
            end
        end
        
        % Add statistical comparison between sessions
        add_session_comparison_stats(type_data, session_nums, x_positions, options);
        
        % Set labels and title
        if isfield(options, 'title') && ~isempty(options.title)
            title([options.title, ' - ', run_type_label, ' Runs ', category_label], 'FontWeight', 'bold');
        else
            title([run_type_label, ' Run Correlations ', category_label], 'FontWeight', 'bold');
        end
        
        xlabel('Session', 'FontSize', 14);
        ylabel('Pearson Correlation (r)', 'FontSize', 14);
        
        % Set x-axis
        session_labels = cell(1, length(session_nums));
        for i = 1:length(session_nums)
            session_labels{i} = ['sess', num2str(session_nums(i))];
        end
        set(gca, 'XTick', x_positions, 'XTickLabel', session_labels, 'FontSize', 14);
        
        % Set y-axis limits
        all_data = [];
        for i = 1:length(type_data)
            if ~isempty(type_data(i).values)
                all_data = [all_data, type_data(i).values];
            end
        end
        
        if ~isempty(all_data)
            y_max = max(max(all_data) + 0.3, 0.4);
            y_min = min(min(all_data) - 0.2, -0.3);
            ylim([-0.8, 1.1]);
        end
        
        % Add grid
        grid off;
        box off;
        hold off;
        
        % Print statistics for this run type
        print_single_direction_stats(type_data, session_nums, run_type_label, options);
    end
end

function add_session_comparison_stats(type_data, session_nums, x_positions, options)
    % Add statistical comparison between sessions for a single run type
    
    if length(session_nums) ~= 2 || type_data(1).count == 0 || type_data(2).count == 0
        return;
    end
    
    % Find common mice
    common_mice = intersect(type_data(1).mice, type_data(2).mice);
    
    if length(common_mice) >= 2
        % Paired t-test
        sess0_values = [];
        sess1_values = [];
        
        for i = 1:length(common_mice)
            mouse_id = common_mice{i};
            idx1 = find(strcmp(type_data(1).mice, mouse_id));
            idx2 = find(strcmp(type_data(2).mice, mouse_id));
            
            if ~isempty(idx1) && ~isempty(idx2)
                sess0_values = [sess0_values, type_data(1).values(idx1)];
                sess1_values = [sess1_values, type_data(2).values(idx2)];
            end
        end
        
        % Apply Fisher z-transformation if enabled
        if options.use_fisher_z
            z0 = 0.5 * log((1 + sess0_values) ./ (1 - sess0_values));
            z1 = 0.5 * log((1 + sess1_values) ./ (1 - sess1_values));
            [~, p_value] = ttest(z0, z1);
        else
            [~, p_value] = ttest(sess0_values, sess1_values);
        end
    else
        % Unpaired t-test
        if options.use_fisher_z
            z0 = 0.5 * log((1 + type_data(1).values) ./ (1 - type_data(1).values));
            z1 = 0.5 * log((1 + type_data(2).values) ./ (1 - type_data(2).values));
            [~, p_value] = ttest2(z0, z1);
        else
            [~, p_value] = ttest2(type_data(1).values, type_data(2).values);
        end
    end
    
    % Add significance indicator if significant
    if p_value < 0.05
        x1 = x_positions(1);
        x2 = x_positions(2);
        y_pos = max(type_data(1).mean + type_data(1).sem, type_data(2).mean + type_data(2).sem) + 0.25;
        
        plot([x1, x2], [y_pos, y_pos], 'k-', 'LineWidth', 1.5);
        
        if p_value < 0.001
            sig_str = '***';
        elseif p_value < 0.01
            sig_str = '**';
        else
            sig_str = '*';
        end
        
        text((x1 + x2)/2, y_pos + 0.06, sig_str, 'HorizontalAlignment', 'center', 'FontSize', 14);
        text((x1 + x2)/2, y_pos + 0.14, sprintf('p = %.3f', p_value), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end

function print_single_direction_stats(type_data, session_nums, run_type_label, options)
    % Print statistical summary for a single run type
    
    fprintf('\n=== %s RUNS STATISTICAL SUMMARY ===\n', upper(run_type_label));
    
    if length(session_nums) == 2 && type_data(1).count > 0 && type_data(2).count > 0
        common_mice = intersect(type_data(1).mice, type_data(2).mice);
        
        fprintf('Sessions: sess%d vs sess%d\n', session_nums(1), session_nums(2));
        fprintf('Sample sizes: n=%d vs n=%d\n', type_data(1).count, type_data(2).count);
        fprintf('Means ± SEM: %.3f ± %.3f vs %.3f ± %.3f\n', ...
            type_data(1).mean, type_data(1).sem, type_data(2).mean, type_data(2).sem);
        
        if length(common_mice) >= 2
            % Paired test
            sess0_values = [];
            sess1_values = [];
            
            for i = 1:length(common_mice)
                mouse_id = common_mice{i};
                idx1 = find(strcmp(type_data(1).mice, mouse_id));
                idx2 = find(strcmp(type_data(2).mice, mouse_id));
                
                if ~isempty(idx1) && ~isempty(idx2)
                    sess0_values = [sess0_values, type_data(1).values(idx1)];
                    sess1_values = [sess1_values, type_data(2).values(idx2)];
                end
            end
            
            % Apply Fisher z-transformation if enabled
            if options.use_fisher_z
                z0 = 0.5 * log((1 + sess0_values) ./ (1 - sess0_values));
                z1 = 0.5 * log((1 + sess1_values) ./ (1 - sess1_values));
                [~, p_value, ~, stats] = ttest(z0, z1);
                test_method = 'paired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest(sess0_values, sess1_values);
                test_method = 'paired t-test on raw correlation values';
            end
            
            fprintf('Common mice: n=%d (%s)\n', length(common_mice), test_method);
            fprintf('t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        else
            % Unpaired test
            if options.use_fisher_z
                z0 = 0.5 * log((1 + type_data(1).values) ./ (1 - type_data(1).values));
                z1 = 0.5 * log((1 + type_data(2).values) ./ (1 - type_data(2).values));
                [~, p_value, ~, stats] = ttest2(z0, z1);
                test_method = 'unpaired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest2(type_data(1).values, type_data(2).values);
                test_method = 'unpaired t-test on raw correlation values';
            end
            
            fprintf('Independent samples (%s)\n', test_method);
            fprintf('t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        end
        
        if p_value < 0.001
            fprintf('Result: *** (p < 0.001)\n');
        elseif p_value < 0.01
            fprintf('Result: ** (p < 0.01)\n');
        elseif p_value < 0.05
            fprintf('Result: * (p < 0.05)\n');
        else
            fprintf('Result: n.s. (p ≥ 0.05)\n');
        end
    else
        fprintf('Insufficient data for comparison\n');
    end
    
    fprintf('=== END %s RUNS SUMMARY ===\n\n', upper(run_type_label));
end

function plot_towards_vs_away(avg_corr, options, session_nums)
    % Plot comparison of towards vs. away runs with grouped bars (original function)
    
    % Exit if no data
    if isempty(avg_corr)
        fprintf('No data available for the selected sessions and run category.\n');
        return;
    end
    
    % Add run category label
    switch options.run_category
        case 'food'
            category_label = '(Food Runs)';
        case 'not_food'
            category_label = '(Non-Food Runs)';
        case 'all'
            category_label = '(All Runs)';
    end
    
    % Create figure
    figure('Position', [100, 100, 500, 600]);
    
    % Use sessions from options, convert to numbers if needed
    if nargin < 3
        if isfield(options, 'sessions') && ~isempty(options.sessions)
            sessions = options.sessions;
            % Convert session names to numbers
            session_nums = [];
            for i = 1:length(sessions)
                session_str = char(sessions{i});
                if contains(session_str, 'sess')
                    session_num = str2double(extractAfter(session_str, 'sess'));
                    if ~isnan(session_num)
                        session_nums = [session_nums, session_num];
                    end
                end
            end
        end
    end
    
    % Sort sessions and ensure we only use sess0 and sess1
    session_nums = sort(session_nums);
    if length(session_nums) > 2
        warning('Only sess0 and sess1 will be plotted');
        session_nums = session_nums(1:2);
    end
    
    % Define colors
    % Towards: Blue (light blue for sess0, dark blue for sess1)
    towards_colors = [0.6, 0.8, 1.0; 0.2, 0.4, 0.8];  % Light blue, Dark blue
    % Away: Red (light red for sess0, dark red for sess1)
    away_colors = [1.0, 0.6, 0.6; 0.8, 0.2, 0.2];     % Light red, Dark red
    
    % Group positions
    group_width = 0.8;
    bar_width = 0.4;  % Reduced bar width
    bar_gap = 0.45;     % Gap between sess0 and sess1 bars within each group
    group_gap = 1.5;
    
    % X positions for groups
    x_towards_group = 1;
    x_away_group = 1 + group_gap;
    
    % Store data for connecting lines
    towards_data = [];
    away_data = [];
    
    % Process each session for towards runs
    for i = 1:length(session_nums)
        sess = session_nums(i);
        
        % Get towards data for this session
        towards_idx = [avg_corr.session] == sess & strcmp({avg_corr.run_type}, 'towards');
        
        if sum(towards_idx) > 0
            towards_mean = mean([avg_corr(towards_idx).avg_r]);
            towards_sem = std([avg_corr(towards_idx).avg_r]) / sqrt(sum(towards_idx));
            towards_count = sum(towards_idx);
            towards_values = [avg_corr(towards_idx).avg_r];
            towards_mice = {avg_corr(towards_idx).mouse_id};
        else
            towards_mean = NaN;
            towards_sem = NaN;
            towards_count = 0;
            towards_values = [];
            towards_mice = {};
        end
        
        % Store data for connecting lines
        towards_data(i).mean = towards_mean;
        towards_data(i).sem = towards_sem;
        towards_data(i).values = towards_values;
        towards_data(i).mice = towards_mice;
        towards_data(i).count = towards_count;
        
        % X position for this bar
        x_pos = x_towards_group + (i-1) * bar_gap - bar_gap/2;
        
        % Plot bar
        hold on;
        bar(x_pos, towards_mean, bar_width, 'FaceColor', towards_colors(i,:));
        
        % Add error bar
        errorbar(x_pos, towards_mean, towards_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        
        % Add individual data points (store positions for connecting lines)
        if ~isempty(towards_values)
            x_jitter = repmat(x_pos, 1, length(towards_values));  % No jittering for connecting lines
            scatter(x_jitter, towards_values, 15, 'k', 'filled');
        end
    end
    
    % Process each session for away runs
    for i = 1:length(session_nums)
        sess = session_nums(i);
        
        % Get away data for this session
        away_idx = [avg_corr.session] == sess & strcmp({avg_corr.run_type}, 'away');
        
        if sum(away_idx) > 0
            away_mean = mean([avg_corr(away_idx).avg_r]);
            away_sem = std([avg_corr(away_idx).avg_r]) / sqrt(sum(away_idx));
            away_count = sum(away_idx);
            away_values = [avg_corr(away_idx).avg_r];
            away_mice = {avg_corr(away_idx).mouse_id};
        else
            away_mean = NaN;
            away_sem = NaN;
            away_count = 0;
            away_values = [];
            away_mice = {};
        end
        
        % Store data for connecting lines
        away_data(i).mean = away_mean;
        away_data(i).sem = away_sem;
        away_data(i).values = away_values;
        away_data(i).mice = away_mice;
        away_data(i).count = away_count;
        
        % X position for this bar
        x_pos = x_away_group + (i-1) * bar_gap - bar_gap/2;
        
        % Plot bar
        bar(x_pos, away_mean, bar_width, 'FaceColor', away_colors(i,:));
        
        % Add error bar
        errorbar(x_pos, away_mean, away_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        
        % Add individual data points (store positions for connecting lines)
        if ~isempty(away_values)
            x_jitter = repmat(x_pos, 1, length(away_values));  % No jittering for connecting lines
            scatter(x_jitter, away_values, 15, 'k', 'filled');
        end
    end
    
    % Add connecting lines for individual mice (towards group)
    if length(session_nums) == 2 && ~isempty(towards_data(1).mice) && ~isempty(towards_data(2).mice)
        common_mice = intersect(towards_data(1).mice, towards_data(2).mice);
        
        for i = 1:length(common_mice)
            mouse_id = common_mice{i};
            
            % Find values for this mouse in both sessions
            idx1 = find(strcmp(towards_data(1).mice, mouse_id));
            idx2 = find(strcmp(towards_data(2).mice, mouse_id));
            
            if ~isempty(idx1) && ~isempty(idx2)
                x1 = x_towards_group - bar_gap/2;
                x2 = x_towards_group + bar_gap/2;
                y1 = towards_data(1).values(idx1);
                y2 = towards_data(2).values(idx2);
                
                plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.7]);
            end
        end
    end
    
    % Add connecting lines for individual mice (away group)
    if length(session_nums) == 2 && ~isempty(away_data(1).mice) && ~isempty(away_data(2).mice)
        common_mice = intersect(away_data(1).mice, away_data(2).mice);
        
        for i = 1:length(common_mice)
            mouse_id = common_mice{i};
            
            % Find values for this mouse in both sessions
            idx1 = find(strcmp(away_data(1).mice, mouse_id));
            idx2 = find(strcmp(away_data(2).mice, mouse_id));
            
            if ~isempty(idx1) && ~isempty(idx2)
                x1 = x_away_group - bar_gap/2;
                x2 = x_away_group + bar_gap/2;
                y1 = away_data(1).values(idx1);
                y2 = away_data(2).values(idx2);
                
                plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.7]);
            end
        end
    end
    
    % Add statistical comparisons within groups (sess0 vs sess1)
    add_within_group_stats(towards_data, session_nums, x_towards_group, bar_gap, 'towards', options);
    add_within_group_stats(away_data, session_nums, x_away_group, bar_gap, 'away', options);
    
    % Add statistical comparison between groups (towards vs away) for each session
    add_between_group_stats(towards_data, away_data, session_nums, x_towards_group, x_away_group, bar_gap, options);
    
    % Print statistical summary to command window
    print_statistical_summary(towards_data, away_data, session_nums, options);
    
    % Set labels and title
    if isfield(options, 'title') && ~isempty(options.title)
        title([options.title, ' ', category_label], 'FontWeight', 'bold');
    else
        title(['Towards vs Away Run Correlations ', category_label], 'FontWeight', 'bold');
    end
    
    xlabel('Run Direction','FontSize', 14);
    ylabel('Pearson Correlation (r)', 'FontSize', 14);
    
    % Set x-axis
    set(gca, 'XTick', [x_towards_group, x_away_group], 'XTickLabel', {'Towards', 'Away'}, 'FontSize', 14);
    
    % Create custom legend
    legend_elements = [];
    legend_labels = {};
    
    % Add legend elements for each session
    for i = 1:length(session_nums)
        % Towards
        h1 = bar(NaN, NaN, 'FaceColor', towards_colors(i,:));
        legend_elements = [legend_elements, h1];
        legend_labels{end+1} = ['Towards sess', num2str(session_nums(i))];
        
        % Away
        h2 = bar(NaN, NaN, 'FaceColor', away_colors(i,:));
        legend_elements = [legend_elements, h2];
        legend_labels{end+1} = ['Away sess', num2str(session_nums(i))];
    end
    
    legend(legend_elements, legend_labels, 'Location', 'best', 'Box', 'off');
    
    % Set y-axis limits
    all_data = [];
    for i = 1:length(towards_data)
        if ~isempty(towards_data(i).values)
            all_data = [all_data, towards_data(i).values];
        end
    end
    for i = 1:length(away_data)
        if ~isempty(away_data(i).values)
            all_data = [all_data, away_data(i).values];
        end
    end
    
    if ~isempty(all_data)
        y_max = max(max(all_data) + 0.3, 0.4);
        y_min = min(min(all_data) - 0.2, -0.3);
        ylim([y_min, y_max]);
    end
    
    % Add grid
    grid off;
    box off;
    hold off;
end

function add_within_group_stats(group_data, session_nums, x_center, bar_gap, group_name, options)
    % Add statistical comparison within a group (sess0 vs sess1)
    
    if length(session_nums) ~= 2 || length(group_data) ~= 2
        return;
    end
    
    % Check if both sessions have data
    if group_data(1).count == 0 || group_data(2).count == 0
        return;
    end
    
    % Find common mice
    common_mice = intersect(group_data(1).mice, group_data(2).mice);
    
    if length(common_mice) >= 2
        % Paired t-test
        sess0_values = [];
        sess1_values = [];
        
        for i = 1:length(common_mice)
            mouse_id = common_mice{i};
            idx1 = find(strcmp(group_data(1).mice, mouse_id));
            idx2 = find(strcmp(group_data(2).mice, mouse_id));
            
            if ~isempty(idx1) && ~isempty(idx2)
                sess0_values = [sess0_values, group_data(1).values(idx1)];
                sess1_values = [sess1_values, group_data(2).values(idx2)];
            end
        end
        
        % Apply Fisher z-transformation if enabled
        if options.use_fisher_z
            z0 = 0.5 * log((1 + sess0_values) ./ (1 - sess0_values));
            z1 = 0.5 * log((1 + sess1_values) ./ (1 - sess1_values));
            [~, p_value] = ttest(z0, z1);
        else
            [~, p_value] = ttest(sess0_values, sess1_values);
        end
        test_type = 'paired';
    else
        % Unpaired t-test
        if options.use_fisher_z
            z0 = 0.5 * log((1 + group_data(1).values) ./ (1 - group_data(1).values));
            z1 = 0.5 * log((1 + group_data(2).values) ./ (1 - group_data(2).values));
            [~, p_value] = ttest2(z0, z1);
        else
            [~, p_value] = ttest2(group_data(1).values, group_data(2).values);
        end
        test_type = 'unpaired';
    end
    
    % Add significance indicator if significant
    if p_value < 0.05
        x1 = x_center - bar_gap/2;
        x2 = x_center + bar_gap/2;
        y_pos = max(group_data(1).mean + group_data(1).sem, group_data(2).mean + group_data(2).sem) + 0.30;  % Moved higher
        
        plot([x1, x2], [y_pos, y_pos], 'k-', 'LineWidth', 1.5);
        
        if p_value < 0.001
            sig_str = '***';
        elseif p_value < 0.01
            sig_str = '**';
        else
            sig_str = '*';
        end
        
        text(x_center, y_pos + 0.06, sig_str, 'HorizontalAlignment', 'center', 'FontSize', 14);  % Bigger font
        text(x_center, y_pos + 0.14, sprintf('p = %.3f', p_value), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);  % Bigger font
    end
end

function add_between_group_stats(towards_data, away_data, session_nums, x_towards, x_away, bar_gap, options)
    % Add statistical comparison between groups (towards vs away) for each session
    
    for i = 1:length(session_nums)
        if towards_data(i).count == 0 || away_data(i).count == 0
            continue;
        end
        
        % Find common mice
        common_mice = intersect(towards_data(i).mice, away_data(i).mice);
        
        if length(common_mice) >= 2
            % Paired t-test
            towards_values = [];
            away_values = [];
            
            for j = 1:length(common_mice)
                mouse_id = common_mice{j};
                idx_towards = find(strcmp(towards_data(i).mice, mouse_id));
                idx_away = find(strcmp(away_data(i).mice, mouse_id));
                
                if ~isempty(idx_towards) && ~isempty(idx_away)
                    towards_values = [towards_values, towards_data(i).values(idx_towards)];
                    away_values = [away_values, away_data(i).values(idx_away)];
                end
            end
            
            % Apply Fisher z-transformation if enabled
            if options.use_fisher_z
                z_towards = 0.5 * log((1 + towards_values) ./ (1 - towards_values));
                z_away = 0.5 * log((1 + away_values) ./ (1 - away_values));
                [~, p_value] = ttest(z_towards, z_away);
            else
                [~, p_value] = ttest(towards_values, away_values);
            end
            test_type = 'paired';
        else
            % Unpaired t-test
            if options.use_fisher_z
                z_towards = 0.5 * log((1 + towards_data(i).values) ./ (1 - towards_data(i).values));
                z_away = 0.5 * log((1 + away_data(i).values) ./ (1 - away_data(i).values));
                [~, p_value] = ttest2(z_towards, z_away);
            else
                [~, p_value] = ttest2(towards_data(i).values, away_data(i).values);
            end
            test_type = 'unpaired';
        end
        
        % Add significance indicator if significant
        if p_value < 0.05
            x1 = x_towards + (i-1) * bar_gap - bar_gap/2;
            x2 = x_away + (i-1) * bar_gap - bar_gap/2;
            y_pos = max(towards_data(i).mean + towards_data(i).sem, away_data(i).mean + away_data(i).sem) + 0.30 + (i-1)*0.12;  % Moved higher
            
            plot([x1, x2], [y_pos, y_pos], 'k-', 'LineWidth', 1.0);
            
            if p_value < 0.001
                sig_str = '***';
            elseif p_value < 0.01
                sig_str = '**';
            else
                sig_str = '*';
            end
            
            text((x1 + x2)/2, y_pos + 0.06, sig_str, 'HorizontalAlignment', 'center', 'FontSize', 12);  % Bigger font
            text((x1 + x2)/2, y_pos + 0.14, sprintf('sess%d: p = %.3f', session_nums(i), p_value), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);  % Bigger font
        end
    end
end

function print_statistical_summary(towards_data, away_data, session_nums, options)
    % Print detailed statistical summary to command window
    
    fprintf('\n=== STATISTICAL SUMMARY (%s) ===\n', upper(options.run_category));
    if options.use_fisher_z
        fprintf('(Using Fisher z-transformation for statistical tests)\n\n');
    else
        fprintf('(Using raw correlation values for statistical tests)\n\n');
    end
    
    % Session comparison for towards runs
    fprintf('TOWARDS RUNS - Session Comparison:\n');
    if length(session_nums) == 2 && towards_data(1).count > 0 && towards_data(2).count > 0
        common_mice = intersect(towards_data(1).mice, towards_data(2).mice);
        
        fprintf('  Sessions: sess%d vs sess%d\n', session_nums(1), session_nums(2));
        fprintf('  Sample sizes: n=%d vs n=%d\n', towards_data(1).count, towards_data(2).count);
        fprintf('  Means ± SEM: %.3f ± %.3f vs %.3f ± %.3f\n', ...
            towards_data(1).mean, towards_data(1).sem, towards_data(2).mean, towards_data(2).sem);
        
        if length(common_mice) >= 2
            % Paired test
            sess0_values = [];
            sess1_values = [];
            
            for i = 1:length(common_mice)
                mouse_id = common_mice{i};
                idx1 = find(strcmp(towards_data(1).mice, mouse_id));
                idx2 = find(strcmp(towards_data(2).mice, mouse_id));
                
                if ~isempty(idx1) && ~isempty(idx2)
                    sess0_values = [sess0_values, towards_data(1).values(idx1)];
                    sess1_values = [sess1_values, towards_data(2).values(idx2)];
                end
            end
            
            if options.use_fisher_z
                z0 = 0.5 * log((1 + sess0_values) ./ (1 - sess0_values));
                z1 = 0.5 * log((1 + sess1_values) ./ (1 - sess1_values));
                [~, p_value, ~, stats] = ttest(z0, z1);
                test_method = 'paired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest(sess0_values, sess1_values);
                test_method = 'paired t-test on raw correlation values';
            end
            
            fprintf('  Common mice: n=%d (%s)\n', length(common_mice), test_method);
            fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        else
            % Unpaired test
            if options.use_fisher_z
                z0 = 0.5 * log((1 + towards_data(1).values) ./ (1 - towards_data(1).values));
                z1 = 0.5 * log((1 + towards_data(2).values) ./ (1 - towards_data(2).values));
                [~, p_value, ~, stats] = ttest2(z0, z1);
                test_method = 'unpaired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest2(towards_data(1).values, towards_data(2).values);
                test_method = 'unpaired t-test on raw correlation values';
            end
            
            fprintf('  Independent samples (%s)\n', test_method);
            fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        end
        
        if p_value < 0.001
            fprintf('  Result: *** (p < 0.001)\n');
        elseif p_value < 0.01
            fprintf('  Result: ** (p < 0.01)\n');
        elseif p_value < 0.05
            fprintf('  Result: * (p < 0.05)\n');
        else
            fprintf('  Result: n.s. (p ≥ 0.05)\n');
        end
    else
        fprintf('  Insufficient data for comparison\n');
    end
    
    fprintf('\n');
    
    % Session comparison for away runs
    fprintf('AWAY RUNS - Session Comparison:\n');
    if length(session_nums) == 2 && away_data(1).count > 0 && away_data(2).count > 0
        common_mice = intersect(away_data(1).mice, away_data(2).mice);
        
        fprintf('  Sessions: sess%d vs sess%d\n', session_nums(1), session_nums(2));
        fprintf('  Sample sizes: n=%d vs n=%d\n', away_data(1).count, away_data(2).count);
        fprintf('  Means ± SEM: %.3f ± %.3f vs %.3f ± %.3f\n', ...
            away_data(1).mean, away_data(1).sem, away_data(2).mean, away_data(2).sem);
        
        if length(common_mice) >= 2
            % Paired test
            sess0_values = [];
            sess1_values = [];
            
            for i = 1:length(common_mice)
                mouse_id = common_mice{i};
                idx1 = find(strcmp(away_data(1).mice, mouse_id));
                idx2 = find(strcmp(away_data(2).mice, mouse_id));
                
                if ~isempty(idx1) && ~isempty(idx2)
                    sess0_values = [sess0_values, away_data(1).values(idx1)];
                    sess1_values = [sess1_values, away_data(2).values(idx2)];
                end
            end
            
            if options.use_fisher_z
                z0 = 0.5 * log((1 + sess0_values) ./ (1 - sess0_values));
                z1 = 0.5 * log((1 + sess1_values) ./ (1 - sess1_values));
                [~, p_value, ~, stats] = ttest(z0, z1);
                test_method = 'paired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest(sess0_values, sess1_values);
                test_method = 'paired t-test on raw correlation values';
            end
            
            fprintf('  Common mice: n=%d (%s)\n', length(common_mice), test_method);
            fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        else
            % Unpaired test
            if options.use_fisher_z
                z0 = 0.5 * log((1 + away_data(1).values) ./ (1 - away_data(1).values));
                z1 = 0.5 * log((1 + away_data(2).values) ./ (1 - away_data(2).values));
                [~, p_value, ~, stats] = ttest2(z0, z1);
                test_method = 'unpaired t-test on Fisher z-transformed values';
            else
                [~, p_value, ~, stats] = ttest2(away_data(1).values, away_data(2).values);
                test_method = 'unpaired t-test on raw correlation values';
            end
            
            fprintf('  Independent samples (%s)\n', test_method);
            fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
        end
        
        if p_value < 0.001
            fprintf('  Result: *** (p < 0.001)\n');
        elseif p_value < 0.01
            fprintf('  Result: ** (p < 0.01)\n');
        elseif p_value < 0.05
            fprintf('  Result: * (p < 0.05)\n');
        else
            fprintf('  Result: n.s. (p ≥ 0.05)\n');
        end
    else
        fprintf('  Insufficient data for comparison\n');
    end
    
    fprintf('\n');
    
    % Towards vs Away comparison for each session
    for i = 1:length(session_nums)
        fprintf('TOWARDS vs AWAY - sess%d:\n', session_nums(i));
        
        if towards_data(i).count > 0 && away_data(i).count > 0
            common_mice = intersect(towards_data(i).mice, away_data(i).mice);
            
            fprintf('  Sample sizes: Towards n=%d, Away n=%d\n', towards_data(i).count, away_data(i).count);
            fprintf('  Means ± SEM: %.3f ± %.3f vs %.3f ± %.3f\n', ...
                towards_data(i).mean, towards_data(i).sem, away_data(i).mean, away_data(i).sem);
            
            if length(common_mice) >= 2
                % Paired test
                towards_values = [];
                away_values = [];
                
                for j = 1:length(common_mice)
                    mouse_id = common_mice{j};
                    idx_towards = find(strcmp(towards_data(i).mice, mouse_id));
                    idx_away = find(strcmp(away_data(i).mice, mouse_id));
                    
                    if ~isempty(idx_towards) && ~isempty(idx_away)
                        towards_values = [towards_values, towards_data(i).values(idx_towards)];
                        away_values = [away_values, away_data(i).values(idx_away)];
                    end
                end
                
                if options.use_fisher_z
                    z_towards = 0.5 * log((1 + towards_values) ./ (1 - towards_values));
                    z_away = 0.5 * log((1 + away_values) ./ (1 - away_values));
                    [~, p_value, ~, stats] = ttest(z_towards, z_away);
                    test_method = 'paired t-test on Fisher z-transformed values';
                else
                    [~, p_value, ~, stats] = ttest(towards_values, away_values);
                    test_method = 'paired t-test on raw correlation values';
                end
                
                fprintf('  Common mice: n=%d (%s)\n', length(common_mice), test_method);
                fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
            else
                % Unpaired test
                if options.use_fisher_z
                    z_towards = 0.5 * log((1 + towards_data(i).values) ./ (1 - towards_data(i).values));
                    z_away = 0.5 * log((1 + away_data(i).values) ./ (1 - away_data(i).values));
                    [~, p_value, ~, stats] = ttest2(z_towards, z_away);
                    test_method = 'unpaired t-test on Fisher z-transformed values';
                else
                    [~, p_value, ~, stats] = ttest2(towards_data(i).values, away_data(i).values);
                    test_method = 'unpaired t-test on raw correlation values';
                end
                
                fprintf('  Independent samples (%s)\n', test_method);
                fprintf('  t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p_value);
            end
            
            if p_value < 0.001
                fprintf('  Result: *** (p < 0.001)\n');
            elseif p_value < 0.01
                fprintf('  Result: ** (p < 0.01)\n');
            elseif p_value < 0.05
                fprintf('  Result: * (p < 0.05)\n');
            else
                fprintf('  Result: n.s. (p ≥ 0.05)\n');
            end
        else
            fprintf('  Insufficient data for comparison\n');
        end
        
        fprintf('\n');
    end
    
    fprintf('=== END STATISTICAL SUMMARY ===\n\n');
end