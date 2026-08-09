function plot_dff_distance_single_sourceDA(runs_fed_food, runs_fasted_food, runs_fed_gel, runs_fasted_gel, options)
    % PLOT_DFF_DISTANCE_SINGLE_SOURCE Creates a figure comparing dF/F vs distance
    % across sessions for a single-source maze, with flexible condition selection
    %
    % Inputs:
    %   runs_fed_food - Run data for fed mice with food source
    %   runs_fasted_food - Run data for fasted mice with food source
    %   runs_fed_gel - Run data for fed mice with gel source
    %   runs_fasted_gel - Run data for fasted mice with gel source
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to include (e.g., {'sess0', 'sess1', 'sess2'})
    %       .conditions - Cell array of conditions to include (e.g., {'fed_food', 'fasted_food'})
    %                      Possible values: 'fed_food', 'fasted_food', 'fed_gel', 'fasted_gel'
    %                      Default: All conditions
    %       .ylim - Y-axis limits [min max], default [-2 2]
    %       .smoothing - Window size for smoothing, default 5
    %       .figure_size - Figure size [width height], default [1200 900]
    %       .title - Main figure title, default 'dF/F vs Distance - Session Comparison'
    %       .plot_sem - Whether to plot SEM shading, default true
    %       .session_colors - Custom colors for sessions (optional)
    
    % Default options
    if nargin < 5
        options = struct();
    end
    
    % Set default options if not provided
    if ~isfield(options, 'sessions')
        % Find all unique session numbers across all run data
        all_sessions = [];
        for data_set = {runs_fed_food, runs_fasted_food, runs_fed_gel, runs_fasted_gel}
            run_data = data_set{1};
            if ~isempty(run_data)
                for m = 1:length(run_data)
                    if ~ismember(run_data(m).session, all_sessions)
                        all_sessions = [all_sessions, run_data(m).session];
                    end
                end
            end
        end
        sessions_to_plot_nums = sort(all_sessions);
        
        % Convert numbers to strings
        options.sessions = cell(1, length(sessions_to_plot_nums));
        for i = 1:length(sessions_to_plot_nums)
            options.sessions{i} = ['sess' num2str(sessions_to_plot_nums(i))];
        end
    end
    
    % Set default conditions if not provided (all conditions)
    if ~isfield(options, 'conditions')
        options.conditions = {'fed_food', 'fasted_food', 'fed_gel', 'fasted_gel'};
    end

    % Other default options
    if ~isfield(options, 'ylim'), options.ylim = [-1.5 1.5]; end
    if ~isfield(options, 'xlim'), options.xlim = [0 210]; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'figure_size'), options.figure_size = [1200 900]; end
    if ~isfield(options, 'title'), options.title = 'dF/F vs Distance - Session Comparison'; end
    if ~isfield(options, 'plot_sem'), options.plot_sem = true; end

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
    
    % Display which sessions and conditions are being plotted
    fprintf('Comparing dF/F across sessions: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf('\nFor conditions: ');
    for i = 1:length(options.conditions)
        fprintf('%s ', options.conditions{i});
    end
    fprintf('\n');
    
    % Create session name mapping for plot labels
    session_labels = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        sess_num = session_numbers(i);
        session_labels{i} = sprintf('Session %d', sess_num);
    end
    
    % Create the figure
    figure('Name', 'dF/F vs Distance - Session Comparison', 'Position', [100, 100, options.figure_size(1), options.figure_size(2)]);
    
    % Subplot titles
    direction_titles = {'Towards', 'Away'};
    
    % Define default colors for sessions
    default_colors = [
        0, 0.4470, 0.7410;  % Blue
        0.8500, 0.3250, 0.0980;  % Orange
        0.9290, 0.6940, 0.1250;  % Yellow
        0.4940, 0.1840, 0.5560;  % Purple
        0.4660, 0.6740, 0.1880;  % Green
        0.3010, 0.7450, 0.9330;  % Light Blue
        0.6350, 0.0780, 0.1840;  % Dark Red
    ];
    
    % Use custom session colors if provided, otherwise generate them
    if ~isfield(options, 'session_colors')
        % Generate session colors
        options.session_colors = cell(1, length(options.sessions));
        for i = 1:length(options.sessions)
            % Use default MATLAB colors or cycle if we have more sessions than colors
            color_idx = mod(i-1, size(default_colors, 1)) + 1;
            options.session_colors{i} = default_colors(color_idx, :);
        end
    end
    
    % Define run types and map conditions to data
    type_names = {'towards', 'away'};
    condition_data_map = containers.Map(...
        {'fed_food', 'fasted_food', 'fed_gel', 'fasted_gel'}, ...
        {runs_fed_food, runs_fasted_food, runs_fed_gel, runs_fasted_gel});
    
    % Find maximum distance across all runs in the sessions we'll plot
    max_distance = 0;
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        if ~isKey(condition_data_map, cond) || isempty(condition_data_map(cond))
            continue;
        end
        
        run_data = condition_data_map(cond);
        for m = 1:length(run_data)
            % Skip if session not in sessions to plot
            sess_str = ['sess' num2str(run_data(m).session)];
            if ~ismember(sess_str, options.sessions)
                continue;
            end
            
            for r = 1:length(run_data(m).runs)
                max_distance = max(max_distance, max(run_data(m).runs(r).distance));
            end
        end
    end
    
    % Create distance bins
    bin_width = 1; % 1 unit bins for better resolution
    dist_bins = 0:bin_width:ceil(max_distance + 5); % Add padding
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialize data structure to collect all z-scored dF/F values by bin
    all_data = struct();
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        for s = 1:length(session_numbers)
            sess = session_numbers(s);
            for t = 1:length(type_names)
                type = type_names{t};
                key = sprintf('%s_sess%d_%s', cond, sess, type);
                all_data.(key) = cell(length(dist_bins)-1, 1);
                for i = 1:length(dist_bins)-1
                    all_data.(key){i} = [];
                end
            end
        end
    end
    
    % Process each condition
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        if ~isKey(condition_data_map, cond) || isempty(condition_data_map(cond))
            warning('No data for condition %s', cond);
            continue;
        end
        
        run_data = condition_data_map(cond);
        
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
                
                % Skip if not a recognized type
                if ~ismember(run.type, type_names)
                    continue;
                end
                
                % Z-score this run's dF/F
                z_scored_dff = (run.dff - dff_mean) / dff_std;
                
                % Bin the z-scored data
                for i = 1:length(dist_bins)-1
                    indices = run.distance >= dist_bins(i) & run.distance < dist_bins(i+1);
                    if any(indices)
                        key = sprintf('%s_sess%d_%s', cond, mouse_sess, run.type);
                        all_data.(key){i} = [all_data.(key){i}; z_scored_dff(indices)];
                    end
                end
            end
        end
    end
    
    % Determine number of rows based on selected conditions
    num_rows = length(options.conditions);
    
    % Create subplot layout (num_conditions rows, 2 columns for directions)
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        
        for t = 1:length(type_names)
            type = type_names{t};
            
            % Calculate subplot index
            sp_idx = (cond_idx-1)*2 + t;
            subplot(num_rows, 2, sp_idx);
            hold on;
            
            % Plot all selected sessions
            for s = 1:length(session_numbers)
                sess = session_numbers(s);
                
                % Create the key for this combination
                key = sprintf('%s_sess%d_%s', cond, sess, type);
                
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
                    
                    % Apply moving average smoothing
                    y_smoothed = movmean(y_valid, options.smoothing);
                    sem_smoothed = movmean(sem_valid, options.smoothing);
                    
                    % Get session label
                    sess_label = session_labels{s};
                    
                    % Plot SEM shading if requested
                    if options.plot_sem
                        % Create shaded error region
                        for i = 1:length(x_valid)-1
                            x_patch = [x_valid(i), x_valid(i+1), x_valid(i+1), x_valid(i)];
                            y_patch = [y_smoothed(i) - sem_smoothed(i), ...
                                      y_smoothed(i+1) - sem_smoothed(i+1), ...
                                      y_smoothed(i+1) + sem_smoothed(i+1), ...
                                      y_smoothed(i) + sem_smoothed(i)];
                            patch(x_patch, y_patch, options.session_colors{s}, 'EdgeColor', 'none', ...
                                'FaceAlpha', 0.2, 'HandleVisibility', 'off');
                        end
                    end
                    
                    % Plot the smoothed line with session label
                    plot(x_valid, y_smoothed, 'Color', options.session_colors{s}, 'LineWidth', 2, 'DisplayName', sess_label);
                end
            end
            
            % Format condition name for display
            cond_parts = strsplit(cond, '_');
            cond_display = [upper(cond_parts{1}(1)) cond_parts{1}(2:end), ' + ', ...
                           upper(cond_parts{2}(1)) cond_parts{2}(2:end)];
            
            % Set subplot title and labels
            title(sprintf('%s - %s', cond_display, direction_titles{t}));
            xlabel('Distance from Food');
            ylabel('Z-scored dF/F');
            
            % Set y-axis limits
            ylim(options.ylim);
            xlim(options.xlim);
            
            % Reverse x-axis for 'towards' runs
            if strcmp(type, 'towards')
                set(gca, 'XDir', 'reverse');
            end
            
            % Add horizontal line at y=0
            line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            % Add legend
            legend('show', 'Location', 'best');
            legend('boxoff');
            grid off;
            box off
        end
    end
    
    % Add overall super title
    sgtitle(options.title, 'FontSize', 14);
end