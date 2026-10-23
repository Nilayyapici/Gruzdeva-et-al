function plot_distance_time_single_source(runs_fed_food, runs_fasted_food, runs_fed_gel, runs_fasted_gel, options)
    % PLOT_DISTANCE_TIME_SINGLE_SOURCE Creates a figure comparing distance vs time
    % across sessions for a single-source maze, with flexible condition selection
    %
    % This function plots how distance to food changes over time for runs.
    %
    % Inputs:
    %   runs_fed_food - Run data for fed mice with food source
    %   runs_fasted_food - Run data for fasted mice with food source
    %   runs_fed_gel - Run data for fed mice with gel source
    %   runs_fasted_gel - Run data for fed mice with gel source
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of session names to include (e.g., {'sess0', 'sess1'})
    %       .conditions - Cell array of conditions to include (e.g., {'fed_food', 'fasted_food'})
    %                      Possible values: 'fed_food', 'fasted_food', 'fed_gel', 'fasted_gel'
    %                      Default: All conditions
    %       .ylim - Y-axis limits [min max], default [0 50]
    %       .smoothing - Window size for smoothing, default 5
    %       .figure_size - Figure size [width height], default [1200 900]
    %       .title - Main figure title, default 'Distance vs Time - Session Comparison'
    %       .plot_sem - Whether to plot SEM shading, default true
    %       .session_colors - Custom colors for sessions (optional)
    %       .normalize_time - Whether to normalize time to [0,1], default false
    
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
    if ~isfield(options, 'ylim'), options.ylim = [0 50]; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'figure_size'), options.figure_size = [1200 900]; end
    if ~isfield(options, 'title'), options.title = 'Distance vs Time - Session Comparison'; end
    if ~isfield(options, 'plot_sem'), options.plot_sem = true; end
    if ~isfield(options, 'normalize_time'), options.normalize_time = false; end
    
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
    fprintf('Comparing distance across sessions: ');
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
    figure('Name', 'Distance vs Time - Session Comparison', 'Position', [100, 100, options.figure_size(1), options.figure_size(2)]);
    
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
    
    % Determine number of rows based on selected conditions
    num_rows = length(options.conditions);
    
    % Initialize data structure to collect all distance values by time bin
    all_data = struct();
    
    % Process each condition and organize data
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        
        if ~isKey(condition_data_map, cond) || isempty(condition_data_map(cond))
            fprintf('No data available for condition: %s\n', cond);
            continue;
        end
        
        run_data = condition_data_map(cond);
        
        % Process each mouse in this condition
        for m = 1:length(run_data)
            mouse_id = run_data(m).mouse_id;
            session = run_data(m).session;
            
            % Skip if session not in sessions to plot
            sess_str = ['sess' num2str(session)];
            if ~ismember(sess_str, options.sessions)
                continue;
            end
            
            runs = run_data(m).runs;
            
            % Skip if no runs
            if isempty(runs)
                continue;
            end
            
            % Process each run
            for r = 1:length(runs)
                run = runs(r);
                
                % Store the run info for later processing
                run_info = struct();
                run_info.distance = run.distance;
                run_info.time = run.time - run.time(1); % Start time at 0
                run_info.type = run.type;
                run_info.duration = run.duration;
                run_info.mouse_id = mouse_id;
                
                % Create the key for this combination
                key = sprintf('%s_sess%d_%s', cond, session, run.type);
                
                % Initialize the array if needed
                if ~isfield(all_data, key)
                    all_data.(key) = [];
                end
                
                % Add this run to the data structure
                all_data.(key) = [all_data.(key); {run_info}];
            end
        end
    end
    
    % Create subplot layout (num_conditions rows, 2 columns for towards/away)
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
                
                % Skip if no data for this combination
                if ~isfield(all_data, key) || isempty(all_data.(key))
                    continue;
                end
                
                % Get all runs for this condition, session, and run type
                runs_list = all_data.(key);
                
                % Time-normalize and bin if needed
                if options.normalize_time
                    % Normalize each run's time to [0,1]
                    for i = 1:length(runs_list)
                        run_info = runs_list{i};
                        if run_info.duration > 0
                            if strcmp(run_info.type, 'towards')
                                % For towards runs, make time negative and end at 0
                                run_info.time = (run_info.time / run_info.duration) - 1;
                            else
                                % For away runs, start at 0
                                run_info.time = run_info.time / run_info.duration;
                            end
                            runs_list{i} = run_info;
                        end
                    end
                    
                    % Common normalized time grid
                    if strcmp(type, 'towards')
                        common_time = linspace(-1, 0, 100); % Negative to 0 for towards
                    else
                        common_time = linspace(0, 1, 100);  % 0 to 1 for away
                    end
                else
                    % Find the maximum duration for time binning
                    max_duration = 0;
                    for i = 1:length(runs_list)
                        max_duration = max(max_duration, runs_list{i}.duration);
                    end
                    
                    % Create common time grid (with a reasonable bin size)
                    bin_size = 0.1; % 100ms bins
                    
                    if strcmp(type, 'towards')
                        % For towards runs, make time negative and end at 0
                        common_time = -ceil(max_duration):(bin_size):0;
                        
                        % Also adjust each run's time to end at 0
                        for i = 1:length(runs_list)
                            run_info = runs_list{i};
                            run_info.time = run_info.time - run_info.duration; % Shift so end time is 0
                            runs_list{i} = run_info;
                        end
                    else
                        % For away runs, start at 0 as before
                        common_time = 0:bin_size:ceil(max_duration);
                    end
                end
                
                % Interpolate each run to the common time grid
                all_dist_interp = nan(length(common_time), length(runs_list));
                
                for i = 1:length(runs_list)
                    run_info = runs_list{i};
                    
                    % Ensure we have enough points to interpolate
                    if length(run_info.time) >= 2
                        % Interpolate this run's distance to the common time grid
                        run_time = run_info.time;
                        run_dist = run_info.distance;
                        
                        % Only interpolate up to the end of this run
                        valid_times = common_time <= max(run_time);
                        if sum(valid_times) >= 2
                            try
                                all_dist_interp(valid_times, i) = interp1(run_time, run_dist, common_time(valid_times), 'linear');
                            catch
                                % If interpolation fails, skip this run
                                warning('Interpolation failed for run with duration %.2f', run_info.duration);
                            end
                        end
                    end
                end
                
                % Calculate mean and SEM for each time bin
                mean_dist = nanmean(all_dist_interp, 2);
                sem_dist = nanstd(all_dist_interp, 0, 2) ./ sqrt(sum(~isnan(all_dist_interp), 2));
                
                % Apply smoothing to the means and SEMs
                if options.smoothing > 1
                    valid = ~isnan(mean_dist);
                    if sum(valid) > options.smoothing
                        mean_dist(valid) = movmean(mean_dist(valid), options.smoothing);
                        sem_dist(valid) = movmean(sem_dist(valid), options.smoothing);
                    end
                end
                
                % Get session label
                sess_label = session_labels{s};
                
                % Plot SEM shading if requested
                if options.plot_sem
                    for i = 1:length(common_time)-1
                        if ~isnan(mean_dist(i)) && ~isnan(mean_dist(i+1)) && ...
                           ~isnan(sem_dist(i)) && ~isnan(sem_dist(i+1))
                            x_patch = [common_time(i), common_time(i+1), common_time(i+1), common_time(i)];
                            y_patch = [mean_dist(i)-sem_dist(i), mean_dist(i+1)-sem_dist(i+1), ...
                                      mean_dist(i+1)+sem_dist(i+1), mean_dist(i)+sem_dist(i)];
                            patch(x_patch, y_patch, options.session_colors{s}, 'EdgeColor', 'none', ...
                                  'FaceAlpha', 0.2, 'HandleVisibility', 'off');
                        end
                    end
                end
                
                % Plot the smoothed line with session label
                valid = ~isnan(mean_dist);
                if sum(valid) > 0
                    plot(common_time(valid), mean_dist(valid), 'Color', options.session_colors{s}, 'LineWidth', 2, 'DisplayName', sess_label);
                end
            end
            
            % Format condition name for display
            cond_parts = strsplit(cond, '_');
            cond_display = [upper(cond_parts{1}(1)) cond_parts{1}(2:end), ' + ', ...
                           upper(cond_parts{2}(1)) cond_parts{2}(2:end)];
            
            % Set subplot title and labels
            title(sprintf('%s - %s', cond_display, direction_titles{t}));
            if options.normalize_time
                if strcmp(type, 'towards')
                    xlabel('Normalized Time to Food (0 = arrival)');
                else
                    xlabel('Normalized Time from Food (0 = departure)');
                end
            else
                if strcmp(type, 'towards')
                    xlabel('Time to Food (s, 0 = arrival)');
                else
                    xlabel('Time from Food (s, 0 = departure)');
                end
            end
            ylabel('Distance to Food');
            
            % Set y-axis limits
            ylim(options.ylim);
            
            % Add legend
            legend('show', 'Location', 'best');
            legend('boxoff');
            grid on;
            hold off;
        end
    end
    
    % Add overall super title
    sgtitle(options.title, 'FontSize', 14);
end