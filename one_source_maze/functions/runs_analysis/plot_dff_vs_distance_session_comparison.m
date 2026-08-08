function plot_dff_vs_distance_session_comparison(run_data, options)
    % Plot dF/F vs distance for food runs, comparing session 0 and session 1
    % - Both towards and away runs on the same subplot
    % - Session 0 runs on top row, Session 1 runs on bottom row
    % - Only for a single food source
    % - Clean, minimalist style with no boxes, axes, or grid
    % - Scale bars instead of traditional axes
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_single_arm_runs_with_speed_check
    %   options - Structure with the following fields (all optional):
    %       .food_source - Food source to analyze (default: 'food')
    %       .smoothing - Window size for smoothing (default: 5)
    %       .max_distance - Maximum distance range to plot (default: 200)
    %       .ylim_dff - Y-axis limits for dF/F (default: automatic)
    %       .distance_points - Number of distance points for interpolation (default: 200)
    %       .min_mice - Minimum number of mice required for a run (default: 1)
    %       .max_runs - Maximum number of runs to plot (default: 10)
    %       .use_zscore - Whether to use z-scored dF/F (default: true)
    %       .gap - Gap between towards and away runs (default: 10)
    %       .scalebar_distance - Length of distance scale bar in units (default: 50)
    %       .scalebar_dff - Height of dF/F scale bar (default: 1)
    %       .scalebar_position - Position of scale bar [x,y] (default: [0.1, 0.1])
    %       .run_category - Category of runs to include: 'food', 'not_food', 'all' (default: 'food')
    
    % Set default options if not provided
    if nargin < 2
        options = struct();
    end
    
    % Extract and validate options, set defaults if not specified
    if ~isfield(options, 'smoothing') || isempty(options.smoothing)
        options.smoothing = 5; % Default smoothing window
    end
    
    if ~isfield(options, 'food_source') || isempty(options.food_source)
        options.food_source = 'food'; % Default food source to analyze
    end
    
    if ~isfield(options, 'max_distance') || isempty(options.max_distance)
        options.max_distance = 200; % Default maximum distance range
    end
    
    if ~isfield(options, 'ylim_dff')
        options.ylim_dff = []; % Default: automatic y-limits for dF/F
    end
    
    if ~isfield(options, 'distance_points') || isempty(options.distance_points)
        options.distance_points = 200; % Default number of distance points for interpolation
    end
    
    if ~isfield(options, 'min_mice') || isempty(options.min_mice)
        options.min_mice = 1; % Default minimum number of mice
    end
    
    if ~isfield(options, 'max_runs') || isempty(options.max_runs)
        options.max_runs = 10; % Default maximum number of runs to plot
    end
    
    if ~isfield(options, 'use_zscore') || isempty(options.use_zscore)
        options.use_zscore = true; % Default: use z-scored dF/F values
    end
    
    if ~isfield(options, 'gap') || isempty(options.gap)
        options.gap = 10; % Default gap between towards and away runs
    end
    
    if ~isfield(options, 'scalebar_distance') || isempty(options.scalebar_distance)
        options.scalebar_distance = 50; % Default length of distance scale bar
    end
    
    if ~isfield(options, 'scalebar_dff') || isempty(options.scalebar_dff)
        options.scalebar_dff = 1; % Default height of dF/F scale bar
    end
    
    if ~isfield(options, 'scalebar_position') || isempty(options.scalebar_position)
        options.scalebar_position = [0.1, 0.1]; % Default position [x,y] as fraction of figure size
    end
    
    if ~isfield(options, 'run_category') || isempty(options.run_category)
        options.run_category = 'food'; % Options: 'food', 'not_food', 'all'
    end
    
    % Determine which run types to include based on run_category option
    switch options.run_category
        case 'food'
            % Only include runs that reached/started from food
            valid_towards_types = {'towards'};
            valid_away_types = {'away'};
            category_str = ' (Food Runs)';
        case 'not_food'
            % Only include runs that didn't reach/start from food
            valid_towards_types = {'not_food_towards'};
            valid_away_types = {'not_food_away'};
            category_str = ' (Non-Food Runs)';
        case 'all'
            % Include all run types
            valid_towards_types = {'towards', 'not_food_towards'};
            valid_away_types = {'away', 'not_food_away'};
            category_str = ' (All Runs)';
        otherwise
            error('Invalid run_category option. Must be ''food'', ''not_food'', or ''all''');
    end
    
    % We'll be comparing sessions 0 and 1
    session_nums = [0, 1];
    
    % Display which sessions and food source are being analyzed
    fprintf('Analyzing runs for sessions 0 and 1, food source: %s%s\n', options.food_source, category_str);
    fprintf('Maximum runs to plot: %d, Minimum mice per run: %d\n', options.max_runs, options.min_mice);
    
    % First, organize mice by session
    session_data = struct('session', {}, 'mice', {});
    for m = 1:length(run_data)
        sess = run_data(m).session;
        
        % Only process sessions 0 and 1
        if ~ismember(sess, session_nums)
            continue;
        end
        
        % Find if session already exists in our structure
        sess_idx = find([session_data.session] == sess, 1);
        if isempty(sess_idx)
            % Create new session entry
            sess_idx = length(session_data) + 1;
            session_data(sess_idx).session = sess;
            session_data(sess_idx).mice = [];
        end
        
        % Add mouse to this session
        session_data(sess_idx).mice = [session_data(sess_idx).mice, m];
    end
    
    % If either session is missing, warn and return
    session_indices = [session_data.session];
    if ~ismember(0, session_indices)
        warning('No data found for session 0');
        return;
    end
    if ~ismember(1, session_indices)
        warning('No data found for session 1');
        return;
    end
    
    % Process data from each session to gather run information
    session_run_data = struct();
    
    for s = 1:length(session_data)
        sess = session_data(s).session;
        mice_indices = session_data(s).mice;
        
        % Store session number
        session_run_data.(sprintf('sess%d', sess)).session = sess;
        
        % Process towards and away runs separately
        for direction_idx = 1:2
            if direction_idx == 1
                valid_types = valid_towards_types;
                type_label = 'towards';
            else
                valid_types = valid_away_types;
                type_label = 'away';
            end
            
            % Create empty structure for this run type
            field_name = sprintf('sess%d_%s', sess, type_label);
            session_run_data.(field_name) = struct();
            
            % First, collect all food runs of this type for each mouse
            mouse_food_runs = cell(length(mice_indices), 1);
            max_runs = 0;
            
            for mi = 1:length(mice_indices)
                m = mice_indices(mi);
                mouse_runs = run_data(m).runs;
                
                % Get indices of food source runs of this type (matching valid_types)
                run_indices = [];
                for r = 1:length(mouse_runs)
                    if any(strcmp(mouse_runs(r).type, valid_types))
                        run_indices = [run_indices, r];
                    end
                end
                
                % Store these runs for this mouse
                mouse_food_runs{mi} = run_indices;
                max_runs = max(max_runs, length(run_indices));
            end
            
            % If no runs found, continue
            if max_runs == 0
                warning('No %s runs found for session %d%s', type_label, sess, category_str);
                session_run_data.(field_name).valid_runs = [];
                session_run_data.(field_name).mouse_food_runs = {};
                continue;
            end
            
            % Count how many mice have each run number
            mice_per_run = zeros(max_runs, 1);
            for mi = 1:length(mice_indices)
                for r = 1:length(mouse_food_runs{mi})
                    if r <= max_runs
                        mice_per_run(r) = mice_per_run(r) + 1;
                    end
                end
            end
            
            % Find which runs meet the minimum mice requirement AND are within max_runs limit
            valid_runs = find(mice_per_run >= options.min_mice);
            valid_runs = valid_runs(valid_runs <= options.max_runs); % Limit to max_runs
            
            % Store run data for this session and direction
            session_run_data.(field_name).valid_runs = valid_runs;
            session_run_data.(field_name).mouse_food_runs = mouse_food_runs;
            session_run_data.(field_name).mice_indices = mice_indices;
            session_run_data.(field_name).max_runs = max_runs;
        end
    end
    
    % Find all valid runs from both sessions and both directions
    all_valid_runs = unique([
        session_run_data.sess0_towards.valid_runs;
        session_run_data.sess0_away.valid_runs;
        session_run_data.sess1_towards.valid_runs;
        session_run_data.sess1_away.valid_runs
    ]);
    
    % If no valid runs, return
    if isempty(all_valid_runs)
        warning('No runs with >= %d mice found for sessions 0 and 1 within first %d runs%s', options.min_mice, options.max_runs, category_str);
        return;
    end
    
    % Sort the valid runs in ascending order and limit to max_runs
    all_valid_runs = sort(all_valid_runs);
    all_valid_runs = all_valid_runs(all_valid_runs <= options.max_runs); % Ensure we don't exceed max_runs
    num_valid_runs = length(all_valid_runs);
    
    fprintf('Plotting %d runs (runs %s)%s\n', num_valid_runs, mat2str(all_valid_runs), category_str);
    
    % Calculate global min/max values for consistent plotting
    [min_dff, max_dff] = calculate_global_dff_limits(run_data, session_run_data, session_nums, ...
                                                      all_valid_runs, num_valid_runs, options, ...
                                                      valid_towards_types, valid_away_types);
    
    % If ylim_dff is specified, use it instead
    if ~isempty(options.ylim_dff)
        min_dff = options.ylim_dff(1);
        max_dff = options.ylim_dff(2);
    end
    
    % Create figure with 2 rows (session 0 and 1) and columns equal to number of valid runs
    figure('Name', sprintf('dF/F vs Distance - Sessions Comparison for %s (First %d Runs)%s', options.food_source, options.max_runs, category_str), ...
           'Position', [100, 100, 250*num_valid_runs, 500], ...
           'Color', 'w'); % White background
    
    % Process each valid run for each session
    first_subplot_handle = plot_all_sessions_and_runs(run_data, session_run_data, session_nums, ...
                                                       all_valid_runs, num_valid_runs, ...
                                                       min_dff, max_dff, options, ...
                                                       valid_towards_types, valid_away_types);
    
    % Add scale bars on the first subplot only
    add_scale_bars(first_subplot_handle, options, min_dff, max_dff);
    
    % Add overall title
    figure_title = sprintf('dF/F vs Distance Comparison - %s (First %d Runs)%s', options.food_source, options.max_runs, category_str);
    if options.use_zscore
        figure_title = [figure_title, ' (Z-scored)'];
    end
    sgtitle(figure_title, 'FontSize', 14);
    
    % Add text annotations for session titles on the left side
    annotation('textbox', [0.01, 0.75, 0.08, 0.1], ...
               'String', 'Session 0', 'FontSize', 14, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
               'EdgeColor', 'none', 'BackgroundColor', 'none', 'Rotation', 90);
    
    annotation('textbox', [0.01, 0.25, 0.08, 0.1], ...
               'String', 'Session 1', 'FontSize', 14, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
               'EdgeColor', 'none', 'BackgroundColor', 'none', 'Rotation', 90);
end

function [min_dff, max_dff] = calculate_global_dff_limits(run_data, session_run_data, session_nums, ...
                                                           all_valid_runs, num_valid_runs, options, ...
                                                           valid_towards_types, valid_away_types)
    % Calculate global min/max dF/F values for consistent plotting across all subplots
    
    min_dff = Inf;
    max_dff = -Inf;
    
    % Process all valid runs for both sessions to get global min/max dF/F
    for sess = session_nums
        for run_idx = 1:num_valid_runs
            run_num = all_valid_runs(run_idx);
            
            % Check data for both towards and away runs
            for direction_idx = 1:2
                if direction_idx == 1
                    valid_types = valid_towards_types;
                    type_label = 'towards';
                else
                    valid_types = valid_away_types;
                    type_label = 'away';
                end
                
                field_name = sprintf('sess%d_%s', sess, type_label);
                
                % Skip if this field doesn't exist or has no valid runs
                if ~isfield(session_run_data, field_name) || ...
                   isempty(session_run_data.(field_name).valid_runs) || ...
                   ~ismember(run_num, session_run_data.(field_name).valid_runs)
                    continue;
                end
                
                % Get averaged dF/F for this run
                avg_dff = process_run_data(run_data, session_run_data.(field_name), run_num, ...
                                          type_label, options, false);
                
                % Update global min/max
                if ~isempty(avg_dff)
                    valid_avg = avg_dff(~isnan(avg_dff));
                    if ~isempty(valid_avg)
                        min_dff = min(min_dff, min(valid_avg));
                        max_dff = max(max_dff, max(valid_avg));
                    end
                end
            end
        end
    end
    
    % Add padding
    padding = (max_dff - min_dff) * 0.1;
    min_dff = min_dff - padding;
    max_dff = max_dff + padding;
end

function first_subplot_handle = plot_all_sessions_and_runs(run_data, session_run_data, session_nums, ...
                                                            all_valid_runs, num_valid_runs, ...
                                                            min_dff, max_dff, options, ...
                                                            valid_towards_types, valid_away_types)
    % Plot all sessions and runs
    
    first_subplot_handle = [];
    
    % Calculate extended x-limits for all plots (to accommodate scale bar)
    xlim_extended = [-options.max_distance * 1.6, options.max_distance];
    
    for sess_idx = 1:2
        sess = session_nums(sess_idx);
        
        for run_idx = 1:num_valid_runs
            run_num = all_valid_runs(run_idx);
            
            % Create subplot
            ax = subplot(2, num_valid_runs, (sess_idx-1)*num_valid_runs + run_idx);
            
            if sess_idx == 1 && run_idx == 1
                first_subplot_handle = ax;
            end
            
            hold on;
            box off;
            set(ax, 'XColor', 'none', 'YColor', 'none');
            
            % Process both directions for this run
            for direction_idx = 1:2
                if direction_idx == 1
                    valid_types = valid_towards_types;
                    type_label = 'towards';
                    color = [0, 0, 0.8]; % Blue
                else
                    valid_types = valid_away_types;
                    type_label = 'away';
                    color = [0.8, 0, 0]; % Red
                end
                
                field_name = sprintf('sess%d_%s', sess, type_label);
                
                % Skip if no valid runs
                if ~isfield(session_run_data, field_name) || ...
                   isempty(session_run_data.(field_name).valid_runs) || ...
                   ~ismember(run_num, session_run_data.(field_name).valid_runs)
                    continue;
                end
                
                % Get processed data for this run
                [distance_grid, avg_dff, sem_dff, mice_with_run] = ...
                    get_processed_run_data(run_data, session_run_data.(field_name), run_num, ...
                                          type_label, options);
                
                if isempty(distance_grid)
                    continue;
                end
                
                % Plot with shaded error bars
                shadedErrorBar(distance_grid, avg_dff, sem_dff, {'Color', color, 'LineWidth', 1.5}, 0.2);
                
                % Add mouse count label
                if strcmp(type_label, 'towards')
                    text_x = -options.max_distance * 0.95;
                else
                    text_x = options.max_distance * 0.05;
                end
                text_y = max_dff - 0.1 * (max_dff - min_dff);
                text(text_x, text_y, sprintf('n=%d', mice_with_run), 'Color', color, ...
                     'FontSize', 8, 'VerticalAlignment', 'top');
            end
            
            % Set SAME limits for ALL subplots (including extended space for scale bar)
            xlim(xlim_extended);
            ylim([min_dff, max_dff]);
            
            % Add reference line (only in the visible data range)
            plot([-options.max_distance, options.max_distance], [0, 0], 'k--', 'LineWidth', 0.5);
            title(sprintf('Run %d', run_num), 'FontSize', 10);
            grid off;
        end
    end
end

function [distance_grid, avg_dff, sem_dff, mice_with_run] = get_processed_run_data(run_data, field_data, run_num, type_label, options)
    % Process run data and return averaged results
    
    % Setup distance grid
    if strcmp(type_label, 'towards')
        distance_grid = linspace(-options.max_distance, -options.gap/2, options.distance_points);
    else
        distance_grid = linspace(options.gap/2, options.max_distance, options.distance_points);
    end
    
    % Collect interpolated data from all mice
    interp_dff = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options);
    
    % Return empty if no valid data
    if isempty(interp_dff) || all(isnan(interp_dff(:)))
        distance_grid = [];
        avg_dff = [];
        sem_dff = [];
        mice_with_run = 0;
        return;
    end
    
    % Calculate statistics
    mice_with_run = size(interp_dff, 1);
    avg_dff = nanmean(interp_dff, 1);
    
    if mice_with_run > 1
        valid_count = sum(~isnan(interp_dff), 1);
        valid_count(valid_count == 0) = 1;
        sem_dff = nanstd(interp_dff, 0, 1) ./ sqrt(valid_count);
    else
        sem_dff = zeros(size(avg_dff));
    end
    
    % Apply smoothing
    if sum(~isnan(avg_dff)) > options.smoothing
        valid_dff = ~isnan(avg_dff);
        if sum(valid_dff) > options.smoothing
            avg_dff(valid_dff) = movmean(avg_dff(valid_dff), options.smoothing);
            sem_dff(valid_dff) = movmean(sem_dff(valid_dff), options.smoothing);
        end
    end
end

function interp_dff = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options)
    % Collect and interpolate data from all mice for a specific run
    
    mouse_food_runs = field_data.mouse_food_runs;
    mice_indices = field_data.mice_indices;
    interp_dff = [];
    
    for mi = 1:length(mice_indices)
        m = mice_indices(mi);
        mouse_run_indices = mouse_food_runs{mi};
        
        if run_num > length(mouse_run_indices)
            continue;
        end
        
        run_idx = mouse_run_indices(run_num);
        run = run_data(m).runs(run_idx);
        
        % Apply z-scoring if enabled
        if options.use_zscore
            all_dff = [];
            for r = 1:length(run_data(m).runs)
                all_dff = [all_dff; run_data(m).runs(r).dff];
            end
            
            dff_mean = mean(all_dff);
            dff_std = std(all_dff);
            
            if dff_std < 1e-10
                continue;
            end
            
            run.dff = (run.dff - dff_mean) / dff_std;
        end
        
        % Process distance and dF/F
        [norm_distance, dff_data] = normalize_distance_data(run, type_label, options);
        
        if length(norm_distance) < 5
            continue;
        end
        
        % Interpolate to common grid
        interp_dff_mouse = interp1(norm_distance, dff_data, distance_grid, 'linear', NaN);
        interp_dff = [interp_dff; interp_dff_mouse];
    end
end

function [norm_distance, dff_data] = normalize_distance_data(run, type_label, options)
    % Normalize distance data based on run type
    
    distance_data = run.distance;
    dff_data = run.dff;
    min_dist = min(distance_data);
    
    if strcmp(type_label, 'towards')
        norm_distance = -(distance_data - min_dist) - options.gap/2;
        
        if min(norm_distance) < -options.max_distance
            valid_indices = norm_distance >= -options.max_distance;
            norm_distance = norm_distance(valid_indices);
            dff_data = dff_data(valid_indices);
        end
        
        [norm_distance, sort_idx] = sort(norm_distance, 'descend');
    else
        norm_distance = (distance_data - min_dist) + options.gap/2;
        
        if max(norm_distance) > options.max_distance
            valid_indices = norm_distance <= options.max_distance;
            norm_distance = norm_distance(valid_indices);
            dff_data = dff_data(valid_indices);
        end
        
        [norm_distance, sort_idx] = sort(norm_distance, 'ascend');
    end
    
    dff_data = dff_data(sort_idx);
    [norm_distance, unique_idx] = unique(norm_distance, 'stable');
    dff_data = dff_data(unique_idx);
end

function avg_dff = process_run_data(run_data, field_data, run_num, type_label, options, return_full_data)
    % Simplified version that only returns avg_dff for limit calculation
    
    if strcmp(type_label, 'towards')
        distance_grid = linspace(-options.max_distance, -options.gap/2, options.distance_points);
    else
        distance_grid = linspace(options.gap/2, options.max_distance, options.distance_points);
    end
    
    interp_dff = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options);
    
    if isempty(interp_dff) || all(isnan(interp_dff(:)))
        avg_dff = [];
    else
        avg_dff = nanmean(interp_dff, 1);
    end
end

function add_scale_bars(first_subplot_handle, options, min_dff, max_dff)
    % Add scale bars to the first subplot
    
    if isempty(first_subplot_handle)
        return;
    end
    
    subplot(first_subplot_handle);
    % Don't change xlim here - it's already set correctly in plot_all_sessions_and_runs
    
    scalebar_x = -options.max_distance * 1.5;
    scalebar_y = min_dff + (max_dff - min_dff) * 0.001;
    
    % Distance scale bar (horizontal)
    plot([scalebar_x, scalebar_x + options.scalebar_distance], [scalebar_y, scalebar_y], 'k-', 'LineWidth', 2);
    
    % dF/F scale bar (vertical)
    plot([scalebar_x, scalebar_x], [scalebar_y, scalebar_y + options.scalebar_dff], 'k-', 'LineWidth', 2);
    
    % Add labels
    if options.use_zscore
        text(scalebar_x, scalebar_y + options.scalebar_dff * 1.05, sprintf('%.1f ZdF/F', options.scalebar_dff), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    else
        text(scalebar_x, scalebar_y + options.scalebar_dff * 1.05, sprintf('%.1f dF/F', options.scalebar_dff), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    
    text(scalebar_x + options.scalebar_distance, scalebar_y * 1.0, sprintf('%d cm', options.scalebar_distance), ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 7);
end

function h = shadedErrorBar(x, y, errBar, lineProps, patchSaturation)
    % Function to create shaded error bars
    
    if nargin < 4 || isempty(lineProps)
        lineProps = {'b-', 'LineWidth', 1.5};
    end
    
    if nargin < 5 || isempty(patchSaturation)
        patchSaturation = 0.5;
    end
    
    % Extract line color
    if isstruct(lineProps)
        lineColor = lineProps.Color;
    elseif iscell(lineProps)
        colorIdx = find(strcmpi('Color', lineProps));
        if ~isempty(colorIdx) && colorIdx < length(lineProps)
            lineColor = lineProps{colorIdx+1};
        else
            lineColor = [0, 0, 1];
        end
    elseif ischar(lineProps)
        if lineProps(1) == 'b', lineColor = [0, 0, 1];
        elseif lineProps(1) == 'r', lineColor = [1, 0, 0];
        elseif lineProps(1) == 'g', lineColor = [0, 1, 0];
        else, lineColor = [0, 0, 1];
        end
    else
        lineColor = [0, 0, 1];
    end
    
    % Handle NaN values
    valid = ~isnan(y) & ~isnan(errBar);
    if sum(valid) < 2
        h = [];
        return;
    end
    
    x_valid = x(valid);
    y_valid = y(valid);
    errBar_valid = errBar(valid);
    
    % Create patch
    patchColor = lineColor * patchSaturation + (1 - patchSaturation);
    
    if length(x_valid) >= 2
        patch([x_valid fliplr(x_valid)], [y_valid+errBar_valid fliplr(y_valid-errBar_valid)], ...
              patchColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    end
    
    % Plot line
    if iscell(lineProps)
        h = plot(x_valid, y_valid, lineProps{:});
    elseif isstruct(lineProps)
        fields = fieldnames(lineProps);
        propPairs = cell(1, 2*length(fields));
        for i = 1:length(fields)
            propPairs{2*i-1} = fields{i};
            propPairs{2*i} = lineProps.(fields{i});
        end
        h = plot(x_valid, y_valid, propPairs{:});
    else
        h = plot(x_valid, y_valid, lineProps, 'LineWidth', 1.5);
    end
end