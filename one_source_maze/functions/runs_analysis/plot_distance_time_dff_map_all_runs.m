function plot_distance_time_dff_map_all_runs(run_data, options)
    % PLOT_DISTANCE_TIME_DFF_MAP_ALL_RUNS Plots distance vs time trajectories with dF/F as color    
    % Default options
    if nargin < 2
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'normalize_method'), options.normalize_method = 'zscore'; end
    if ~isfield(options, 'cmap'), options.cmap = 'jet'; end
    if ~isfield(options, 'clim'), options.clim = [-2 2]; end
    if ~isfield(options, 'time_lim_towards'), options.time_lim_towards = [-20 0]; end % Towards: negative to 0
    if ~isfield(options, 'time_lim_away'), options.time_lim_away = [0 20]; end % Away: 0 to positive
    if ~isfield(options, 'dist_lim'), options.dist_lim = [0 200]; end
    if ~isfield(options, 'grid_size'), options.grid_size = [100, 100]; end
    if ~isfield(options, 'sampling_rate'), options.sampling_rate = 10; end
    if ~isfield(options, 'smooth_factor'), options.smooth_factor = 3; end
    if ~isfield(options, 'run_category'), options.run_category = 'food'; end % Options: 'food', 'not_food', 'all'
    
    % Print data structure info
    fprintf('Input run_data contains %d mice\n', length(run_data));
    
    % Filter by session if specified
    if isfield(options, 'session')
        fprintf('Filtering for session %d...\n', options.session);
        run_data = filter_by_session(run_data, options.session);
        session_str = sprintf('Session %d: ', options.session);
    else
        session_str = 'All Sessions: ';
    end
    
    fprintf('After filtering: %d mice\n', length(run_data));
    
    % Process data for towards and away runs
    towards_data = struct('time', {}, 'distance', {}, 'dff', {});
    away_data = struct('time', {}, 'distance', {}, 'dff', {});
    
    % Count valid runs
    total_runs = 0;
    valid_runs = 0;
    towards_count = 0;
    away_count = 0;
    
    % Determine which run types to include based on run_category option
    switch options.run_category
        case 'food'
            % Only include runs that reached/started from food
            valid_towards_types = {'towards'};
            valid_away_types = {'away'};
            category_str = ' (Food Runs Only)';
        case 'not_food'
            % Only include runs that didn't reach/start from food
            valid_towards_types = {'not_food_towards'};
            valid_away_types = {'not_food_away'};
            category_str = ' (Non-Food Runs Only)';
        case 'all'
            % Include all run types
            valid_towards_types = {'towards', 'not_food_towards'};
            valid_away_types = {'away', 'not_food_away'};
            category_str = ' (All Runs)';
        otherwise
            error('Invalid run_category option. Must be ''food'', ''not_food'', or ''all''');
    end
    
    % Extract run data
    for m = 1:length(run_data)
        for r = 1:length(run_data(m).runs)
            run = run_data(m).runs(r);
            total_runs = total_runs + 1;
            
            % Skip if missing data
            if isempty(run.distance) || isempty(run.dff)
                continue;
            end
            
            % Determine which speed field to use
            if isfield(run, 'calculated_speed') && ~isempty(run.calculated_speed)
                run_speed = run.calculated_speed;
            else
                run_speed = run.speed;
            end
            
            if isempty(run_speed)
                continue;
            end
            
            % Check for NaN or Inf values
            if any(isnan(run.distance)) || any(isnan(run.dff)) || any(isnan(run_speed)) || ...
               any(isinf(run.distance)) || any(isinf(run.dff)) || any(isinf(run_speed))
                continue;
            end
            
            % Calculate time for the run
            num_frames = length(run.dff);
            
            % Create time vector
            if isfield(run, 'time') && ~isempty(run.time)
                raw_time = run.time;
            else
                raw_time = (0:(num_frames-1))' / options.sampling_rate;
            end
            
            % Normalize dF/F
            dff = run.dff;
            switch options.normalize_method
                case 'zscore'
                    if std(dff) > 0
                        dff_norm = (dff - mean(dff)) / std(dff);
                    else
                        continue;
                    end
                case 'minmax'
                    range_dff = max(dff) - min(dff);
                    if range_dff > 0
                        dff_norm = (dff - min(dff)) / range_dff;
                    else
                        continue;
                    end
                case 'none'
                    dff_norm = dff;
                otherwise
                    dff_norm = dff;
            end
            
            % Apply smoothing if requested
            if options.smooth_factor > 0
                dff_norm = smoothdata(dff_norm, 'gaussian', options.smooth_factor);
                dist_smooth = smoothdata(run.distance, 'gaussian', options.smooth_factor);
            else
                dist_smooth = run.distance;
            end
            
            % Check if this run type should be included
            is_towards_run = any(strcmp(run.type, valid_towards_types));
            is_away_run = any(strcmp(run.type, valid_away_types));
            
            if ~is_towards_run && ~is_away_run
                continue; % Skip this run if it doesn't match the category filter
            end
            
            % Calculate time - SIMPLE MIRRORING
            if is_towards_run
                % For towards runs: flip the time so max becomes 0, and 0 becomes -max
                normalized_time = raw_time - min(raw_time);  % Start from 0
                time_data = max(normalized_time) - normalized_time;  % Flip: max->0, 0->max
                time_data = -time_data;  % Make it negative: 0->0, max->-max
            else % away runs
                % For away runs: keep normal time progression from 0
                time_data = raw_time - min(raw_time);     % Start from 0, go positive
            end
            
            % Store data in appropriate category
            run_struct = struct('time', time_data, 'distance', dist_smooth, 'dff', dff_norm);
            
            valid_runs = valid_runs + 1;
            
            if is_towards_run
                towards_data(end+1) = run_struct;
                towards_count = towards_count + 1;
            else % away runs
                away_data(end+1) = run_struct;
                away_count = away_count + 1;
            end
        end
    end

    % Print summary
    fprintf('Total runs processed: %d, Valid runs: %d (%.1f%%)%s\n', total_runs, valid_runs, 100*valid_runs/total_runs, category_str);
    fprintf('Towards runs: %d, Away runs: %d\n', towards_count, away_count);


    figure('Position', [100, 100, 1100, 400], 'Name', [session_str 'Distance vs Time with dF/F Color Map' category_str]);

    % Plot data for towards runs
    subplot(1, 2, 1);
    plot_distance_time_heatmap(towards_data, 'Towards', options, options.time_lim_towards);

    % Plot data for away runs
    subplot(1, 2, 2);
    plot_distance_time_heatmap(away_data, 'Away', options, options.time_lim_away);

    % Add colorbar
    colormap(options.cmap);
    caxis(options.clim);
    cb = colorbar();
    cb.Position = [0.93, 0.2, 0.02, 0.6];
    ylabel(cb, 'Normalized dF/F', 'FontSize', 14);

    % Add overall title
    sgtitle([session_str 'Distance vs Time Trajectories with dF/F Color' category_str], 'FontSize', 14);
end

function filtered_data = filter_by_session(run_data, session_number)
    % Filter run_data to include only the specified session
    filtered_data = struct('mouse_id', {}, 'session', {}, 'runs', {});
    
    for m = 1:length(run_data)
        if run_data(m).session == session_number
            filtered_data = [filtered_data; run_data(m)];
        end
    end
    
    % Check if any data was found for the session
    if isempty(filtered_data)
        warning('No data found for session %d', session_number);
    else
        fprintf('Found %d mice for session %d\n', length(filtered_data), session_number);
    end
end

function plot_distance_time_heatmap(run_data, plot_title, options, time_lim)
    % Create a heatmap of distance vs time with dF/F as color
    
    if isempty(run_data)
        text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(plot_title, 'FontSize', 14);
        xlabel('Time (s)');
        ylabel('Distance (cm)');
        return;
    end
    
    % Initialize grid for the average heatmap
    time_grid = linspace(time_lim(1), time_lim(2), options.grid_size(2));
    dist_grid = linspace(options.dist_lim(1), options.dist_lim(2), options.grid_size(1));
    
    % Create meshgrid for interpolation
    [T, D] = meshgrid(time_grid, dist_grid);
    
    % Initialize count and sum grids for averaging
    dff_sum = zeros(size(T));
    count_grid = zeros(size(T));
    
    % Process each run
    hold on;
    
    % First pass: collect data for average heatmap
    for i = 1:length(run_data)
        run = run_data(i);
        
        % Extract data
        time_data = run.time;
        dist_data = run.distance;
        dff_data = run.dff;
        
        % Only use data within limits
        valid_indices = time_data >= time_lim(1) & time_data <= time_lim(2) & ...
                        dist_data >= options.dist_lim(1) & dist_data <= options.dist_lim(2);
        
        if ~any(valid_indices)
            continue;
        end
        
        time_data = time_data(valid_indices);
        dist_data = dist_data(valid_indices);
        dff_data = dff_data(valid_indices);
        
        % Interpolate this run's dF/F values onto the grid
        for j = 1:length(time_data)
            % Find the nearest grid points
            [~, t_idx] = min(abs(time_grid - time_data(j)));
            [~, d_idx] = min(abs(dist_grid - dist_data(j)));
            
            % Add to sum and count grids
            dff_sum(d_idx, t_idx) = dff_sum(d_idx, t_idx) + dff_data(j);
            count_grid(d_idx, t_idx) = count_grid(d_idx, t_idx) + 1;
        end
    end
    
    % Calculate average dF/F (avoiding division by zero)
    avg_dff = zeros(size(dff_sum));
    valid_grid = count_grid > 0;
    avg_dff(valid_grid) = dff_sum(valid_grid) ./ count_grid(valid_grid);
    
    % Apply smoothing to the heatmap for better visualization
    avg_dff = smoothdata(avg_dff, 1, 'gaussian', 5);
    avg_dff = smoothdata(avg_dff, 2, 'gaussian', 5);
    
    % Apply colormap limits
    avg_dff = max(min(avg_dff, options.clim(2)), options.clim(1));

    % Plot the heatmap
    pcolor(T, D, avg_dff);
    shading interp;

    % Set axis limits and labels
    xlim(time_lim);
    ylim(options.dist_lim);
    
    % Label axes based on run type
    if contains(plot_title, 'Towards')
        xlabel('Time to Food (s)', 'FontSize', 14);
    else
        xlabel('Time from Food (s)', 'FontSize', 14);
    end
    
    ylabel('Distance (cm)', 'FontSize', 14);
    title([plot_title ' (n=' num2str(length(run_data)) ' runs)']);
    
    % Add grid and box
    grid on;
    box on;
    
    % Set consistent colormap
    colormap(options.cmap);
    clim(options.clim);
end