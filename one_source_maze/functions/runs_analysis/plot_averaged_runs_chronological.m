function plot_averaged_runs_chronological(run_data, session_number, options)
    % PLOT_AVERAGED_RUNS_CHRONOLOGICAL Plot averaged runs across all mice chronologically
    % Creates one figure with averaged runs in chronological order
    % Towards runs on top row, away runs on bottom row
    % Each subplot shows the average across all mice for that run number
    %
    % Inputs:
    %   run_data - Structure from analyze_single_arm_runs_with_speed_check
    %   session_number - Which session to plot (required)
    %   options - Structure with optional fields:
    %       .normalize_dff - Whether to z-score normalize dF/F per run (default: false)
    %       .dist_lim - Distance axis limits [min max] (default: [0 200])
    %       .dff_lim - dF/F axis limits [min max] (default: automatic)
    %       .max_runs - Maximum number of runs to plot (default: 10)
    %       .distance_points - Number of distance points for interpolation (default: 100)
    %       .run_category - Category of runs to include: 'food', 'not_food', 'all' (default: 'food')
    
    % Check inputs
    if nargin < 2
        error('Session number is required. Usage: plot_averaged_runs_chronological(run_data, session_number, options)');
    end
    
    % Default options
    if nargin < 3
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'normalize_dff')
        options.normalize_dff = true;  % Default to z-score normalization
    end
    if ~isfield(options, 'dist_lim')
        options.dist_lim = [0 200];
    end
    if ~isfield(options, 'dff_lim')
        options.dff_lim = [];
    end
    if ~isfield(options, 'max_runs')
        options.max_runs = 10;
    end
    if ~isfield(options, 'distance_points')
        options.distance_points = 100;
    end
    if ~isfield(options, 'run_category')
        options.run_category = 'food';  % Options: 'food', 'not_food', 'all'
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
    
    % Filter by session
    fprintf('Filtering for session %d%s...\n', session_number, category_str);
    session_data = run_data([run_data.session] == session_number);
    
    if isempty(session_data)
        fprintf('No data found for session %d\n', session_number);
        return;
    end
    
    fprintf('Processing %d mice for session %d%s...\n', length(session_data), session_number, category_str);
    
    % Collect all runs from all mice, organized by run number
    towards_runs_by_number = cell(options.max_runs, 1);
    away_runs_by_number = cell(options.max_runs, 1);
    
    % Loop through each mouse and collect runs by their chronological order
    for m = 1:length(session_data)
        mouse_info = session_data(m);
        runs = mouse_info.runs;
        
        if isempty(runs)
            continue;
        end
        
        % Sort runs by start time (chronological order)
        start_times = [runs.start_time];
        [~, sort_idx] = sort(start_times);
        runs = runs(sort_idx);
        
        % Separate towards and away runs while preserving order, filtering by category
        towards_runs = [];
        away_runs = [];
        
        for r = 1:length(runs)
            if any(strcmp(runs(r).type, valid_towards_types))
                towards_runs = [towards_runs; runs(r)];
            elseif any(strcmp(runs(r).type, valid_away_types))
                away_runs = [away_runs; runs(r)];
            end
        end
        
        % Add towards runs to collection (by run number)
        for r = 1:min(length(towards_runs), options.max_runs)
            if isempty(towards_runs_by_number{r})
                towards_runs_by_number{r} = {};
            end
            towards_runs_by_number{r}{end+1} = towards_runs(r);
        end
        
        % Add away runs to collection (by run number)
        for r = 1:min(length(away_runs), options.max_runs)
            if isempty(away_runs_by_number{r})
                away_runs_by_number{r} = {};
            end
            away_runs_by_number{r}{end+1} = away_runs(r);
        end
    end
    
    % Find maximum number of runs that have data
    max_towards = 0;
    max_away = 0;
    for r = 1:options.max_runs
        if ~isempty(towards_runs_by_number{r})
            max_towards = r;
        end
        if ~isempty(away_runs_by_number{r})
            max_away = r;
        end
    end
    
    max_runs_to_plot = max(max_towards, max_away);
    
    if max_runs_to_plot == 0
        fprintf('No runs found for session %d%s\n', session_number, category_str);
        return;
    end
    
    % Create figure
    fig_title = sprintf('Session %d - Averaged Runs Across All Mice%s', session_number, category_str);
    figure('Position', [100, 100, 200*max_runs_to_plot, 600], 'Name', fig_title);
    
    % Plot averaged towards runs (top row)
    for r = 1:max_towards
        runs_for_this_number = towards_runs_by_number{r};
        
        if isempty(runs_for_this_number)
            continue;
        end
        
        % Create subplot (top row)
        subplot(2, max_runs_to_plot, r);
        
        % Average the runs for this run number
        [avg_distance, avg_dff, sem_dff, num_mice] = average_runs(runs_for_this_number, 'towards', options);
        
        if isempty(avg_distance)
            continue;
        end
        
        % Plot averaged run with SEM shading
        plot_with_sem(avg_distance, avg_dff, sem_dff, [0.2 0.6 0.8]);
        hold on;
        
        % Mark start and end points
        plot(avg_distance(1), avg_dff(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
        plot(avg_distance(end), avg_dff(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
        
        % Set axis limits
        if ~isempty(options.dff_lim)
            ylim(options.dff_lim);
        end
        
        % Labels and title
        if r == 1
            ylabel('dF/F');
        end
        title(sprintf('T%d (n=%d)', r, num_mice), 'FontSize', 10);
        grid on;
        
        hold off;
    end
    
    % Plot averaged away runs (bottom row)
    for r = 1:max_away
        runs_for_this_number = away_runs_by_number{r};
        
        if isempty(runs_for_this_number)
            continue;
        end
        
        % Create subplot (bottom row)
        subplot(2, max_runs_to_plot, max_runs_to_plot + r);
        
        % Average the runs for this run number
        [avg_distance, avg_dff, sem_dff, num_mice] = average_runs(runs_for_this_number, 'away', options);
        
        if isempty(avg_distance)
            continue;
        end
        
        % Plot averaged run with SEM shading
        plot_with_sem(avg_distance, avg_dff, sem_dff, [0.8 0.4 0.2]);
        hold on;
        
        % Mark start and end points
        plot(avg_distance(1), avg_dff(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
        plot(avg_distance(end), avg_dff(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
        
        % Set axis limits
        if ~isempty(options.dist_lim)
            xlim(options.dist_lim);
        end
        
        if ~isempty(options.dff_lim)
            ylim(options.dff_lim);
        end
        
        % Labels and title
        xlabel('Distance (cm)');
        if r == 1
            ylabel('dF/F');
        end
        title(sprintf('A%d (n=%d)', r, num_mice), 'FontSize', 10);
        grid on;
        
        hold off;
    end
    
    % Add overall title
    sgtitle(sprintf('Session %d - Averaged Runs (%dT, %dA max runs)%s', ...
            session_number, max_towards, max_away, category_str), 'FontSize', 12);
    
    fprintf('Created averaged plot for session %d (%d towards, %d away max runs)%s\n', ...
            session_number, max_towards, max_away, category_str);
end

function [avg_distance, avg_dff, sem_dff, num_mice] = average_runs(runs, run_type, options)
    % Average multiple runs by using all original data points
    
    if isempty(runs)
        avg_distance = [];
        avg_dff = [];
        sem_dff = [];
        num_mice = 0;
        return;
    end
    
    num_mice = length(runs);
    
    % Sample points along the distance range
    % First determine the range from all runs
    all_distances_for_range = [];
    for i = 1:length(runs)
        run = runs{i};
        distance = run.distance;
        
        % Transform distance for towards runs
        if strcmp(run_type, 'towards')
            end_distance = distance(end);
            distance = -(distance - end_distance);
        end
        
        all_distances_for_range = [all_distances_for_range; distance];
    end
    
    min_dist = min(all_distances_for_range);
    max_dist = max(all_distances_for_range);
    sample_distances = linspace(min_dist, max_dist, options.distance_points);
    
    % For each sample distance, collect values from all runs
    avg_distance = sample_distances;
    avg_dff = zeros(size(sample_distances));
    sem_dff = zeros(size(sample_distances));
    
    for i = 1:length(sample_distances)
        target_dist = sample_distances(i);
        values_at_this_distance = [];
        
        % Collect values from each run at this distance
        for r = 1:length(runs)
            run = runs{r};
            distance = run.distance;
            dff = run.dff;
            
            % Transform distance for towards runs
            if strcmp(run_type, 'towards')
                end_distance = distance(end);
                distance = -(distance - end_distance);
            end
            
            % Normalize dF/F if requested
            if options.normalize_dff
                if std(dff) > 0
                    dff = (dff - mean(dff)) / std(dff);
                end
            end
            
            % Find closest point to target distance
            if ~isempty(distance)
                [~, closest_idx] = min(abs(distance - target_dist));
                values_at_this_distance(end+1) = dff(closest_idx);
            end
        end
        
        % Calculate mean and SEM
        if ~isempty(values_at_this_distance)
            avg_dff(i) = mean(values_at_this_distance);
            if length(values_at_this_distance) > 1
                sem_dff(i) = std(values_at_this_distance) / sqrt(length(values_at_this_distance));
            else
                sem_dff(i) = 0;
            end
        else
            avg_dff(i) = NaN;
            sem_dff(i) = NaN;
        end
    end
    
    % Remove NaN points
    valid_points = ~isnan(avg_dff);
    avg_distance = avg_distance(valid_points);
    avg_dff = avg_dff(valid_points);
    sem_dff = sem_dff(valid_points);
end

function plot_with_sem(x, y, sem, color)
    % Plot line with SEM shading
    hold on;
    
    % Create shaded area for SEM
    if length(x) > 1 && length(y) == length(sem) && length(y) == length(x)
        % Remove any NaN values
        valid = ~isnan(y) & ~isnan(sem);
        if any(valid)
            x_valid = x(valid);
            y_valid = y(valid);
            sem_valid = sem(valid);
            
            % Create patch for shaded area
            x_patch = [x_valid, fliplr(x_valid)];
            y_patch = [y_valid + sem_valid, fliplr(y_valid - sem_valid)];
            
            % Plot shaded area
            patch_color = color + (1-color)*0.7; % Lighter version of the line color
            fill(x_patch, y_patch, patch_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            
            % Plot main line
            plot(x_valid, y_valid, '-', 'Color', color, 'LineWidth', 2);
        end
    end
end