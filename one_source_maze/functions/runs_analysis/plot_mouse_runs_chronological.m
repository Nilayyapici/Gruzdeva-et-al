function plot_mouse_runs_chronological(run_data, session_number, options)
    % PLOT_MOUSE_RUNS_CHRONOLOGICAL Plot runs chronologically in subplots
    % Creates one figure per mouse with runs in chronological order
    % Towards runs on top row, away runs on bottom row
    % Each run is plotted as dF/F (y-axis) vs distance (x-axis)
    %
    % For towards runs: distance is transformed so end of run = 0, start = negative
    %
    % Inputs:
    %   run_data - Structure from analyze_single_arm_runs_with_speed_check
    %   session_number - Which session to plot (required)
    %   options - Structure with optional fields:
    %       .normalize_dff - Whether to z-score normalize dF/F per run (default: false)
    %       .dist_lim - Distance axis limits [min max] (default: [0 200])
    %       .dff_lim - dF/F axis limits [min max] (default: automatic)
    %       .mouse_filter - Specific mouse to plot (optional, if empty plots all mice)
    
    % Check inputs
    if nargin < 2
        error('Session number is required. Usage: plot_mouse_runs_chronological(run_data, session_number, options)');
    end
    
    % Default options
    if nargin < 3
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'normalize_dff')
        options.normalize_dff = false;
    end
    if ~isfield(options, 'dist_lim')
        options.dist_lim = [0 200];
    end
    if ~isfield(options, 'dff_lim')
        options.dff_lim = [];
    end
    if ~isfield(options, 'mouse_filter')
        options.mouse_filter = '';
    end
    
    % Filter by session
    fprintf('Filtering for session %d...\n', session_number);
    session_data = run_data([run_data.session] == session_number);
    
    if isempty(session_data)
        fprintf('No data found for session %d\n', session_number);
        return;
    end
    
    % Filter by specific mouse if requested
    if ~isempty(options.mouse_filter)
        session_data = session_data(strcmp({session_data.mouse_id}, options.mouse_filter));
        if isempty(session_data)
            fprintf('No data found for mouse %s in session %d\n', options.mouse_filter, session_number);
            return;
        end
    end
    
    fprintf('Creating chronological plots for %d mice in session %d...\n', length(session_data), session_number);
    
    % Loop through each mouse
    for m = 1:length(session_data)
        mouse_info = session_data(m);
        mouse_id = mouse_info.mouse_id;
        runs = mouse_info.runs;
        
        if isempty(runs)
            fprintf('No runs found for mouse %s\n', mouse_id);
            continue;
        end
        
        % Sort runs by start time (chronological order)
        start_times = [runs.start_time];
        [~, sort_idx] = sort(start_times);
        runs = runs(sort_idx);
        
        % Separate towards and away runs while preserving order
        towards_runs = runs(strcmp({runs.type}, 'towards'));
        away_runs = runs(strcmp({runs.type}, 'away'));
        
        % Skip if no runs
        if isempty(towards_runs) && isempty(away_runs)
            fprintf('No towards or away runs found for mouse %s\n', mouse_id);
            continue;
        end
        
        % Determine subplot layout
        max_runs = max(length(towards_runs), length(away_runs));
        if max_runs == 0
            continue;
        end
        
        % Create figure for this mouse
        fig_title = sprintf('Mouse %s - Session %d - Chronological Runs', mouse_id, session_number);
        figure('Position', [100 + 50, 100 + 50, 200*max_runs, 600], ...
               'Name', fig_title);
        
        % Plot towards runs (top row)
        for r = 1:length(towards_runs)
            run = towards_runs(r);
            
            % Create subplot (top row)
            subplot(2, max_runs, r);
            
            % Get data
            distance = run.distance;
            dff = run.dff;
            
            % Transform distance for towards runs:
            % End becomes 0, start becomes -start_distance  
            end_distance = distance(end);
            distance = -(distance - end_distance);
            
            % Normalize dF/F if requested
            if options.normalize_dff
                if std(dff) > 0
                    dff = (dff - mean(dff)) / std(dff);
                end
            end
            
            % Plot this run
            plot(distance, dff, '-', 'Color', [0.2 0.6 0.8], 'LineWidth', 1.5);
            hold on;
            
            % Mark start and end points
            plot(distance(1), dff(1), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'green');
            plot(distance(end), dff(end), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
            
            % Set axis limits (adjusted for negative distances)
            if ~isempty(options.dist_lim)
                xlim([-options.dist_lim(2), 0]);
            end
            
            if ~isempty(options.dff_lim)
                ylim(options.dff_lim);
            end
            
            % Labels and title
            if r == 1
                ylabel('dF/F');
            end
            title(sprintf('T%d (%.1fs)', run.id, run.duration), 'FontSize', 10);
            grid on;
            
            hold off;
        end
        
        % Plot away runs (bottom row)
        for r = 1:length(away_runs)
            run = away_runs(r);
            
            % Create subplot (bottom row)
            subplot(2, max_runs, max_runs + r);
            
            % Get data
            distance = run.distance;
            dff = run.dff;
            
            % For away runs, keep distance as is (starting from food, going away)
            
            % Normalize dF/F if requested
            if options.normalize_dff
                if std(dff) > 0
                    dff = (dff - mean(dff)) / std(dff);
                end
            end
            
            % Plot this run
            plot(distance, dff, '-', 'Color', [0.8 0.4 0.2], 'LineWidth', 1.5);
            hold on;
            
            % Mark start and end points
            plot(distance(1), dff(1), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'green');
            plot(distance(end), dff(end), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
            
            % Set axis limits (normal positive distances)
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
            title(sprintf('A%d (%.1fs)', run.id, run.duration), 'FontSize', 10);
            grid on;
            
            hold off;
        end
        
        % Add overall title
        sgtitle(sprintf('%s - Session %d (%dT, %dA runs)', ...
                mouse_id, session_number, length(towards_runs), length(away_runs)), ...
                'FontSize', 12);
        
        fprintf('Created chronological plot for mouse %s (%d towards, %d away runs)\n', ...
                mouse_id, length(towards_runs), length(away_runs));
    end
    
    fprintf('Completed chronological plotting for session %d.\n', session_number);
end