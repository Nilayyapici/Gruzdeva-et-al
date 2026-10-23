function plot_data = plot_runs_by_behavior(run_data, options)
    % PLOT_RUNS_BY_BEHAVIOR Plot dF/F vs distance for runs classified by behavior
    %
    % Inputs:
    %   run_data - Output from analyze_runs_with_behavior_classification
    %   options - Structure with optional parameters:
    %       .sessions - Cell array of sessions (default: all available)
    %       .ylim - Y-axis limits, default [-1.5 1.5]
    %       .xlim_towards - X-axis for towards, default [0 210]
    %       .xlim_away - X-axis for away, default [0 210]
    %       .smoothing - Smoothing window, default 5
    %       .plot_sem - Show SEM, default true
    %       .title - Custom title (optional)
    %       .save - Save the plotted trace data to Excel, default true
    %       .save_dir - Folder for the Excel file, default pwd
    %       .save_name - File name, default
    %                    'runs_by_behavior_<sessions>_<yyyymmdd>.xlsx'
    %
    % Output:
    %   plot_data - table with every plotted point (same as the Excel
    %               'traces' sheet)
    %
    % SAVED EXCEL FILE
    %   Sheet 'traces' (long format), one row per direction x behavior x
    %   session x distance bin:
    %     Direction, Behavior, Session, Distance_cm (bin center),
    %     Mean_raw, SEM_raw, N_frames  (before smoothing)
    %     Mean_smoothed, SEM_smoothed  (the plotted line and shading;
    %                                   NaN where the line was not drawn)
    %     Plotted (1 if the line was drawn)
    %   Sheet 'counts': number of runs and mice per direction x behavior x
    %     session.
    %   NOTE: SEM is computed over all z-scored frames in a bin (pooled
    %   across runs and mice), exactly as plotted -- it is NOT an SEM across
    %   mice, so it is not suitable for statistics between conditions.
    
    if nargin < 2
        options = struct();
    end
    
    % Find available sessions
    if ~isfield(options, 'sessions')
        all_sessions = [];
        for m = 1:length(run_data)
            if ~ismember(run_data(m).session, all_sessions)
                all_sessions = [all_sessions, run_data(m).session];
            end
        end
        session_numbers = sort(all_sessions);
        options.sessions = cell(1, length(session_numbers));
        for i = 1:length(session_numbers)
            options.sessions{i} = ['sess' num2str(session_numbers(i))];
        end
    else
        session_numbers = [];
        for i = 1:length(options.sessions)
            session_numbers = [session_numbers, str2double(options.sessions{i}(5:end))];
        end
    end
    
    % Set defaults
    if ~isfield(options, 'ylim'), options.ylim = [-1.5 1.5]; end
    if ~isfield(options, 'xlim_towards'), options.xlim_towards = [0 210]; end
    if ~isfield(options, 'xlim_away'), options.xlim_away = [0 210]; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'plot_sem'), options.plot_sem = true; end
    if ~isfield(options, 'save'), options.save = true; end
    if ~isfield(options, 'save_dir'), options.save_dir = pwd; end
    if ~isfield(options, 'save_name') || isempty(options.save_name)
        options.save_name = sprintf('runs_by_behavior_%s_%s.xlsx', ...
            strjoin(options.sessions, '-'), datestr(now, 'yyyymmdd'));
    end
    
    fprintf('Plotting runs by behavior...\n');
    fprintf('Sessions: ');
    for i = 1:length(options.sessions)
        fprintf('%s ', options.sessions{i});
    end
    fprintf('\n');
    
    % Create figure
    if isfield(options, 'title')
        fig_name = options.title;
    else
        fig_name = 'dF/F by Behavior Classification';
    end
    
    figure('Name', fig_name, 'Position', [100, 100, 1400, 600]);
    
    % Plot layout: 1 row x 2 columns (towards, away)
    direction_types = {'towards', 'away'};
    direction_labels = {'TOWARDS', 'AWAY'};
    
    tables = cell(1, length(direction_types));
    counts = cell(1, length(direction_types));
    for d = 1:length(direction_types)
        subplot(1, 2, d);
        [tables{d}, counts{d}] = plot_direction_with_behaviors(run_data, direction_types{d}, ...
                                     session_numbers, options);
        title(direction_labels{d}, 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Add overall title
    if isfield(options, 'title')
        sgtitle(options.title, 'FontSize', 16, 'FontWeight', 'bold');
    else
        sgtitle('dF/F vs Distance by Behavior', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    plot_data  = vertcat(tables{:});
    count_data = vertcat(counts{:});
    
    % Save plotted data
    if options.save
        if isempty(plot_data)
            fprintf('Nothing was plotted -- no data saved.\n');
            return;
        end
        try
            if ~exist(options.save_dir, 'dir'), mkdir(options.save_dir); end
            xl_path = fullfile(options.save_dir, options.save_name);
            writetable(plot_data, xl_path, 'Sheet', 'traces');
            if ~isempty(count_data)
                writetable(count_data, xl_path, 'Sheet', 'counts');
            end
            fprintf('Plotted trace data saved to: %s\n', xl_path);
        catch ME
            fprintf('Could not save the trace data: %s\n', ME.message);
        end
    end
end

function [T, C] = plot_direction_with_behaviors(run_data, direction, session_numbers, options)
    % Plot a single direction with eating and visit runs in different colors.
    % Returns the plotted data (T) and run / mouse counts (C).
    
    hold on;
    
    T = table();
    C = table();
    
    behavior_types = {'eating', 'visit'};
    
    % Process each behavior type
    plotted_any = false;
    for beh_idx = 1:length(behavior_types)
        behavior = behavior_types{beh_idx};
        
        % Collect all runs for this direction and behavior
        all_runs = [];
        for m = 1:length(run_data)
            for r = 1:length(run_data(m).runs)
                run = run_data(m).runs(r);
                if strcmp(run.type, direction) && strcmp(run.behavior, behavior)
                    % Add mouse_id and session info to run
                    run.mouse_id = run_data(m).mouse_id;
                    run.session_num = run_data(m).session;
                    all_runs = [all_runs; run];
                end
            end
        end
        
        if isempty(all_runs)
            continue;
        end
        
        % Find max distance
        max_distance = 0;
        for i = 1:length(all_runs)
            max_distance = max(max_distance, max(all_runs(i).distance));
        end
        
        % Create distance bins
        bin_width = 1;
        dist_bins = 0:bin_width:ceil(max_distance + 5);
        bin_centers = dist_bins(1:end-1) + bin_width/2;
        
        % Initialize data collection by session
        all_data = struct();
        for s = 1:length(session_numbers)
            sess = session_numbers(s);
            key = sprintf('sess%d', sess);
            all_data.(key) = cell(length(dist_bins)-1, 1);
            for i = 1:length(dist_bins)-1
                all_data.(key){i} = [];
            end
        end
        
        % Group runs by mouse and session for z-scoring
        for sess_idx = 1:length(session_numbers)
            sess = session_numbers(sess_idx);
            
            % Get runs for this session
            session_runs = all_runs([all_runs.session_num] == sess);
            
            if isempty(session_runs)
                continue;
            end
            
            % Group by mouse
            unique_mice = unique({session_runs.mouse_id});
            
            % run / mouse counts for the export
            C = [C; table({direction}, {behavior}, sess, numel(session_runs), ...
                numel(unique_mice), 'VariableNames', ...
                {'Direction', 'Behavior', 'Session', 'N_runs', 'N_mice'})]; %#ok<AGROW>
            
            for mouse_idx = 1:length(unique_mice)
                mouse_id = unique_mice{mouse_idx};
                mouse_runs = session_runs(strcmp({session_runs.mouse_id}, mouse_id));
                
                % Collect all dF/F for this mouse
                all_dff = [];
                for r = 1:length(mouse_runs)
                    all_dff = [all_dff; mouse_runs(r).dff];
                end
                
                % Z-score parameters
                dff_mean = mean(all_dff);
                dff_std = std(all_dff);
                
                if dff_std < 1e-10
                    continue;
                end
                
                % Z-score and bin each run
                key = sprintf('sess%d', sess);
                for r = 1:length(mouse_runs)
                    run = mouse_runs(r);
                    z_dff = (run.dff - dff_mean) / dff_std;
                    
                    % Bin the data
                    for i = 1:length(dist_bins)-1
                        indices = run.distance >= dist_bins(i) & run.distance < dist_bins(i+1);
                        if any(indices)
                            all_data.(key){i} = [all_data.(key){i}; z_dff(indices)];
                        end
                    end
                end
            end
        end
        
        % Plot each session for this behavior
        for sess_idx = 1:length(session_numbers)
            sess = session_numbers(sess_idx);
            key = sprintf('sess%d', sess);
            
            if ~isfield(all_data, key)
                continue;
            end
            
            % Calculate mean and SEM
            means = nan(length(bin_centers), 1);
            sems = nan(length(bin_centers), 1);
            n_frames = zeros(length(bin_centers), 1);
            
            for i = 1:length(dist_bins)-1
                bin_data = all_data.(key){i};
                n_frames(i) = numel(bin_data);
                if ~isempty(bin_data)
                    means(i) = mean(bin_data);
                    sems(i) = std(bin_data) / sqrt(length(bin_data));
                end
            end
            
            valid = ~isnan(means);
            if ~any(valid)
                continue;
            end
            x_valid   = bin_centers(valid);
            y_valid   = means(valid);
            sem_valid = sems(valid);
            y_smoothed   = NaN(size(y_valid));
            sem_smoothed = NaN(size(sem_valid));
            is_plotted   = false;
            
            % Apply smoothing
            if sum(valid) > options.smoothing
                y_smoothed = movmean(y_valid, options.smoothing);
                sem_smoothed = movmean(sem_valid, options.smoothing);
                
                % Define dark color gradients for each behavior across sessions
                if strcmp(behavior, 'eating')
                    session_colors_eating = {
                        [0.6, 0.0, 0.0],  % Sess 0: Dark maroon
                        [0.8, 0.1, 0.1],  % Sess 1: Dark red
                        [1.0, 0.2, 0.2],  % Sess 2: Red
                        [0.9, 0.3, 0.1]   % Sess 3: Red-orange
                    };
                    session_color = session_colors_eating{min(sess_idx, 4)};
                else % visit
                    session_colors_visit = {
                        [0.0, 0.0, 0.6],  % Sess 0: Dark navy
                        [0.1, 0.2, 0.8],  % Sess 1: Dark blue
                        [0.2, 0.4, 1.0],  % Sess 2: Blue
                        [0.1, 0.5, 0.9]   % Sess 3: Sky blue
                    };
                    session_color = session_colors_visit{min(sess_idx, 4)};
                end
                
                % Plot SEM shading
                if options.plot_sem
                    for i = 1:length(x_valid)-1
                        x_patch = [x_valid(i), x_valid(i+1), x_valid(i+1), x_valid(i)];
                        y_patch = [y_smoothed(i) - sem_smoothed(i), ...
                                  y_smoothed(i+1) - sem_smoothed(i+1), ...
                                  y_smoothed(i+1) + sem_smoothed(i+1), ...
                                  y_smoothed(i) + sem_smoothed(i)];
                        
                        patch(x_patch, y_patch, session_color, ...
                            'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
                    end
                end
                
                % Plot line with dark, distinct colors
                behavior_label = behavior;
                behavior_label(1) = upper(behavior_label(1));
                plot(x_valid, y_smoothed, 'Color', session_color, ...
                    'LineWidth', 2.5, 'DisplayName', sprintf('%s - Sess%d', behavior_label, sess));
                plotted_any = true;
                is_plotted  = true;
            end
            
            % ---- collect the plotted data for export ----
            nv = numel(x_valid);
            T = [T; table(repmat({direction}, nv, 1), repmat({behavior}, nv, 1), ...
                repmat(sess, nv, 1), x_valid(:), y_valid(:), sem_valid(:), ...
                n_frames(valid), y_smoothed(:), sem_smoothed(:), ...
                repmat(double(is_plotted), nv, 1), ...
                'VariableNames', {'Direction', 'Behavior', 'Session', 'Distance_cm', ...
                'Mean_raw', 'SEM_raw', 'N_frames', 'Mean_smoothed', 'SEM_smoothed', ...
                'Plotted'})]; %#ok<AGROW>
        end
    end
    
    % Format plot
    xlabel('Distance from Food (cm)', 'FontSize', 12);
    ylabel('Z-scored dF/F', 'FontSize', 12);
    ylim(options.ylim);
    
    if strcmp(direction, 'towards')
        xlim(options.xlim_towards);
        set(gca, 'XDir', 'reverse');
        xticks_current = get(gca, 'XTick');
        set(gca, 'XTickLabel', arrayfun(@(x) sprintf('-%g', x), xticks_current, 'UniformOutput', false));
    else
        xlim(options.xlim_away);
    end
    
    % Add horizontal line at y=0
    line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    if plotted_any
        legend('show', 'Location', 'best');
        legend('boxoff');
    else
        text(0.5, 0.5, sprintf('No %s runs with eating or visit', direction), ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10);
    end
    
    grid off;
    box off;
    hold off;
end