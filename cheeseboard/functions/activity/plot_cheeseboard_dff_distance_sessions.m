function plot_cheeseboard_dff_distance_sessions(run_data, options)
% PLOT_CHEESEBOARD_DFF_DISTANCE_SESSIONS - Plot z-scored dF/F vs distance for towards and away runs
%
% Creates a 1x2 subplot comparing pre and test sessions for towards (left) and away (right) runs
% Towards runs shown in blue colors, away runs in red colors
% X-axis is inverted for towards runs (approaching food = negative distance)
%
% Usage:
%   options = struct();
%   options.max_distance = 35; % Plot up to 35 cm from food
%   options.ylim = [-0.5 0.5]; % Y-axis limits for z-scored dF/F
%   options.smoothing = 10; % Smoothing window size
%   options.run_stats = true; % Include statistical tests
%   plot_cheeseboard_dff_distance_sessions(run_data, options);
%
% Inputs:
%   run_data - Structure from analyze_cheeseboard_runs
%   options - Structure with optional parameters:
%       .max_distance - Maximum distance to plot in cm (default: 15)
%       .ylim - Y-axis limits [min max] (default: [-2 2])
%       .smoothing - Window size for smoothing (default: 5)
%       .figure_width - Width of the figure in pixels (default: 1200)
%       .figure_height - Height of the figure in pixels (default: 600)
%       .plot_sem - Whether to plot SEM shading (default: true)
%       .run_stats - Whether to perform statistical tests (default: true)
%       .stat_alpha - Alpha level for significance (default: 0.05)
%       .min_mice_per_bin - Minimum mice per bin for stats (default: 3)
%       .bin_width - Width of distance bins in cm (default: 1)
%       .title - Main figure title (optional)

% Default options
if nargin < 2
    options = struct();
end

% Set default options if not provided
if ~isfield(options, 'max_distance'), options.max_distance = 15; end
if ~isfield(options, 'ylim'), options.ylim = [-2 2]; end
if ~isfield(options, 'smoothing'), options.smoothing = 5; end
if ~isfield(options, 'figure_width'), options.figure_width = 1200; end
if ~isfield(options, 'figure_height'), options.figure_height = 600; end
if ~isfield(options, 'plot_sem'), options.plot_sem = true; end
if ~isfield(options, 'run_stats'), options.run_stats = true; end
if ~isfield(options, 'stat_alpha'), options.stat_alpha = 0.05; end
if ~isfield(options, 'min_mice_per_bin'), options.min_mice_per_bin = 3; end
if ~isfield(options, 'bin_width'), options.bin_width = 1; end
if ~isfield(options, 'title'), options.title = 'Z-scored dF/F vs Distance - Cheeseboard Runs'; end

fprintf('=== CHEESEBOARD dF/F vs DISTANCE PLOT ===\n');
fprintf('Max distance: %.1f cm\n', options.max_distance);
fprintf('Smoothing window: %d\n', options.smoothing);
fprintf('Y-axis limits: [%.2f, %.2f]\n', options.ylim(1), options.ylim(2));

% Define sessions and run types
sessions = {'pre', 'test'};
run_types = {'towards', 'away'};

% Define colors for sessions and run types
% Towards: blue colors, Away: red colors
session_colors = struct();
session_colors.towards = struct('pre', [0.4 0.6 1.0], 'test', [0.1 0.3 0.8]); % Light blue, dark blue
session_colors.away = struct('pre', [1.0 0.5 0.5], 'test', [0.8 0.1 0.1]);     % Light red, dark red

% Create distance bins
max_distance = options.max_distance;
bin_width = options.bin_width;
dist_bins = 0:bin_width:max_distance;
bin_centers = dist_bins(1:end-1) + bin_width/2;

fprintf('Distance bins: %d bins from 0 to %.1f cm\n', length(bin_centers), max_distance);

% Initialize data structure
all_data = struct();
mouse_bin_data = struct(); % For statistical testing

for s = 1:length(sessions)
    for t = 1:length(run_types)
        session = sessions{s};
        run_type = run_types{t};
        key = sprintf('%s_%s', session, run_type);
        
        all_data.(key) = cell(length(bin_centers), 1);
        mouse_bin_data.(key) = cell(length(bin_centers), 1);
        
        for i = 1:length(bin_centers)
            all_data.(key){i} = [];
            mouse_bin_data.(key){i} = containers.Map();
        end
    end
end

% Process each mouse
processed_mice = 0;
for m = 1:length(run_data)
    mouse_id = run_data(m).mouse_id;
    mouse_session = run_data(m).session;
    runs = run_data(m).runs;
    
    if isempty(runs)
        continue;
    end
    
    % Convert session number to session name
    session_name = '';
    
    if isnumeric(mouse_session)
        if any(mouse_session == 0)
            session_name = 'pre';
        elseif any(mouse_session == 1) || any(mouse_session == 2)
            session_name = 'test';
        end
    elseif ischar(mouse_session) || isstring(mouse_session)
        mouse_session_str = char(mouse_session);
        if strcmp(mouse_session_str, 'pre')
            session_name = 'pre';
        elseif strcmp(mouse_session_str, 'test')
            session_name = 'test';
        end
    end
    
    if isempty(session_name)
        continue; % Skip unknown sessions
    end
    
    % Collect all dF/F values for this mouse to calculate z-score parameters
    all_dff_values = [];
    for r = 1:length(runs)
        if isfield(runs(r), 'dff') && ~isempty(runs(r).dff)
            all_dff_values = [all_dff_values; runs(r).dff(:)];
        end
    end
    
    if isempty(all_dff_values)
        continue;
    end
    
    % Calculate mean and std for z-scoring
    dff_mean = mean(all_dff_values);
    dff_std = std(all_dff_values);
    
    % Skip mouse if std is zero or very small
    if dff_std < 1e-10
        warning('Mouse %s has near-zero dF/F standard deviation, skipping.', mouse_id);
        continue;
    end
    
    % Initialize temporary storage for this mouse's binned data
    mouse_temp_bins = struct();
    for t = 1:length(run_types)
        run_type = run_types{t};
        key = sprintf('%s_%s', session_name, run_type);
        mouse_temp_bins.(key) = cell(length(bin_centers), 1);
        for i = 1:length(bin_centers)
            mouse_temp_bins.(key){i} = [];
        end
    end
    
    % Process runs for this mouse
    for r = 1:length(runs)
        run = runs(r);
        
        % Check if this run type is recognized
        if ~ismember(run.type, run_types)
            continue;
        end
        
        % Apply distance limit
        valid_distance_idx = run.distance <= max_distance;
        if ~any(valid_distance_idx)
            continue;
        end
        
        % Filter the run data
        filtered_distance = run.distance(valid_distance_idx);
        filtered_dff = run.dff(valid_distance_idx);
        
        % Z-score this run's dF/F
        z_scored_dff = (filtered_dff - dff_mean) / dff_std;
        
        % Bin the z-scored data
        for i = 1:length(bin_centers)
            indices = filtered_distance >= dist_bins(i) & filtered_distance < dist_bins(i+1);
            if any(indices)
                key = sprintf('%s_%s', session_name, run.type);
                all_data.(key){i} = [all_data.(key){i}; z_scored_dff(indices)];
                mouse_temp_bins.(key){i} = [mouse_temp_bins.(key){i}; z_scored_dff(indices)];
            end
        end
    end
    
    % Calculate per-mouse averages for statistical testing
    for t = 1:length(run_types)
        run_type = run_types{t};
        key = sprintf('%s_%s', session_name, run_type);
        
        for i = 1:length(bin_centers)
            if ~isempty(mouse_temp_bins.(key){i})
                mouse_bin_avg = mean(mouse_temp_bins.(key){i});
                
                if mouse_bin_data.(key){i}.isKey(mouse_id)
                    existing_avg = mouse_bin_data.(key){i}(mouse_id);
                    mouse_bin_data.(key){i}(mouse_id) = (existing_avg + mouse_bin_avg) / 2;
                else
                    mouse_bin_data.(key){i}(mouse_id) = mouse_bin_avg;
                end
            end
        end
    end
    
    processed_mice = processed_mice + 1;
end

fprintf('Processed %d mice\n', processed_mice);

% Create figure
figure('Name', 'Cheeseboard dF/F vs Distance', ...
       'Position', [100, 100, options.figure_width, options.figure_height]);

% Plot towards and away runs
for t = 1:length(run_types)
    run_type = run_types{t};
    
    subplot(1, 2, t);
    hold on;
    
    % Plot both sessions
    for s = 1:length(sessions)
        session = sessions{s};
        key = sprintf('%s_%s', session, run_type);
        
        % Calculate mean and SEM for each distance bin
        means = nan(length(bin_centers), 1);
        sems = nan(length(bin_centers), 1);
        
        for i = 1:length(bin_centers)
            bin_data = all_data.(key){i};
            if ~isempty(bin_data)
                means(i) = mean(bin_data);
                sems(i) = std(bin_data) / sqrt(length(bin_data));
            end
        end
        
        % Apply smoothing
        valid = ~isnan(means);
        if sum(valid) > options.smoothing
            x_valid = bin_centers(valid);
            y_valid = means(valid);
            sem_valid = sems(valid);
            
            % For towards runs, make x-coordinates negative
            if strcmp(run_type, 'towards')
                x_valid = -x_valid;
            end
            
            % Apply moving average smoothing
            y_smoothed = movmean(y_valid, options.smoothing);
            sem_smoothed = movmean(sem_valid, options.smoothing);
            
            % Get session color
            line_color = session_colors.(run_type).(session);
            
            % Plot SEM shading if requested
            if options.plot_sem && length(x_valid) > 1 && length(y_smoothed) == length(sem_smoothed) && length(x_valid) == length(y_smoothed)
                % Ensure all vectors are the same size
                x_fill = [x_valid(:)', fliplr(x_valid(:)')];
                y_fill = [y_smoothed(:)' - sem_smoothed(:)', fliplr(y_smoothed(:)' + sem_smoothed(:)')];
                
                % Remove any NaN values that could cause fill to fail
                valid_fill = ~isnan(x_fill) & ~isnan(y_fill);
                if sum(valid_fill) > 2  % Need at least 3 points for a polygon
                    fill(x_fill(valid_fill), y_fill(valid_fill), line_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
            end
            
            % Plot the main line
            if strcmp(session, 'pre')
                session_label = 'Before';
            else
                session_label = 'Test';
            end
            plot(x_valid, y_smoothed, 'Color', line_color, 'LineWidth', 3, 'DisplayName', session_label);
            
            fprintf('  %s %s: %d valid bins\n', session, run_type, sum(valid));
        end
    end
    
    % Add statistical comparisons between sessions if requested
    if options.run_stats
        addSessionSignificanceBars(run_type, mouse_bin_data, bin_centers, options);
    end
    
    % Formatting
    if strcmp(run_type, 'towards')
        title('Towards Runs', 'FontWeight', 'bold', 'FontSize', 14, 'Color', [0.1 0.3 0.8]);
        xlabel('Distance to Food (cm)', 'FontSize', 12);
        xlim([-max_distance, 0]);
    else
        title('Away Runs', 'FontWeight', 'bold', 'FontSize', 14, 'Color', [0.8 0.1 0.1]);
        xlabel('Distance from Food (cm)', 'FontSize', 12);
        xlim([0, max_distance]);
    end
    
    ylabel('Z-scored dF/F', 'FontSize', 12);
    ylim(options.ylim);
    
    % Add zero lines
    plot(xlim, [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Add legend
    legend('show', 'Location', 'best');
    legend('boxoff');
    
    grid off;
    box off;
end

% Add overall title
sgtitle(options.title, 'FontSize', 16, 'FontWeight', 'bold');

%% Save trace data to xlsx
try
    xlsx_filename = sprintf('dff_distance_traces_%s.xlsx', strrep(lower(options.title), ' ', '_'));
    xlsx_filename = regexprep(xlsx_filename, '[^a-zA-Z0-9_.]', '');
    
    for t = 1:length(run_types)
        run_type = run_types{t};
        
        for s = 1:length(sessions)
            session = sessions{s};
            key = sprintf('%s_%s', session, run_type);
            
            % Recalculate means and SEMs (same logic as plotting)
            means = nan(length(bin_centers), 1);
            sems  = nan(length(bin_centers), 1);
            n_pts = nan(length(bin_centers), 1);
            
            for i = 1:length(bin_centers)
                bin_data = all_data.(key){i};
                if ~isempty(bin_data)
                    means(i) = mean(bin_data);
                    sems(i)  = std(bin_data) / sqrt(length(bin_data));
                    n_pts(i) = length(bin_data);
                end
            end
            
            % Apply smoothing (same as plot)
            valid = ~isnan(means);
            x_out    = nan(length(bin_centers), 1);
            y_out    = nan(length(bin_centers), 1);
            sem_out  = nan(length(bin_centers), 1);
            
            if sum(valid) > options.smoothing
                x_valid   = bin_centers(valid)';
                y_valid   = means(valid);
                sem_valid = sems(valid);
                
                if strcmp(run_type, 'towards')
                    x_valid = -x_valid;
                end
                
                y_smoothed   = movmean(y_valid,   options.smoothing);
                sem_smoothed = movmean(sem_valid, options.smoothing);
                
                x_out(valid)   = x_valid;
                y_out(valid)   = y_smoothed;
                sem_out(valid) = sem_smoothed;
            end
            
            % Per-mouse bin averages for this session/run_type
            mouse_ids_list = {};
            mouse_means_mat = nan(length(bin_centers), 0);
            
            for i = 1:length(bin_centers)
                m_map = mouse_bin_data.(key){i};
                if m_map.Count > 0
                    ids = keys(m_map);
                    for k = 1:length(ids)
                        if ~ismember(ids{k}, mouse_ids_list)
                            mouse_ids_list{end+1} = ids{k};
                        end
                    end
                end
            end
            
            if ~isempty(mouse_ids_list)
                mouse_means_mat = nan(length(bin_centers), length(mouse_ids_list));
                for i = 1:length(bin_centers)
                    m_map = mouse_bin_data.(key){i};
                    for k = 1:length(mouse_ids_list)
                        if m_map.isKey(mouse_ids_list{k})
                            mouse_means_mat(i, k) = m_map(mouse_ids_list{k});
                        end
                    end
                end
            end
            
            % Build table
            sheet_name = sprintf('%s_%s', session, run_type);  % e.g. pre_towards
            
            T = table(x_out, y_out, sem_out, n_pts, ...
                'VariableNames', {'Distance_cm', 'Mean_Zscored_dFF', 'SEM', 'N_datapoints'});
            
            if ~isempty(mouse_ids_list)
                mouse_table = array2table(mouse_means_mat, 'VariableNames', ...
                    cellfun(@(x) ['Mouse_' strrep(x,'-','_')], mouse_ids_list, 'UniformOutput', false));
                T = [T, mouse_table];
            end
            
            writetable(T, xlsx_filename, 'Sheet', sheet_name);
        end
    end
    
    fprintf('Trace data saved to: %s\n', xlsx_filename);
catch err
    warning('Could not save xlsx: %s', err.message);
end


fprintf('\nPlot completed!\n');
end

function addSessionSignificanceBars(run_type, mouse_bin_data, bin_centers, options)
    % Add significance bars comparing pre vs test sessions
    
    % Get keys for pre and test sessions
    pre_key = sprintf('pre_%s', run_type);
    test_key = sprintf('test_%s', run_type);
    
    % Check if both sessions have data
    if ~isfield(mouse_bin_data, pre_key) || ~isfield(mouse_bin_data, test_key)
        return;
    end
    
    % Get current y-limits
    yl = ylim;
    y_range = yl(2) - yl(1);
    bar_y = yl(2) - 0.05 * y_range;
    
    % Find significant bins
    sig_bins = [];
    sig_x_coords = [];
    
    for i = 1:length(bin_centers)
        % Get data for this bin
        pre_data_map = mouse_bin_data.(pre_key){i};
        test_data_map = mouse_bin_data.(test_key){i};
        
        % Check if we have enough data
        if pre_data_map.Count >= options.min_mice_per_bin && test_data_map.Count >= options.min_mice_per_bin
            pre_data = cell2mat(values(pre_data_map));
            test_data = cell2mat(values(test_data_map));
            
            % Perform t-test
            try
                [~, p_val] = ttest2(pre_data, test_data);
                
                if p_val < options.stat_alpha
                    sig_bins = [sig_bins, i];
                    
                    % Calculate x coordinate
                    x_coord = bin_centers(i);
                    if strcmp(run_type, 'towards')
                        x_coord = -x_coord;
                    end
                    sig_x_coords = [sig_x_coords, x_coord];
                end
            catch
                % Skip if t-test fails
                continue;
            end
        end
    end
    
    % Draw significance bars
    if ~isempty(sig_bins)
        % Group consecutive significant bins
        if length(sig_bins) == 1
            run_starts = 1;
            run_ends = 1;
        else
            diff_bins = diff(sig_bins);
            break_points = find(diff_bins > 1);
            
            if isempty(break_points)
                run_starts = 1;
                run_ends = length(sig_bins);
            else
                run_starts = [1, break_points + 1];
                run_ends = [break_points, length(sig_bins)];
            end
        end
        
        % Draw bars for consecutive runs
        for run_idx = 1:length(run_starts)
            start_idx = run_starts(run_idx);
            end_idx = run_ends(run_idx);
            
            if start_idx <= length(sig_x_coords) && end_idx <= length(sig_x_coords)
                x_start = sig_x_coords(start_idx);
                x_end = sig_x_coords(end_idx);
                
                % Extend bar slightly
                bin_width = abs(bin_centers(2) - bin_centers(1));
                if strcmp(run_type, 'towards')
                    bin_width = -bin_width;
                end
                x_start = x_start - bin_width/2;
                x_end = x_end + bin_width/2;
                
                % Draw significance bar
                plot([x_start, x_end], [bar_y, bar_y], 'k-', 'LineWidth', 3, 'HandleVisibility', 'off');
            end
        end
        
        fprintf('  %s: Before vs Test significant in %d bins\n', run_type, length(sig_bins));
        
        % Adjust y-limits to accommodate significance bar
        new_ylim = [yl(1), yl(2) + 0.1 * y_range];
        ylim(new_ylim);
    end
end