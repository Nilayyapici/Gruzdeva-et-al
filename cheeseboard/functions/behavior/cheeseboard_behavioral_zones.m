function results = cheeseboard_behavioral_zones(mice, options)
% CHEESEBOARD_BEHAVIORAL_ANALYSIS - Analyze mouse behavior in cheeseboard maze
%
% Usage:
%   options.group = 'saline';           % Analyze saline group only
%   options.group = 'CNO';              % Analyze CNO group only  
%   options.group = {'saline', 'CNO'};  % Compare both groups
%   options.time_limit = 10;            % Analyze only first 10 minutes (optional)
%   options.plot_trajectories = true;   % Plot spatial trajectories for zone validation (optional)
%   
%   results = cheeseboard_behavioral_analysis(mice, options);
%
% Input:
%   mice - cell array with structure:
%     mice{i,1} - mouse name
%     mice{i,2} - [x,y] food position
%     mice{i,3} - pre-test data (Nx10 double)
%     mice{i,4} - test data (Nx10 double)
%     mice{i,5} - [x_center,y_center] maze center (calculated if empty)
%     mice{i,6} - latencies [pre;test]
%
%   options - structure with fields:
%     .group - 'saline', 'CNO', or {'saline', 'CNO'}
%     .time_limit - (optional) analyze only first X minutes of each session
%     .plot_trajectories - (optional) true/false to plot x,y trajectories colored by zone
%
% Data columns:
%   1-time, 2,3-x,y position, 4-465nm, 5-405nm, 6-dff, 7-speed, 
%   8-zones (0-outside, 1-food, 2-area2, 3-area3), 9-distance, 10-grooming

%% Input validation
if nargin < 2
    error('Both mice data and options are required');
end

if ~isfield(options, 'group')
    error('options.group must be specified');
end

% Check if time limit is specified
if isfield(options, 'time_limit') && ~isempty(options.time_limit)
    time_limit_min = options.time_limit;
    time_limit_sec = time_limit_min * 60; % Convert to seconds
    use_time_limit = true;
    fprintf('Using time limit: %.1f minutes (%.0f seconds)\n', time_limit_min, time_limit_sec);
else
    use_time_limit = false;
    fprintf('Analyzing full session duration\n');
end

% Check if trajectory plotting is requested
if isfield(options, 'plot_trajectories') && options.plot_trajectories
    plot_trajectories = true;
    fprintf('Will generate trajectory plots for zone validation\n');
else
    plot_trajectories = false;
end

%% Extract mouse information
mouse_names = mice(:,1);
saline_idx = contains(mouse_names, 'saline');
cno_idx = contains(mouse_names, 'CNO');

%% Determine which mice to analyze
if ischar(options.group)
    % Single group analysis
    if strcmpi(options.group, 'saline')
        selected_mice = find(saline_idx);
        group_names = {'Saline'};
        compare_groups = false;
    elseif strcmpi(options.group, 'CNO')
        selected_mice = find(cno_idx);
        group_names = {'CNO'};
        compare_groups = false;
    else
        error('Invalid group option. Use "saline", "CNO", or {"saline", "CNO"}');
    end
elseif iscell(options.group) && length(options.group) == 2
    % Group comparison
    if any(strcmpi(options.group, 'saline')) && any(strcmpi(options.group, 'CNO'))
        saline_selected = find(saline_idx);
        cno_selected = find(cno_idx);
        selected_mice = [saline_selected; cno_selected];
        group_names = {'Saline', 'CNO'};
        compare_groups = true;
    else
        error('For group comparison, use {"saline", "CNO"}');
    end
else
    error('Invalid group option. Use "saline", "CNO", or {"saline", "CNO"}');
end

n_mice = length(selected_mice);
fprintf('Analyzing %d mice\n', n_mice);

if compare_groups
    fprintf('  - Saline: %d mice\n', length(saline_selected));
    fprintf('  - CNO: %d mice\n', length(cno_selected));
end

%% Calculate maze center for mice that don't have it
fprintf('\n=== CALCULATING MAZE CENTERS ===\n');
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    
    % Check if maze center exists
    if isempty(mice{mouse_idx, 5}) || length(mice{mouse_idx, 5}) < 2
        fprintf('Calculating maze center for mouse: %s\n', mice{mouse_idx, 1});
        
        % Get all position data from both pre and test sessions
        pre_data = mice{mouse_idx, 3};
        test_data = mice{mouse_idx, 4};
        
        % Combine x and y positions from both sessions
        all_x = [pre_data(:, 2); test_data(:, 2)];
        all_y = [pre_data(:, 3); test_data(:, 3)];
        
        % Calculate center as midpoint of ranges
        x_center = min(all_x) + (max(all_x) - min(all_x))/2;
        y_center = min(all_y) + (max(all_y) - min(all_y))/2;
        
        % Store the calculated center
        mice{mouse_idx, 5} = [x_center, y_center];
        
        fprintf('  Center: [%.2f, %.2f]\n', x_center, y_center);
    else
        fprintf('Mouse %s already has maze center: [%.2f, %.2f]\n', ...
                mice{mouse_idx, 1}, mice{mouse_idx, 5}(1), mice{mouse_idx, 5}(2));
    end
end

%% TRAJECTORY PLOTTING (Zone Validation)
if plot_trajectories
    fprintf('\n=== PLOTTING TRAJECTORIES FOR ZONE VALIDATION ===\n');
    
    % Define zone colors
    zone_colors = [0.7 0.7 0.7;    % Zone 0 (Outside) - gray
                   1.0 0.2 0.2;     % Zone 1 (Food) - red
                   0.2 0.6 1.0;     % Zone 2 (Area 2) - blue
                   0.2 0.8 0.2];    % Zone 3 (Area 3) - green
    zone_names_plot = {'Outside', 'Food', 'Area 2', 'Area 3'};
    
    % Determine layout based on number of mice
    if n_mice <= 4
        n_rows = 2;
        n_cols = n_mice;
    elseif n_mice <= 8
        n_rows = 2;
        n_cols = ceil(n_mice/2);
    else
        n_rows = ceil(sqrt(n_mice*2));
        n_cols = ceil(n_mice*2/n_rows);
    end
    
    % Create figure for trajectories
    fig_traj = figure('Name', 'Trajectory Plots - Zone Validation', ...
                      'Position', [50, 50, 300*n_cols, 300*n_rows]);
    
    for i = 1:n_mice
        mouse_idx = selected_mice(i);
        mouse_name = mice{mouse_idx, 1};
        food_pos = mice{mouse_idx, 2};
        
        % Get pre and test data
        pre_data = mice{mouse_idx, 3};
        test_data = mice{mouse_idx, 4};
        
        % Apply time limit if specified
        if use_time_limit
            pre_time = pre_data(:, 1);
            test_time = test_data(:, 1);
            pre_start_time = min(pre_time);
            test_start_time = min(test_time);
            pre_time_mask = (pre_time - pre_start_time) <= time_limit_sec;
            test_time_mask = (test_time - test_start_time) <= time_limit_sec;
            pre_data_filtered = pre_data(pre_time_mask, :);
            test_data_filtered = test_data(test_time_mask, :);
        else
            pre_data_filtered = pre_data;
            test_data_filtered = test_data;
        end
        
        % Extract position and zone data
        pre_x = pre_data_filtered(:, 2);
        pre_y = pre_data_filtered(:, 3);
        pre_zones = pre_data_filtered(:, 8);
        
        test_x = test_data_filtered(:, 2);
        test_y = test_data_filtered(:, 3);
        test_zones = test_data_filtered(:, 8);
        
        % Plot Pre session
        subplot(n_rows, n_cols, i);
        hold on;
        
        % Plot trajectory points colored by zone
        for zone = 0:3
            zone_idx = pre_zones == zone;
            if any(zone_idx)
                scatter(pre_x(zone_idx), pre_y(zone_idx), 20, zone_colors(zone+1,:), 'filled', 'MarkerFaceAlpha', 0.6);
            end
        end
        
        % Plot food location
        scatter(food_pos(1), food_pos(2), 150, 'k', 'p', 'filled', 'MarkerEdgeColor', 'y', 'LineWidth', 2);
        
        % Add maze center
        maze_center = mice{mouse_idx, 5};
        scatter(maze_center(1), maze_center(2), 100, 'k', 'o', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        
        title(sprintf('%s - Pre', strrep(mouse_name, '_', ' ')), 'FontSize', 9, 'Interpreter', 'none');
        xlabel('X Position');
        ylabel('Y Position');
        axis equal;
        grid on;
        set(gca, 'GridAlpha', 0.3);
        
        % Add legend only for first plot
        if i == 1
            legend_handles = [];
            legend_labels = {};
            for zone = 0:3
                h = scatter(NaN, NaN, 20, zone_colors(zone+1,:), 'filled');
                legend_handles = [legend_handles, h];
                legend_labels{end+1} = sprintf('Zone %d: %s', zone, zone_names_plot{zone+1});
            end
            h_food = scatter(NaN, NaN, 150, 'k', 'p', 'filled', 'MarkerEdgeColor', 'y', 'LineWidth', 2);
            h_center = scatter(NaN, NaN, 100, 'k', 'o', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
            legend_handles = [legend_handles, h_food, h_center];
            legend_labels{end+1} = 'Food';
            legend_labels{end+1} = 'Maze Center';
            legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 7);
        end
        
        % Plot Test session
        subplot(n_rows, n_cols, i + n_cols);
        hold on;
        
        % Plot trajectory points colored by zone
        for zone = 0:3
            zone_idx = test_zones == zone;
            if any(zone_idx)
                scatter(test_x(zone_idx), test_y(zone_idx), 20, zone_colors(zone+1,:), 'filled', 'MarkerFaceAlpha', 0.6);
            end
        end
        
        % Plot food location
        scatter(food_pos(1), food_pos(2), 150, 'k', 'p', 'filled', 'MarkerEdgeColor', 'y', 'LineWidth', 2);
        
        % Add maze center
        scatter(maze_center(1), maze_center(2), 100, 'k', 'o', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
        
        title(sprintf('%s - Test', strrep(mouse_name, '_', ' ')), 'FontSize', 9, 'Interpreter', 'none');
        xlabel('X Position');
        ylabel('Y Position');
        axis equal;
        grid on;
        set(gca, 'GridAlpha', 0.3);
    end
    
    % Add overall title
    if use_time_limit
        sgtitle(sprintf('Spatial Trajectories Colored by Zone (First %.1f min)', time_limit_min), 'FontSize', 12, 'FontWeight', 'bold');
    else
        sgtitle('Spatial Trajectories Colored by Zone (Full Session)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    fprintf('Trajectory plots created. Review to validate zone definitions.\n');
end

%% 1. ZONE ANALYSIS
fprintf('\n=== ZONE ANALYSIS ===\n');

% Initialize storage
zone_percent_pre = zeros(n_mice, 4);  % zones 0-3
zone_percent_test = zeros(n_mice, 4);
zone_names = {'Outside', 'Food', 'Area 2', 'Area 3'};

% Calculate zone percentages for each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    fprintf('Processing mouse %d/%d: %s\n', i, n_mice, mice{mouse_idx,1});
    
    % Get pre and test data
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};
    
    % Apply time limit if specified
    if use_time_limit
        % Get time columns
        pre_time = pre_data(:, 1);
        test_time = test_data(:, 1);
        
        % Find indices within time limit (from start of session)
        pre_start_time = min(pre_time);
        test_start_time = min(test_time);
        
        pre_time_mask = (pre_time - pre_start_time) <= time_limit_sec;
        test_time_mask = (test_time - test_start_time) <= time_limit_sec;
        
        % Filter data
        pre_data_filtered = pre_data(pre_time_mask, :);
        test_data_filtered = test_data(test_time_mask, :);
        
        % Check if we have enough data
        if size(pre_data_filtered, 1) < 10
            warning('Mouse %s: Less than 10 data points in pre session within time limit', mice{mouse_idx,1});
        end
        if size(test_data_filtered, 1) < 10
            warning('Mouse %s: Less than 10 data points in test session within time limit', mice{mouse_idx,1});
        end
        
        fprintf('  Pre: %d/%d points (%.1f%% of session)\n', size(pre_data_filtered,1), size(pre_data,1), size(pre_data_filtered,1)/size(pre_data,1)*100);
        fprintf('  Test: %d/%d points (%.1f%% of session)\n', size(test_data_filtered,1), size(test_data,1), size(test_data_filtered,1)/size(test_data,1)*100);
    else
        % Use full data
        pre_data_filtered = pre_data;
        test_data_filtered = test_data;
    end
    
    % Extract zone information
    pre_zones = pre_data_filtered(:, 8);
    test_zones = test_data_filtered(:, 8);
    
    % Calculate percentage time in each zone (using actual time, not just data points)
    for zone = 0:3
        % Method 1: Using data points (assumes regular sampling)
        zone_percent_pre_points = sum(pre_zones == zone) / length(pre_zones) * 100;
        zone_percent_test_points = sum(test_zones == zone) / length(test_zones) * 100;
        
        % Method 2: Using actual time durations
        if length(pre_data_filtered) > 1 && length(test_data_filtered) > 1
            pre_times = pre_data_filtered(:, 1);
            test_times = test_data_filtered(:, 1);
            
            % Calculate total time duration
            total_pre_duration = max(pre_times) - min(pre_times);
            total_test_duration = max(test_times) - min(test_times);
            
            % Calculate time spent in this zone
            pre_zone_indices = find(pre_zones == zone);
            test_zone_indices = find(test_zones == zone);
            
            if ~isempty(pre_zone_indices) && total_pre_duration > 0
                % Estimate time in zone (assuming each point represents equal time intervals)
                zone_percent_pre(i, zone+1) = length(pre_zone_indices) / length(pre_zones) * 100;
            else
                zone_percent_pre(i, zone+1) = 0;
            end
            
            if ~isempty(test_zone_indices) && total_test_duration > 0
                % Estimate time in zone (assuming each point represents equal time intervals)
                zone_percent_test(i, zone+1) = length(test_zone_indices) / length(test_zones) * 100;
            else
                zone_percent_test(i, zone+1) = 0;
            end
        else
            zone_percent_pre(i, zone+1) = 0;
            zone_percent_test(i, zone+1) = 0;
        end
    end
    
    % Report what the analysis represents
    if i == 1  % Only print once
        if use_time_limit
            fprintf('  Zone percentages represent %% of time within first %.1f minutes\n', time_limit_min);
        else
            fprintf('  Zone percentages represent %% of time within full session\n');
        end
        fprintf('  (Assuming regular sampling intervals)\n');
    end
end

%% Plot zone analysis
if compare_groups
    % Group comparison plots - single row with individual points
    saline_mask = ismember(selected_mice, saline_selected);
    cno_mask = ismember(selected_mice, cno_selected);
    
    % Create comparison figure
    figure('Name', 'Zone Analysis: Saline vs CNO', 'Position', [100, 100, 1400, 400]);
    
    for zone = 1:4
        % Get data for each group
        saline_pre = zone_percent_pre(saline_mask, zone);
        saline_test = zone_percent_test(saline_mask, zone);
        cno_pre = zone_percent_pre(cno_mask, zone);
        cno_test = zone_percent_test(cno_mask, zone);
        
        % Plot for this zone
        subplot(1, 4, zone);
        
        % Prepare data for grouped bar plot
        group_means = [mean(saline_pre), mean(saline_test); 
                       mean(cno_pre), mean(cno_test)];
        group_errors = [std(saline_pre)/sqrt(sum(saline_mask)), std(saline_test)/sqrt(sum(saline_mask));
                        std(cno_pre)/sqrt(sum(cno_mask)), std(cno_test)/sqrt(sum(cno_mask))];
        
        % Create grouped bar plot
        b = bar(group_means);
        b(1).FaceColor = [0.7 0.7 0.7];  % Pre (light gray)
        b(2).FaceColor = [0.3 0.3 0.3];  % Test (dark gray)
        b(1).FaceAlpha = 0.7;
        b(2).FaceAlpha = 0.7;
        
        hold on;
        
        % Add error bars
        ngroups = size(group_means, 1);
        nbars = size(group_means, 2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for j = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
            errorbar(x, group_means(:,j), group_errors(:,j), 'k.', 'LineWidth', 1.5);
        end
        
        % Add individual mouse points
        % Saline group points
        if sum(saline_mask) > 0
            x_pre_saline = 1 - groupwidth/2 + groupwidth/(2*nbars);
            x_test_saline = 1 - groupwidth/2 + 3*groupwidth/(2*nbars);
            
            % Add some jitter for visibility
            jitter_pre = x_pre_saline + (rand(sum(saline_mask),1)-0.5)*0.1;
            jitter_test = x_test_saline + (rand(sum(saline_mask),1)-0.5)*0.1;
            
            % Use different colors for each mouse
            saline_mouse_names = mouse_names(selected_mice(saline_mask));
            colors_saline = lines(sum(saline_mask)); % Generate distinct colors
            
            for m = 1:sum(saline_mask)
                scatter(jitter_pre(m), saline_pre(m), 40, colors_saline(m,:), 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                scatter(jitter_test(m), saline_test(m), 40, colors_saline(m,:), 'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                plot([jitter_pre(m), jitter_test(m)], [saline_pre(m), saline_test(m)], '-', 'Color', colors_saline(m,:), 'LineWidth', 1);
            end
        end
        
        % CNO group points
        if sum(cno_mask) > 0
            x_pre_cno = 2 - groupwidth/2 + groupwidth/(2*nbars);
            x_test_cno = 2 - groupwidth/2 + 3*groupwidth/(2*nbars);
            
            % Add some jitter for visibility
            jitter_pre = x_pre_cno + (rand(sum(cno_mask),1)-0.5)*0.1;
            jitter_test = x_test_cno + (rand(sum(cno_mask),1)-0.5)*0.1;
            
            % Use different colors for each mouse (starting from where saline left off)
            cno_mouse_names = mouse_names(selected_mice(cno_mask));
            colors_cno = lines(sum(cno_mask) + sum(saline_mask)); % Generate colors
            colors_cno = colors_cno(sum(saline_mask)+1:end,:); % Take second half
            
            for m = 1:sum(cno_mask)
                scatter(jitter_pre(m), cno_pre(m), 40, colors_cno(m,:), 'filled', 's', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                scatter(jitter_test(m), cno_test(m), 40, colors_cno(m,:), 'filled', 's', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                plot([jitter_pre(m), jitter_test(m)], [cno_pre(m), cno_test(m)], '-', 'Color', colors_cno(m,:), 'LineWidth', 1);
            end
        end
        
        % Add legend for mouse identification (only for food zone)
        if zone == 2 && (sum(saline_mask) > 0 || sum(cno_mask) > 0)
            legend_entries = {};
            legend_handles = [];
            
            % Add saline mice to legend
            if sum(saline_mask) > 0
                for m = 1:sum(saline_mask)
                    short_name = strrep(saline_mouse_names{m}, '_saline', '');
                    h = plot(NaN, NaN, 'o-', 'Color', colors_saline(m,:), 'MarkerFaceColor', colors_saline(m,:), ...
                            'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerSize', 6);
                    legend_handles = [legend_handles, h];
                    legend_entries{end+1} = [short_name ' (S)'];
                end
            end
            
            % Add CNO mice to legend
            if sum(cno_mask) > 0
                for m = 1:sum(cno_mask)
                    short_name = strrep(cno_mouse_names{m}, '_CNO', '');
                    h = plot(NaN, NaN, 's-', 'Color', colors_cno(m,:), 'MarkerFaceColor', colors_cno(m,:), ...
                            'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerSize', 6);
                    legend_handles = [legend_handles, h];
                    legend_entries{end+1} = [short_name ' (C)'];
                end
            end
            
            % Create legend outside the plot
            legend(legend_handles, legend_entries, 'Location', 'eastoutside', 'FontSize', 8);
        end
        
        % Statistical tests
        [~, p_saline] = ttest(saline_pre, saline_test);
        [~, p_cno] = ttest(cno_pre, cno_test);
        [~, p_between_pre] = ttest2(saline_pre, cno_pre);
        [~, p_between_test] = ttest2(saline_test, cno_test);
        
        % Add significance indicators
        max_y = max([saline_pre; saline_test; cno_pre; cno_test]) * 1.1;
        
        % Within-group significance
        if p_saline < 0.05
            text(1, max_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'red', 'FontWeight', 'bold');
        end
        if p_cno < 0.05
            text(2, max_y, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'red', 'FontWeight', 'bold');
        end
        
        % Between-group significance (if significant for either timepoint)
        if p_between_pre < 0.05 || p_between_test < 0.05
            text(1.5, max_y*1.1, '†', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'blue', 'FontWeight', 'bold');
        end
        
        title(sprintf('%s Zone', zone_names{zone}), 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Time (%)', 'FontSize', 11);
        set(gca, 'XTickLabel', {'Saline', 'CNO'});
        
        if zone == 1
            legend([b(1), b(2)], {'Pre', 'Test'}, 'Location', 'northeast');
        end
        
        ylim([0, max_y*1.2]);
        grid on;
        set(gca, 'GridAlpha', 0.3);
        
        % Print statistics
        fprintf('%s Zone:\n', zone_names{zone});
        fprintf('  Saline Pre vs Test: p=%.3f\n', p_saline);
        fprintf('  CNO Pre vs Test: p=%.3f\n', p_cno);
        fprintf('  Between groups (Pre): p=%.3f\n', p_between_pre);
        fprintf('  Between groups (Test): p=%.3f\n', p_between_test);
    end
    
    % Add overall figure annotations
    if use_time_limit
        sgtitle(sprintf('Zone Analysis: Saline vs CNO Comparison (First %.1f min)', time_limit_min), 'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle('Zone Analysis: Saline vs CNO Comparison (Full Session)', 'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Add legend for significance
    annotation('textbox', [0.02, 0.95, 0.3, 0.05], 'String', '* = within-group difference, † = between-group difference', ...
               'FontSize', 10, 'EdgeColor', 'none', 'BackgroundColor', 'white');
    
else
    % Single group analysis
    if use_time_limit
        figure('Name', sprintf('Zone Analysis: %s Group (First %.1f min)', group_names{1}, time_limit_min), 'Position', [100, 100, 1200, 400]);
    else
        figure('Name', sprintf('Zone Analysis: %s Group (Full Session)', group_names{1}), 'Position', [100, 100, 1200, 400]);
    end
    
    for zone = 1:4
        subplot(1, 4, zone);
        
        pre_vals = zone_percent_pre(:, zone);
        test_vals = zone_percent_test(:, zone);
        
        % Plot individual mouse trajectories with different colors
        colors_all = lines(n_mice); % Generate distinct colors for each mouse
        
        for m = 1:n_mice
            plot([1, 2], [pre_vals(m), test_vals(m)], 'o-', 'Color', colors_all(m,:), ...
                'LineWidth', 1.5, 'MarkerFaceColor', colors_all(m,:), 'MarkerSize', 6, 'MarkerEdgeColor', 'k');
            hold on;
        end
        
        % Plot group means
        bar_data = [mean(pre_vals), mean(test_vals)];
        bar_err = [std(pre_vals)/sqrt(n_mice), std(test_vals)/sqrt(n_mice)];
        
        bar(bar_data, 'FaceAlpha', 0.7);
        errorbar(1:2, bar_data, bar_err, 'k.', 'LineWidth', 2);
        
        % Add legend for mouse identification (only for food zone)
        if zone == 2
            legend_entries = {};
            legend_handles = [];
            
            for m = 1:n_mice
                mouse_name = mouse_names{selected_mice(m)};
                short_name = strrep(strrep(mouse_name, '_saline', ''), '_CNO', '');
                h = plot(NaN, NaN, 'o-', 'Color', colors_all(m,:), 'MarkerFaceColor', colors_all(m,:), ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 6);
                legend_handles = [legend_handles, h];
                legend_entries{end+1} = short_name;
            end
            
            % Create legend outside the plot
            legend(legend_handles, legend_entries, 'Location', 'eastoutside', 'FontSize', 8);
        end
        
        % Statistical test
        [~, p_val] = ttest(pre_vals, test_vals);
        max_y = max(bar_data) + max(bar_err);
        
        if p_val < 0.05
            text(1.5, max_y*1.1, sprintf('p=%.3f*', p_val), 'HorizontalAlignment', 'center');
        else
            text(1.5, max_y*1.1, sprintf('p=%.3f', p_val), 'HorizontalAlignment', 'center');
        end
        
        title(sprintf('%s Zone', zone_names{zone}));
        ylabel('Time (%)');
        set(gca, 'XTick', 1:2, 'XTickLabel', {'Pre', 'Test'});
        ylim([0, max_y*1.3]);
        
        fprintf('%s Zone: Pre vs Test p=%.3f\n', zone_names{zone}, p_val);
    end
end

%% Store results
results = struct();
results.mice = mice;  % Return updated mice array with calculated centers
results.zone_percent_pre = zone_percent_pre;
results.zone_percent_test = zone_percent_test;
results.selected_mice = selected_mice;
results.mouse_names = mouse_names(selected_mice);
results.group_names = group_names;
results.compare_groups = compare_groups;
results.time_limit_used = use_time_limit;
results.plot_trajectories = plot_trajectories;
if use_time_limit
    results.time_limit_min = time_limit_min;
    results.time_limit_sec = time_limit_sec;
end

if compare_groups
    results.saline_indices = find(ismember(selected_mice, saline_selected));
    results.cno_indices = find(ismember(selected_mice, cno_selected));
end

%% Summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
if use_time_limit
    fprintf('Analysis period: First %.1f minutes of each session\n', time_limit_min);
else
    fprintf('Analysis period: Full session duration\n');
end

if compare_groups
    saline_mask = ismember(selected_mice, saline_selected);
    cno_mask = ismember(selected_mice, cno_selected);
    
    fprintf('\nSALINE GROUP (n=%d):\n', sum(saline_mask));
    for zone = 1:4
        fprintf('  %s Zone - Pre: %.1f±%.1f%%, Test: %.1f±%.1f%%\n', ...
                zone_names{zone}, ...
                mean(zone_percent_pre(saline_mask,zone)), std(zone_percent_pre(saline_mask,zone)), ...
                mean(zone_percent_test(saline_mask,zone)), std(zone_percent_test(saline_mask,zone)));
    end
    
    fprintf('\nCNO GROUP (n=%d):\n', sum(cno_mask));
    for zone = 1:4
        fprintf('  %s Zone - Pre: %.1f±%.1f%%, Test: %.1f±%.1f%%\n', ...
                zone_names{zone}, ...
                mean(zone_percent_pre(cno_mask,zone)), std(zone_percent_pre(cno_mask,zone)), ...
                mean(zone_percent_test(cno_mask,zone)), std(zone_percent_test(cno_mask,zone)));
    end
else
    fprintf('\n%s GROUP (n=%d):\n', upper(group_names{1}), n_mice);
    for zone = 1:4
        fprintf('  %s Zone - Pre: %.1f±%.1f%%, Test: %.1f±%.1f%%\n', ...
                zone_names{zone}, ...
                mean(zone_percent_pre(:,zone)), std(zone_percent_pre(:,zone)), ...
                mean(zone_percent_test(:,zone)), std(zone_percent_test(:,zone)));
    end
end

fprintf('\nAnalysis completed!\n');
fprintf('Updated mice array (with calculated centers) available in results.mice\n');

end