function plot_spatial_difference_heatmaps(mice_all, condition, colormap_name, options)
    % PLOT_SPATIAL_DIFFERENCE_HEATMAPS Creates spatial difference heatmaps (Learning-Before, Test-Before)
    %
    % Inputs:
    %   mice_all - Cell array containing mouse data (organized as described)
    %   condition - String specifying the condition to analyze (e.g., 'control', 'CNO', 'chrimson')
    %   colormap_name - String specifying the colormap (e.g., 'RdBu', 'coolwarm', 'seismic')
    %   options - Structure with optional parameters
    %
    % Options:
    %   .dist_lim - Distance threshold for food zone detection (default: 30)
    %   .bin_size - Spatial bin size for heatmap (default: 4)
    %   .x_range - X coordinate range [min max] (default: auto-detected)
    %   .y_range - Y coordinate range [min max] (default: auto-detected)
    %   .smooth_factor - Gaussian smoothing factor (default: 1.5)
    %   .min_time_threshold - Minimum time spent in arm to include (default: 1 second)
    %   .diff_clim - Custom difference color limits [min max] (default: auto-calculated)
    %   .exclude_grooming - Exclude data points during grooming behavior (default: false)
    
    % Set default options
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 0; end
    if ~isfield(options, 'bin_size'), options.bin_size = 4; end
    if ~isfield(options, 'smooth_factor'), options.smooth_factor = 1.5; end
    if ~isfield(options, 'min_time_threshold'), options.min_time_threshold = 1; end
    if ~isfield(options, 'exclude_grooming'), options.exclude_grooming = false; end
    
    % Set default colormap for difference plots
    if nargin < 3 || isempty(colormap_name)
        colormap_name = 'jet'; % Use standard MATLAB colormap as default
    end
    
    % Session information
    sessions = {'sess0', 'sess1', 'sess2'};
    session_names = {'Before', 'Learning', 'Test'};
    
    % Get unique food arm locations for this condition
    condition_indices = strcmp(mice_all(:,2), condition) & strcmp(mice_all(:,5), 'ok_signal');
    if sum(condition_indices) == 0
        error('No data found for condition "%s" with good signal quality', condition);
    end
    
    food_arms = unique(mice_all(condition_indices, 3));
    fprintf('Found food arm locations for condition "%s": %s\n', condition, strjoin(food_arms, ', '));
    
    % Initialize data collection structures for each food arm location
    arm_data = struct();
    all_coordinates = [];
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        arm_data.(food_arm_loc) = struct();
        arm_data.(food_arm_loc).all_data = cell(3, 1);
    end
    
    % First pass: collect all coordinates and determine ranges
    mice_count_by_food_arm = containers.Map();
    
    for fa = 1:length(food_arms)
        mice_count_by_food_arm(food_arms{fa}) = 0;
    end
    
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        mouse_condition = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        signal_quality = mice_all{i, 5};
        
        % Skip if not the desired condition or bad signal
        if ~strcmp(mouse_condition, condition) || strcmp(signal_quality, 'bad_signal')
            continue;
        end
        
        % Apply grooming filter if requested
        if options.exclude_grooming && size(data, 2) >= 13
            non_grooming_indices = data(:, 13) == 0;
            data = data(non_grooming_indices, :);
        end
        
        % Collect all x,y coordinates for range determination
        all_coordinates = [all_coordinates; data(:, 2:3)];
    end
    
    % Set coordinate ranges (ensure consistency across all plots)
    if ~isfield(options, 'x_range')
        options.x_range = [min(all_coordinates(:,1)) - 10, max(all_coordinates(:,1)) + 10];
    end
    if ~isfield(options, 'y_range')
        options.y_range = [min(all_coordinates(:,2)) - 10, max(all_coordinates(:,2)) + 10];
    end
    
    % Create spatial bins (same for all plots)
    x_edges = options.x_range(1):options.bin_size:options.x_range(2);
    y_edges = options.y_range(1):options.bin_size:options.y_range(2);
    x_centers = x_edges(1:end-1) + options.bin_size/2;
    y_centers = y_edges(1:end-1) + options.bin_size/2;
    
    fprintf('Spatial grid: %d x %d bins\n', length(x_edges)-1, length(y_edges)-1);
    
    % Second pass: process data and collect coordinates
    unique_mice_by_food_arm = containers.Map();
    for fa = 1:length(food_arms)
        unique_mice_by_food_arm(food_arms{fa}) = {};
    end
    
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        mouse_condition = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        signal_quality = mice_all{i, 5};
        
        % Skip if not the desired condition or bad signal
        if ~strcmp(mouse_condition, condition) || strcmp(signal_quality, 'bad_signal')
            continue;
        end
        
        % Determine session number
        sess_idx = 0;
        for s = 1:length(sessions)
            if contains(session_info, sessions{s})
                sess_idx = s;
                break;
            end
        end
        
        if sess_idx == 0
            continue;
        end
        
        % Apply grooming filter if requested
        if options.exclude_grooming && size(data, 2) >= 13
            non_grooming_indices = data(:, 13) == 0;
            data = data(non_grooming_indices, :);
        end
        
        % Extract mouse ID for counting
        mouse_id = '';
        if contains(session_info, '_')
            parts = strsplit(session_info, '_');
            mouse_id = parts{1};
        else
            mouse_id = session_info;
        end
        
        % Track unique mice
        current_mice = unique_mice_by_food_arm(food_arm);
        if ~any(strcmp(current_mice, mouse_id))
            current_mice{end+1} = mouse_id;
            unique_mice_by_food_arm(food_arm) = current_mice;
        end
        
        % Store all data for this session
        sampling_rate = 10; % Assuming 10 Hz sampling
        min_points = options.min_time_threshold * sampling_rate;
        
        if size(data, 1) > min_points
            arm_data.(food_arm).all_data{sess_idx} = [arm_data.(food_arm).all_data{sess_idx}; data];
        end
    end
    
    % Update mice counts
    for fa = 1:length(food_arms)
        mice_count_by_food_arm(food_arms{fa}) = length(unique_mice_by_food_arm(food_arms{fa}));
    end
    
    % Pre-calculate all heatmaps for each session and compute differences
    all_differences = [];
    heatmap_data = struct();
    difference_data = struct();
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        heatmap_data.(food_arm_loc) = cell(3, 1);
        heatmap_data.([food_arm_loc '_mask']) = cell(3, 1);
        difference_data.(food_arm_loc) = cell(2, 1); % Learning-Before, Test-Before
        difference_data.([food_arm_loc '_mask']) = cell(2, 1); % Masks for difference plots
        
        % Calculate individual session heatmaps
        for sess_idx = 1:3
            current_data = arm_data.(food_arm_loc).all_data{sess_idx};
            
            if ~isempty(current_data)
                % Create 2D histogram of all positions
                [counts, ~, ~] = histcounts2(current_data(:, 2), current_data(:, 3), ...
                                           x_edges, y_edges);
                
                % CREATE MASK: Track which bins have data (non-zero counts)
                data_mask = counts > 0;
                
                % Normalize to occupancy percentage
                counts = counts / sum(counts(:)) * 100;
                
                % Apply smoothing
                if options.smooth_factor > 0
                    counts = imgaussfilt(counts, options.smooth_factor);
                end
                
                % Store both the heatmap and the mask
                heatmap_data.(food_arm_loc){sess_idx} = counts;
                heatmap_data.([food_arm_loc '_mask']){sess_idx} = data_mask;
            else
                heatmap_data.(food_arm_loc){sess_idx} = zeros(length(x_centers), length(y_centers));
                heatmap_data.([food_arm_loc '_mask']){sess_idx} = false(length(x_centers), length(y_centers));
            end
        end
        
        % Calculate differences: Learning-Before and Test-Before
        if ~isempty(heatmap_data.(food_arm_loc){1}) % Before session exists
            % Learning - Before
            if ~isempty(heatmap_data.(food_arm_loc){2}) % Learning session exists
                diff_learning = heatmap_data.(food_arm_loc){2} - heatmap_data.(food_arm_loc){1};
                % Combine masks: only show where EITHER session had data
                combined_mask = heatmap_data.([food_arm_loc '_mask']){1} | heatmap_data.([food_arm_loc '_mask']){2};
                difference_data.(food_arm_loc){1} = diff_learning;
                difference_data.([food_arm_loc '_mask']){1} = combined_mask;
                all_differences = [all_differences; diff_learning(:)];
            end
            
            % Test - Before
            if ~isempty(heatmap_data.(food_arm_loc){3}) % Test session exists
                diff_test = heatmap_data.(food_arm_loc){3} - heatmap_data.(food_arm_loc){1};
                % Combine masks: only show where EITHER session had data
                combined_mask = heatmap_data.([food_arm_loc '_mask']){1} | heatmap_data.([food_arm_loc '_mask']){3};
                difference_data.(food_arm_loc){2} = diff_test;
                difference_data.([food_arm_loc '_mask']){2} = combined_mask;
                all_differences = [all_differences; diff_test(:)];
            end
        end
    end
    
    % Calculate global difference color limits (symmetric around 0)
    if isfield(options, 'diff_clim') && ~isempty(options.diff_clim)
        diff_clim_max = max(abs(options.diff_clim));
        diff_clim_min = -diff_clim_max;
        fprintf('Using custom difference color limits: [%.2f, %.2f]%% occupancy change\n', diff_clim_min, diff_clim_max);
    else
        if ~isempty(all_differences)
            diff_clim_max = max(abs(prctile(all_differences, [5, 95]))); % Use 5th-95th percentile range
            diff_clim_min = -diff_clim_max;
        else
            diff_clim_max = 1;
            diff_clim_min = -1;
        end
        fprintf('Auto-calculated difference color limits: [%.2f, %.2f]%% occupancy change\n', diff_clim_min, diff_clim_max);
    end
    
    % Create plots for each food arm location
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        
        % Add grooming status to figure name
        grooming_status = '';
        if options.exclude_grooming
            grooming_status = ' - No Grooming';
        end
        
        figure('Name', sprintf('Spatial Difference Heatmaps - %s - Food Arm: %s - %s%s', condition, food_arm_loc, colormap_name, grooming_status), ...
               'Position', [100 + (fa-1)*50, 100 + (fa-1)*50, 800, 400]);
        
        % Plot difference heatmaps
        diff_titles = {'Learning - Before', 'Test - Before'};
        
        for diff_idx = 1:2
            subplot(1, 2, diff_idx);
            
            diff_counts = difference_data.(food_arm_loc){diff_idx};
            
            if isempty(diff_counts) || all(diff_counts(:) == 0)
                title(sprintf('%s\n(No Data)', diff_titles{diff_idx}));
                set(gca, 'XLim', options.x_range, 'YLim', options.y_range);
                axis off;
                continue;
            end
            
            % Apply rotation and reflection for specific chrimson conditions
            if contains(condition, 'chrimson_control') || contains(condition, 'chrimson_opto') || contains(condition, 'chrimson_optotest')
                diff_display = flipud(fliplr(diff_counts));
                data_mask = flipud(fliplr(difference_data.([food_arm_loc '_mask']){diff_idx}));
            else
                diff_display = diff_counts;
                data_mask = difference_data.([food_arm_loc '_mask']){diff_idx};
            end
            
            % Create the heatmap with transparency for areas with no data
            imagesc(x_centers, y_centers, diff_display', 'AlphaData', data_mask');
            axis xy;
            colormap(gca, colormap_name);
            set(gca, 'CLim', [diff_clim_min, diff_clim_max]);
            set(gca, 'Color', 'white'); % Background is white where there's no data
            hold on;
            
            title(sprintf('%s\n(Δ%% Occupancy)', diff_titles{diff_idx}), 'FontSize', 14);
            xlabel('X Position');
            ylabel('Y Position');
            
            % Add colorbar
            cb = colorbar;
            cb.Label.String = 'Occupancy Change (%)';
            
            axis off;
            box off;
        end
        
        % Add overall title
        main_title = sprintf('Spatial Difference Analysis - %s - Food Arm: %s (n=%d mice)', ...
                           condition, food_arm_loc, mice_count_by_food_arm(food_arm_loc));
        if options.exclude_grooming
            main_title = [main_title, ' - Grooming Excluded'];
        end
        
        sgtitle(main_title, 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    % Print summary statistics
    fprintf('\n=== SPATIAL DIFFERENCE ANALYSIS SUMMARY ===\n');
    fprintf('Condition: %s\n', condition);
    fprintf('Analysis: Difference heatmaps (Learning-Before, Test-Before)\n');
    if options.exclude_grooming
        fprintf('Grooming exclusion: ENABLED\n');
    else
        fprintf('Grooming exclusion: DISABLED\n');
    end
    
    fprintf('Difference color limits: [%.2f, %.2f]%% occupancy change\n', diff_clim_min, diff_clim_max);
    fprintf('Coordinate ranges: X [%.1f, %.1f], Y [%.1f, %.1f]\n', ...
            options.x_range(1), options.x_range(2), options.y_range(1), options.y_range(2));
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        fprintf('\nFood Arm Location: %s (%d mice)\n', food_arm_loc, mice_count_by_food_arm(food_arm_loc));
        
        % Calculate and report change statistics
        diff_titles = {'Learning - Before', 'Test - Before'};
        for diff_idx = 1:2
            diff_name = diff_titles{diff_idx};
            diff_counts = difference_data.(food_arm_loc){diff_idx};
            
            if ~isempty(diff_counts)
                pos_change = sum(diff_counts(:) > 1); % Areas with >1% increase
                neg_change = sum(diff_counts(:) < -1); % Areas with >1% decrease
                total_bins = numel(diff_counts);
                
                fprintf('  %s: %d bins increased (%.1f%%), %d bins decreased (%.1f%%)\n', ...
                        diff_name, pos_change, pos_change/total_bins*100, ...
                        neg_change, neg_change/total_bins*100);
                
                % Report extreme changes
                max_increase = max(diff_counts(:));
                max_decrease = min(diff_counts(:));
                fprintf('    Max increase: %.2f%%, Max decrease: %.2f%%\n', max_increase, max_decrease);
            else
                fprintf('  %s: No data available\n', diff_name);
            end
        end
    end
end