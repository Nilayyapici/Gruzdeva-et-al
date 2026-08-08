function plot_activity_difference_heatmaps_smoothM3_fixed(mice_all, condition, colormap_name, options)
    % PLOT_ACTIVITY_DIFFERENCE_HEATMAPS_SMOOTH_FIXED
    % Fixed version that averages per-mouse maps (like the non-smooth version)
    % instead of pooling all data points together
    %
    % Key fix: Compute heatmaps for each mouse separately, then average them
    % This gives equal weight to each mouse regardless of data quantity
    
    % Column definitions
    COL_TIME = 1;
    COL_X = 2;
    COL_Y = 3;
    COL_DFF = 6;
    COL_SPEED = 7;
    COL_DIST_FOOD = 9;
    COL_DOOR = 12;
    COL_GROOM = 13;
    
    % Set default options
    if nargin < 4
        options = struct();
    end
    
    if ~isfield(options, 'bin_size'), options.bin_size = 4; end
    if ~isfield(options, 'smooth_factor'), options.smooth_factor = 2.0; end
    if ~isfield(options, 'min_time_threshold'), options.min_time_threshold = 1; end
    if ~isfield(options, 'exclude_grooming'), options.exclude_grooming = false; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 0; end
    if ~isfield(options, 'interpolate_display'), options.interpolate_display = true; end
    if ~isfield(options, 'interp_factor'), options.interp_factor = 4; end
    
    % Set default colormap
    if nargin < 3 || isempty(colormap_name)
        colormap_name = 'RdBu';
    end
    
    % Session information
    sessions = {'sess0', 'sess1', 'sess2'};
    session_names = {'Before', 'Learning', 'Test'};
    
    % Get unique food arm locations
    condition_indices = strcmp(mice_all(:,2), condition) & strcmp(mice_all(:,5), 'ok_signal');
    if sum(condition_indices) == 0
        error('No data found for condition "%s" with good signal quality', condition);
    end
    
    food_arms = unique(mice_all(condition_indices, 3));
    fprintf('Found food arm locations for condition "%s": %s\n', condition, strjoin(food_arms, ', '));
    
    % Initialize structures to store per-mouse data
    mouse_data_store = struct();
    all_coordinates = [];
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        mouse_data_store.(food_arm_loc) = struct();
    end
    
    % Track statistics
    total_points_before = 0;
    total_points_after = 0;
    unique_mice_by_food_arm = containers.Map();
    for fa = 1:length(food_arms)
        unique_mice_by_food_arm(food_arms{fa}) = {};
    end
    
    % First pass: organize data by mouse and collect coordinates
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        mouse_condition = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        signal_quality = mice_all{i, 5};
        
        if ~strcmp(mouse_condition, condition) || strcmp(signal_quality, 'bad_signal')
            continue;
        end
        
        % Determine session
        sess_idx = 0;
        for s = 1:length(sessions)
            if contains(session_info, sessions{s})
                sess_idx = s;
                break;
            end
        end
        if sess_idx == 0, continue; end
        
        % Extract mouse ID
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
        
        % Filter data
        total_points_before = total_points_before + size(data, 1);
        
        valid_idx = true(size(data, 1), 1);
        if options.exclude_grooming && size(data, 2) >= COL_GROOM
            valid_idx = valid_idx & (data(:, COL_GROOM) == 0);
        end
        if options.speed_threshold > 0 && size(data, 2) >= COL_SPEED
            valid_idx = valid_idx & (data(:, COL_SPEED) >= options.speed_threshold);
        end
        
        data_filtered = data(valid_idx, :);
        total_points_after = total_points_after + size(data_filtered, 1);
        
        if ~isempty(data_filtered)
            all_coordinates = [all_coordinates; data_filtered(:, [COL_X, COL_Y])];
            
            % Store data organized by mouse and session
            if ~isfield(mouse_data_store.(food_arm), mouse_id)
                mouse_data_store.(food_arm).(mouse_id) = struct();
                for s = 1:3
                    sess_key = sprintf('sess%d', s);
                    mouse_data_store.(food_arm).(mouse_id).(sess_key) = struct();
                    mouse_data_store.(food_arm).(mouse_id).(sess_key).x = [];
                    mouse_data_store.(food_arm).(mouse_id).(sess_key).y = [];
                    mouse_data_store.(food_arm).(mouse_id).(sess_key).dff = [];
                end
            end
            
            % Z-score dF/F for this mouse in this session
            dff_values = data_filtered(:, COL_DFF);
            if ~isempty(dff_values) && std(dff_values) > 0
                dff_zscore = (dff_values - mean(dff_values)) / std(dff_values);
            else
                dff_zscore = dff_values;
            end
            
            % Store this mouse's session data
            sess_key = sprintf('sess%d', sess_idx);
            mouse_data_store.(food_arm).(mouse_id).(sess_key).x = data_filtered(:, COL_X);
            mouse_data_store.(food_arm).(mouse_id).(sess_key).y = data_filtered(:, COL_Y);
            mouse_data_store.(food_arm).(mouse_id).(sess_key).dff = dff_zscore;
        end
    end
    
    % Set coordinate ranges
    if ~isfield(options, 'x_range')
        options.x_range = [min(all_coordinates(:,1)) - 10, max(all_coordinates(:,1)) + 10];
    end
    if ~isfield(options, 'y_range')
        options.y_range = [min(all_coordinates(:,2)) - 10, max(all_coordinates(:,2)) + 10];
    end
    
    % Create spatial bins
    x_edges = options.x_range(1):options.bin_size:options.x_range(2);
    y_edges = options.y_range(1):options.bin_size:options.y_range(2);
    x_centers = x_edges(1:end-1) + options.bin_size/2;
    y_centers = y_edges(1:end-1) + options.bin_size/2;
    
    fprintf('\nSpatial grid: %d x %d bins (%.1f cm bins)\n', length(x_centers), length(y_centers), options.bin_size);
    fprintf('Smoothing factor: %.1f\n', options.smooth_factor);
    fprintf('Data retention: %.1f%%\n', total_points_after/total_points_before*100);
    
    % Now process each mouse separately to create per-mouse heatmaps
    per_mouse_difference_maps = struct();
    all_differences = [];
    mice_count_by_food_arm = containers.Map();
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        per_mouse_difference_maps.(food_arm_loc) = struct();
        per_mouse_difference_maps.(food_arm_loc).diff_learning_maps = {};
        per_mouse_difference_maps.(food_arm_loc).diff_test_maps = {};
        
        % Get list of mice for this food arm
        mouse_ids = fieldnames(mouse_data_store.(food_arm_loc));
        mice_count_by_food_arm(food_arm_loc) = length(mouse_ids);
        
        % Process each mouse
        for m = 1:length(mouse_ids)
            mouse_id = mouse_ids{m};
            mouse_sessions = mouse_data_store.(food_arm_loc).(mouse_id);
            
            % Create heatmaps for each session of this mouse
            session_maps = cell(3, 1);
            for sess_idx = 1:3
                sess_key = sprintf('sess%d', sess_idx);
                session_maps{sess_idx} = createSmoothedHeatmap(...
                    mouse_sessions.(sess_key).x, ...
                    mouse_sessions.(sess_key).y, ...
                    mouse_sessions.(sess_key).dff, ...
                    x_edges, y_edges, x_centers, y_centers, ...
                    options.smooth_factor);
            end
            
            % Compute difference maps for this mouse
            before_map = session_maps{1};
            learning_map = session_maps{2};
            test_map = session_maps{3};
            
            % Learning - Before
            if ~all(isnan(before_map(:))) && ~all(isnan(learning_map(:)))
                diff_learning = learning_map - before_map;
                per_mouse_difference_maps.(food_arm_loc).diff_learning_maps{end+1} = diff_learning;
                all_differences = [all_differences; diff_learning(:)];
            end
            
            % Test - Before
            if ~all(isnan(before_map(:))) && ~all(isnan(test_map(:)))
                diff_test = test_map - before_map;
                per_mouse_difference_maps.(food_arm_loc).diff_test_maps{end+1} = diff_test;
                all_differences = [all_differences; diff_test(:)];
            end
        end
    end
    
    % Average the difference maps across mice
    averaged_difference_maps = struct();
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        averaged_difference_maps.(food_arm_loc) = struct();
        
        % Average Learning - Before maps
        if ~isempty(per_mouse_difference_maps.(food_arm_loc).diff_learning_maps)
            map_stack = cat(3, per_mouse_difference_maps.(food_arm_loc).diff_learning_maps{:});
            averaged_difference_maps.(food_arm_loc).diff_learning = nanmean(map_stack, 3);
        else
            averaged_difference_maps.(food_arm_loc).diff_learning = nan(length(x_centers), length(y_centers));
        end
        
        % Average Test - Before maps
        if ~isempty(per_mouse_difference_maps.(food_arm_loc).diff_test_maps)
            map_stack = cat(3, per_mouse_difference_maps.(food_arm_loc).diff_test_maps{:});
            averaged_difference_maps.(food_arm_loc).diff_test = nanmean(map_stack, 3);
        else
            averaged_difference_maps.(food_arm_loc).diff_test = nan(length(x_centers), length(y_centers));
        end
    end
    
    % Calculate color limits
    if isfield(options, 'diff_clim') && ~isempty(options.diff_clim)
        diff_clim_max = max(abs(options.diff_clim));
        diff_clim_min = -diff_clim_max;
    else
        all_diff_clean = all_differences(~isnan(all_differences));
        if ~isempty(all_diff_clean)
            diff_clim_max = max(abs(prctile(all_diff_clean, [5, 95])));
            diff_clim_min = -diff_clim_max;
        else
            diff_clim_max = 1;
            diff_clim_min = -1;
        end
    end
    
    fprintf('Difference color limits: [%.2f, %.2f] z-scored dF/F\n', diff_clim_min, diff_clim_max);
    
    % Create plots for each food arm
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        
        figure('Name', sprintf('Smooth Activity Difference (Fixed) - %s - %s', condition, food_arm_loc), ...
               'Position', [100 + (fa-1)*50, 100 + (fa-1)*50, 1000, 400]);
        
        diff_maps = {averaged_difference_maps.(food_arm_loc).diff_learning, ...
                     averaged_difference_maps.(food_arm_loc).diff_test};
        diff_titles = {'Learning - Before', 'Test - Before'};
        
        for diff_idx = 1:2
            subplot(1, 2, diff_idx);
            
            diff_map = diff_maps{diff_idx};
            
            if all(isnan(diff_map(:)))
                title(sprintf('%s\n(No Data)', diff_titles{diff_idx}), 'FontSize', 14);
                axis off;
                continue;
            end
            
            % Optional: Interpolate to finer grid for ultra-smooth display
            if options.interpolate_display
                [X, Y] = meshgrid(x_centers, y_centers);
                x_fine = linspace(x_centers(1), x_centers(end), length(x_centers) * options.interp_factor);
                y_fine = linspace(y_centers(1), y_centers(end), length(y_centers) * options.interp_factor);
                [X_fine, Y_fine] = meshgrid(x_fine, y_fine);
                
                % Interpolate using cubic method
                diff_map_fine = interp2(X, Y, diff_map', X_fine, Y_fine, 'cubic');
                
                imagesc(x_fine, y_fine, diff_map_fine, 'AlphaData', ~isnan(diff_map_fine));
            else
                imagesc(x_centers, y_centers, diff_map', 'AlphaData', ~isnan(diff_map'));
            end
            
            axis xy equal tight;
            colormap(gca, colormap_name);
            set(gca, 'CLim', [diff_clim_min, diff_clim_max]);
            set(gca, 'Color', 'white');
            
            cb = colorbar;
            cb.Label.String = '\Delta(Z-scored \Delta F/F)';
            
            title(sprintf('%s\n(Δ Z-scored ΔF/F)', diff_titles{diff_idx}), 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('X Position', 'FontSize', 12);
            ylabel('Y Position', 'FontSize', 12);
            
            axis off;
            box off;
        end
        
        sgtitle(sprintf('Smooth Activity Difference (Equal Mouse Weighting) - %s - %s (n=%d mice)', ...
                       condition, food_arm_loc, mice_count_by_food_arm(food_arm_loc)), ...
                'FontSize', 16, 'FontWeight', 'bold');
    end
    
    % Summary
    fprintf('\n=== FIXED SMOOTH ACTIVITY DIFFERENCE ANALYSIS ===\n');
    fprintf('Condition: %s\n', condition);
    fprintf('Method: Per-mouse heatmaps averaged (equal mouse weighting)\n');
    fprintf('Bin size: %.1f cm\n', options.bin_size);
    fprintf('Smoothing factor: %.1f\n', options.smooth_factor);
    
    for fa = 1:length(food_arms)
        food_arm_loc = food_arms{fa};
        fprintf('\n%s (%d mice)\n', food_arm_loc, mice_count_by_food_arm(food_arm_loc));
        
        for diff_idx = 1:2
            if diff_idx == 1
                diff_map = averaged_difference_maps.(food_arm_loc).diff_learning;
                diff_name = 'Learning - Before';
            else
                diff_map = averaged_difference_maps.(food_arm_loc).diff_test;
                diff_name = 'Test - Before';
            end
            
            valid_vals = diff_map(~isnan(diff_map));
            if ~isempty(valid_vals)
                fprintf('  %s: Mean=%.3f, Max=%.3f, Min=%.3f\n', ...
                        diff_name, mean(valid_vals), max(valid_vals), min(valid_vals));
            end
        end
    end
end

function heatmap = createSmoothedHeatmap(x_data, y_data, dff_data, x_edges, y_edges, x_centers, y_centers, smooth_factor)
    % Create a smoothed heatmap for a single dataset
    
    if isempty(x_data)
        heatmap = nan(length(x_centers), length(y_centers));
        return;
    end
    
    % Create 2D histogram with weighted values
    dff_grid = zeros(length(x_centers), length(y_centers));
    count_grid = zeros(length(x_centers), length(y_centers));
    
    % Bin the data
    x_bin = discretize(x_data, x_edges);
    y_bin = discretize(y_data, y_edges);
    
    % Accumulate dF/F values
    for i = 1:length(x_data)
        if ~isnan(x_bin(i)) && ~isnan(y_bin(i))
            dff_grid(x_bin(i), y_bin(i)) = dff_grid(x_bin(i), y_bin(i)) + dff_data(i);
            count_grid(x_bin(i), y_bin(i)) = count_grid(x_bin(i), y_bin(i)) + 1;
        end
    end
    
    % Calculate mean dF/F per bin
    mean_dff = zeros(size(dff_grid));
    valid_bins = count_grid > 0;
    mean_dff(valid_bins) = dff_grid(valid_bins) ./ count_grid(valid_bins);
    mean_dff(~valid_bins) = NaN;
    
    % Apply Gaussian smoothing
    if smooth_factor > 0
        smooth_input = mean_dff;
        smooth_input(isnan(smooth_input)) = 0;
        
        smoothed = imgaussfilt(smooth_input, smooth_factor);
        mask_smooth = imgaussfilt(double(valid_bins), smooth_factor);
        
        smoothed(mask_smooth < 0.1) = NaN;
        heatmap = smoothed;
    else
        heatmap = mean_dff;
    end
end