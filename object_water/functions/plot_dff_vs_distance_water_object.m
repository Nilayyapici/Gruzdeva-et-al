function plot_dff_vs_distance_water_object(run_data, smooth_window, x_lim, y_lim)
    % Plot Z-scored dF/F as a function of distance to water/object for towards and away runs
    % with shaded SEM areas
    %
    % Inputs:
    %   run_data - Structure containing run information (output from analyze_mouse_runs_water_object)
    %   smooth_window - Window size for smoothing (default: 5)
    %   x_lim - Distance limits for x-axis [min max] (default: [-10 5])
    %   y_lim - Z-score dF/F limits for y-axis [min max] (default: [-1 3])
    
    % Default values
    if nargin < 2
        smooth_window = 5;
    end
    
    if nargin < 3
        x_lim = [-10 5]; % Default distance range to plot
    end
    
    if nargin < 4
        y_lim = [-1 3]; % Default z-score range
    end
    
    % Create figure
    figure('Position', [100, 100, 1200, 800]);
    
    % Define colors for runs
    towards_color = [0, 0, 0.8]; % Blue
    away_color = [0.8, 0, 0]; % red
    
    % Area names for labeling
    area_names = {'water', 'object'};
    
    % Z-score all dF/F data across all runs
    all_dff = [];
    
    % Collect all dF/F values from all runs
    for i = 1:length(run_data)
        for j = 1:length(run_data(i).runs)
            all_dff = [all_dff; run_data(i).runs(j).dff];
        end
    end
    
    % Calculate mean and std for z-scoring
    dff_mean = mean(all_dff, 'omitnan');
    dff_std = std(all_dff, 'omitnan');
    
    % Create dummy lines for legend
    dummy_axes = axes('Position', [0 0 0.001 0.001], 'Visible', 'off');
    towards_line = line([0 1], [0 1], 'Color', towards_color, 'LineWidth', 2, 'Parent', dummy_axes);
    away_line = line([0 1], [0 1], 'Color', away_color, 'LineWidth', 2, 'Parent', dummy_axes);
    
    % First create the top row - Before sessions (2 plots)
    for area_idx = 1:2
        area_name = area_names{area_idx};
        subplot(2, 2, area_idx);
        
        % Extract data for this area from before sessions
        before_sessions = find(strcmp({run_data.session}, 'before'));
        
        plot_dff_vs_distance_panel_zscore(run_data, before_sessions, area_name, ...
            towards_color, away_color, smooth_window, x_lim, y_lim, dff_mean, dff_std);
        
        if area_idx == 1
            title({'Before Session', 'Z-scored dF/F vs Distance to Water'});
        else
            title({'Before Session', 'Z-scored dF/F vs Distance to Object'});
        end
    end
    
    % Second row - water and object sessions
    subplot(2, 2, 3);
    % Water sessions, water area
    water_sessions = find(strcmp({run_data.session}, 'water'));
    plot_dff_vs_distance_panel_zscore(run_data, water_sessions, 'water', ...
        towards_color, away_color, smooth_window, x_lim, y_lim, dff_mean, dff_std);
    title({'Water Session', 'Z-scored dF/F vs Distance to Water'});
    
    subplot(2, 2, 4);
    % Object sessions, object area
    object_sessions = find(strcmp({run_data.session}, 'object'));
    plot_dff_vs_distance_panel_zscore(run_data, object_sessions, 'object', ...
        towards_color, away_color, smooth_window, x_lim, y_lim, dff_mean, dff_std);
    title({'Object Session', 'Z-scored dF/F vs Distance to Object'});
    
    % Add overall title
    sgtitle('Z-scored dF/F vs Distance to Water/Object', 'FontSize', 14);
    
    % Add legend using dummy lines with correct colors
    legend([towards_line, away_line], {'Towards Runs', 'Away Runs'}, ...
           'Position', [0.4, 0.05, 0.2, 0.02], 'Units', 'normalized', 'Orientation', 'horizontal');
end

function plot_dff_vs_distance_panel_zscore(run_data, session_indices, area_name, ...
    towards_color, away_color, smooth_window, x_lim, y_lim, dff_mean, dff_std)
    % Helper function to plot a single panel with z-scored data and SEM shading
    
    % Check if we have valid session indices
    if isempty(session_indices)
        text(0.5, 0.5, 'No data available', 'HorizontalAlignment', 'center');
        set(gca, 'XLim', x_lim, 'YLim', y_lim, 'Box', 'on');
        xlabel('Distance (cm)');
        ylabel('Z-scored dF/F');
        return;
    end
    
    % Number of bins for distance grouping
    num_bins = 100;
    
    % Find unique mice for averaging
    all_mouse_ids = {};
    for i = session_indices
        if i <= length(run_data)
            all_mouse_ids{end+1} = run_data(i).mouse_id;
        end
    end
    unique_mouse_ids = unique(all_mouse_ids);
    
    % Initialize containers for per-mouse binned data
    mouse_towards_binned = cell(length(unique_mouse_ids), num_bins);
    mouse_away_binned = cell(length(unique_mouse_ids), num_bins);
    
    % Initialize distance bins based on x_lim
    % For towards runs, we use the negative range of x_lim (reversed)
    % Note: we're adjusting the binning to ensure proper range coverage
    towards_dist_bins = linspace(0, -x_lim(1), num_bins);
    away_dist_bins = linspace(0, x_lim(2), num_bins);
    
    % Second pass: bin data for each mouse
    for mouse_idx = 1:length(unique_mouse_ids)
        mouse_id = unique_mouse_ids{mouse_idx};
        
        % Find all runs for this mouse and area
        mouse_towards_distances = [];
        mouse_towards_dff = [];
        mouse_away_distances = [];
        mouse_away_dff = [];
        
        for sess_idx = 1:length(session_indices)
            session_idx = session_indices(sess_idx);
            
            if session_idx > length(run_data) || ~strcmp(run_data(session_idx).mouse_id, mouse_id)
                continue;
            end
            
            % Get runs for this mouse/session
            runs = run_data(session_idx).runs;
            
            for run_idx = 1:length(runs)
                run = runs(run_idx);
                
                % Only use runs for the specified area
                if ~strcmp(run.area, area_name)
                    continue;
                end
                
                % Z-score dF/F values
                z_dff = (run.dff - dff_mean) / dff_std;
                
                % Collect data based on run type
                if strcmp(run.type, 'towards')
                    mouse_towards_distances = [mouse_towards_distances; run.distance];
                    mouse_towards_dff = [mouse_towards_dff; z_dff];
                elseif strcmp(run.type, 'away')
                    mouse_away_distances = [mouse_away_distances; run.distance];
                    mouse_away_dff = [mouse_away_dff; z_dff];
                end
            end
        end
        
        % Bin towards run data for this mouse
        if ~isempty(towards_dist_bins) && ~isempty(mouse_towards_distances)
            for bin_idx = 1:length(towards_dist_bins)-1
                bin_start = towards_dist_bins(bin_idx);
                bin_end = towards_dist_bins(bin_idx+1);
                bin_mask = mouse_towards_distances >= bin_start & mouse_towards_distances < bin_end;
                
                if any(bin_mask)
                    % Average dF/F in this bin for this mouse
                    mouse_towards_binned{mouse_idx, bin_idx} = mean(mouse_towards_dff(bin_mask), 'omitnan');
                else
                    mouse_towards_binned{mouse_idx, bin_idx} = NaN;
                end
            end
        end
        
        % Bin away run data for this mouse
        if ~isempty(away_dist_bins) && ~isempty(mouse_away_distances)
            for bin_idx = 1:length(away_dist_bins)-1
                bin_start = away_dist_bins(bin_idx);
                bin_end = away_dist_bins(bin_idx+1);
                bin_mask = mouse_away_distances >= bin_start & mouse_away_distances < bin_end;
                
                if any(bin_mask)
                    % Average dF/F in this bin for this mouse
                    mouse_away_binned{mouse_idx, bin_idx} = mean(mouse_away_dff(bin_mask), 'omitnan');
                else
                    mouse_away_binned{mouse_idx, bin_idx} = NaN;
                end
            end
        end
    end
    
    % Calculate mean and SEM across mice for each bin
    towards_mean = nan(num_bins-1, 1);
    towards_sem = nan(num_bins-1, 1);
    away_mean = nan(num_bins-1, 1);
    away_sem = nan(num_bins-1, 1);
    
    for bin_idx = 1:num_bins-1
        % Extract bin data across all mice for towards runs
        bin_data = cellfun(@(x) x, mouse_towards_binned(:, bin_idx), 'UniformOutput', false);
        bin_data = bin_data(~cellfun('isempty', bin_data));
        bin_data_array = zeros(length(bin_data), 1);
        
        for i = 1:length(bin_data)
            if ~isnan(bin_data{i})
                bin_data_array(i) = bin_data{i};
            else
                bin_data_array(i) = NaN;
            end
        end
        
        if ~isempty(bin_data_array) && any(~isnan(bin_data_array))
            towards_mean(bin_idx) = mean(bin_data_array, 'omitnan');
            towards_sem(bin_idx) = std(bin_data_array, 'omitnan') / sqrt(sum(~isnan(bin_data_array)));
        end
        
        % Extract bin data across all mice for away runs
        bin_data = cellfun(@(x) x, mouse_away_binned(:, bin_idx), 'UniformOutput', false);
        bin_data = bin_data(~cellfun('isempty', bin_data));
        bin_data_array = zeros(length(bin_data), 1);
        
        for i = 1:length(bin_data)
            if ~isnan(bin_data{i})
                bin_data_array(i) = bin_data{i};
            else
                bin_data_array(i) = NaN;
            end
        end
        
        if ~isempty(bin_data_array) && any(~isnan(bin_data_array))
            away_mean(bin_idx) = mean(bin_data_array, 'omitnan');
            away_sem(bin_idx) = std(bin_data_array, 'omitnan') / sqrt(sum(~isnan(bin_data_array)));
        end
    end
    
    % Smooth the data
    valid_indices = ~isnan(towards_mean);
    if sum(valid_indices) > smooth_window
        towards_mean(valid_indices) = smoothdata(towards_mean(valid_indices), 'gaussian', smooth_window);
        towards_sem(valid_indices) = smoothdata(towards_sem(valid_indices), 'gaussian', smooth_window);
    end
    
    valid_indices = ~isnan(away_mean);
    if sum(valid_indices) > smooth_window
        away_mean(valid_indices) = smoothdata(away_mean(valid_indices), 'gaussian', smooth_window);
        away_sem(valid_indices) = smoothdata(away_sem(valid_indices), 'gaussian', smooth_window);
    end
    
    % Plot the data
    hold on;
    
    % For towards runs (negate the x-axis values)
    if ~isempty(towards_dist_bins) && length(towards_dist_bins) > 1
        % Calculate the center points for the bins
        bin_centers = -(towards_dist_bins(1:end-1) + diff(towards_dist_bins)/2);
        
        % Only plot where we have valid data
        valid_idx = find(~isnan(towards_mean));
        
        if ~isempty(valid_idx)
            % Plot SEM shading
            for i = 1:length(valid_idx)-1
                % Use polygon sections for shading
                current_idx = valid_idx(i);
                next_idx = valid_idx(i+1);
                
                % Only add shading if indices are adjacent (to avoid crossing gaps)
                if next_idx == current_idx + 1
                    x = [bin_centers(current_idx), bin_centers(next_idx), bin_centers(next_idx), bin_centers(current_idx)];
                    y = [towards_mean(current_idx) - towards_sem(current_idx), ...
                         towards_mean(next_idx) - towards_sem(next_idx), ...
                         towards_mean(next_idx) + towards_sem(next_idx), ...
                         towards_mean(current_idx) + towards_sem(current_idx)];
                    
                    fill(x, y, towards_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                end
            end
            
            % Plot mean line on top
            plot(bin_centers(valid_idx), towards_mean(valid_idx), 'Color', towards_color, 'LineWidth', 2);
        end
    end
    
    % For away runs (use the actual x-axis values)
    if ~isempty(away_dist_bins) && length(away_dist_bins) > 1
        % Calculate the center points for the bins
        bin_centers = away_dist_bins(1:end-1) + diff(away_dist_bins)/2;
        
        % Only plot where we have valid data
        valid_idx = find(~isnan(away_mean));
        
        if ~isempty(valid_idx)
            % Plot SEM shading
            for i = 1:length(valid_idx)-1
                % Use polygon sections for shading
                current_idx = valid_idx(i);
                next_idx = valid_idx(i+1);
                
                % Only add shading if indices are adjacent (to avoid crossing gaps)
                if next_idx == current_idx + 1
                    x = [bin_centers(current_idx), bin_centers(next_idx), bin_centers(next_idx), bin_centers(current_idx)];
                    y = [away_mean(current_idx) - away_sem(current_idx), ...
                         away_mean(next_idx) - away_sem(next_idx), ...
                         away_mean(next_idx) + away_sem(next_idx), ...
                         away_mean(current_idx) + away_sem(current_idx)];
                    
                    fill(x, y, away_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                end
            end
            
            % Plot mean line on top
            plot(bin_centers(valid_idx), away_mean(valid_idx), 'Color', away_color, 'LineWidth', 2);
        end
    end
    
    % Set axis limits and labels
    set(gca, 'XLim', x_lim, 'YLim', y_lim, 'Box', 'on');
    xlabel('Distance (cm)');
    ylabel('Z-scored dF/F');
    grid on;
    
    % Add a vertical line at x=0 for reference
    xline(0, 'k--');
    yline(0, 'k--')
    
    hold off;
end