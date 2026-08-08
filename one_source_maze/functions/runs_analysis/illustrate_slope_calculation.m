function illustrate_slope_calculation(run_data, options)
    % illustrate_slope_calculation - Shows step-by-step how slopes are calculated
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_single_arm_runs_with_speed_check
    %   options - Structure with optional parameters:
    %       .example_mouse - Mouse ID to use for example (default: auto-select)
    %       .example_session - Session number (default: 1)
    %       .example_direction - 'towards' or 'away' (default: 'towards')
    %       .example_run - Run number within mouse/session/direction (default: 1)
    %       .bin_width - Width of distance bins (default: 1)
    %       .slope_range - Distance range for slope fitting (default: 50)
    %       .smoothing - Smoothing window size (default: 5)
    %       .max_distance - Maximum distance to analyze (default: 200)
    %       .browse_mode - Show multiple examples to choose from (default: false)
    %       .show_mouse_runs - Show all runs from specific mouse (default: false)
    %       .n_examples - Number of examples to show in browse mode (default: 9)

    if nargin < 2
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'example_mouse'), options.example_mouse = []; end
    if ~isfield(options, 'example_session'), options.example_session = 1; end
    if ~isfield(options, 'example_direction'), options.example_direction = 'towards'; end
    if ~isfield(options, 'example_run'), options.example_run = 1; end
    if ~isfield(options, 'bin_width'), options.bin_width = 1; end
    if ~isfield(options, 'slope_range'), options.slope_range = 50; end
    if ~isfield(options, 'smoothing'), options.smoothing = 5; end
    if ~isfield(options, 'max_distance'), options.max_distance = 200; end
    if ~isfield(options, 'browse_mode'), options.browse_mode = false; end
    if ~isfield(options, 'show_mouse_runs'), options.show_mouse_runs = false; end
    if ~isfield(options, 'n_examples'), options.n_examples = 9; end
    
    % Find suitable examples
    suitable_runs = find_suitable_examples(run_data, options);
    
    if isempty(suitable_runs)
        fprintf('No suitable runs found for illustration.\n');
        return;
    end
    
    % Choose mode
    if options.browse_mode
        show_multiple_examples(suitable_runs, options);
    elseif options.show_mouse_runs
        show_mouse_runs(suitable_runs, options);
    else
        show_single_example(suitable_runs, options);
    end
end

function suitable_runs = find_suitable_examples(run_data, options)
    % Find runs that are good for illustration
    
    suitable_runs = [];
    
    for m = 1:length(run_data)
        mouse_id = run_data(m).mouse_id;
        session = run_data(m).session;
        
        % Filter by session if specified
        if session ~= options.example_session
            continue;
        end
        
        runs = run_data(m).runs;
        
        for r = 1:length(runs)
            run = runs(r);
            
            % Filter by direction
            if ~strcmp(run.type, options.example_direction)
                continue;
            end
            
            % Filter by distance limit
            valid_distance_idx = run.distance <= options.max_distance;
            if sum(valid_distance_idx) < 20
                continue;
            end
            
            % Calculate data quality metrics
            filtered_distance = run.distance(valid_distance_idx);
            filtered_dff = run.dff(valid_distance_idx);
            
            % Check for reasonable data range
            distance_range = max(filtered_distance) - min(filtered_distance);
            if distance_range < options.slope_range / 2
                continue;
            end
            
            % Store suitable run info
            run_info = struct();
            run_info.mouse_id = mouse_id;
            run_info.session = session;
            run_info.run_idx = r;
            run_info.direction = run.type;
            run_info.distance = filtered_distance;
            run_info.dff = filtered_dff;
            run_info.time = run.time(valid_distance_idx);
            run_info.n_points = length(filtered_distance);
            run_info.distance_range = distance_range;
            run_info.dff_range = max(filtered_dff) - min(filtered_dff);
            
            suitable_runs = [suitable_runs; run_info];
        end
    end
    
    % Sort by data quality
    if ~isempty(suitable_runs)
        quality_scores = [suitable_runs.n_points] .* [suitable_runs.distance_range] .* [suitable_runs.dff_range];
        [~, sort_idx] = sort(quality_scores, 'descend');
        suitable_runs = suitable_runs(sort_idx);
    end
end

function show_multiple_examples(suitable_runs, options)
    % Show multiple examples in a grid
    
    n_examples = min(options.n_examples, length(suitable_runs));
    
    % Calculate grid layout
    n_cols = min(4, ceil(sqrt(n_examples)));
    n_rows = ceil(n_examples / n_cols);
    
    figure('Name', 'Choose Best Example for Slope Illustration', ...
           'Position', [50, 50, 350*n_cols, 300*n_rows]);
    
    for i = 1:n_examples
        subplot(n_rows, n_cols, i);
        
        run_info = suitable_runs(i);
        [slope_stats, bin_data] = calculate_slope_for_plot(run_info, options);
        
        % Apply coordinate transformation
        if strcmp(run_info.direction, 'towards')
            x_distance = -run_info.distance;
            x_centers = -bin_data.centers;
        else
            x_distance = run_info.distance;
            x_centers = bin_data.centers;
        end
        
        hold on;
        scatter(x_distance, run_info.dff, 10, [0.7, 0.7, 0.7], 'filled');
        plot(x_centers, bin_data.means, 'b.-', 'LineWidth', 2);
        
        % Plot fitted line only within the actual fitting range used for slope calculation
        if slope_stats.is_significant
            % Get the actual fitting range used in slope calculation
            if strcmp(run_info.direction, 'towards')
                range_min = -options.slope_range;
                range_max = 0;
                x_fit_range = x_centers(x_centers >= range_min & x_centers <= range_max);
            else
                range_min = 0;
                range_max = options.slope_range;
                x_fit_range = x_centers(x_centers >= range_min & x_centers <= range_max);
            end
            
            if ~isempty(x_fit_range)
                % Calculate the fitted line using the slope from the actual fitting
                if strcmp(run_info.direction, 'towards')
                    % For towards: slope was calculated on negative coordinates
                    y_fit = slope_stats.intercept + slope_stats.slope * x_fit_range;
                else
                    % For away: slope was calculated on positive coordinates  
                    y_fit = slope_stats.intercept + slope_stats.slope * x_fit_range;
                end
                plot(x_fit_range, y_fit, 'r-', 'LineWidth', 2);
            end
        end
        
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        
        title(sprintf('%s (n=%d)\nSlope=%.4f, p=%.3f', ...
              run_info.mouse_id, run_info.n_points, slope_stats.slope, slope_stats.p_value));
        xlabel('Distance');
        ylabel('dF/F');
        grid off;
    end
    
    sgtitle(sprintf('%s Direction Examples (Session %d)', ...
            options.example_direction, options.example_session));
    
    % Print mouse list
    fprintf('\n=== Available Examples ===\n');
    for i = 1:n_examples
        fprintf('%d. Mouse: %s, Run: %d\n', i, suitable_runs(i).mouse_id, suitable_runs(i).run_idx);
    end
    fprintf('Set options.example_mouse = ''MouseID'' to select\n');
end

function show_mouse_runs(suitable_runs, options)
    % Show all runs from specific mouse
    
    if isempty(options.example_mouse)
        fprintf('Error: show_mouse_runs requires options.example_mouse\n');
        return;
    end
    
    % Filter for this mouse
    mouse_indices = strcmp({suitable_runs.mouse_id}, options.example_mouse);
    mouse_runs = suitable_runs(mouse_indices);
    
    if isempty(mouse_runs)
        fprintf('No runs found for mouse %s\n', options.example_mouse);
        return;
    end
    
    n_runs = length(mouse_runs);
    n_cols = min(4, ceil(sqrt(n_runs)));
    n_rows = ceil(n_runs / n_cols);
    
    figure('Name', sprintf('All Runs for Mouse %s', options.example_mouse), ...
           'Position', [50, 50, 350*n_cols, 300*n_rows]);
    
    for i = 1:n_runs
        subplot(n_rows, n_cols, i);
        
        run_info = mouse_runs(i);
        [slope_stats, bin_data] = calculate_slope_for_plot(run_info, options);
        
        % Apply coordinate transformation
        if strcmp(run_info.direction, 'towards')
            x_distance = -run_info.distance;
            x_centers = -bin_data.centers;
        else
            x_distance = run_info.distance;
            x_centers = bin_data.centers;
        end
        
        hold on;
        scatter(x_distance, run_info.dff, 10, [0.7, 0.7, 0.7], 'filled');
        plot(x_centers, bin_data.means, 'b.-', 'LineWidth', 2);
        
        % Plot fitted line only within the actual fitting range used for slope calculation
        if slope_stats.is_significant
            % Get the actual fitting range used in slope calculation
            if strcmp(run_info.direction, 'towards')
                range_min = -options.slope_range;
                range_max = 0;
                x_fit_range = x_centers(x_centers >= range_min & x_centers <= range_max);
            else
                range_min = 0;
                range_max = options.slope_range;
                x_fit_range = x_centers(x_centers >= range_min & x_centers <= range_max);
            end
            
            if ~isempty(x_fit_range)
                % Calculate the fitted line using the slope from the actual fitting
                if strcmp(run_info.direction, 'towards')
                    % For towards: slope was calculated on negative coordinates
                    y_fit = slope_stats.intercept + slope_stats.slope * x_fit_range;
                else
                    % For away: slope was calculated on positive coordinates  
                    y_fit = slope_stats.intercept + slope_stats.slope * x_fit_range;
                end
                plot(x_fit_range, y_fit, 'r-', 'LineWidth', 2);
            end
        end
        
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
        
        title(sprintf('Run %d (n=%d)\nSlope=%.4f, p=%.3f', ...
              run_info.run_idx, run_info.n_points, slope_stats.slope, slope_stats.p_value));
        xlabel('Distance');
        ylabel('dF/F');
        grid off;
        
        % Highlight run number
        text(0.05, 0.95, sprintf('RUN %d', run_info.run_idx), ...
             'Units', 'normalized', 'BackgroundColor', 'yellow', 'FontWeight', 'bold');
    end
    
    sgtitle(sprintf('Mouse %s - %s Runs (Session %d)', ...
            options.example_mouse, options.example_direction, options.example_session));
    
    % Print run stats
    fprintf('\n=== Runs for Mouse %s ===\n', options.example_mouse);
    for i = 1:n_runs
        run_info = mouse_runs(i);
        [slope_stats, ~] = calculate_slope_for_plot(run_info, options);
        fprintf('Run %d: Slope=%.4f, p=%.3f', run_info.run_idx, slope_stats.slope, slope_stats.p_value);
        if slope_stats.is_significant
            fprintf(' *SIG*');
        end
        fprintf('\n');
    end
    fprintf('Set options.example_run = X for detailed view\n');
end

function show_single_example(suitable_runs, options)
    % Show detailed single example
    
    % Find the requested run
    if isempty(options.example_mouse)
        example_run = suitable_runs(1);
    else
        mouse_indices = strcmp({suitable_runs.mouse_id}, options.example_mouse);
        mouse_runs = suitable_runs(mouse_indices);
        
        if isempty(mouse_runs)
            fprintf('Mouse %s not found, using %s\n', options.example_mouse, suitable_runs(1).mouse_id);
            example_run = suitable_runs(1);
        else
            if options.example_run <= length(mouse_runs)
                example_run = mouse_runs(options.example_run);
            else
                fprintf('Run %d not found, using run 1\n', options.example_run);
                example_run = mouse_runs(1);
            end
        end
    end
    
    create_detailed_illustration(example_run, options);
end

function create_detailed_illustration(run_info, options)
    % Create detailed 6-panel illustration
    
    figure('Name', 'Slope Calculation Illustration', 'Position', [50, 50, 1400, 1000]);
    
    % Calculate all steps
    [slope_stats, bin_data, steps] = calculate_detailed_slope(run_info, options);
    
    % Apply coordinate transformation
    if strcmp(run_info.direction, 'towards')
        x_distance = -run_info.distance;
        x_label = 'Distance to Food (towards)';
        fit_color = [0.2, 0.4, 0.8]; % Blue for towards
    else
        x_distance = run_info.distance;
        x_label = 'Distance from Food (away)';
        fit_color = [0.8, 0.2, 0.2]; % Red for away
    end
    
    % Panel 1: Raw data
    subplot(2, 3, 1);
    scatter(x_distance, run_info.dff, 30, run_info.time, 'filled');
    try
        colormap(gca, 'parula');
    catch
        colormap(gca, 'jet');
    end
    cb1 = colorbar;
    cb1.Label.String = 'Time (s)';
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    title('Step 1: Raw Data');
    xlabel(x_label);
    ylabel('Raw dF/F');
    grid off;
    
    % Panel 2: Z-scored data
    subplot(2, 3, 2);
    scatter(x_distance, steps.z_scored, 30, run_info.time, 'filled');
    try
        colormap(gca, 'parula');
    catch
        colormap(gca, 'jet');
    end
    cb2 = colorbar;
    cb2.Label.String = 'Time (s)';
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    title('Step 2: Z-scored');
    xlabel(x_label);
    ylabel('Z-scored dF/F');
    grid off;
    
    % Panel 3: Binned data
    subplot(2, 3, 3);
    hold on;
    scatter(x_distance, steps.z_scored, 15, [0.7, 0.7, 0.7], 'filled', 'MarkerFaceAlpha', 0.3);
    
    if strcmp(run_info.direction, 'towards')
        x_bins = -bin_data.centers;
    else
        x_bins = bin_data.centers;
    end
    
    plot(x_bins, bin_data.means, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(x_bins, bin_data.means, bin_data.sems, 'b', 'LineStyle', 'none');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    title(sprintf('Step 3: Binned (%.2f units)', options.bin_width));
    xlabel(x_label);
    ylabel('Z-scored dF/F');
    grid off;
    
    % Panel 4: Smoothed data
    subplot(2, 3, 4);
    hold on;
    plot(x_bins, bin_data.means, 'b.-', 'Color', [0.5, 0.5, 1]);
    plot(steps.x_smooth, steps.y_smooth, 'r.-', 'LineWidth', 2);
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    title(sprintf('Step 4: Smoothed (window=%d)', options.smoothing));
    xlabel(x_label);
    ylabel('Z-scored dF/F');
    legend('Binned', 'Smoothed', 'Food', 'Location', 'best');
    legend('boxoff')
    grid off;
    
    % Panel 5: Fitted range
    subplot(2, 3, 5);
    hold on;
    plot(steps.x_smooth, steps.y_smooth, 'k.-', 'Color', [0.7, 0.7, 0.7]);
    plot(steps.x_fit, steps.y_fit, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    
    if ~isempty(steps.x_line)
        plot(steps.x_line, steps.y_line, 'b-', 'LineWidth', 3);
    end
    
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    title('Step 5: Slope Fitting');
    xlabel(x_label);
    ylabel('Z-scored dF/F');
    legend('All Data', 'Fit Range', 'Fitted Line', 'Food', 'Location', 'best');
    legend('boxoff')
    grid off;
    
    % Panel 6: Summary Plot with requested styling
    subplot(2, 3, 6);
    hold on;
    
    % Raw data in gray
    scatter(x_distance, steps.z_scored, 40, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Fit range data points in direction-specific color
    if ~isempty(steps.x_fit) && ~isempty(steps.y_fit)
        scatter(steps.x_fit, steps.y_fit, 30, fit_color, 'o', 'filled', 'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
    end
    
    % Fitted line in magenta
    if ~isempty(steps.x_line) && ~isempty(steps.y_line)
        plot(steps.x_line, steps.y_line, 'Color', [1, 0, 1], 'LineWidth', 4); % Magenta
    end
    % 
    % % Food location
    % line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    
    title(sprintf('SUMMARY: %s Run %d', run_info.mouse_id, run_info.run_idx));
    xlabel(x_label);
    ylabel('Z-scored dF/F');
    
    % Create legend
    if strcmp(run_info.direction, 'towards')
        legend('Raw Data', 'Fit Range (Towards)', 'Fitted Line', 'Food Location', 'Location', 'best');
    else
        legend('Raw Data', 'Fit Range (Away)', 'Fitted Line', 'Food Location', 'Location', 'best');
    end
    legend('boxoff')
    grid off;
    axis off
    
    % Add slope info as text
    if slope_stats.is_significant
        sig_text = sprintf('Slope: %.4f (p=%.3f) *SIG*', slope_stats.slope, slope_stats.p_value);
        text_color = [0, 0.7, 0];
    else
        sig_text = sprintf('Slope: %.4f (p=%.3f)', slope_stats.slope, slope_stats.p_value);
        text_color = [0.7, 0, 0];
    end
    
    text(0.05, 0.95, sig_text, 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Color', text_color, 'BackgroundColor', 'white', 'EdgeColor', 'none');
    
    sgtitle(sprintf('Slope Analysis: %s - %s Direction', run_info.mouse_id, run_info.direction), 'FontSize', 16);
end

function [slope_stats, bin_data] = calculate_slope_for_plot(run_info, options)
    % Simple slope calculation for plotting
    
    % Z-score
    z_scored = (run_info.dff - mean(run_info.dff)) / std(run_info.dff);
    
    % Bin data
    bin_edges = 0:options.bin_width:options.max_distance;
    bin_centers = bin_edges(1:end-1) + options.bin_width/2;
    
    bin_means = nan(length(bin_centers), 1);
    bin_sems = nan(length(bin_centers), 1);
    
    for i = 1:length(bin_centers)
        idx = run_info.distance >= bin_edges(i) & run_info.distance < bin_edges(i+1);
        if any(idx)
            values = z_scored(idx);
            bin_means(i) = mean(values);
            bin_sems(i) = std(values) / sqrt(length(values));
        end
    end
    
    % Remove empty bins
    valid = ~isnan(bin_means);
    bin_centers = bin_centers(valid);
    bin_means = bin_means(valid);
    bin_sems = bin_sems(valid);
    
    % Smooth if needed
    if length(bin_means) > options.smoothing
        bin_means = movmean(bin_means, options.smoothing);
    end
    
    % Transform coordinates
    if strcmp(run_info.direction, 'towards')
        x_data = -bin_centers;
        range_min = -options.slope_range;
        range_max = 0;
    else
        x_data = bin_centers;
        range_min = 0;
        range_max = options.slope_range;
    end
    
    % Filter to range
    range_idx = x_data >= range_min & x_data <= range_max;
    x_fit = x_data(range_idx);
    y_fit = bin_means(range_idx);
    
    % Linear regression
    if length(x_fit) >= 3
        X = [ones(length(x_fit), 1), x_fit(:)];
        beta = X \ y_fit(:);
        
        y_pred = X * beta;
        residuals = y_fit(:) - y_pred;
        mse = sum(residuals.^2) / (length(x_fit) - 2);
        se_slope = sqrt(mse * inv(X' * X));
        se_slope = se_slope(2, 2);
        
        t_stat = beta(2) / se_slope;
        p_value = 2 * (1 - tcdf(abs(t_stat), length(x_fit) - 2));
        
        ss_tot = sum((y_fit(:) - mean(y_fit)).^2);
        ss_res = sum(residuals.^2);
        r_squared = 1 - (ss_res / ss_tot);
        
        slope_stats = struct('slope', beta(2), 'intercept', beta(1), ...
                           'se_slope', se_slope, 't_stat', t_stat, ...
                           'p_value', p_value, 'r_squared', r_squared, ...
                           'is_significant', p_value < 0.05);
    else
        slope_stats = struct('slope', NaN, 'intercept', NaN, 'se_slope', NaN, ...
                           't_stat', NaN, 'p_value', NaN, 'r_squared', NaN, ...
                           'is_significant', false);
    end
    
    bin_data = struct('centers', bin_centers, 'means', bin_means, 'sems', bin_sems);
end

function [slope_stats, bin_data, steps] = calculate_detailed_slope(run_info, options)
    % Detailed slope calculation with all intermediate steps
    
    % Step 1: Z-score
    dff_mean = mean(run_info.dff);
    dff_std = std(run_info.dff);
    z_scored = (run_info.dff - dff_mean) / dff_std;
    
    % Step 2: Bin data
    bin_edges = 0:options.bin_width:options.max_distance;
    bin_centers = bin_edges(1:end-1) + options.bin_width/2;
    
    bin_means = nan(length(bin_centers), 1);
    bin_sems = nan(length(bin_centers), 1);
    
    for i = 1:length(bin_centers)
        idx = run_info.distance >= bin_edges(i) & run_info.distance < bin_edges(i+1);
        if any(idx)
            values = z_scored(idx);
            bin_means(i) = mean(values);
            bin_sems(i) = std(values) / sqrt(length(values));
        end
    end
    
    % Remove empty bins
    valid = ~isnan(bin_means);
    bin_centers = bin_centers(valid);
    bin_means = bin_means(valid);
    bin_sems = bin_sems(valid);
    
    % Step 3: Smooth
    if length(bin_means) > options.smoothing
        y_smoothed = movmean(bin_means, options.smoothing);
    else
        y_smoothed = bin_means;
    end
    
    % Step 4: Transform coordinates
    if strcmp(run_info.direction, 'towards')
        x_smooth = -bin_centers;
        range_min = -options.slope_range;
        range_max = 0;
    else
        x_smooth = bin_centers;
        range_min = 0;
        range_max = options.slope_range;
    end
    
    % Step 5: Filter to fitting range
    range_idx = x_smooth >= range_min & x_smooth <= range_max;
    x_fit = x_smooth(range_idx);
    y_fit = y_smoothed(range_idx);
    
    % Step 6: Linear regression
    if length(x_fit) >= 3
        X = [ones(length(x_fit), 1), x_fit(:)];
        beta = X \ y_fit(:);
        
        y_pred = X * beta;
        residuals = y_fit(:) - y_pred;
        mse = sum(residuals.^2) / (length(x_fit) - 2);
        cov_matrix = mse * inv(X' * X);
        se_slope = sqrt(cov_matrix(2, 2));
        
        t_stat = beta(2) / se_slope;
        p_value = 2 * (1 - tcdf(abs(t_stat), length(x_fit) - 2));
        
        ss_tot = sum((y_fit(:) - mean(y_fit)).^2);
        ss_res = sum(residuals.^2);
        r_squared = 1 - (ss_res / ss_tot);
        
        % Create fitted line
        x_line = linspace(min(x_fit), max(x_fit), 100);
        y_line = beta(1) + beta(2) * x_line;
        
        slope_stats = struct('slope', beta(2), 'intercept', beta(1), ...
                           'se_slope', se_slope, 't_stat', t_stat, ...
                           'p_value', p_value, 'r_squared', r_squared, ...
                           'is_significant', p_value < 0.05);
    else
        slope_stats = struct('slope', NaN, 'intercept', NaN, 'se_slope', NaN, ...
                           't_stat', NaN, 'p_value', NaN, 'r_squared', NaN, ...
                           'is_significant', false);
        x_line = [];
        y_line = [];
    end
    
    % Store intermediate steps
    steps = struct('z_scored', z_scored, 'x_smooth', x_smooth, 'y_smooth', y_smoothed, ...
                  'x_fit', x_fit, 'y_fit', y_fit, 'x_line', x_line, 'y_line', y_line);
    
    bin_data = struct('centers', bin_centers, 'means', bin_means, 'sems', bin_sems);
end