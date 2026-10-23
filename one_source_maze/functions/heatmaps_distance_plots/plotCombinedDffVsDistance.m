function plotCombinedDffVsDistance(mice_all, options)
    % Plot combined df/f vs distance for all mice from same state and source
    % before and after food discovery with various normalization and plotting options
    %
    % Parameters:
    %   mice_all: cell array with mouse data
    %   options: struct with options for filtering
    %     - state: 'fasted' or 'fed' 
    %     - source: 'food' or 'gel'
    %     - dist_limit: minimum distance threshold (default: 5)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - zscore_method: 'pooled', 'within_mouse' (default), or 'both'
    %     - show_individual_corr: boolean, show individual mouse correlations (default: false)
    %     - bin_data: boolean, whether to bin distance data (default: false)
    %     - bin_size: size of distance bins in cm (default: 1)
    %     - plot_style: 'scatter' (default), 'binned_scatter', 'line', or 'errorbar'
    %     - subsample_factor: if > 1, take every Nth point (default: 1, no subsampling)
    
    % Define constants
    COL_DIST = 5;     % Distance to food 
    COL_DOOR = 7;     % Door status
    COL_GROOM = 10;   % Grooming
    COL_EATING = 9;   % Eating
    COL_DFF = 11;     % DFF data
    
    % Define color scheme
    color_before = [0.4, 0.7, 0.9]; % Light blue for before
    color_after = [0.2, 0.4, 0.7];  % Dark blue for after
    color_before_line = [0.2, 0.5, 0.7]; % Slightly darker for regression line
    color_after_line = [0.1, 0.2, 0.5];  % Slightly darker for regression line
    
    % Set default options if not provided
    if ~isfield(options, 'dist_limit')
        options.dist_limit = 5;
    end
    
    if ~isfield(options, 'remove_grooming')
        options.remove_grooming = true;
    end
    
    if ~isfield(options, 'zscore_method')
        options.zscore_method = 'within_mouse';
    end
    
    if ~isfield(options, 'show_individual_corr')
        options.show_individual_corr = false;
    end
    
    if ~isfield(options, 'bin_data')
        options.bin_data = false;
    end
    
    if ~isfield(options, 'bin_size')
        options.bin_size = 1; % 1 cm bins
    end
    
    if ~isfield(options, 'plot_style')
        options.plot_style = 'scatter';
    end
    
    if ~isfield(options, 'subsample_factor')
        options.subsample_factor = 1;
    end
    
    if ~isfield(options, 'state')
        error('Please specify state: "fasted" or "fed"');
    end
    
    if ~isfield(options, 'source')
        error('Please specify source: "food" or "gel"');
    end
    
    % Initialize arrays to collect data
    dff_before_all = [];
    dist_before_all = [];
    dff_after_all = [];
    dist_after_all = [];
    mouse_ids = {};
    individual_corr_before = [];
    individual_corr_after = [];
    
    % Loop through all mice and filter by state and source
    valid_mice_count = 0;
    
    for i = 1:size(mice_all, 1)
        % Extract mouse info
        mouse_id = mice_all{i, 1};
        condition = mice_all{i, 2};
        stimulus = mice_all{i, 3};
        data = mice_all{i, 4};
        
        % Check if discovery frame exists
        if size(mice_all, 2) >= 6 && ~isempty(mice_all{i, 6})
            discovery = mice_all{i, 6};
        else
            fprintf('Warning: No discovery frame for mouse %s, skipping\n', mouse_id);
            continue;
        end
        
        % Filter by state and source
        if ~strcmp(condition, options.state) || ~strcmp(stimulus, options.source)
            continue;
        end
        
        % Check if we have enough data
        if isempty(data) || discovery <= 1 || discovery >= size(data, 1)
            fprintf('Warning: Invalid data or discovery frame for mouse %s, skipping\n', mouse_id);
            continue;
        end
        
        valid_mice_count = valid_mice_count + 1;
        mouse_ids{end+1} = mouse_id;
        
        % Find end frame (second closed door or end of data)
        if length(data) > 11000
            closed_indices = find(data(discovery:end, COL_DOOR) < 1);
            if ~isempty(closed_indices)
                end_frame = discovery + closed_indices(1) - 1;
            else
                end_frame = length(data);
            end
        else
            end_frame = length(data);
        end
        
        % Process data before discovery
        if options.remove_grooming
            valid_before = data(1:discovery, COL_DIST) > options.dist_limit & ...
                          data(1:discovery, COL_GROOM) == 0 & ...
                          data(1:discovery, COL_EATING) == 0;
        else
            valid_before = data(1:discovery, COL_DIST) > options.dist_limit & ...
                          data(1:discovery, COL_EATING) == 0;
        end
        
        % Process data after discovery
        if options.remove_grooming
            valid_after = data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                         data(discovery:end_frame, COL_GROOM) == 0 & ...
                         data(discovery:end_frame, COL_EATING) == 0;
        else
            valid_after = data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                         data(discovery:end_frame, COL_EATING) == 0;
        end
        
        % Get df/f and distance values for this mouse
        if any(valid_before)
            dff_before_mouse = data(1:discovery, COL_DFF);
            dist_before_mouse = data(1:discovery, COL_DIST);
            dff_before_mouse = dff_before_mouse(valid_before);
            dist_before_mouse = dist_before_mouse(valid_before);
        else
            dff_before_mouse = [];
            dist_before_mouse = [];
        end
        
        if any(valid_after)
            dff_after_mouse = data(discovery:end_frame, COL_DFF);
            dist_after_mouse = data(discovery:end_frame, COL_DIST);
            dff_after_mouse = dff_after_mouse(valid_after);
            dist_after_mouse = dist_after_mouse(valid_after);
        else
            dff_after_mouse = [];
            dist_after_mouse = [];
        end
        
        % Calculate individual mouse correlations
        if length(dff_before_mouse) > 3
            [r_before, ~] = corr(dff_before_mouse, dist_before_mouse);
            individual_corr_before(end+1) = r_before;
        else
            individual_corr_before(end+1) = NaN;
        end
        
        if length(dff_after_mouse) > 3
            [r_after, ~] = corr(dff_after_mouse, dist_after_mouse);
            individual_corr_after(end+1) = r_after;
        else
            individual_corr_after(end+1) = NaN;
        end
        
        % Apply within-mouse z-scoring if needed
        if strcmp(options.zscore_method, 'within_mouse') || strcmp(options.zscore_method, 'both')
            % Z-score within each mouse (normalize each mouse's data separately)
            all_mouse_dff = [dff_before_mouse; dff_after_mouse];
            if length(all_mouse_dff) > 1
                mouse_mean = mean(all_mouse_dff);
                mouse_std = std(all_mouse_dff);
                if mouse_std > 0
                    dff_before_mouse_z = (dff_before_mouse - mouse_mean) / mouse_std;
                    dff_after_mouse_z = (dff_after_mouse - mouse_mean) / mouse_std;
                else
                    dff_before_mouse_z = zeros(size(dff_before_mouse));
                    dff_after_mouse_z = zeros(size(dff_after_mouse));
                end
            else
                dff_before_mouse_z = dff_before_mouse;
                dff_after_mouse_z = dff_after_mouse;
            end
        else
            % Keep raw values for pooled z-scoring
            dff_before_mouse_z = dff_before_mouse;
            dff_after_mouse_z = dff_after_mouse;
        end
        
        % Collect data
        dff_before_all = [dff_before_all; dff_before_mouse_z];
        dist_before_all = [dist_before_all; dist_before_mouse];
        dff_after_all = [dff_after_all; dff_after_mouse_z];
        dist_after_all = [dist_after_all; dist_after_mouse];
    end
    
    if valid_mice_count == 0
        error('No mice found matching state "%s" and source "%s"', options.state, options.source);
    end
    
    fprintf('Found %d mice matching state "%s" and source "%s"\n', valid_mice_count, options.state, options.source);
    
    % Apply pooled z-scoring if needed
    if strcmp(options.zscore_method, 'pooled')
        all_dff = [dff_before_all; dff_after_all];
        dff_mean = mean(all_dff);
        dff_std = std(all_dff);
        dff_before_all = (dff_before_all - dff_mean) / dff_std;
        dff_after_all = (dff_after_all - dff_mean) / dff_std;
    end
    
    % Apply subsampling if requested
    if options.subsample_factor > 1
        % Subsample before data
        indices_before = 1:options.subsample_factor:length(dff_before_all);
        dff_before_all = dff_before_all(indices_before);
        dist_before_all = dist_before_all(indices_before);
        
        % Subsample after data
        indices_after = 1:options.subsample_factor:length(dff_after_all);
        dff_after_all = dff_after_all(indices_after);
        dist_after_all = dist_after_all(indices_after);
    end
    
    % Apply binning if requested
    if options.bin_data
        [dff_before_all, dist_before_all, sem_before_all] = bin_data_by_distance(dff_before_all, dist_before_all, options.bin_size);
        [dff_after_all, dist_after_all, sem_after_all] = bin_data_by_distance(dff_after_all, dist_after_all, options.bin_size);
    end
    
    % Calculate correlations
    [rho_before, pval_before] = corr(dff_before_all, dist_before_all, 'Type', 'Pearson');
    [rho_after, pval_after] = corr(dff_after_all, dist_after_all, 'Type', 'Pearson');
    
    % Determine figure layout
    if options.show_individual_corr
        figure('Position', [100, 100, 1800, 500]);
        subplot_layout = [1, 3];
    else
        figure('Position', [100, 100, 1200, 500]);
        subplot_layout = [1, 2];
    end
    
    % Before discovery plot
    subplot(subplot_layout(1), subplot_layout(2), 1);
    
    % Create plot based on style
    plot_data(dist_before_all, dff_before_all, options, color_before);
    if options.bin_data && exist('sem_before_all', 'var') && strcmp(options.plot_style, 'errorbar')
        hold on;
        errorbar(dist_before_all, dff_before_all, sem_before_all, 'k.', 'LineWidth', 1);
    end
    
    % Add regression line
    if length(dist_before_all) > 1
        hold on;
        p_before = polyfit(dist_before_all, dff_before_all, 1);
        x_range = linspace(min(dist_before_all), max(dist_before_all), 100);
        y_fit = polyval(p_before, x_range);
        plot(x_range, y_fit, 'LineWidth', 3, 'Color', color_before_line);
    end
    
    % Add correlation info
    add_correlation_text(rho_before, pval_before, 0.05, 0.95);
    
    % Format plot
    xlabel('Distance to Food (cm)', 'FontSize', 14);
    ylabel('Z-scored \Delta F/F', 'FontSize', 14);
    title('Before Food Discovery', 'FontSize', 16, 'FontWeight', 'bold');
    grid off; box off;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    
    % After discovery plot
    subplot(subplot_layout(1), subplot_layout(2), 2);
    
    % Create plot based on style  
    plot_data(dist_after_all, dff_after_all, options, color_after);
    if options.bin_data && exist('sem_after_all', 'var') && strcmp(options.plot_style, 'errorbar')
        hold on;
        errorbar(dist_after_all, dff_after_all, sem_after_all, 'k.', 'LineWidth', 1);
    end
    
    % Add regression line
    if length(dist_after_all) > 1
        hold on;
        p_after = polyfit(dist_after_all, dff_after_all, 1);
        x_range = linspace(min(dist_after_all), max(dist_after_all), 100);
        y_fit = polyval(p_after, x_range);
        plot(x_range, y_fit, 'LineWidth', 3, 'Color', color_after_line);
    end
    
    % Add correlation info
    add_correlation_text(rho_after, pval_after, 0.05, 0.95);
    
    % Format plot
    xlabel('Distance to Food (cm)', 'FontSize', 14);
    ylabel('Z-scored \Delta F/F', 'FontSize', 14);
    title('After Food Discovery', 'FontSize', 16, 'FontWeight', 'bold');
    grid off; box off;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    
    % Individual correlations plot
    if options.show_individual_corr
        subplot(subplot_layout(1), subplot_layout(2), 3);
        
        % Create paired data for individual correlations
        valid_pairs = ~isnan(individual_corr_before) & ~isnan(individual_corr_after);
        corr_before_valid = individual_corr_before(valid_pairs);
        corr_after_valid = individual_corr_after(valid_pairs);
        
        if ~isempty(corr_before_valid) && ~isempty(corr_after_valid)
            % Bar plot of means
            means = [mean(corr_before_valid), mean(corr_after_valid)];
            sems = [std(corr_before_valid)/sqrt(length(corr_before_valid)), ...
                    std(corr_after_valid)/sqrt(length(corr_after_valid))];
            
            b = bar(1:2, means, 'FaceColor', 'flat');
            b.CData(1,:) = color_before;
            b.CData(2,:) = color_after;
            hold on;
            
            % Error bars
            errorbar(1:2, means, sems, 'k.', 'LineWidth', 1.5);
            
            % Individual points
            for j = 1:length(corr_before_valid)
                plot([1, 2], [corr_before_valid(j), corr_after_valid(j)], 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.7]);
            end
            scatter(ones(size(corr_before_valid))+randn(size(corr_before_valid))*0.05, corr_before_valid, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
            scatter(2*ones(size(corr_after_valid))+randn(size(corr_after_valid))*0.05, corr_after_valid, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
            
            % Add horizontal line at zero
            line([0.5, 2.5], [0, 0], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
            
            % Statistical test
            if length(corr_before_valid) > 1 && length(corr_after_valid) > 1
                [~, p_paired] = ttest(corr_before_valid, corr_after_valid);
                add_significance_bar([1, 2], p_paired, max([corr_before_valid, corr_after_valid]) + 0.1);
            end
        end
        
        % Format plot
        set(gca, 'XTick', 1:2, 'XTickLabel', {'Before', 'After'});
        ylabel('Individual Mouse Correlations', 'FontSize', 14);
        title(sprintf('Individual Correlations\n(n=%d mice)', sum(valid_pairs)), 'FontSize', 14, 'FontWeight', 'bold');
        xlim([0.5, 2.5]);
        grid off; box off;
        set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    end
    
    % Synchronize y-axis limits for main plots
    sync_axis_limits(subplot_layout);
    
    % Add main title
    create_main_title(options, valid_mice_count, length(dff_before_all), length(dff_after_all));
    
    % Print summary
    print_summary(options, valid_mice_count, rho_before, pval_before, rho_after, pval_after, ...
                  individual_corr_before, individual_corr_after, mouse_ids);
end

function add_correlation_text(rho, pval, x_pos, y_pos)
    if ~isnan(rho)
        if pval < 0.001
            p_text = 'p < 0.001***';
        elseif pval < 0.01
            p_text = 'p < 0.01**';
        elseif pval < 0.05
            p_text = 'p < 0.05*';
        else
            p_text = sprintf('p = %.3f', pval);
        end
        
        text(x_pos, y_pos, {sprintf('r = %.3f', rho), p_text}, ...
             'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'none');
    end
end

function add_significance_bar(positions, p_value, height)
    if isnan(p_value)
        return;
    end
    
    plot(positions, [height, height], 'k-', 'LineWidth', 1.5);
    
    if p_value < 0.001
        p_text = '***';
    elseif p_value < 0.01
        p_text = '**';
    elseif p_value < 0.05
        p_text = '*';
    else
        p_text = 'ns';
    end
    
    text(mean(positions), height*1.02, p_text, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

function sync_axis_limits(subplot_layout)
    % Get all subplot handles for the first two plots
    ax1 = subplot(subplot_layout(1), subplot_layout(2), 1);
    ax2 = subplot(subplot_layout(1), subplot_layout(2), 2);
    
    % Sync y-limits
    y1 = ylim(ax1);
    y2 = ylim(ax2);
    y_common = [min(y1(1), y2(1)), max(y1(2), y2(2))];
    ylim(ax1, y_common);
    ylim(ax2, y_common);
    
    % Sync x-limits
    x1 = xlim(ax1);
    x2 = xlim(ax2);
    x_common = [min(x1(1), x2(1)), max(x1(2), x2(2))];
    xlim(ax1, x_common);
    xlim(ax2, x_common);
end

function create_main_title(options, valid_mice_count, n_before, n_after)
    state_title = options.state;
    state_title(1) = upper(state_title(1));
    source_title = options.source;
    source_title(1) = upper(source_title(1));
    
    method_text = '';
    if strcmp(options.zscore_method, 'within_mouse')
        method_text = ' (Within-mouse Z-scoring)';
    elseif strcmp(options.zscore_method, 'pooled')
        method_text = ' (Pooled Z-scoring)';
    end
    
    sgtitle({
        sprintf('%s - %s: ΔF/F vs Distance to Food%s', state_title, source_title, method_text), ...
        sprintf('Combined data from %d mice (%d points before, %d points after)', ...
               valid_mice_count, n_before, n_after)
    }, 'FontSize', 16, 'FontWeight', 'bold');
end

function print_summary(options, valid_mice_count, rho_before, pval_before, rho_after, pval_after, ...
                      individual_corr_before, individual_corr_after, mouse_ids)
    fprintf('\n=== Summary Statistics ===\n');
    fprintf('State: %s, Source: %s\n', options.state, options.source);
    fprintf('Z-scoring method: %s\n', options.zscore_method);
    fprintf('Plot style: %s\n', options.plot_style);
    if options.bin_data
        fprintf('Data binned: %g cm bins\n', options.bin_size);
    end
    if options.subsample_factor > 1
        fprintf('Subsampled: every %d points\n', options.subsample_factor);
    end
    fprintf('Number of mice: %d\n', valid_mice_count);
    fprintf('Combined correlation before: r = %.4f, p = %.4f\n', rho_before, pval_before);
    fprintf('Combined correlation after: r = %.4f, p = %.4f\n', rho_after, pval_after);
    
    % Individual mouse statistics
    valid_before = ~isnan(individual_corr_before);
    valid_after = ~isnan(individual_corr_after);
    
    if any(valid_before)
        fprintf('Individual correlations before: mean = %.4f ± %.4f (n=%d)\n', ...
                mean(individual_corr_before(valid_before)), std(individual_corr_before(valid_before)), sum(valid_before));
    end
    
    if any(valid_after)
        fprintf('Individual correlations after: mean = %.4f ± %.4f (n=%d)\n', ...
                mean(individual_corr_after(valid_after)), std(individual_corr_after(valid_after)), sum(valid_after));
    end
    
    % Test if individual correlations are significantly different from zero
    if any(valid_before)
        [~, p_before_vs_zero] = ttest(individual_corr_before(valid_before));
        fprintf('Before correlations vs zero: p = %.4f\n', p_before_vs_zero);
    end
    
    if any(valid_after)
        [~, p_after_vs_zero] = ttest(individual_corr_after(valid_after));
        fprintf('After correlations vs zero: p = %.4f\n', p_after_vs_zero);
    end
    
    % Paired comparison
    valid_pairs = valid_before & valid_after;
    if sum(valid_pairs) > 1
        [~, p_paired] = ttest(individual_corr_before(valid_pairs), individual_corr_after(valid_pairs));
        fprintf('Paired comparison (before vs after): p = %.4f\n', p_paired);
    end
    
    fprintf('\nMice included:\n');
    for i = 1:length(mouse_ids)
        if i <= length(individual_corr_before) && i <= length(individual_corr_after)
            fprintf('  %s: r_before=%.3f, r_after=%.3f\n', mouse_ids{i}, individual_corr_before(i), individual_corr_after(i));
        else
            fprintf('  %s\n', mouse_ids{i});
        end
    end
end

function plot_data(dist_data, dff_data, options, default_color)
    % Plot data according to the specified style
    
    switch options.plot_style
        case 'scatter'
            % Standard scatter plot
            scatter(dist_data, dff_data, 20, 'filled', 'MarkerFaceColor', default_color, 'MarkerFaceAlpha', 0.6);
            
        case 'binned_scatter'
            % Larger points for binned data
            scatter(dist_data, dff_data, 60, 'filled', 'MarkerFaceColor', default_color, 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'k');
            
        case 'line'
            % Line plot
            plot(dist_data, dff_data, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', default_color, 'MarkerFaceColor', default_color);
            
        case 'errorbar'
            % Line plot (error bars added separately)
            plot(dist_data, dff_data, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', default_color, 'MarkerFaceColor', default_color);
    end
end

function [dff_binned, dist_binned, sem_binned] = bin_data_by_distance(dff_data, dist_data, bin_size)
    % Bin data by distance and calculate means and SEMs
    
    if isempty(dff_data) || isempty(dist_data)
        dff_binned = [];
        dist_binned = [];
        sem_binned = [];
        return;
    end
    
    % Define bin edges
    min_dist = floor(min(dist_data) / bin_size) * bin_size;
    max_dist = ceil(max(dist_data) / bin_size) * bin_size;
    bin_edges = min_dist:bin_size:max_dist;
    
    % Initialize output arrays
    dff_binned = [];
    dist_binned = [];
    sem_binned = [];
    
    % Process each bin
    for i = 1:(length(bin_edges)-1)
        % Find data points in this bin
        in_bin = dist_data >= bin_edges(i) & dist_data < bin_edges(i+1);
        
        if sum(in_bin) > 0
            % Calculate statistics for this bin
            bin_dff = dff_data(in_bin);
            bin_dist = dist_data(in_bin);
            
            % Store results
            dff_binned(end+1) = mean(bin_dff);
            dist_binned(end+1) = mean(bin_dist);  % Use mean distance in bin
            sem_binned(end+1) = std(bin_dff) / sqrt(length(bin_dff));
        end
    end
    
    % Convert to column vectors
    dff_binned = dff_binned(:);
    dist_binned = dist_binned(:);
    sem_binned = sem_binned(:);
end