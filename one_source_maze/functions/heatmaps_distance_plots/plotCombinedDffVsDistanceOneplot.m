function plotCombinedDffVsDistanceOneplot(mice_all, options)
    % Plot combined df/f vs distance for all mice from same state and source
    % before and after food discovery with various normalization and plotting options
    % NOW PLOTS BOTH BEFORE AND AFTER ON THE SAME PLOT WITH TWO REGRESSION LINES
    %
    % Parameters:
    %   mice_all: cell array with mouse data
    %   options: struct with options for filtering
    %     - state: 'fasted' or 'fed' 
    %     - source: 'food' or 'gel'
    %     - dist_limit: minimum distance threshold (default: 5)
    %     - speed_threshold: minimum speed threshold in cm/s (default: 0.5)
    %     - remove_grooming: boolean, whether to remove grooming periods (default: true)
    %     - zscore_method: 'pooled', 'within_mouse' (default), or 'both'
    %     - show_individual_corr: boolean, show individual mouse correlations (default: false)
    %     - bin_data: boolean, whether to bin distance data (default: false)
    %     - bin_size: size of distance bins in cm (default: 1)
    %     - plot_style: 'scatter' (default), 'binned_scatter', 'line', or 'errorbar'
    %     - subsample_factor: if > 1, take every Nth point (default: 1, no subsampling)
    %     - show_confidence_intervals: boolean, show 95% CI on regression lines (default: false)
    
    % Define constants
    COL_SPEED = 4;    % Speed
    COL_DIST = 5;     % Distance to food
    COL_DOOR = 7;     % Door status
    COL_FOOD_INT = 8; % Food interaction
    COL_EATING = 9;   % Eating
    COL_GROOM = 10;   % Grooming
    COL_DFF = 11;     % DFF data
    
    % Define color scheme
    color_before = [0.4, 0.7, 0.9]; % Light blue for before
    color_after = [0.2, 0.4, 0.7];  % Dark blue for after
    color_before_line = [0.2, 0.5, 0.7]; % Slightly darker for regression line
    color_after_line = [0.1, 0.2, 0.5];  % Darker for regression line
    
    % Set default options if not provided
    if ~isfield(options, 'dist_limit')
        options.dist_limit = 5;
    end
    
    if ~isfield(options, 'speed_threshold')
        options.speed_threshold = 0.5; % Default minimum speed threshold in cm/s
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
    
    if ~isfield(options, 'show_confidence_intervals')
        options.show_confidence_intervals = false; % Default: don't show CI
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
            valid_before = data(1:discovery, COL_SPEED) > options.speed_threshold & ...
                          data(1:discovery, COL_DIST) > options.dist_limit & ...
                          data(1:discovery, COL_GROOM) == 0 & ...
                          data(1:discovery, COL_EATING) == 0 & ...
                          data(1:discovery, COL_FOOD_INT) == 0;
        else
            valid_before = data(1:discovery, COL_SPEED) > options.speed_threshold & ...
                          data(1:discovery, COL_DIST) > options.dist_limit & ...
                          data(1:discovery, COL_EATING) == 0 & ...
                          data(1:discovery, COL_FOOD_INT) == 0;
        end
        
        % Process data after discovery
        if options.remove_grooming
            valid_after = data(discovery:end_frame, COL_SPEED) > options.speed_threshold & ...
                         data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                         data(discovery:end_frame, COL_GROOM) == 0 & ...
                         data(discovery:end_frame, COL_EATING) == 0 & ...
                         data(discovery:end_frame, COL_FOOD_INT) == 0;
        else
            valid_after = data(discovery:end_frame, COL_SPEED) > options.speed_threshold & ...
                         data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
                         data(discovery:end_frame, COL_EATING) == 0 & ...
                         data(discovery:end_frame, COL_FOOD_INT) == 0;
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
    fprintf('Applied speed filtering: speed ≥ %.2f cm/s\n', options.speed_threshold);
    
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
        figure('Position', [100, 100, 1400, 500]);
        subplot_layout = [1, 2];
        main_plot_subplot = 1;
        individual_plot_subplot = 2;
    else
        figure('Position', [100, 100, 800, 600]);
        subplot_layout = [1, 1];
        main_plot_subplot = 1;
    end
    
    % Main combined plot (before and after on same axes)
    if options.show_individual_corr
        subplot(subplot_layout(1), subplot_layout(2), main_plot_subplot);
    end
    
    % Plot before discovery data
    plot_data(dist_before_all, dff_before_all, options, color_before, 'before');
    if options.bin_data && exist('sem_before_all', 'var') && strcmp(options.plot_style, 'errorbar')
        hold on;
        errorbar(dist_before_all, dff_before_all, sem_before_all, 'Color', color_before, 'LineWidth', 1);
    end
    
    hold on;
    
    % Plot after discovery data
    plot_data(dist_after_all, dff_after_all, options, color_after, 'after');
    if options.bin_data && exist('sem_after_all', 'var') && strcmp(options.plot_style, 'errorbar')
        errorbar(dist_after_all, dff_after_all, sem_after_all, 'Color', color_after, 'LineWidth', 1);
    end
    
    % Add regression lines
    if length(dist_before_all) > 1
        p_before = polyfit(dist_before_all, dff_before_all, 1);
        x_range_before = linspace(min(dist_before_all), max(dist_before_all), 100);
        y_fit_before = polyval(p_before, x_range_before);
        plot(x_range_before, y_fit_before, 'LineWidth', 3, 'Color', color_before_line, 'LineStyle', '-');
        
        % Add confidence intervals if requested
        if options.show_confidence_intervals
            add_confidence_interval(dist_before_all, dff_before_all, x_range_before, color_before_line);
        end
    end
    
    if length(dist_after_all) > 1
        p_after = polyfit(dist_after_all, dff_after_all, 1);
        x_range_after = linspace(min(dist_after_all), max(dist_after_all), 100);
        y_fit_after = polyval(p_after, x_range_after);
        plot(x_range_after, y_fit_after, 'LineWidth', 3, 'Color', color_after_line, 'LineStyle', '-');
        
        % Add confidence intervals if requested
        if options.show_confidence_intervals
            add_confidence_interval(dist_after_all, dff_after_all, x_range_after, color_after_line);
        end
    end

    % Add correlation text for both periods (REMOVED - text box disabled)
    % add_correlation_text_combined(rho_before, pval_before, rho_after, pval_after, length(dff_before_all), length(dff_after_all));

    % Add legend
    legend_entries = {'Before Discovery', 'After Discovery', 'Before Fit', 'After Fit'};

    % Create legend manually with proper handles
    legend_handles = [];
    legend_labels = {};
    
    % Add data point handles
    h_before = scatter(NaN, NaN, 50, 'filled', 'MarkerFaceColor', color_before, 'MarkerFaceAlpha', 0.6);
    legend_handles(end+1) = h_before;
    legend_labels{end+1} = 'Before Discovery';
    
    h_after = scatter(NaN, NaN, 50, 'filled', 'MarkerFaceColor', color_after, 'MarkerFaceAlpha', 0.6);
    legend_handles(end+1) = h_after;
    legend_labels{end+1} = 'After Discovery';
    
    % Add regression line handles
    if length(dist_before_all) > 1
        h_fit_before = plot(NaN, NaN, 'LineWidth', 3, 'Color', color_before_line);
        legend_handles(end+1) = h_fit_before;
        legend_labels{end+1} = sprintf('Before Fit (r=%.3f)', rho_before);
    end
    
    if length(dist_after_all) > 1
        h_fit_after = plot(NaN, NaN, 'LineWidth', 3, 'Color', color_after_line);
        legend_handles(end+1) = h_fit_after;
        legend_labels{end+1} = sprintf('After Fit (r=%.3f)', rho_after);
    end
    
    legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 10);
    legend('boxoff')
    
    % Format plot
    xlabel('Distance to Food (cm)', 'FontSize', 14);
    ylabel('Z-scored \Delta F/F', 'FontSize', 14);
    title('ΔF/F vs Distance: Before and After Food Discovery', 'FontSize', 16, 'FontWeight', 'bold');
    grid off; box off;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    
    % Individual correlations plot
    if options.show_individual_corr
        subplot(subplot_layout(1), subplot_layout(2), individual_plot_subplot);
        
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
    
    % Add main title
    create_main_title(options, valid_mice_count, length(dff_before_all), length(dff_after_all));
    
    % Print summary
    print_summary(options, valid_mice_count, rho_before, pval_before, rho_after, pval_after, ...
                  individual_corr_before, individual_corr_after, mouse_ids, ...
                  length(dff_before_all), length(dff_after_all));
end

function add_correlation_text_combined(rho_before, pval_before, rho_after, pval_after, n_before, n_after)
    % Add correlation text for both before and after on the same plot
    % Now includes practical significance assessment
    
    % Define practical significance threshold
    practical_threshold = 0.1; % |r| > 0.1 considered practically meaningful
    
    % Format correlation info with practical significance
    function [text_str] = format_correlation_info(rho, pval, n, period_name)
        % Calculate R-squared (% variance explained)
        r_squared = rho^2 * 100;
        
        % Statistical significance
        if ~isnan(pval)
            if pval < 0.001
                stat_sig = '***';
            elseif pval < 0.01
                stat_sig = '**';
            elseif pval < 0.05
                stat_sig = '*';
            else
                stat_sig = 'ns';
            end
            p_text = sprintf('p = %.4f%s', pval, stat_sig);
        else
            p_text = 'p = NaN';
            stat_sig = '';
        end
        
        % Practical significance
        if abs(rho) >= practical_threshold
            practical_sig = ' [MEANINGFUL]';
            practical_color = '';
        else
            practical_sig = ' [weak effect]';
            practical_color = '';
        end
        
        % Create text
        text_str = {
            sprintf('%s:', period_name),
            sprintf('r = %.4f%s', rho, practical_sig),
            sprintf('R² = %.3f%% var. explained', r_squared),
            sprintf('%s (n = %d)', p_text, n)
        };
    end
    
    % Format both periods
    before_text = format_correlation_info(rho_before, pval_before, n_before, 'Before Discovery');
    after_text = format_correlation_info(rho_after, pval_after, n_after, 'After Discovery');
    
    % Combine text
    text_str = [before_text; {''}; after_text; {''}; ...
               {sprintf('Practical threshold: |r| > %.1f', practical_threshold)}];
    
    % Add warning for large sample sizes with weak correlations
    max_n = max(n_before, n_after);
    if max_n > 1000 && (abs(rho_before) < practical_threshold || abs(rho_after) < practical_threshold)
        text_str{end+1} = '';
        text_str{end+1} = 'Note: Large n can make weak';
        text_str{end+1} = 'correlations statistically significant';
    end
    
    text(0.02, 0.98, text_str, ...
         'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1);
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
    
end

function print_summary(options, valid_mice_count, rho_before, pval_before, rho_after, pval_after, ...
                      individual_corr_before, individual_corr_after, mouse_ids, n_before, n_after)
    fprintf('\n=== Summary Statistics ===\n');
    fprintf('State: %s, Source: %s\n', options.state, options.source);
    fprintf('Speed threshold: %.2f cm/s\n', options.speed_threshold);
    fprintf('Distance threshold: %.1f cm\n', options.dist_limit);
    fprintf('Z-scoring method: %s\n', options.zscore_method);
    fprintf('Plot style: %s\n', options.plot_style);
    if options.bin_data
        fprintf('Data binned: %g cm bins\n', options.bin_size);
    end
    if options.subsample_factor > 1
        fprintf('Subsampled: every %d points\n', options.subsample_factor);
    end
    fprintf('Number of mice: %d\n', valid_mice_count);
    
    % Add warning about pooled correlation analysis
    max_n = max(n_before, n_after);
    if max_n > 1000
        fprintf('\n*** IMPORTANT NOTE ***\n');
        fprintf('Large sample size (n > 1000) can make very weak correlations\n');
        fprintf('statistically significant. Focus on effect size (r-value) and R².\n');
    end
    
    fprintf('\nCombined correlation before: r = %.4f, p = %.4f (R² = %.3f%%, n = %d)\n', ...
            rho_before, pval_before, rho_before^2*100, n_before);
    fprintf('Combined correlation after: r = %.4f, p = %.4f (R² = %.3f%%, n = %d)\n', ...
            rho_after, pval_after, rho_after^2*100, n_after);
    
    % Add practical significance assessment
    practical_threshold = 0.1;
    fprintf('\nPractical significance assessment (|r| > %.1f):\n', practical_threshold);
    if abs(rho_before) >= practical_threshold
        fprintf('  Before: MEANINGFUL effect (r = %.4f)\n', rho_before);
    else
        fprintf('  Before: Weak effect (r = %.4f) - may not be practically meaningful\n', rho_before);
    end
    
    if abs(rho_after) >= practical_threshold
        fprintf('  After: MEANINGFUL effect (r = %.4f)\n', rho_after);
    else
        fprintf('  After: Weak effect (r = %.4f) - may not be practically meaningful\n', rho_after);
    end
    
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
    
    % Final interpretation guidance
    fprintf('\n=== INTERPRETATION GUIDANCE ===\n');
    fprintf('1. Statistical significance (p-value) can be misleading with large samples\n');
    fprintf('2. Focus on effect size: |r| > 0.1 = small, |r| > 0.3 = medium, |r| > 0.5 = large effect\n');
    fprintf('3. R² shows percentage of variance explained by the relationship\n');
    fprintf('4. Individual mouse correlations may be more meaningful than pooled analysis\n');
end

function plot_data(dist_data, dff_data, options, default_color, period)
    % Plot data according to the specified style
    % period parameter helps distinguish between 'before' and 'after' for styling
    
    % Adjust marker properties based on period
    if strcmp(period, 'before')
        marker_alpha = 0.4;  % More transparent for before
        marker_size = 15;
    else
        marker_alpha = 0.6;  % Less transparent for after
        marker_size = 18;
    end
    
    switch options.plot_style
        case 'scatter'
            % Standard scatter plot
            scatter(dist_data, dff_data, marker_size, 'filled', 'MarkerFaceColor', default_color, 'MarkerFaceAlpha', marker_alpha);
            
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

function add_confidence_interval(x_data, y_data, x_range, line_color)
    % Add 95% confidence interval to regression line
    % This helps visualize the uncertainty in the fit
    
    if length(x_data) < 3
        return; % Need at least 3 points for meaningful CI
    end
    
    % Fit linear model for confidence intervals
    n = length(x_data);
    
    % Calculate regression statistics
    x_mean = mean(x_data);
    y_mean = mean(y_data);
    
    % Calculate slope and intercept
    b1 = sum((x_data - x_mean) .* (y_data - y_mean)) / sum((x_data - x_mean).^2);
    b0 = y_mean - b1 * x_mean;
    
    % Calculate residual standard error
    y_pred = b0 + b1 * x_data;
    residuals = y_data - y_pred;
    mse = sum(residuals.^2) / (n - 2); % Mean squared error
    se = sqrt(mse); % Standard error of residuals
    
    % Calculate confidence intervals for the regression line
    sxx = sum((x_data - x_mean).^2);
    
    % For each point in x_range, calculate CI
    ci_upper = zeros(size(x_range));
    ci_lower = zeros(size(x_range));
    
    t_crit = 1.96; % Approximate 95% CI (for large n)
    if n < 30
        % Use t-distribution for small samples
        t_crit = 2.5; % Conservative estimate
    end
    
    for i = 1:length(x_range)
        x_val = x_range(i);
        y_pred_val = b0 + b1 * x_val;
        
        % Standard error for this prediction
        se_pred = se * sqrt(1/n + (x_val - x_mean)^2 / sxx);
        
        % Confidence interval
        margin = t_crit * se_pred;
        ci_upper(i) = y_pred_val + margin;
        ci_lower(i) = y_pred_val - margin;
    end
    
    % Plot confidence interval as filled area
    ci_color = line_color + (1 - line_color) * 0.7; % Lighter version of line color
    fill([x_range, fliplr(x_range)], [ci_upper, fliplr(ci_lower)], ...
         ci_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end