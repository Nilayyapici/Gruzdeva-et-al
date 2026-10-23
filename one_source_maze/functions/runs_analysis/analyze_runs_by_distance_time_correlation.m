function [results] = analyze_runs_by_distance_time_correlation(run_data, options)
    % Separate runs based on distance-time correlation and analyze dF/F relationships
    % 
    % This function:
    % 1. Calculates correlation between distance and time for each run
    % 2. Separates runs into high vs low distance-time correlation groups
    % 3. Analyzes how distance and time correlate with dF/F in each group
    %
    % Inputs:
    %   run_data - Output from analyze_food_runs_by_visits
    %   options - Structure with fields:
    %       .corr_threshold - Threshold for separating runs (default: 0.5)
    %       .session - Session to analyze (0, 1, 2, 3, or 'all')
    %       .group - Group filter: 'fasted', 'fed', or 'all' (default: 'all')
    %       .source - Source filter: 'food', 'gel', or 'all' (default: 'all')
    %       .min_run_length - Minimum points per run (default: 5)
    %       .plot_results - Generate plots (default: true)
    %
    % Output:
    %   results - Structure with separated run data and correlations
    
    % Set default options
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'corr_threshold')
        options.corr_threshold = 0.5;
    end
    if ~isfield(options, 'session')
        options.session = 'all';
    end
    if ~isfield(options, 'group')
        options.group = 'all';
    end
    if ~isfield(options, 'source')
        options.source = 'all';
    end
    if ~isfield(options, 'min_run_length')
        options.min_run_length = 5;
    end
    if ~isfield(options, 'plot_results')
        options.plot_results = true;
    end
    
    % Initialize results structure
    results = struct();
    directions = {'towards', 'away'};
    
    for d = 1:length(directions)
        direction = directions{d};
        
        % Initialize storage for this direction
        results.(direction).high_corr = struct(...
            'runs', [], ...
            'dist_time_corr', [], ...
            'dist_dff_corr', [], ...
            'time_dff_corr', [], ...
            'dist_dff_partial', [], ...
            'time_dff_partial', []);
        
        results.(direction).low_corr = struct(...
            'runs', [], ...
            'dist_time_corr', [], ...
            'dist_dff_corr', [], ...
            'time_dff_corr', [], ...
            'dist_dff_partial', [], ...
            'time_dff_partial', []);
    end
    
    % Process each mouse/session
    for i = 1:length(run_data)
        mouse_id = run_data(i).mouse_id;
        session = run_data(i).session;
        group = run_data(i).group;
        source = run_data(i).source;
        
        % Apply filters
        if ~strcmp(options.group, 'all') && ~strcmp(group, options.group)
            continue;
        end
        if ~strcmp(options.source, 'all') && ~strcmp(source, options.source)
            continue;
        end
        if ~strcmp(options.session, 'all') && options.session ~= session
            continue;
        end
        
        % Process each run
        for j = 1:length(run_data(i).runs)
            run = run_data(i).runs(j);
            
            % Skip runs that are too short
            if length(run.time) < options.min_run_length
                continue;
            end
            
            % Get run direction
            direction = run.type;
            
            % Calculate time from beginning of run
            time_relative = run.time - run.time(1);
            
            % Get data
            distance = run.distance;
            dff = run.dff;
            
            % Remove NaN values
            valid_idx = ~isnan(distance) & ~isnan(time_relative) & ~isnan(dff);
            if sum(valid_idx) < options.min_run_length
                continue;
            end
            
            distance = distance(valid_idx);
            time_relative = time_relative(valid_idx);
            dff = dff(valid_idx);
            
            % === Step 1: Calculate distance-time correlation ===
            [r_dist_time, p_dist_time] = corr(distance, time_relative, 'Type', 'Pearson');
            
            % === Step 2: Calculate distance-dF/F and time-dF/F correlations ===
            [r_dist_dff, p_dist_dff] = corr(distance, dff, 'Type', 'Pearson');
            [r_time_dff, p_time_dff] = corr(time_relative, dff, 'Type', 'Pearson');
            
            % === Step 3: Calculate partial correlations ===
            % Distance-dF/F controlling for time
            [r_dist_dff_partial, p_dist_dff_partial] = partialcorr(distance, dff, time_relative);
            
            % Time-dF/F controlling for distance
            [r_time_dff_partial, p_time_dff_partial] = partialcorr(time_relative, dff, distance);
            
            % === Step 4: Categorize run based on distance-time correlation ===
            if abs(r_dist_time) >= options.corr_threshold
                % High correlation between distance and time
                category = 'high_corr';
            else
                % Low correlation between distance and time
                category = 'low_corr';
            end
            
            % === Step 5: Store results ===
            run_info = struct(...
                'mouse_id', mouse_id, ...
                'session', session, ...
                'group', group, ...
                'source', source, ...
                'run_id', run.id, ...
                'n_points', length(distance), ...
                'duration', run.duration, ...
                'dist_change', run.dist_change, ...
                'r_dist_time', r_dist_time, ...
                'p_dist_time', p_dist_time, ...
                'r_dist_dff', r_dist_dff, ...
                'p_dist_dff', p_dist_dff, ...
                'r_time_dff', r_time_dff, ...
                'p_time_dff', p_time_dff, ...
                'r_dist_dff_partial', r_dist_dff_partial, ...
                'p_dist_dff_partial', p_dist_dff_partial, ...
                'r_time_dff_partial', r_time_dff_partial, ...
                'p_time_dff_partial', p_time_dff_partial);
            
            % Append to appropriate category
            results.(direction).(category).runs = [results.(direction).(category).runs; run_info];
            results.(direction).(category).dist_time_corr = [results.(direction).(category).dist_time_corr; r_dist_time];
            results.(direction).(category).dist_dff_corr = [results.(direction).(category).dist_dff_corr; r_dist_dff];
            results.(direction).(category).time_dff_corr = [results.(direction).(category).time_dff_corr; r_time_dff];
            results.(direction).(category).dist_dff_partial = [results.(direction).(category).dist_dff_partial; r_dist_dff_partial];
            results.(direction).(category).time_dff_partial = [results.(direction).(category).time_dff_partial; r_time_dff_partial];
        end
    end
    
    % Generate plots
    if options.plot_results
        plot_separated_runs_analysis(results, options);
    end
    
    % Print summary
    print_separated_runs_summary(results, options);
end

function plot_separated_runs_analysis(results, options)
    % Create comprehensive plots for the separated runs analysis
    
    directions = {'towards', 'away'};
    
    % Create main figure
    figure('Position', [100, 100, 1600, 900], ...
        'Name', sprintf('Distance-Time Correlation Analysis (threshold=%.2f, session=%s)', ...
        options.corr_threshold, num2str(options.session)));
    
    for d = 1:length(directions)
        direction = directions{d};
        
        % Get data for high and low correlation groups
        high = results.(direction).high_corr;
        low = results.(direction).low_corr;
        
        if isempty(high.runs) && isempty(low.runs)
            continue;
        end
        
        % Base subplot index for this direction
        base_idx = (d-1) * 4;
        
        % === Subplot 1: Distribution of distance-time correlations ===
        subplot(2, 4, base_idx + 1);
        
        all_dist_time = [high.dist_time_corr; low.dist_time_corr];
        histogram(all_dist_time, 20, 'FaceColor', [0.5 0.5 0.5]);
        hold on;
        
        % Add vertical line at threshold
        yl = ylim;
        plot([options.corr_threshold, options.corr_threshold], yl, 'r--', 'LineWidth', 2);
        plot([-options.corr_threshold, -options.corr_threshold], yl, 'r--', 'LineWidth', 2);
        
        xlabel('Distance-Time Correlation (r)');
        ylabel('Number of Runs');
        title(sprintf('%s: Distance-Time Correlation', [upper(direction(1)), direction(2:end)]));
        grid on;
        
        % Add text with counts
        n_high = length(high.dist_time_corr);
        n_low = length(low.dist_time_corr);
        text(0.05, 0.95, sprintf('High: n=%d\nLow: n=%d', n_high, n_low), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'EdgeColor', 'k');
        
        % === Subplot 2: Distance-dF/F correlation by group ===
        subplot(2, 4, base_idx + 2);
        
        if ~isempty(high.dist_dff_corr) && ~isempty(low.dist_dff_corr)
            % Calculate means and SEMs
            high_dist = high.dist_dff_corr(:);
            low_dist = low.dist_dff_corr(:);
            
            means = [mean(high_dist), mean(low_dist)];
            sems = [std(high_dist)/sqrt(length(high_dist)), std(low_dist)/sqrt(length(low_dist))];
            
            % Create bar plot
            b = bar(1:2, means, 'FaceColor', 'flat');
            b.CData(1,:) = [0.2, 0.6, 0.8];
            b.CData(2,:) = [0.8, 0.4, 0.2];
            hold on;
            
            % Add error bars
            errorbar(1:2, means, sems, 'k.', 'LineWidth', 1.5);
            
            % Add individual data points with jitter
            jitter_amount = 0.15;
            x1 = 1 + (rand(size(high_dist)) - 0.5) * jitter_amount;
            x2 = 2 + (rand(size(low_dist)) - 0.5) * jitter_amount;
            
            scatter(x1, high_dist, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
            scatter(x2, low_dist, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
            
            set(gca, 'XTick', 1:2, 'XTickLabel', {'High Dist-Time', 'Low Dist-Time'});
            ylabel('Distance-dF/F Correlation (r)');
            title(sprintf('%s: Distance Effect on dF/F', [upper(direction(1)), direction(2:end)]));
            grid on;
            plot(xlim, [0 0], 'k--');
            xlim([0.5 2.5]);
            
            % Add significance test
            [~, p_val] = ttest2(high_dist, low_dist);
            add_significance_star([1, 2], p_val);
        end
        
        % === Subplot 3: Time-dF/F correlation by group ===
        subplot(2, 4, base_idx + 3);
        
        if ~isempty(high.time_dff_corr) && ~isempty(low.time_dff_corr)
            % Calculate means and SEMs
            high_time = high.time_dff_corr(:);
            low_time = low.time_dff_corr(:);
            
            means = [mean(high_time), mean(low_time)];
            sems = [std(high_time)/sqrt(length(high_time)), std(low_time)/sqrt(length(low_time))];
            
            % Create bar plot
            b = bar(1:2, means, 'FaceColor', 'flat');
            b.CData(1,:) = [0.2, 0.6, 0.8];
            b.CData(2,:) = [0.8, 0.4, 0.2];
            hold on;
            
            % Add error bars
            errorbar(1:2, means, sems, 'k.', 'LineWidth', 1.5);
            
            % Add individual data points with jitter
            jitter_amount = 0.15;
            x1 = 1 + (rand(size(high_time)) - 0.5) * jitter_amount;
            x2 = 2 + (rand(size(low_time)) - 0.5) * jitter_amount;
            
            scatter(x1, high_time, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
            scatter(x2, low_time, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
            
            set(gca, 'XTick', 1:2, 'XTickLabel', {'High Dist-Time', 'Low Dist-Time'});
            ylabel('Time-dF/F Correlation (r)');
            title(sprintf('%s: Time Effect on dF/F', [upper(direction(1)), direction(2:end)]));
            grid on;
            plot(xlim, [0 0], 'k--');
            xlim([0.5 2.5]);
            
            % Add significance test
            [~, p_val] = ttest2(high_time, low_time);
            add_significance_star([1, 2], p_val);
        end
        
        % === Subplot 4: Partial correlations comparison ===
        subplot(2, 4, base_idx + 4);
        
        if ~isempty(high.runs) && ~isempty(low.runs)
            % Prepare data for grouped bar plot
            means_high = [mean(high.dist_dff_partial), mean(high.time_dff_partial)];
            sems_high = [std(high.dist_dff_partial)/sqrt(length(high.dist_dff_partial)), ...
                        std(high.time_dff_partial)/sqrt(length(high.time_dff_partial))];
            
            means_low = [mean(low.dist_dff_partial), mean(low.time_dff_partial)];
            sems_low = [std(low.dist_dff_partial)/sqrt(length(low.dist_dff_partial)), ...
                       std(low.time_dff_partial)/sqrt(length(low.time_dff_partial))];
            
            % Create grouped bar plot
            x = 1:2;
            width = 0.35;
            bar(x - width/2, means_high, width, 'FaceColor', [0.2 0.6 0.8]);
            hold on;
            bar(x + width/2, means_low, width, 'FaceColor', [0.8 0.4 0.2]);
            
            % Add error bars
            errorbar(x - width/2, means_high, sems_high, 'k.', 'LineWidth', 1.5);
            errorbar(x + width/2, means_low, sems_low, 'k.', 'LineWidth', 1.5);
            
            set(gca, 'XTick', x, 'XTickLabel', {'Distance', 'Time'});
            ylabel('Partial Correlation with dF/F');
            title(sprintf('%s: Partial Correlations', [upper(direction(1)), direction(2:end)]));
            legend({'High Dist-Time Corr', 'Low Dist-Time Corr'}, 'Location', 'best');
            grid on;
            hold on;
            plot(xlim, [0 0], 'k--');
        end
    end
    
    % Add overall title
    sgtitle(sprintf('Runs Separated by Distance-Time Correlation (Group=%s, Source=%s, Session=%s)', ...
        options.group, options.source, num2str(options.session)), 'FontWeight', 'bold');
    
    % Create second figure: Scatter plots
    figure('Position', [150, 150, 1400, 600], ...
        'Name', 'Distance vs Time Effects on dF/F');
    
    for d = 1:length(directions)
        direction = directions{d};
        
        high = results.(direction).high_corr;
        low = results.(direction).low_corr;
        
        if isempty(high.runs) && isempty(low.runs)
            continue;
        end
        
        % Subplot for high correlation runs
        subplot(2, 2, (d-1)*2 + 1);
        if ~isempty(high.runs)
            scatter(abs(high.dist_dff_corr), abs(high.time_dff_corr), 50, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            xlabel('|Distance-dF/F Correlation|');
            ylabel('|Time-dF/F Correlation|');
            title(sprintf('%s: High Dist-Time Corr (n=%d)', ...
                [upper(direction(1)), direction(2:end)], length(high.runs)));
            hold on;
            plot([0 1], [0 1], 'k--');
            axis equal;
            max_val = max([abs(high.dist_dff_corr); abs(high.time_dff_corr)]);
            xlim([0 min(1, max_val*1.1)]);
            ylim([0 min(1, max_val*1.1)]);
            grid on;
        end
        
        % Subplot for low correlation runs
        subplot(2, 2, (d-1)*2 + 2);
        if ~isempty(low.runs)
            scatter(abs(low.dist_dff_corr), abs(low.time_dff_corr), 50, 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
            xlabel('|Distance-dF/F Correlation|');
            ylabel('|Time-dF/F Correlation|');
            title(sprintf('%s: Low Dist-Time Corr (n=%d)', ...
                [upper(direction(1)), direction(2:end)], length(low.runs)));
            hold on;
            plot([0 1], [0 1], 'k--');
            axis equal;
            max_val = max([abs(low.dist_dff_corr); abs(low.time_dff_corr)]);
            xlim([0 min(1, max_val*1.1)]);
            ylim([0 min(1, max_val*1.1)]);
            grid on;
        end
    end
    
    sgtitle('Distance vs Time Effects: Separated by Distance-Time Correlation', ...
        'FontWeight', 'bold');
end

function add_significance_star(positions, p_val)
    % Add significance annotation to plot
    
    if isnan(p_val) || p_val >= 0.05
        return;
    end
    
    % Get y-axis limits
    yl = ylim;
    y_pos = yl(2) * 0.9;
    
    % Draw line
    plot(positions, [y_pos, y_pos], 'k-', 'LineWidth', 1.5);
    
    % Add stars
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    else
        stars = '*';
    end
    
    text(mean(positions), y_pos*1.05, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
end

function print_separated_runs_summary(results, options)
    % Print summary statistics
    
    fprintf('\n=== Runs Separated by Distance-Time Correlation ===\n');
    fprintf('Threshold: |r| >= %.2f\n', options.corr_threshold);
    fprintf('Session: %s\n', num2str(options.session));
    fprintf('Group: %s, Source: %s\n\n', options.group, options.source);
    
    directions = {'towards', 'away'};
    
    for d = 1:length(directions)
        direction = directions{d};
        
        fprintf('--- %s RUNS ---\n', upper(direction));
        
        high = results.(direction).high_corr;
        low = results.(direction).low_corr;
        
        % High correlation group
        fprintf('\n  HIGH Distance-Time Correlation (|r| >= %.2f):\n', options.corr_threshold);
        if ~isempty(high.runs)
            fprintf('    N runs: %d\n', length(high.runs));
            fprintf('    Distance-Time r: %.3f ± %.3f\n', ...
                mean(high.dist_time_corr), std(high.dist_time_corr));
            fprintf('    Distance-dF/F r: %.3f ± %.3f\n', ...
                mean(high.dist_dff_corr), std(high.dist_dff_corr));
            fprintf('    Time-dF/F r: %.3f ± %.3f\n', ...
                mean(high.time_dff_corr), std(high.time_dff_corr));
            fprintf('    Distance-dF/F partial r: %.3f ± %.3f\n', ...
                mean(high.dist_dff_partial), std(high.dist_dff_partial));
            fprintf('    Time-dF/F partial r: %.3f ± %.3f\n', ...
                mean(high.time_dff_partial), std(high.time_dff_partial));
        else
            fprintf('    No runs in this category\n');
        end
        
        % Low correlation group
        fprintf('\n  LOW Distance-Time Correlation (|r| < %.2f):\n', options.corr_threshold);
        if ~isempty(low.runs)
            fprintf('    N runs: %d\n', length(low.runs));
            fprintf('    Distance-Time r: %.3f ± %.3f\n', ...
                mean(low.dist_time_corr), std(low.dist_time_corr));
            fprintf('    Distance-dF/F r: %.3f ± %.3f\n', ...
                mean(low.dist_dff_corr), std(low.dist_dff_corr));
            fprintf('    Time-dF/F r: %.3f ± %.3f\n', ...
                mean(low.time_dff_corr), std(low.time_dff_corr));
            fprintf('    Distance-dF/F partial r: %.3f ± %.3f\n', ...
                mean(low.dist_dff_partial), std(low.dist_dff_partial));
            fprintf('    Time-dF/F partial r: %.3f ± %.3f\n', ...
                mean(low.time_dff_partial), std(low.time_dff_partial));
        else
            fprintf('    No runs in this category\n');
        end
        
        % Statistical comparison between groups
        if ~isempty(high.runs) && ~isempty(low.runs)
            fprintf('\n  Comparison (High vs Low):\n');
            
            [~, p_dist] = ttest2(high.dist_dff_corr, low.dist_dff_corr);
            [~, p_time] = ttest2(high.time_dff_corr, low.time_dff_corr);
            [~, p_dist_partial] = ttest2(high.dist_dff_partial, low.dist_dff_partial);
            [~, p_time_partial] = ttest2(high.time_dff_partial, low.time_dff_partial);
            
            fprintf('    Distance-dF/F: p = %.4f %s\n', p_dist, get_sig_stars(p_dist));
            fprintf('    Time-dF/F: p = %.4f %s\n', p_time, get_sig_stars(p_time));
            fprintf('    Distance-dF/F (partial): p = %.4f %s\n', p_dist_partial, get_sig_stars(p_dist_partial));
            fprintf('    Time-dF/F (partial): p = %.4f %s\n', p_time_partial, get_sig_stars(p_time_partial));
        end
        
        fprintf('\n');
    end
    
    fprintf('====================================================\n\n');
end

function sig_stars = get_sig_stars(p)
    % Get significance stars
    if p < 0.001
        sig_stars = '***';
    elseif p < 0.01
        sig_stars = '**';
    elseif p < 0.05
        sig_stars = '*';
    else
        sig_stars = 'ns';
    end
end