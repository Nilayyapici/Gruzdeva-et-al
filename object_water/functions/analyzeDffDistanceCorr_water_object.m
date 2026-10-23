function [corr_data] = analyzeDffDistanceCorr_water_object(mice_all, options)
    % ANALYZEDFF_DISTANCE_CORRELATIONS_WATER_OBJECT Analyzes correlations between dF/F and 
    % distances to water and object across different sessions, especially focusing
    % on post-discovery periods
    %
    % Parameters:
    %   mice_all - Cell array containing mouse data
    %   options - Struct with analysis parameters:
    %     - filter_order_dist: Order of Butterworth filter for distance (default: 1)
    %     - cutoff_freq_dist: Cutoff frequency for distance filter (default: 0.08)
    %     - filter_order_dff: Order of Butterworth filter for dF/F (default: 1)
    %     - cutoff_freq_dff: Cutoff frequency for dF/F filter (default: 0.05)
    %     - dist_lim: Upper limit for distance values (default: 150)
    %     - dist_too_close: Lower limit for distance values (default: 0)
    %     - remove_grooming: Whether to exclude grooming periods (default: true)
    %     - remove_climbing: Whether to exclude climbing periods (default: true)
    %     - title: Title for the plot (default: 'dF/F and Distance Correlations')
    
    % Default options
    if nargin < 2
        options = struct();
    end
    
    % Basic options
    if ~isfield(options, 'title'), options.title = 'dF/F and Distance Correlations'; end
    if ~isfield(options, 'ylimit'), options.ylimit = [-0.4, 0.4]; end
    
    % Distance limits
    if ~isfield(options, 'dist_lim'), options.dist_lim = 150; end  % for positions as far as this limit
    if ~isfield(options, 'dist_too_close'), options.dist_too_close = 0; end  % for removing too close positions
    
    % Filter parameters
    if ~isfield(options, 'filter_order_dist'), options.filter_order_dist = 1; end
    if ~isfield(options, 'cutoff_freq_dist'), options.cutoff_freq_dist = 0.08; end
    if ~isfield(options, 'filter_order_dff'), options.filter_order_dff = 1; end
    if ~isfield(options, 'cutoff_freq_dff'), options.cutoff_freq_dff = 0.05; end
    
    % Filtering options
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'remove_climbing'), options.remove_climbing = true; end
    
    % Column indices
    col_time = 1;
    col_dff = 6;
    col_dist_water = 9;
    col_dist_object = 10;
    col_grooming = 13;
    col_climbing = 14;
    
    % Initialize correlation data structure
    corr_data = struct();
    session_types = {'before', 'water', 'object'};
    
    % Initialize session data
    for s = 1:length(session_types)
        session = session_types{s};
        corr_data.(session) = struct();
        corr_data.(session).water_corr = [];
        corr_data.(session).object_corr = [];
        corr_data.(session).water_pval = [];
        corr_data.(session).object_pval = [];
        corr_data.(session).mouse_ids = {};
    end
    
    % Process each session
    for i = 1:size(mice_all, 1)
        % Get session info
        session_id = mice_all{i, 1};
        
        % Skip if not a character
        if ~ischar(session_id)
            continue;
        end
        
        % Determine session type
        if contains(session_id, 'before')
            session_type = 'before';
        elseif contains(session_id, 'water')
            session_type = 'water';
        elseif contains(session_id, 'object')
            session_type = 'object';
        else
            continue; % Skip unknown session types
        end
        
        % Get data and discovery info
        data = mice_all{i, 4};
        if size(mice_all, 2) >= 6 && ~isempty(mice_all{i, 6})
            discovery = mice_all{i, 6};
        else
            discovery = [NaN, NaN];
        end
        
        % Extract mouse ID (everything before first underscore)
        parts = strsplit(session_id, '_');
        mouse_id = parts{1};
        
        % Skip if empty data
        if isempty(data)
            continue;
        end
        
        % Filter data based on discovery time for water and object sessions
        if strcmp(session_type, 'water') && length(discovery) >= 1 && ~isnan(discovery(1))
            discovery_idx = discovery(1);
            discovery_time = data(discovery_idx, col_time);
            data = data(data(:, col_time) >= discovery_time, :);
            fprintf('Session %s: Using data after water discovery at time %.2f\n', session_id, discovery_time);
        elseif strcmp(session_type, 'object') && length(discovery) >= 2 && ~isnan(discovery(2))
            discovery_idx = discovery(2);
            discovery_time = data(discovery_idx, col_time);
            data = data(data(:, col_time) >= discovery_time, :);
            fprintf('Session %s: Using data after object discovery at time %.2f\n', session_id, discovery_time);
        end
        
        % Skip if no data after filtering by discovery time
        if isempty(data)
            continue;
        end
        
        % Apply behavioral filtering (remove grooming and climbing)
        mask = true(size(data, 1), 1);
        if options.remove_grooming && size(data, 2) >= col_grooming
            mask = mask & (data(:, col_grooming) == 0);
        end
        if options.remove_climbing && size(data, 2) >= col_climbing
            mask = mask & (data(:, col_climbing) == 0);
        end
        data = data(mask, :);
        
        % Skip if no data after behavioral filtering
        if isempty(data)
            continue;
        end
        
        % Apply distance limits
        water_mask = (data(:, col_dist_water) > options.dist_too_close) & ...
                     (data(:, col_dist_water) <= options.dist_lim);
        object_mask = (data(:, col_dist_object) > options.dist_too_close) & ...
                      (data(:, col_dist_object) <= options.dist_lim);
        
        % Calculate correlations
        if sum(water_mask) > 5
            [water_corr, water_pval] = corr(data(water_mask, col_dist_water), ...
                                           data(water_mask, col_dff), 'Type', 'Pearson');
        else
            water_corr = NaN;
            water_pval = NaN;
        end
        
        if sum(object_mask) > 5
            [object_corr, object_pval] = corr(data(object_mask, col_dist_object), ...
                                             data(object_mask, col_dff), 'Type', 'Pearson');
        else
            object_corr = NaN;
            object_pval = NaN;
        end
        
        % Store results
        corr_data.(session_type).mouse_ids{end+1} = mouse_id;
        corr_data.(session_type).water_corr(end+1) = water_corr;
        corr_data.(session_type).object_corr(end+1) = object_corr;
        corr_data.(session_type).water_pval(end+1) = water_pval;
        corr_data.(session_type).object_pval(end+1) = object_pval;
        
        % Print summary
        fprintf('Mouse %s, session %s: Water corr=%.3f (p=%.3f), Object corr=%.3f (p=%.3f)\n', ...
            mouse_id, session_type, water_corr, water_pval, object_corr, object_pval);
    end
    
    % Plot results
    plotDistanceCorrelations(corr_data, options);
    
    % Return results
    return;
end

function z = fisher_z(r)
    % Fisher's z-transformation for correlation coefficients
    z = NaN(size(r));
    valid = ~isnan(r) & abs(r) < 1; % Avoid edge cases where |r| = 1di
    z(valid) = 0.5 * log((1 + r(valid)) ./ (1 - r(valid)));
end

function plotDistanceCorrelations(corr_data, options)
    % Plot correlation results with individual mouse data points
    
 % Define session types and colors
    session_types = {'before', 'water', 'object'};
    water_color = [0.3, 0.7, 0.9]; % Light blue
    object_color = [0.8, 0.6, 0]; % Brown-yellow
    
    % Create figure
    figure('Position', [100, 100, 1000, 600]);
    hold on;
    
    % Bar configuration
    bar_width = 0.2;      % Width of each bar
    group_spacing = 1.5;   % Spacing between groups
    bar_gap = 0.15;         % Gap between bars within a group
    
    % Calculate group and bar positions
    group_centers = zeros(1, length(session_types));
    for i = 1:length(session_types)
        group_centers(i) = (i-1) * group_spacing + 1;
    end
    
    % Water and object bar centers
    water_centers = group_centers - bar_width/2 - bar_gap/2;
    object_centers = group_centers + bar_width/2 + bar_gap/2;
    
    % Compute means for each group
    water_means = zeros(1, length(session_types));
    object_means = zeros(1, length(session_types));
    water_sems = zeros(1, length(session_types));
    object_sems = zeros(1, length(session_types));
    water_n = zeros(1, length(session_types));
    object_n = zeros(1, length(session_types));
    
    for i = 1:length(session_types)
        type = session_types{i};
        
        % Get correlations for this session type
        water_corrs = corr_data.(type).water_corr;
        object_corrs = corr_data.(type).object_corr;
        
        % Apply Fisher z-transformation for averaging
        water_z = fisher_z(water_corrs);
        object_z = fisher_z(object_corrs);
        
        % Compute means and SEMs in z-space
        water_mean_z = nanmean(water_z);
        object_mean_z = nanmean(object_z);
        water_sem_z = nanstd(water_z) / sqrt(sum(~isnan(water_z)));
        object_sem_z = nanstd(object_z) / sqrt(sum(~isnan(object_z)));
        
        % Convert back to correlation space for plotting
        water_mean = tanh(water_mean_z);
        object_mean = tanh(object_mean_z);
        water_sem_upper = tanh(water_mean_z + water_sem_z) - water_mean;
        water_sem_lower = water_mean - tanh(water_mean_z - water_sem_z);
        object_sem_upper = tanh(object_mean_z + object_sem_z) - object_mean;
        object_sem_lower = object_mean - tanh(object_mean_z - object_sem_z);
        
        % Use average of upper and lower SEM for error bars
        water_sem = (water_sem_upper + water_sem_lower) / 2;
        object_sem = (object_sem_upper + object_sem_lower) / 2;
        
        % Store for plotting
        water_means(i) = water_mean;
        object_means(i) = object_mean;
        water_sems(i) = water_sem;
        object_sems(i) = object_sem;
        water_n(i) = sum(~isnan(water_corrs));
        object_n(i) = sum(~isnan(object_corrs));
    end
    
    % Plot water bars
    water_bars = bar(water_centers, water_means, bar_width, 'FaceColor', water_color);
    % Plot object bars
    object_bars = bar(object_centers, object_means, bar_width, 'FaceColor', object_color);
    
    % Add error bars
    errorbar(water_centers, water_means, water_sems, 'k.', 'LineWidth', 1);
    errorbar(object_centers, object_means, object_sems, 'k.', 'LineWidth', 1);
    
    % Add individual points per mouse
    for i = 1:length(session_types)
        type = session_types{i};
        
        % Get data for this session type
        water_corrs = corr_data.(type).water_corr;
        object_corrs = corr_data.(type).object_corr;
        
        % Plot individual points (jittered)
        for j = 1:length(water_corrs)
            if ~isnan(water_corrs(j))
                jitter = (rand - 0.5) * bar_width * 0.2;
                plot(water_centers(i) + jitter, water_corrs(j), 'o', ...
                    'MarkerFaceColor', water_color * 0.8, ...
                    'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            end
        end
        
        for j = 1:length(object_corrs)
            if ~isnan(object_corrs(j))
                jitter = (rand - 0.5) * bar_width * 0.2;
                plot(object_centers(i) + jitter, object_corrs(j), 'o', ...
                    'MarkerFaceColor', object_color * 0.8, ...
                    'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            end
        end
    end
    
    % Add a horizontal line at y=0
    plot([0.5, max(group_centers) + 0.5], [0, 0], 'k--', 'LineWidth', 1);
    
    % Format plot
    % xlabel('Session Type');
    ylabel('Correlation (r) with dF/F');
    title(options.title);
    
    % Set x-ticks at group centers
    xticks(group_centers);
    xticklabels(session_types);
    
    % Remove grid and set axis limits
    grid off;
    ylim(options.ylimit);
    xlim([group_centers(1) - 1, group_centers(end) + 1]);
    
    % Add legend
    legend([water_bars, object_bars], {'to Water', 'to Object'});
    legend('boxoff')
    
    % Add significance markers
    for i = 1:length(session_types)
        type = session_types{i};
        water_corrs = corr_data.(type).water_corr;
        object_corrs = corr_data.(type).object_corr;
        
        % T-test against zero
        if sum(~isnan(water_corrs)) >= 3
            [~, p_water] = ttest(water_corrs);
            add_significance_marker(water_centers(i), water_means(i), water_sems(i), p_water);
        end
        
        if sum(~isnan(object_corrs)) >= 3
            [~, p_object] = ttest(object_corrs);
            add_significance_marker(object_centers(i), object_means(i), object_sems(i), p_object);
        end
    end

    
    % Print statistics
    fprintf('\nCorrelation Statistics:\n');
    for i = 1:length(session_types)
        type = session_types{i};
        fprintf('%s session:\n', upper(type(1)) + type(2:end));
        fprintf('  Water correlations: Mean = %.3f, SEM = %.3f, n = %d\n', ...
            water_means(i), water_sems(i), water_n(i));
        fprintf('  Object correlations: Mean = %.3f, SEM = %.3f, n = %d\n', ...
            object_means(i), object_sems(i), object_n(i));
        
        % Print statistical tests
        water_corrs = corr_data.(type).water_corr;
        object_corrs = corr_data.(type).object_corr;
        
        if sum(~isnan(water_corrs)) >= 3
            [~, p_water] = ttest(water_corrs);
            fprintf('  Water correlation vs. zero: p = %.4f\n', p_water);
        end
        
        if sum(~isnan(object_corrs)) >= 3
            [~, p_object] = ttest(object_corrs);
            fprintf('  Object correlation vs. zero: p = %.4f\n', p_object);
        end
        
        % Compare water vs object if enough paired data
        valid_idx = ~isnan(water_corrs) & ~isnan(object_corrs);
        if sum(valid_idx) >= 3
            [~, p_diff] = ttest(water_corrs(valid_idx), object_corrs(valid_idx));
            fprintf('  Water vs. Object correlation: p = %.4f\n', p_diff);
        end
    end
    
    % Add ANOVA comparisons between session types
    fprintf('\n---------- ANOVA Comparisons ----------\n');
    
    % Prepare data for ANOVA tests
    % For water correlations
    water_corrs_before = corr_data.before.water_corr;
    water_corrs_water = corr_data.water.water_corr;
    water_corrs_object = corr_data.object.water_corr;
    
    % For object correlations
    object_corrs_before = corr_data.before.object_corr;
    object_corrs_water = corr_data.water.object_corr;
    object_corrs_object = corr_data.object.object_corr;
    
    % Find mice that have data in all three sessions
    all_mice = unique([corr_data.before.mouse_ids, corr_data.water.mouse_ids, corr_data.object.mouse_ids]);
    
    % 1. Water Correlations ANOVA across sessions
    fprintf('ANOVA: Water Correlations across Sessions\n');
    water_data = [];
    water_groups = [];
    
    for i = 1:length(water_corrs_before)
        if ~isnan(water_corrs_before(i))
            water_data = [water_data; water_corrs_before(i)];
            water_groups = [water_groups; 1]; % 1 = before
        end
    end
    
    for i = 1:length(water_corrs_water)
        if ~isnan(water_corrs_water(i))
            water_data = [water_data; water_corrs_water(i)];
            water_groups = [water_groups; 2]; % 2 = water
        end
    end
    
    for i = 1:length(water_corrs_object)
        if ~isnan(water_corrs_object(i))
            water_data = [water_data; water_corrs_object(i)];
            water_groups = [water_groups; 3]; % 3 = object
        end
    end
    
    % Run ANOVA on water correlations
    if length(unique(water_groups)) > 1 && length(water_data) > 5
        [p_water_anova, tbl_water, stats_water] = anova1(water_data, water_groups, 'off');
        fprintf('  One-way ANOVA p-value: %.4f\n', p_water_anova);
        
        % Post-hoc tests (if ANOVA is significant)
        if p_water_anova < 0.05
            fprintf('  Post-hoc multiple comparisons:\n');
            [c_water, m_water, h_water, gnames_water] = multcompare(stats_water, 'Display', 'off');
            comparisons = {'Before vs Water', 'Before vs Object', 'Water vs Object'};
            for j = 1:size(c_water, 1)
                fprintf('    %s: p = %.4f\n', comparisons{j}, c_water(j, 6));
            end
        end
    else
        fprintf('  Not enough data for ANOVA on water correlations\n');
    end
    
    % 2. Object Correlations ANOVA across sessions
    fprintf('\nANOVA: Object Correlations across Sessions\n');
    object_data = [];
    object_groups = [];
    
    for i = 1:length(object_corrs_before)
        if ~isnan(object_corrs_before(i))
            object_data = [object_data; object_corrs_before(i)];
            object_groups = [object_groups; 1]; % 1 = before
        end
    end
    
    for i = 1:length(object_corrs_water)
        if ~isnan(object_corrs_water(i))
            object_data = [object_data; object_corrs_water(i)];
            object_groups = [object_groups; 2]; % 2 = water
        end
    end
    
    for i = 1:length(object_corrs_object)
        if ~isnan(object_corrs_object(i))
            object_data = [object_data; object_corrs_object(i)];
            object_groups = [object_groups; 3]; % 3 = object
        end
    end
    
    % Run ANOVA on object correlations
    if length(unique(object_groups)) > 1 && length(object_data) > 5
        [p_object_anova, tbl_object, stats_object] = anova1(object_data, object_groups, 'off');
        fprintf('  One-way ANOVA p-value: %.4f\n', p_object_anova);
        
        % Post-hoc tests (if ANOVA is significant)
        if p_object_anova < 0.05
            fprintf('  Post-hoc multiple comparisons:\n');
            [c_object, m_object, h_object, gnames_object] = multcompare(stats_object, 'Display', 'off');
            comparisons = {'Before vs Water', 'Before vs Object', 'Water vs Object'};
            for j = 1:size(c_object, 1)
                fprintf('    %s: p = %.4f\n', comparisons{j}, c_object(j, 6));
            end
        end
    else
        fprintf('  Not enough data for ANOVA on object correlations\n');
    end
    
    % 3. Water vs Object within each session
    fprintf('\nWater vs Object Correlations Within Each Session\n');
    
    % Before session
    water_object_before = [];
    for i = 1:length(water_corrs_before)
        if i <= length(object_corrs_before) && ~isnan(water_corrs_before(i)) && ~isnan(object_corrs_before(i))
            water_object_before = [water_object_before; water_corrs_before(i), object_corrs_before(i)];
        end
    end
    
    if size(water_object_before, 1) >= 3
        [~, p_before] = ttest(water_object_before(:,1), water_object_before(:,2));
        fprintf('  Before session (Water vs Object): p = %.4f\n', p_before);
    else
        fprintf('  Not enough paired data for Before session comparison\n');
    end
    
    % Water session
    water_object_water = [];
    for i = 1:length(water_corrs_water)
        if i <= length(object_corrs_water) && ~isnan(water_corrs_water(i)) && ~isnan(object_corrs_water(i))
            water_object_water = [water_object_water; water_corrs_water(i), object_corrs_water(i)];
        end
    end
    
    if size(water_object_water, 1) >= 3
        [~, p_water] = ttest(water_object_water(:,1), water_object_water(:,2));
        fprintf('  Water session (Water vs Object): p = %.4f\n', p_water);
    else
        fprintf('  Not enough paired data for Water session comparison\n');
    end
    
    % Object session
    water_object_object = [];
    for i = 1:length(water_corrs_object)
        if i <= length(object_corrs_object) && ~isnan(water_corrs_object(i)) && ~isnan(object_corrs_object(i))
            water_object_object = [water_object_object; water_corrs_object(i), object_corrs_object(i)];
        end
    end
    
    if size(water_object_object, 1) >= 3
        [~, p_object] = ttest(water_object_object(:,1), water_object_object(:,2));
        fprintf('  Object session (Water vs Object): p = %.4f\n', p_object);
    else
        fprintf('  Not enough paired data for Object session comparison\n');
    end
end

function add_significance_marker(x, y, sem, p)
    % Add significance markers above bar
    offset = 1.8 * sem;
    if sem < 0.01
        offset = 0.05;
    end
    
    if p < 0.001
        text(x, y + offset, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
    elseif p < 0.01
        text(x, y + offset, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
    elseif p < 0.05
        text(x, y + offset, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end