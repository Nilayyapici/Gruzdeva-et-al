function [data_table, stats_results] = plotSessionComparisonWithStats(mice_all, options)
    % Creates plots with session comparison bars and saves individual mouse data
    % 
    % OUTPUTS:
    %   data_table: Table with individual mouse values (one row per mouse per arm per session)
    %   stats_results: Statistical test results
    %
    % BARS SHOWN:
    %   BLACK bars = Session comparisons (e.g., Before vs Test within Food Arm)
    
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'plot_type'), options.plot_type = 'percent'; end
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'group'), options.group = 'all'; end
    if ~isfield(options, 'title'), options.title = 'Session Comparison Across Arms'; end
    if ~isfield(options, 'use_limits'), options.use_limits = false; end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 50; end
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end
    if ~isfield(options, 'horizontal'), options.horizontal = false; end
    if ~isfield(options, 'connect_points'), options.connect_points = true; end
    if ~isfield(options, 'memory_tracking'), options.memory_tracking = false; end
    if ~isfield(options, 'alpha'), options.alpha = 0.05; end
    if ~isfield(options, 'show_stats'), options.show_stats = true; end
    if ~isfield(options, 'colors'), options.colors = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.8 0.2]}; end
    if ~isfield(options, 'legend_location'), options.legend_location = 'northeast'; end
    if ~isfield(options, 'figure_size'), options.figure_size = [100, 100, 1000, 600]; end
    if ~isfield(options, 'font_size'), options.font_size = 12; end
    if ~isfield(options, 'bar_width'), options.bar_width = 0.8; end
    if ~isfield(options, 'ylim'), options.ylim = []; end
    if ~isfield(options, 'combine_nonfood_arms'), options.combine_nonfood_arms = false; end
    if ~isfield(options, 'save_table'), options.save_table = true; end
    if ~isfield(options, 'table_filename'), options.table_filename = 'individual_mouse_data.csv'; end
    
    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});
    
    % Analyze the data
    [zone_data, ~] = analyzeZoneTimes(mice_all, options);
    
    % Select data based on plot_type
    if strcmp(options.plot_type, 'percent')
        data_to_plot = zone_data.times_percent;
        y_label = 'Time (%)';
        data_units = 'percent';
    else
        data_to_plot = zone_data.times;
        y_label = 'Time (s)';
        data_units = 'seconds';
    end
    
    valid_sessions = intersect(fieldnames(data_to_plot)', options.sessions);
    if isempty(valid_sessions)
        error('No valid sessions found for plotting.');
    end
    
    % Process data to combine non-food arms if requested
    if options.combine_nonfood_arms
        data_to_plot = combineNonfoodArms(data_to_plot, valid_sessions);
        zone_data.mouse_ids = combineNonfoodArmMouseIds(zone_data.mouse_ids, valid_sessions);
    end
    
    % Determine which arms to include
    if options.combine_nonfood_arms
        if options.horizontal
            arms = {'Food Arm', 'Non-Food Arms', 'Horizontal'};
        else
            arms = {'Food Arm', 'Non-Food Arms'};
        end
    else
        if options.horizontal
            arms = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2', 'Horizontal'};
        else
            arms = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
        end
    end
    numArms = length(arms);
    
    % Get data for each session
    sessionData = cell(1, length(valid_sessions));
    sessionMeans = cell(1, length(valid_sessions));
    sessionSEMs = cell(1, length(valid_sessions));
    mouseIds = cell(1, length(valid_sessions));
    
    for i = 1:length(valid_sessions)
        sessionData{i} = data_to_plot.(valid_sessions{i});
        mouseIds{i} = zone_data.mouse_ids.(valid_sessions{i});
        if ~isempty(sessionData{i})
            sessionMeans{i} = mean(sessionData{i}, 1);
            sessionSEMs{i} = std(sessionData{i}, [], 1) / sqrt(size(sessionData{i}, 1));
        else
            sessionMeans{i} = zeros(1, numArms);
            sessionSEMs{i} = zeros(1, numArms);
        end
    end
    
    all_mice = unique([mouseIds{:}]);
    
    % CREATE DATA TABLE with individual mouse values
    data_table = createDataTable(sessionData, mouseIds, valid_sessions, arms, session_labels_map, data_units);
    
    % SAVE TABLE if requested
    if options.save_table
        writetable(data_table, options.table_filename);
        fprintf('\n=== DATA TABLE SAVED ===\n');
        fprintf('File: %s\n', options.table_filename);
        fprintf('Rows: %d (individual mouse measurements)\n', height(data_table));
        fprintf('Columns: %s\n', strjoin(data_table.Properties.VariableNames, ', '));
        fprintf('========================\n\n');
    end
    
    % Prepare data for grouped bar plot
    groupedData = zeros(numArms, length(valid_sessions));
    errorData = zeros(numArms, length(valid_sessions));
    
    for i = 1:numArms
        for j = 1:length(valid_sessions)
            if i <= length(sessionMeans{j})
                groupedData(i, j) = sessionMeans{j}(i);
                errorData(i, j) = sessionSEMs{j}(i);
            end
        end
    end
    
    % Perform statistical tests if requested
    stats_results = struct();
    if options.show_stats
        stats_results = performStatisticalTests(data_to_plot, valid_sessions, arms, options);
    end
    
    % Create figure
    figure('Position', options.figure_size);
    
    % Create grouped bar plot
    h = bar(groupedData, options.bar_width);
    
    % Set colors for each session bar
    for i = 1:length(valid_sessions)
        h(i).FaceColor = options.colors{mod(i-1, length(options.colors))+1};
        h(i).EdgeColor = 'none';
        h(i).FaceAlpha = 0.8;
    end
    
    hold on;
    
    % Add error bars and calculate bar positions
    numBars = length(valid_sessions);
    groupwidth = min(options.bar_width, numBars/(numBars+1.5));
    barPositions = zeros(numArms, numBars);
    
    for i = 1:numArms
        for j = 1:numBars
            if numBars == 1
                x = i;
            else
                x = i + (j - (numBars+1)/2) * groupwidth/numBars;
            end
            barPositions(i, j) = x;
            
            % Add error bars
            errorbar(x, groupedData(i, j), errorData(i, j), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
            
            % Add individual data points (these match the data_table values)
            if ~isempty(sessionData{j}) && i <= size(sessionData{j}, 2)
                indPoints = sessionData{j}(:, i);
                scatter(repmat(x, length(indPoints), 1), indPoints, 25, 'ko', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
            end
        end
    end
    
    % Connect individual mouse data points across sessions if requested
    if options.connect_points && length(valid_sessions) > 1
        connectMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, valid_sessions, numArms);
    end
    
    % Format plot BEFORE adding significance bars
    formatPlot(arms, valid_sessions, session_labels_map, options, y_label);
    
    % Add SESSION comparison bars (BLACK bars only)
    if options.show_stats && ~isempty(fieldnames(stats_results))
        addSessionSignificanceBars(stats_results, barPositions, groupedData, errorData, options, valid_sessions, arms);
    end
    
    % Display statistics summary if available
    if options.show_stats && isfield(stats_results, 'summary')
        displayStatsSummary(stats_results);
    end
    
    hold off;
end

function data_table = createDataTable(sessionData, mouseIds, valid_sessions, arms, session_labels_map, data_units)
    % Create a table with individual mouse data
    all_mouse_ids = {};
    all_sessions = {};
    all_arms = {};
    all_values = [];
    
    for sess_idx = 1:length(valid_sessions)
        session = valid_sessions{sess_idx};
        session_label = getSessionLabel(session, session_labels_map);
        
        if ~isempty(sessionData{sess_idx})
            for mouse_idx = 1:length(mouseIds{sess_idx})
                mouse_id = mouseIds{sess_idx}{mouse_idx};
                
                for arm_idx = 1:min(size(sessionData{sess_idx}, 2), length(arms))
                    all_mouse_ids{end+1} = mouse_id;
                    all_sessions{end+1} = session_label;
                    all_arms{end+1} = arms{arm_idx};
                    all_values(end+1) = sessionData{sess_idx}(mouse_idx, arm_idx);
                end
            end
        end
    end
    
    % Create table
    data_table = table(all_mouse_ids', all_sessions', all_arms', all_values', ...
        'VariableNames', {'MouseID', 'Session', 'Arm', ['Time_' data_units]});
    
    % Sort for readability
    data_table = sortrows(data_table, {'MouseID', 'Session', 'Arm'});
end

function combined_data = combineNonfoodArms(data_to_plot, valid_sessions)
    combined_data = struct();
    
    for i = 1:length(valid_sessions)
        session = valid_sessions{i};
        original_data = data_to_plot.(session);
        
        if ~isempty(original_data) && size(original_data, 2) >= 3
            combined_matrix = zeros(size(original_data, 1), size(original_data, 2) - 1);
            combined_matrix(:, 1) = original_data(:, 1);
            combined_matrix(:, 2) = (original_data(:, 2) + original_data(:, 3)) / 2;
            
            if size(original_data, 2) > 3
                combined_matrix(:, 3:end) = original_data(:, 4:end);
            end
            
            combined_data.(session) = combined_matrix;
        else
            combined_data.(session) = original_data;
        end
    end
end

function combined_mouse_ids = combineNonfoodArmMouseIds(mouse_ids, valid_sessions)
    combined_mouse_ids = mouse_ids;
end

function stats_results = performStatisticalTests(data_to_plot, valid_sessions, arms, options)
    % Perform statistical tests for both arms and sessions
    stats_results = struct();
    stats_results.test_type = 'Separate One-Way Tests';
    stats_results.arm_comparisons = struct();
    stats_results.session_comparisons = struct();
    stats_results.summary = struct();
    
    numArms = length(arms);
    numSessions = length(valid_sessions);
    
    fprintf('\n========== STATISTICAL ANALYSIS ==========\n');
    
    % 1. Compare arms within EACH session
    fprintf('\n--- ARM COMPARISONS (within sessions) ---\n');
    for s = 1:numSessions
        session = valid_sessions{s};
        if ~isempty(data_to_plot.(session)) && size(data_to_plot.(session), 2) >= numArms
            session_data = data_to_plot.(session)(:, 1:numArms);
            values = session_data(:);
            groups = repmat((1:numArms)', size(session_data, 1), 1);
            
            [p, ~, stats] = anova1(values, groups, 'off');
            
            stats_results.arm_comparisons.(session).p_value = p;
            stats_results.arm_comparisons.(session).significant = p < options.alpha;
            
            fprintf('%s: p = %.4f %s\n', session, p, getSignificanceLabel(p < options.alpha));
            
            if p < options.alpha
                [c, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                stats_results.arm_comparisons.(session).posthoc = c;
                
                for i = 1:size(c, 1)
                    if c(i, 6) < options.alpha
                        fprintf('  %s vs %s: p = %.4f %s\n', ...
                            arms{c(i,1)}, arms{c(i,2)}, c(i,6), getStarSymbol(c(i,6)));
                    end
                end
            else
                stats_results.arm_comparisons.(session).posthoc = [];
            end
        end
    end
    
    % 2. Compare sessions within EACH arm (these get bars)
    fprintf('\n--- SESSION COMPARISONS (within arms) ---\n');
    for a = 1:numArms
        arm_name = arms{a};
        arm_field_name = matlab.lang.makeValidName(arm_name);
        arm_data = [];
        session_groups = [];
        
        for s = 1:numSessions
            session = valid_sessions{s};
            if ~isempty(data_to_plot.(session)) && size(data_to_plot.(session), 2) >= a
                arm_values = data_to_plot.(session)(:, a);
                arm_data = [arm_data; arm_values];
                session_groups = [session_groups; s * ones(size(arm_values))];
            end
        end
        
        if length(unique(session_groups)) >= 2
            if length(unique(session_groups)) == 2
                group1_data = arm_data(session_groups == 1);
                group2_data = arm_data(session_groups == 2);
                [~, p] = ttest2(group1_data, group2_data);
                
                stats_results.session_comparisons.(arm_field_name).p_value = p;
                stats_results.session_comparisons.(arm_field_name).significant = p < options.alpha;
                stats_results.session_comparisons.(arm_field_name).test_type = 't-test';
                stats_results.session_comparisons.(arm_field_name).arm_name = arm_name;
                
                fprintf('%s: p = %.4f %s\n', arm_name, p, getSignificanceLabel(p < options.alpha));
                if p < options.alpha
                    fprintf('  %s vs %s: p = %.4f %s\n', ...
                        valid_sessions{1}, valid_sessions{2}, p, getStarSymbol(p));
                end
            else
                [p, ~, stats] = anova1(arm_data, session_groups, 'off');
                
                stats_results.session_comparisons.(arm_field_name).p_value = p;
                stats_results.session_comparisons.(arm_field_name).significant = p < options.alpha;
                stats_results.session_comparisons.(arm_field_name).test_type = 'ANOVA';
                stats_results.session_comparisons.(arm_field_name).arm_name = arm_name;
                
                fprintf('%s: p = %.4f %s\n', arm_name, p, getSignificanceLabel(p < options.alpha));
                
                if p < options.alpha
                    [c, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                    stats_results.session_comparisons.(arm_field_name).posthoc = c;
                    
                    for i = 1:size(c, 1)
                        if c(i, 6) < options.alpha
                            fprintf('  %s vs %s: p = %.4f %s\n', ...
                                valid_sessions{c(i,1)}, valid_sessions{c(i,2)}, c(i,6), getStarSymbol(c(i,6)));
                        end
                    end
                end
            end
        end
    end
    
    fprintf('\n==========================================\n');
    
    % Count significant comparisons
    stats_results.summary.significant_arm_comparisons = 0;
    stats_results.summary.significant_session_comparisons = 0;
    
    session_fields = fieldnames(stats_results.arm_comparisons);
    for i = 1:length(session_fields)
        if stats_results.arm_comparisons.(session_fields{i}).significant
            stats_results.summary.significant_arm_comparisons = stats_results.summary.significant_arm_comparisons + 1;
        end
    end
    
    arm_fields = fieldnames(stats_results.session_comparisons);
    for i = 1:length(arm_fields)
        if stats_results.session_comparisons.(arm_fields{i}).significant
            stats_results.summary.significant_session_comparisons = stats_results.summary.significant_session_comparisons + 1;
        end
    end
end

function addSessionSignificanceBars(stats_results, barPositions, groupedData, errorData, options, valid_sessions, arms)
    % Add BLACK bars for SESSION comparisons only
    
    current_ylim = ylim;
    data_range = current_ylim(2) - current_ylim(1);
    star_height_offset = data_range * 0.04;
    
    if ~isfield(stats_results, 'session_comparisons')
        return;
    end
    
    arm_fields = fieldnames(stats_results.session_comparisons);
    
    for a = 1:length(arm_fields)
        arm_field = arm_fields{a};
        
        if stats_results.session_comparisons.(arm_field).significant
            arm_name = stats_results.session_comparisons.(arm_field).arm_name;
            arm_idx = find(strcmp(arms, arm_name));
            
            if ~isempty(arm_idx)
                % Check if we have post-hoc results with significant comparisons
                if isfield(stats_results.session_comparisons.(arm_field), 'posthoc') && ...
                   ~isempty(stats_results.session_comparisons.(arm_field).posthoc)
                    
                    posthoc = stats_results.session_comparisons.(arm_field).posthoc;
                    sig_comparisons = posthoc(posthoc(:, 6) < options.alpha, :);
                    
                    if ~isempty(sig_comparisons)
                        % Draw bar for each significant comparison
                        for i = 1:size(sig_comparisons, 1)
                            sess1_idx = sig_comparisons(i, 1);
                            sess2_idx = sig_comparisons(i, 2);
                            p_val = sig_comparisons(i, 6);
                            
                            x1 = barPositions(arm_idx, sess1_idx);
                            x2 = barPositions(arm_idx, sess2_idx);
                            
                            max_y = max(groupedData(arm_idx, [sess1_idx, sess2_idx]) + ...
                                       errorData(arm_idx, [sess1_idx, sess2_idx]));
                            
                            drawSignificanceBar(x1, x2, max_y, p_val, star_height_offset * (i + 0.5), ...
                                [0 0 0], options.font_size);
                        end
                    end
                elseif strcmp(stats_results.session_comparisons.(arm_field).test_type, 't-test')
                    % For t-test (2 sessions only)
                    p_val = stats_results.session_comparisons.(arm_field).p_value;
                    
                    x1 = barPositions(arm_idx, 1);
                    x2 = barPositions(arm_idx, 2);
                    max_y = max(groupedData(arm_idx, :) + errorData(arm_idx, :));
                    
                    drawSignificanceBar(x1, x2, max_y, p_val, star_height_offset * 2, ...
                        [0 0 0], options.font_size);
                end
            end
        end
    end
end

function drawSignificanceBar(x1, x2, base_height, p_val, offset, color, font_size)
    % Draw horizontal bar with star
    star_symbol = getStarSymbol(p_val);
    line_y = base_height + offset;
    
    % Horizontal line
    plot([x1, x2], [line_y, line_y], 'Color', color, 'LineWidth', 2, 'LineStyle', '-');
    
    % Vertical caps
    cap_height = offset * 0.2;
    plot([x1, x1], [line_y - cap_height, line_y], 'Color', color, 'LineWidth', 2);
    plot([x2, x2], [line_y - cap_height, line_y], 'Color', color, 'LineWidth', 2);
    
    % Star
    text(mean([x1, x2]), line_y + offset*0.4, star_symbol, ...
        'HorizontalAlignment', 'center', 'FontSize', font_size + 2, ...
        'FontWeight', 'bold', 'Color', color);
end

function connectMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, valid_sessions, numArms)
    for arm_idx = 1:numArms
        for m = 1:length(all_mice)
            current_mouse = all_mice{m};
            
            mouse_x = [];
            mouse_y = [];
            
            for sess_idx = 1:length(valid_sessions)
                mouse_in_session = find(strcmp(mouseIds{sess_idx}, current_mouse));
                
                if ~isempty(mouse_in_session) && arm_idx <= size(sessionData{sess_idx}, 2)
                    x_pos = barPositions(arm_idx, sess_idx);
                    y_val = sessionData{sess_idx}(mouse_in_session, arm_idx);
                    
                    mouse_x = [mouse_x, x_pos];
                    mouse_y = [mouse_y, y_val];
                end
            end
            
            if length(mouse_x) >= 2
                plot(mouse_x, mouse_y, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5 0.3]);
            end
        end
    end
end

function session_label = getSessionLabel(session, session_labels_map)
    if nargin < 2
        session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});
    end
    
    if isKey(session_labels_map, session)
        session_label = session_labels_map(session);
    else
        session_label = strrep(session, 'sess', 'Session ');
    end
end

function star_symbol = getStarSymbol(p_val)
    if p_val < 0.001
        star_symbol = '***';
    elseif p_val < 0.01
        star_symbol = '**';
    elseif p_val < 0.05
        star_symbol = '*';
    else
        star_symbol = 'ns';
    end
end

function label = getSignificanceLabel(is_significant)
    if is_significant
        label = '(significant)';
    else
        label = '(not significant)';
    end
end

function formatPlot(arms, valid_sessions, session_labels_map, options, y_label)
    ylabel(y_label, 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    formatted_group_name = formatGroupName(options.group);
    
    plot_title = options.title;
    if options.combine_nonfood_arms
        plot_title = [plot_title, ' (Non-Food Combined)'];
    end
    title([plot_title, ' - ', formatted_group_name], 'FontSize', options.font_size + 4, 'FontWeight', 'bold');
    
    xticks(1:length(arms));
    xticklabels(arms);
    set(gca, 'FontSize', options.font_size);
    
    sessionLabels = cell(1, length(valid_sessions));
    for i = 1:length(valid_sessions)
        if isKey(session_labels_map, valid_sessions{i})
            sessionLabels{i} = session_labels_map(valid_sessions{i});
        else
            sessionLabels{i} = strrep(valid_sessions{i}, 'sess', 'Session ');
        end
    end
    
    legend(sessionLabels, 'Location', options.legend_location, 'FontSize', options.font_size);
    legend('boxoff');

    if ~isempty(options.ylim)
        ylim(options.ylim);
    elseif strcmp(options.plot_type, 'percent')
        max_val = max(groupedData(:) + errorData(:));
        if max_val > 75
            ylim([0, 120]);  % Extra space for bars
        else
            ylim([0, 100]);
        end
    else
        current_ylim = ylim;
        ylim([current_ylim(1), current_ylim(2) * 1.4]);  % 40% extra for bars
    end

    box off;
    grid off;
    set(gca, 'LineWidth', 1.2);
    set(gcf, 'Color', 'white');
    set(gca, 'Color', 'white');
end

function formatted_name = formatGroupName(group_name)
    if isempty(group_name)
        formatted_name = '';
        return;
    end
    
    formatted_name = strrep(group_name, '_', ' ');
    words = strsplit(formatted_name, ' ');
    for i = 1:length(words)
        if ~isempty(words{i})
            words{i} = [upper(words{i}(1)), lower(words{i}(2:end))];
        end
    end
    formatted_name = strjoin(words, ' ');
end

function displayStatsSummary(stats_results)
    fprintf('\n===== Statistical Summary =====\n');
    fprintf('Test Type: %s\n', stats_results.test_type);
    fprintf('Significant arm comparisons: %d\n', stats_results.summary.significant_arm_comparisons);
    fprintf('Significant session comparisons: %d (shown as BLACK bars)\n', stats_results.summary.significant_session_comparisons);
    fprintf('===============================\n');
end