function plotConditionComparisonWithStats(mice_all, options)
    % Function to compare different conditions across sessions
    % Creates separate figures for each arm (Food Arm, Non-Food Arm 1, Non-Food Arm 2)
    % NOW: Each figure shows CONDITIONS on x-axis with bars for different SESSIONS
    
    % Set defaults for options
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'plot_type'), options.plot_type = 'percent'; end
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'conditions'), error('Must specify conditions to compare in options.conditions'); end
    if ~isfield(options, 'use_limits'), options.use_limits = false; end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 50; end
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end
    if ~isfield(options, 'horizontal'), options.horizontal = false; end
    if ~isfield(options, 'connect_points'), options.connect_points = true; end
    if ~isfield(options, 'alpha'), options.alpha = 0.05; end
    if ~isfield(options, 'show_stats'), options.show_stats = true; end
    if ~isfield(options, 'session_colors')
        % Default colors for sessions (now used instead of condition colors)
        options.session_colors = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.8 0.2], [0.8 0.6 0.2]};
    end
    if ~isfield(options, 'legend_location'), options.legend_location = 'northeast'; end
    if ~isfield(options, 'figure_size'), options.figure_size = [100, 100, 800, 600]; end
    if ~isfield(options, 'font_size'), options.font_size = 12; end
    if ~isfield(options, 'bar_width'), options.bar_width = 0.9; end
    if ~isfield(options, 'ylim'), options.ylim = []; end
    if ~isfield(options, 'combine_nonfood_arms'), options.combine_nonfood_arms = false; end
    if ~isfield(options, 'memory_tracking'), options.memory_tracking = 0; end
    
    % Session labels mapping
    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});
    
    % Define arm names based on whether we're combining non-food arms
    if options.combine_nonfood_arms
        if options.horizontal
            arm_names = {'Food Arm', 'Non-Food Arms Combined', 'Horizontal'};
        else
            arm_names = {'Food Arm', 'Non-Food Arms Combined'};
        end
    else
        if options.horizontal
            arm_names = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2', 'Horizontal'};
        else
            arm_names = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
        end
    end
    
    % Get data for each condition
    condition_data = struct();
    for c = 1:length(options.conditions)
        condition_name = options.conditions{c};
        temp_options = options;
        temp_options.group = condition_name;
        temp_options.memory_tracking = 0;  % Explicitly set memory_tracking to 0
        
        [zone_data, ~] = analyzeZoneTimes(mice_all, temp_options);
        
        % Select data based on plot_type
        if strcmp(options.plot_type, 'percent')
            condition_data.(condition_name).data = zone_data.times_percent;
            y_label = 'Time (%)';
        else
            condition_data.(condition_name).data = zone_data.times;
            y_label = 'Time (s)';
        end
        
        condition_data.(condition_name).mouse_ids = zone_data.mouse_ids;
        
        % Combine non-food arms if requested
        if options.combine_nonfood_arms
            condition_data.(condition_name).data = combineNonfoodArms(condition_data.(condition_name).data, options.sessions);
        end
    end
    
    % Ensure we have data for the requested sessions
    valid_sessions = options.sessions;
    
    % Create a separate figure for each arm
    for arm_idx = 1:length(arm_names)
        createArmFigure(arm_idx, arm_names{arm_idx}, condition_data, options, ...
                       valid_sessions, session_labels_map, y_label);
    end
end

function createArmFigure(arm_idx, arm_name, condition_data, options, valid_sessions, session_labels_map, y_label)
    % Create a figure for a specific arm comparing conditions across sessions
    % NOW: x-axis = conditions, colored bars = sessions
    
    num_conditions = length(options.conditions);
    num_sessions = length(valid_sessions);
    
    % Prepare data matrices - rows = conditions, columns = sessions
    groupedData = zeros(num_conditions, num_sessions);
    errorData = zeros(num_conditions, num_sessions);
    sessionData = cell(num_conditions, num_sessions);
    mouseIds = cell(num_conditions, num_sessions);
    
    % Collect data for this arm across all conditions and sessions
    for c = 1:num_conditions
        condition_name = options.conditions{c};
        
        for s = 1:num_sessions
            session = valid_sessions{s};
            
            if isfield(condition_data.(condition_name).data, session)
                data_matrix = condition_data.(condition_name).data.(session);
                
                if ~isempty(data_matrix) && size(data_matrix, 2) >= arm_idx
                    arm_data = data_matrix(:, arm_idx);
                    sessionData{c, s} = arm_data;
                    mouseIds{c, s} = condition_data.(condition_name).mouse_ids.(session);
                    
                    groupedData(c, s) = mean(arm_data);
                    errorData(c, s) = std(arm_data) / sqrt(length(arm_data));
                else
                    sessionData{c, s} = [];
                    mouseIds{c, s} = {};
                end
            else
                sessionData{c, s} = [];
                mouseIds{c, s} = {};
            end
        end
    end
    
    % Perform statistical tests
    stats_results = struct();
    if options.show_stats
        stats_results = performConditionStatisticalTests(sessionData, valid_sessions, options);
    end
    
    % Create figure
    figure('Position', options.figure_size);
    
    % Create grouped bar plot
    h = bar(groupedData, options.bar_width);
    
    % Set colors for each session
    for s = 1:num_sessions
        h(s).FaceColor = options.session_colors{mod(s-1, length(options.session_colors))+1};
        h(s).EdgeColor = 'none';
        h(s).FaceAlpha = 0.8;
    end
    
    hold on;
    
    % Add error bars and individual data points
    numBars = num_sessions;
    groupwidth = min(options.bar_width, numBars/(numBars+1.5));
    barPositions = zeros(num_conditions, numBars);
    
    for c = 1:num_conditions
        for s = 1:numBars
            % Calculate bar position
            if numBars == 1
                x = c;
            else
                x = c + (s - (numBars+1)/2) * groupwidth/numBars;
            end
            barPositions(c, s) = x;
            
            % Add error bars
            errorbar(x, groupedData(c, s), errorData(c, s), 'k', ...
                    'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
            
            % Add individual data points
            if ~isempty(sessionData{c, s})
                indPoints = sessionData{c, s};
                scatter(repmat(x, length(indPoints), 1), indPoints, 25, 'ko', ...
                       'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
            end
        end
    end
    
    % Connect individual mouse data points if requested
    if options.connect_points && num_sessions > 1
        connectMouseDataPointsAcrossConditions(sessionData, mouseIds, barPositions, options);
    end
    
    % Add significance indicators
    if options.show_stats
        addConditionSignificanceIndicators(stats_results, barPositions, groupedData, ...
                                          errorData, options, valid_sessions);
    end
    
    % Format the plot
    formatConditionPlot(arm_name, valid_sessions, session_labels_map, options, y_label);
    
    % Display statistics summary
    if options.show_stats
        displayConditionStatsSummary(stats_results, arm_name);
    end
    
    hold off;
end

function stats_results = performConditionStatisticalTests(sessionData, valid_sessions, options)
    % Perform statistical tests comparing conditions and sessions
    % sessionData is now: rows = conditions, columns = sessions
    stats_results = struct();
    stats_results.session_comparisons = struct();   % Compare sessions within each condition
    stats_results.condition_comparisons = struct(); % Compare conditions within each session
    
    num_conditions = size(sessionData, 1);
    num_sessions = size(sessionData, 2);
    
    % 1. Compare sessions within each condition (comparing across bars in one group)
    for c = 1:num_conditions
        condition_name = options.conditions{c};
        condition_field = matlab.lang.makeValidName(condition_name);
        
        all_data = [];
        session_groups = [];
        
        for s = 1:num_sessions
            if ~isempty(sessionData{c, s})
                all_data = [all_data; sessionData{c, s}];
                session_groups = [session_groups; s * ones(size(sessionData{c, s}))];
            end
        end
        
        if length(unique(session_groups)) >= 2
            if length(unique(session_groups)) == 2
                % t-test for two sessions
                sess1_data = all_data(session_groups == 1);
                sess2_data = all_data(session_groups == 2);
                [~, p] = ttest2(sess1_data, sess2_data);
                
                stats_results.session_comparisons.(condition_field).p_value = p;
                stats_results.session_comparisons.(condition_field).significant = p < options.alpha;
                stats_results.session_comparisons.(condition_field).test_type = 't-test';
                stats_results.session_comparisons.(condition_field).condition_name = condition_name;
            else
                % ANOVA for multiple sessions
                [p, ~, stats] = anova1(all_data, session_groups, 'off');
                
                stats_results.session_comparisons.(condition_field).p_value = p;
                stats_results.session_comparisons.(condition_field).significant = p < options.alpha;
                stats_results.session_comparisons.(condition_field).test_type = 'ANOVA';
                stats_results.session_comparisons.(condition_field).condition_name = condition_name;
                
                if p < options.alpha
                    [c_comp, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                    stats_results.session_comparisons.(condition_field).posthoc = c_comp;
                end
            end
        end
    end
    
    % 2. Compare conditions within each session (comparing across groups for one bar color)
    for s = 1:num_sessions
        session = valid_sessions{s};
        session_field = matlab.lang.makeValidName(session);
        
        all_data = [];
        condition_groups = [];
        
        for c = 1:num_conditions
            if ~isempty(sessionData{c, s})
                all_data = [all_data; sessionData{c, s}];
                condition_groups = [condition_groups; c * ones(size(sessionData{c, s}))];
            end
        end
        
        if length(unique(condition_groups)) >= 2
            if length(unique(condition_groups)) == 2
                % t-test for two conditions
                cond1_data = all_data(condition_groups == 1);
                cond2_data = all_data(condition_groups == 2);
                [~, p] = ttest2(cond1_data, cond2_data);
                
                stats_results.condition_comparisons.(session_field).p_value = p;
                stats_results.condition_comparisons.(session_field).significant = p < options.alpha;
                stats_results.condition_comparisons.(session_field).test_type = 't-test';
            else
                % ANOVA for multiple conditions
                [p, ~, stats] = anova1(all_data, condition_groups, 'off');
                
                stats_results.condition_comparisons.(session_field).p_value = p;
                stats_results.condition_comparisons.(session_field).significant = p < options.alpha;
                stats_results.condition_comparisons.(session_field).test_type = 'ANOVA';
                
                if p < options.alpha
                    [c_comp, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                    stats_results.condition_comparisons.(session_field).posthoc = c_comp;
                end
            end
        end
    end
end

function connectMouseDataPointsAcrossConditions(sessionData, mouseIds, barPositions, options)
    % Connect individual mouse data points across sessions within each condition
    % sessionData is now: rows = conditions, columns = sessions
    num_conditions = size(sessionData, 1);
    num_sessions = size(sessionData, 2);
    
    for c = 1:num_conditions
        % Collect all unique mice for this condition
        all_mice = {};
        for s = 1:num_sessions
            if ~isempty(mouseIds{c, s})
                all_mice = [all_mice, mouseIds{c, s}];
            end
        end
        all_mice = unique(all_mice);
        
        % Connect points for each mouse across sessions
        for m = 1:length(all_mice)
            current_mouse = all_mice{m};
            mouse_x = [];
            mouse_y = [];
            
            for s = 1:num_sessions
                if ~isempty(mouseIds{c, s})
                    mouse_idx = find(strcmp(mouseIds{c, s}, current_mouse));
                    
                    if ~isempty(mouse_idx)
                        x_pos = barPositions(c, s);
                        y_val = sessionData{c, s}(mouse_idx);
                        
                        mouse_x = [mouse_x, x_pos];
                        mouse_y = [mouse_y, y_val];
                    end
                end
            end
            
            % Connect the points if we have at least 2
            if length(mouse_x) >= 2
                plot(mouse_x, mouse_y, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5 0.3]);
            end
        end
    end
end

function addConditionSignificanceIndicators(stats_results, barPositions, groupedData, errorData, options, valid_sessions)
    % Add significance indicators for both session and condition comparisons
    star_height_offset = max(groupedData(:) + errorData(:)) * 0.03;
    
    % 1. Add stars for session comparisons within conditions (horizontal comparisons within groups)
    if isfield(stats_results, 'session_comparisons')
        condition_fields = fieldnames(stats_results.session_comparisons);
        
        for c = 1:length(condition_fields)
            condition_field = condition_fields{c};
            
            % Find the condition index
            condition_name = stats_results.session_comparisons.(condition_field).condition_name;
            condition_idx = find(strcmp(options.conditions, condition_name));
            
            if ~isempty(condition_idx) && stats_results.session_comparisons.(condition_field).significant
                % Get star symbol
                p_val = stats_results.session_comparisons.(condition_field).p_value;
                star_symbol = getStarSymbol(p_val);
                
                % Position for overall significance
                max_height = max(groupedData(condition_idx, :) + errorData(condition_idx, :));
                condition_x = mean(barPositions(condition_idx, :));
                star_y = max_height + star_height_offset * 3;
                
                text(condition_x, star_y, star_symbol, 'HorizontalAlignment', 'center', ...
                    'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
                
                % Add pairwise comparisons if available
                if isfield(stats_results.session_comparisons.(condition_field), 'posthoc') && ...
                   ~isempty(stats_results.session_comparisons.(condition_field).posthoc)
                    addSessionPairwiseLines(stats_results.session_comparisons.(condition_field).posthoc, ...
                                          condition_idx, barPositions, groupedData, errorData, ...
                                          star_height_offset, options);
                end
            end
        end
    end
    
    % 2. Add stars for condition comparisons within sessions (vertical comparisons across groups)
    if isfield(stats_results, 'condition_comparisons')
        session_fields = fieldnames(stats_results.condition_comparisons);
        
        for s = 1:length(session_fields)
            session_field = session_fields{s};
            
            % Find the session index
            session_idx = find(strcmp(cellfun(@(x) matlab.lang.makeValidName(x), valid_sessions, 'UniformOutput', false), session_field));
            
            if ~isempty(session_idx) && stats_results.condition_comparisons.(session_field).significant
                p_val = stats_results.condition_comparisons.(session_field).p_value;
                star_symbol = getStarSymbol(p_val);
                
                % Add pairwise comparison lines if available
                if isfield(stats_results.condition_comparisons.(session_field), 'posthoc') && ...
                   ~isempty(stats_results.condition_comparisons.(session_field).posthoc)
                    addConditionPairwiseLinesForSession(stats_results.condition_comparisons.(session_field).posthoc, ...
                                                      session_idx, barPositions, groupedData, errorData, ...
                                                      star_height_offset, options);
                else
                    % Show overall significance
                    max_height_sess = max(groupedData(:, session_idx) + errorData(:, session_idx));
                    sess_x = mean(barPositions(:, session_idx));
                    star_y = max_height_sess + star_height_offset * 6;
                    
                    text(sess_x, star_y, star_symbol, 'HorizontalAlignment', 'center', ...
                        'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', [0 0 0]);
                    
                    % Add horizontal line
                    if size(barPositions, 1) > 1
                        x_top = min(barPositions(:, session_idx));
                        x_bottom = max(barPositions(:, session_idx));
                        line_y = star_y - star_height_offset * 1.5;
                        plot([x_top, x_bottom], [line_y, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
                    end
                end
            end
        end
    end
end

function addSessionPairwiseLines(sig_comparisons, condition_idx, barPositions, groupedData, errorData, star_height_offset, options)
    % Add pairwise comparison lines between sessions within a condition
    base_height = max(groupedData(condition_idx, :) + errorData(condition_idx, :));
    line_spacing = star_height_offset * 1;
    
    for i = 1:size(sig_comparisons, 1)
        comp = sig_comparisons(i, :);
        sess1_idx = comp(1);
        sess2_idx = comp(2);
        p_val = comp(6);
        
        star_symbol = getStarSymbol(p_val);
        
        x1 = barPositions(condition_idx, sess1_idx);
        x2 = barPositions(condition_idx, sess2_idx);
        line_y = base_height + line_spacing * (i + 0.5);
        
        % Draw significance line
        plot([x1, x2], [line_y, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x1, x1], [line_y - line_spacing*0.1, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x2, x2], [line_y - line_spacing*0.1, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        
        text(mean([x1, x2]), line_y + line_spacing*0.15, star_symbol, ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size, 'FontWeight', 'bold', 'Color', [0 0 0]);
    end
end

function addConditionPairwiseLinesForSession(sig_comparisons, session_idx, barPositions, groupedData, errorData, star_height_offset, options)
    % Add pairwise comparison lines between conditions within a session
    base_height = max(groupedData(:, session_idx) + errorData(:, session_idx));
    line_spacing = star_height_offset * 1;
    
    for i = 1:size(sig_comparisons, 1)
        comp = sig_comparisons(i, :);
        cond1_idx = comp(1);
        cond2_idx = comp(2);
        p_val = comp(6);
        
        star_symbol = getStarSymbol(p_val);
        
        x1 = barPositions(cond1_idx, session_idx);
        x2 = barPositions(cond2_idx, session_idx);
        line_y = base_height + line_spacing * (i + 0.5);
        
        % Draw significance line
        plot([x1, x2], [line_y, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x1, x1], [line_y - line_spacing*0.1, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x2, x2], [line_y - line_spacing*0.1, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        
        text(mean([x1, x2]), line_y + line_spacing*0.15, star_symbol, ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size, 'FontWeight', 'bold', 'Color', [0 0 0]);
    end
end

function formatConditionPlot(arm_name, valid_sessions, session_labels_map, options, y_label)
    % Format the plot
    ylabel(y_label, 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    % Create title
    plot_title = ['Condition Comparison - ', arm_name];
    if options.combine_nonfood_arms && contains(arm_name, 'Non-Food')
        plot_title = [plot_title, ' (Combined)'];
    end
    title(plot_title, 'FontSize', options.font_size + 4, 'FontWeight', 'bold');
    
    % Create session labels for legend
    sessionLabels = cell(1, length(valid_sessions));
    for i = 1:length(valid_sessions)
        if isKey(session_labels_map, valid_sessions{i})
            sessionLabels{i} = session_labels_map(valid_sessions{i});
        else
            sessionLabels{i} = strrep(valid_sessions{i}, 'sess', 'Session ');
        end
    end
    
    % Format condition names for x-axis
    conditionLabels = cell(1, length(options.conditions));
    for i = 1:length(options.conditions)
        conditionLabels{i} = formatGroupName(options.conditions{i});
    end
    
    % Set x-axis to conditions
    xticks(1:length(options.conditions));
    xticklabels(conditionLabels);
    set(gca, 'FontSize', options.font_size);
    
    % Add legend with session labels
    legend(sessionLabels, 'Location', options.legend_location, 'FontSize', options.font_size);
    legend('boxoff');
    
    % Set y-axis limits
    if ~isempty(options.ylim)
        ylim(options.ylim);
    elseif strcmp(options.plot_type, 'percent')
        ylim([0, 100]);
    else
        current_ylim = ylim;
        ylim([current_ylim(1), current_ylim(2) * 1.3]);
    end
    
    box off;
    grid off;
    set(gca, 'LineWidth', 1.2);
    set(gcf, 'Color', 'white');
    set(gca, 'Color', 'white');
end

function displayConditionStatsSummary(stats_results, arm_name)
    % Display statistics summary
    fprintf('\n===== Statistical Summary for %s =====\n', arm_name);
    
    % Session comparisons within conditions
    if isfield(stats_results, 'session_comparisons')
        fprintf('\nSession comparisons within conditions:\n');
        condition_fields = fieldnames(stats_results.session_comparisons);
        for i = 1:length(condition_fields)
            condition_field = condition_fields{i};
            if stats_results.session_comparisons.(condition_field).significant
                condition_name = stats_results.session_comparisons.(condition_field).condition_name;
                fprintf('  %s: p = %.4f (%s, significant)\n', ...
                    formatGroupName(condition_name), ...
                    stats_results.session_comparisons.(condition_field).p_value, ...
                    stats_results.session_comparisons.(condition_field).test_type);
            end
        end
    end
    
    % Condition comparisons within sessions
    if isfield(stats_results, 'condition_comparisons')
        fprintf('\nCondition comparisons within sessions:\n');
        session_fields = fieldnames(stats_results.condition_comparisons);
        for i = 1:length(session_fields)
            session_field = session_fields{i};
            if stats_results.condition_comparisons.(session_field).significant
                fprintf('  %s: p = %.4f (%s, significant)\n', ...
                    session_field, ...
                    stats_results.condition_comparisons.(session_field).p_value, ...
                    stats_results.condition_comparisons.(session_field).test_type);
            end
        end
    end
    
    fprintf('=====================================\n');
end

function combined_data = combineNonfoodArms(data_to_plot, valid_sessions)
    % Combine non-food arms (columns 2 and 3) by summing them
    combined_data = struct();
    
    for i = 1:length(valid_sessions)
        session = valid_sessions{i};
        original_data = data_to_plot.(session);
        
        if ~isempty(original_data) && size(original_data, 2) >= 3
            combined_matrix = zeros(size(original_data, 1), size(original_data, 2) - 1);
            combined_matrix(:, 1) = original_data(:, 1); % Food arm
            combined_matrix(:, 2) = original_data(:, 2) + original_data(:, 3); % Combined non-food arms
            
            % Add any additional columns (like horizontal)
            if size(original_data, 2) > 3
                combined_matrix(:, 3:end) = original_data(:, 4:end);
            end
            
            combined_data.(session) = combined_matrix;
        else
            combined_data.(session) = original_data;
        end
    end
end

function star_symbol = getStarSymbol(p_val)
    % Get star symbol based on p-value
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

function formatted_name = formatGroupName(group_name)
    % Format group names for display
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