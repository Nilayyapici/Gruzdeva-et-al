function plotCorrelationConditionComparison(mice_all, options)
    % Compare dF/F-distance correlations between different conditions across sessions
    % Creates separate figures for each arm (Food Arm, Non-Food Arm 1, Non-Food Arm 2)
    % NOW: Each figure shows CONDITIONS on x-axis with bars for different SESSIONS
    
    % Set defaults for options
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'conditions'), error('Must specify conditions to compare in options.conditions'); end
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 200; end
    if ~isfield(options, 'dist_too_close'), options.dist_too_close = 0; end
    if ~isfield(options, 'filter_order_dist'), options.filter_order_dist = 1; end
    if ~isfield(options, 'cutoff_freq_dist'), options.cutoff_freq_dist = 0.08; end
    if ~isfield(options, 'filter_order_dff'), options.filter_order_dff = 1; end
    if ~isfield(options, 'cutoff_freq_dff'), options.cutoff_freq_dff = 0.05; end
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 0; end
    if ~isfield(options, 'connect_points'), options.connect_points = true; end
    if ~isfield(options, 'alpha'), options.alpha = 0.05; end
    if ~isfield(options, 'show_stats'), options.show_stats = true; end
    if ~isfield(options, 'ylimit'), options.ylimit = [-0.7, 0.7]; end
    if ~isfield(options, 'session_colors')
        % Now using session colors instead of condition colors
        options.session_colors = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.8 0.2], [0.8 0.6 0.2]};
    end
    if ~isfield(options, 'figure_size'), options.figure_size = [100, 100, 800, 600]; end
    if ~isfield(options, 'font_size'), options.font_size = 12; end
    if ~isfield(options, 'bar_width'), options.bar_width = 0.8; end
    if ~isfield(options, 'group_spacing'), options.group_spacing = 1.5; end
    if ~isfield(options, 'stats_color'), options.stats_color = [0 0 0]; end
    if ~isfield(options, 'show_zero_line'), options.show_zero_line = true; end
    if ~isfield(options, 'use_fisher_z'), options.use_fisher_z = true; end
    if ~isfield(options, 'legend_location'), options.legend_location = 'northeast'; end
    
    % Session labels mapping
    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});
    
    % Arms are always the same for correlations
    arm_names = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
    
    % Get correlation data for each condition
    condition_data = struct();
    for c = 1:length(options.conditions)
        condition_name = options.conditions{c};
        temp_options = options;
        temp_options.group = condition_name;
        temp_options.suppress_figures = true;
        temp_options.run_stats = false;
        
        corr_data = analyzeDffDistanceCorrelations(mice_all, temp_options);
        
        condition_data.(condition_name).corr = corr_data.corr;  % Raw correlations for display
        condition_data.(condition_name).mouse_ids = corr_data.mouse_ids;
        
        % Store Fisher z-transformed data if available
        if isfield(corr_data, 'z_corr')
            condition_data.(condition_name).z_corr = corr_data.z_corr;
        end
    end
    
    % Create a separate figure for each arm
    for arm_idx = 1:length(arm_names)
        createCorrelationArmFigure(arm_idx, arm_names{arm_idx}, condition_data, options, ...
                                   session_labels_map);
    end
end

function createCorrelationArmFigure(arm_idx, arm_name, condition_data, options, session_labels_map)
    % Create a figure for a specific arm comparing conditions across sessions
    % NOW: x-axis = conditions, colored bars = sessions
    
    num_conditions = length(options.conditions);
    num_sessions = length(options.sessions);
    
    % Prepare data matrices - rows = conditions, columns = sessions
    groupedData = zeros(num_conditions, num_sessions);
    errorData = zeros(num_conditions, num_sessions);
    sessionData = cell(num_conditions, num_sessions);  % Raw correlations for display
    sessionDataStats = cell(num_conditions, num_sessions);  % For statistics
    mouseIds = cell(num_conditions, num_sessions);
    
    % Collect data for this arm across all conditions and sessions
    for c = 1:num_conditions
        condition_name = options.conditions{c};
        
        for s = 1:num_sessions
            session = options.sessions{s};
            
            if isfield(condition_data.(condition_name).corr, session)
                corr_matrix = condition_data.(condition_name).corr.(session);
                
                if ~isempty(corr_matrix) && size(corr_matrix, 2) >= arm_idx
                    arm_corr = corr_matrix(:, arm_idx);
                    sessionData{c, s} = arm_corr;  % Raw correlations for display
                    mouseIds{c, s} = condition_data.(condition_name).mouse_ids.(session);
                    
                    % Use Fisher z for statistics if available and requested
                    if options.use_fisher_z && isfield(condition_data.(condition_name), 'z_corr') && ...
                       isfield(condition_data.(condition_name).z_corr, session)
                        z_matrix = condition_data.(condition_name).z_corr.(session);
                        if ~isempty(z_matrix) && size(z_matrix, 2) >= arm_idx
                            sessionDataStats{c, s} = z_matrix(:, arm_idx);
                        else
                            sessionDataStats{c, s} = arm_corr;
                        end
                    else
                        sessionDataStats{c, s} = arm_corr;
                    end
                    
                    % Calculate mean and SEM (always display raw correlations)
                    valid_corr = arm_corr(~isnan(arm_corr));
                    if ~isempty(valid_corr)
                        groupedData(c, s) = nanmean(valid_corr);
                        errorData(c, s) = nanstd(valid_corr) / sqrt(length(valid_corr));
                    else
                        groupedData(c, s) = 0;
                        errorData(c, s) = 0;
                    end
                else
                    sessionData{c, s} = [];
                    sessionDataStats{c, s} = [];
                    mouseIds{c, s} = {};
                end
            else
                sessionData{c, s} = [];
                sessionDataStats{c, s} = [];
                mouseIds{c, s} = {};
            end
        end
    end
    
    % Perform statistical tests
    stats_results = struct();
    if options.show_stats
        stats_results = performCorrelationConditionStatisticalTests(sessionData, sessionDataStats, ...
                                                                     options.sessions, options);
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
    
    % % Add zero reference line
    % if options.show_zero_line
    %     h_zero = plot([0.5, num_conditions + 0.5], [0, 0], 'k--', 'LineWidth', 1.5);
    %     try
    %         h_zero.Color(4) = 0.7;
    %     catch
    %         h_zero.Color = [0.5 0.5 0.5];
    %     end
    % end
    
    % Add error bars and individual data points
    numBars = num_sessions;
    groupwidth = min(options.bar_width, numBars/(numBars+options.group_spacing));
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
            
            % Add individual data points (raw correlations)
            if ~isempty(sessionData{c, s})
                indPoints = sessionData{c, s};
                validPoints = indPoints(~isnan(indPoints));
                if ~isempty(validPoints)
                    scatter(repmat(x, length(validPoints), 1), validPoints, 25, 'ko', ...
                           'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
                end
            end
        end
    end
    
    % Connect individual mouse data points if requested
    if options.connect_points && num_sessions > 1
        connectCorrelationMouseDataPointsAcrossSessions(sessionData, mouseIds, barPositions);
    end
    
    % Add significance indicators
    if options.show_stats
        addCorrelationConditionSignificanceIndicators(stats_results, barPositions, groupedData, ...
                                                      errorData, options, options.sessions);
    end
    
    % Format the plot
    formatCorrelationConditionPlot(arm_name, options.sessions, session_labels_map, options);
    
    % Display statistics summary
    if options.show_stats
        displayCorrelationConditionStatsSummary(stats_results, arm_name);
    end
    
    hold off;
end

function stats_results = performCorrelationConditionStatisticalTests(sessionData, sessionDataStats, valid_sessions, options)
    % Perform statistical tests comparing conditions and sessions for correlations
    % sessionData is now: rows = conditions, columns = sessions
    stats_results = struct();
    stats_results.session_comparisons = struct();   % Compare sessions within each condition
    stats_results.condition_comparisons = struct(); % Compare conditions within each session
    stats_results.fisher_z_used = options.use_fisher_z;
    
    num_conditions = size(sessionData, 1);
    num_sessions = size(sessionData, 2);
    
    % 1. Compare sessions within each condition (comparing across bars in one group)
    for c = 1:num_conditions
        condition_name = options.conditions{c};
        condition_field = matlab.lang.makeValidName(condition_name);
        
        all_data = [];
        session_groups = [];
        
        for s = 1:num_sessions
            if ~isempty(sessionDataStats{c, s})
                valid_values = sessionDataStats{c, s}(~isnan(sessionDataStats{c, s}));
                if ~isempty(valid_values)
                    all_data = [all_data; valid_values];
                    session_groups = [session_groups; s * ones(size(valid_values))];
                end
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
            if ~isempty(sessionDataStats{c, s})
                valid_values = sessionDataStats{c, s}(~isnan(sessionDataStats{c, s}));
                if ~isempty(valid_values)
                    all_data = [all_data; valid_values];
                    condition_groups = [condition_groups; c * ones(size(valid_values))];
                end
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

function connectCorrelationMouseDataPointsAcrossSessions(sessionData, mouseIds, barPositions)
    % Connect individual mouse correlation data points across sessions within each condition
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
                    
                    if ~isempty(mouse_idx) && ~isempty(sessionData{c, s})
                        correlation_val = sessionData{c, s}(mouse_idx);
                        
                        if ~isnan(correlation_val)
                            x_pos = barPositions(c, s);
                            mouse_x = [mouse_x, x_pos];
                            mouse_y = [mouse_y, correlation_val];
                        end
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

function addCorrelationConditionSignificanceIndicators(stats_results, barPositions, groupedData, errorData, options, valid_sessions)
    % Add significance indicators for correlation condition comparisons
    star_height_offset = (max(options.ylimit) - min(options.ylimit)) * 0.025;
    
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
                star_y = max_height + star_height_offset * 6;
                
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
                    star_y = max_height_sess + star_height_offset * 8;
                    
                    text(sess_x, star_y, star_symbol, 'HorizontalAlignment', 'center', ...
                        'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', options.stats_color);
                    
                    % Add horizontal line
                    if size(barPositions, 1) > 1
                        x_top = min(barPositions(:, session_idx));
                        x_bottom = max(barPositions(:, session_idx));
                        line_y = star_y - star_height_offset * 1.5;
                        plot([x_top, x_bottom], [line_y, line_y], 'Color', options.stats_color, 'LineWidth', 2);
                    end
                end
            end
        end
    end
end

function addSessionPairwiseLines(sig_comparisons, condition_idx, barPositions, groupedData, errorData, star_height_offset, options)
    % Add pairwise comparison lines between sessions within a condition
    base_height = max(groupedData(condition_idx, :) + errorData(condition_idx, :));
    line_spacing = star_height_offset * 2;
    
    for i = 1:size(sig_comparisons, 1)
        comp = sig_comparisons(i, :);
        sess1_idx = comp(1);
        sess2_idx = comp(2);
        p_val = comp(6);
        
        star_symbol = getStarSymbol(p_val);
        
        x1 = barPositions(condition_idx, sess1_idx);
        x2 = barPositions(condition_idx, sess2_idx);
        line_y = base_height + line_spacing * (i + 1.5);
        
        % Draw significance line
        plot([x1, x2], [line_y, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x1, x1], [line_y - line_spacing*0.15, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        plot([x2, x2], [line_y - line_spacing*0.15, line_y], 'Color', [0 0 0], 'LineWidth', 1.5);
        
        text(mean([x1, x2]), line_y + line_spacing*0.1, star_symbol, ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size, 'FontWeight', 'bold', 'Color', [0 0 0]);
    end
end

function addConditionPairwiseLinesForSession(sig_comparisons, session_idx, barPositions, groupedData, errorData, star_height_offset, options)
    % Add pairwise comparison lines between conditions within a session
    base_height = max(groupedData(:, session_idx) + errorData(:, session_idx));
    line_spacing = star_height_offset * 2;
    
    for i = 1:size(sig_comparisons, 1)
        comp = sig_comparisons(i, :);
        cond1_idx = comp(1);
        cond2_idx = comp(2);
        p_val = comp(6);
        
        star_symbol = getStarSymbol(p_val);
        
        x1 = barPositions(cond1_idx, session_idx);
        x2 = barPositions(cond2_idx, session_idx);
        line_y = base_height + line_spacing * (i + 1.5);
        
        % Draw significance line
        plot([x1, x2], [line_y, line_y], 'Color', options.stats_color, 'LineWidth', 2);
        plot([x1, x1], [line_y - line_spacing*0.15, line_y], 'Color', options.stats_color, 'LineWidth', 1.5);
        plot([x2, x2], [line_y - line_spacing*0.15, line_y], 'Color', options.stats_color, 'LineWidth', 1.5);
        
        text(mean([x1, x2]), line_y + line_spacing*0.1, star_symbol, ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size, 'FontWeight', 'bold', 'Color', options.stats_color);
    end
end

function formatCorrelationConditionPlot(arm_name, valid_sessions, session_labels_map, options)
    % Format the plot
    ylabel('dF/F-Distance Correlation (r)', 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    % Create title
    title(['Correlation Condition Comparison - ', arm_name], 'FontSize', options.font_size + 4, 'FontWeight', 'bold');
    
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
    ylim(options.ylimit);
    
    box off;
    grid off;
    set(gca, 'LineWidth', 1.2);
    set(gcf, 'Color', 'white');
    set(gca, 'Color', 'white');
end

function displayCorrelationConditionStatsSummary(stats_results, arm_name)
    % Display statistics summary
    fprintf('\n===== Correlation Condition Comparison for %s =====\n', arm_name);
    
    if isfield(stats_results, 'fisher_z_used')
        if stats_results.fisher_z_used
            fprintf('Statistical tests performed using Fisher z-transformation\n');
        else
            fprintf('Statistical tests performed using raw correlations\n');
        end
    end
    
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
    
    fprintf('====================================================\n');
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