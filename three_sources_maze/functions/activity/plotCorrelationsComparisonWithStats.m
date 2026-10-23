function plotCorrelationsComparisonWithStats(mice_all, options)
    % Enhanced plotting function for dF/F-distance correlations with statistical annotations
    % Creates publication-ready plots with significance indicators for correlation comparisons
    
    % Set defaults for options
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'group'), options.group = 'all'; end
    if ~isfield(options, 'title'), options.title = 'dF/F-Distance Correlations Across Sessions'; end
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 200; end
    if ~isfield(options, 'dist_too_close'), options.dist_too_close = 5; end
    if ~isfield(options, 'filter_order_dist'), options.filter_order_dist = 1; end
    if ~isfield(options, 'cutoff_freq_dist'), options.cutoff_freq_dist = 0.08; end
    if ~isfield(options, 'filter_order_dff'), options.filter_order_dff = 1; end
    if ~isfield(options, 'cutoff_freq_dff'), options.cutoff_freq_dff = 0.05; end
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 0; end
    if ~isfield(options, 'connect_points'), options.connect_points = true; end
    if ~isfield(options, 'alpha'), options.alpha = 0.05; end
    if ~isfield(options, 'show_stats'), options.show_stats = true; end
    if ~isfield(options, 'ylimit'), options.ylimit = [-0.8, 0.7]; end
    if ~isfield(options, 'colors'), options.colors = {[0.6 0.6 0.6], [0.2 0.4 0.8], [0.1 0.7 0.9]}; end
    if ~isfield(options, 'figure_size'), options.figure_size = [100, 100, 1000, 600]; end
    if ~isfield(options, 'font_size'), options.font_size = 12; end
    if ~isfield(options, 'bar_width'), options.bar_width = 0.8; end
    if ~isfield(options, 'stats_color'), options.stats_color = [0 0.6 0]; end
    if ~isfield(options, 'show_zero_line'), options.show_zero_line = true; end
    if ~isfield(options, 'suppress_figures'), options.suppress_figures = false; end
    if ~isfield(options, 'use_fisher_z'), options.use_fisher_z = true; end  % Use Fisher z-transformation for statistics
    
    % Session labels mapping
    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});
    
    % Analyze correlation data (suppress internal figures)
    options_internal = options;
    options_internal.suppress_figures = true;
    options_internal.run_stats = false; % We'll handle stats ourselves
    corr_data = analyzeDffDistanceCorrelations(mice_all, options_internal);
    
    if options.suppress_figures
        return; % Exit early if figures are suppressed
    end
    
    % Ensure we have data for the requested sessions
    valid_sessions = intersect(fieldnames(corr_data.corr)', options.sessions);
    if isempty(valid_sessions)
        error('No valid sessions found for plotting correlations.');
    end
    
    % Arms are always the same for correlations
    arms = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
    numArms = length(arms);
    
    % Get data for each session and calculate stats
    sessionData = cell(1, length(valid_sessions));  % Raw correlations for display
    sessionDataStats = cell(1, length(valid_sessions));  % Data for statistics (raw or Fisher z)
    sessionMeans = cell(1, length(valid_sessions));
    sessionSEMs = cell(1, length(valid_sessions));
    mouseIds = cell(1, length(valid_sessions));
    
    for i = 1:length(valid_sessions)
        sessionData{i} = corr_data.corr.(valid_sessions{i});  % Always use raw for display
        
        % Choose data for statistics based on option
        if options.use_fisher_z && isfield(corr_data, 'z_corr') && isfield(corr_data.z_corr, valid_sessions{i})
            sessionDataStats{i} = corr_data.z_corr.(valid_sessions{i});  % Use Fisher z for statistics
        else
            sessionDataStats{i} = corr_data.corr.(valid_sessions{i});    % Use raw correlations
        end
        
        mouseIds{i} = corr_data.mouse_ids.(valid_sessions{i});
        
        if ~isempty(sessionData{i})
            % Always display raw correlations (not Fisher z)
            sessionMeans{i} = nanmean(sessionData{i}, 1);
            valid_counts = sum(~isnan(sessionData{i}), 1);
            sessionSEMs{i} = nanstd(sessionData{i}, [], 1) ./ sqrt(valid_counts);
            % Replace NaN SEMs with 0
            sessionSEMs{i}(isnan(sessionSEMs{i})) = 0;
        else
            sessionMeans{i} = zeros(1, numArms);
            sessionSEMs{i} = zeros(1, numArms);
        end
    end
    
    % Collect all unique mouse IDs
    all_mice = unique([mouseIds{:}]);
    
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
        stats_results = performCorrelationStatisticalTests(corr_data, sessionDataStats, valid_sessions, arms, options);
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
    
    % Add zero reference line if requested
    if options.show_zero_line
        h_zero = plot([0.5, numArms + 0.5], [0, 0], 'k--', 'LineWidth', 1.5);
        % Set alpha using compatible method
        try
            h_zero.Color(4) = 0.7;  % Set alpha channel
        catch
            % Fallback for older MATLAB versions - use lighter gray color
            h_zero.Color = [0.5 0.5 0.5];
        end
    end
    
    % Add error bars and calculate bar positions
    numBars = length(valid_sessions);
    groupwidth = min(options.bar_width, numBars/(numBars+1.5));
    barPositions = zeros(numArms, numBars);
    
    for i = 1:numArms
        for j = 1:numBars
            % Use the same calculation as MATLAB's grouped bar positioning
            if numBars == 1
                x = i;
            else
                x = i + (j - (numBars+1)/2) * groupwidth/numBars;
            end
            barPositions(i, j) = x;
            
            % Add error bars
            errorbar(x, groupedData(i, j), errorData(i, j), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4);
            
            % Add individual data points - perfectly centered on bars
            if ~isempty(sessionData{j}) && i <= size(sessionData{j}, 2)
                indPoints = sessionData{j}(:, i);
                % Remove NaN values for plotting
                validPoints = indPoints(~isnan(indPoints));
                if ~isempty(validPoints)
                    scatter(repmat(x, length(validPoints), 1), validPoints, 25, 'ko', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
                end
            end
        end
    end
    
    % Connect individual mouse data points across sessions if requested
    if options.connect_points && length(valid_sessions) > 1
        connectCorrelationMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, valid_sessions, numArms);
    end
    
    % Add significance indicators if statistics were performed
    if options.show_stats && (isfield(stats_results, 'arm_comparisons') || isfield(stats_results, 'session_comparisons'))
        addCorrelationSignificanceIndicators(stats_results, barPositions, groupedData, errorData, options, valid_sessions, arms);
    end
    
    % Format the plot
    formatCorrelationPlot(arms, valid_sessions, session_labels_map, options);
    
    % Display statistics summary if available
    if options.show_stats && isfield(stats_results, 'summary')
        displayCorrelationStatsSummary(stats_results);
    end
    
    hold off;
end

function stats_results = performCorrelationStatisticalTests(corr_data, sessionDataStats, valid_sessions, arms, options)
    % Perform comprehensive statistical tests for correlation data
    % Uses sessionDataStats which can be either raw correlations or Fisher z-transformed
    stats_results = struct();
    stats_results.arm_comparisons = struct();
    stats_results.session_comparisons = struct();
    stats_results.summary = struct();
    stats_results.zero_tests = struct(); % Test if correlations are significantly different from zero
    stats_results.fisher_z_used = options.use_fisher_z; % Record which method was used
    
    numArms = length(arms);
    numSessions = length(valid_sessions);
    
    % 1. Test if correlations are significantly different from zero
    for s = 1:numSessions
        session = valid_sessions{s};
        if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s})
            session_stats_data = sessionDataStats{s};  % Use stats data (raw or Fisher z)
            
            for a = 1:numArms
                arm_values = session_stats_data(:, a);
                valid_values = arm_values(~isnan(arm_values));
                
                if length(valid_values) > 3
                    % One-sample t-test against zero
                    [~, p] = ttest(valid_values, 0);
                    stats_results.zero_tests.(session).p_values(a) = p;
                    stats_results.zero_tests.(session).significant(a) = p < options.alpha;
                    
                    % Store the display mean (always raw correlation for interpretation)
                    if isfield(corr_data, 'corr') && isfield(corr_data.corr, session)
                        raw_values = corr_data.corr.(session)(:, a);
                        valid_raw = raw_values(~isnan(raw_values));
                        if ~isempty(valid_raw)
                            stats_results.zero_tests.(session).means(a) = mean(valid_raw);
                        else
                            stats_results.zero_tests.(session).means(a) = NaN;
                        end
                    else
                        stats_results.zero_tests.(session).means(a) = NaN;
                    end
                else
                    stats_results.zero_tests.(session).p_values(a) = NaN;
                    stats_results.zero_tests.(session).significant(a) = false;
                    stats_results.zero_tests.(session).means(a) = NaN;
                end
            end
        end
    end
    
    % 2. Compare arms within each session (one-way ANOVA)
    for s = 1:numSessions
        session = valid_sessions{s};
        if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s})
            session_data = sessionDataStats{s};  % Use stats data (raw or Fisher z)
            
            % Remove rows with any NaN values for ANOVA
            valid_rows = ~any(isnan(session_data), 2);
            if sum(valid_rows) > 3 && size(session_data, 2) >= numArms
                clean_data = session_data(valid_rows, 1:numArms);
                
                % Reshape data for ANOVA
                values = clean_data(:);
                groups = repmat((1:numArms)', size(clean_data, 1), 1);
                
                % Perform ANOVA
                [p, ~, stats] = anova1(values, groups, 'off');
                
                stats_results.arm_comparisons.(session).p_value = p;
                stats_results.arm_comparisons.(session).significant = p < options.alpha;
                
                % If significant, perform post-hoc tests
                if p < options.alpha
                    [c, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                    stats_results.arm_comparisons.(session).posthoc = c;
                else
                    stats_results.arm_comparisons.(session).posthoc = [];
                end
            end
        end
    end
    
    % 3. Compare sessions within each arm (one-way ANOVA or t-test)
    for a = 1:numArms
        arm_name = arms{a};
        % Create valid field name
        arm_field_name = matlab.lang.makeValidName(arm_name);
        arm_data = [];
        session_groups = [];
        
        % Collect data for this arm across sessions (using stats data)
        for s = 1:numSessions
            if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s}) && size(sessionDataStats{s}, 2) >= a
                arm_values = sessionDataStats{s}(:, a);  % Use stats data
                valid_values = arm_values(~isnan(arm_values));
                
                if ~isempty(valid_values)
                    arm_data = [arm_data; valid_values];
                    session_groups = [session_groups; s * ones(size(valid_values))];
                end
            end
        end
        
        if length(unique(session_groups)) >= 2
            if length(unique(session_groups)) == 2
                % Perform t-test for two sessions
                group1_data = arm_data(session_groups == 1);
                group2_data = arm_data(session_groups == 2);
                [~, p] = ttest2(group1_data, group2_data);
                
                stats_results.session_comparisons.(arm_field_name).p_value = p;
                stats_results.session_comparisons.(arm_field_name).significant = p < options.alpha;
                stats_results.session_comparisons.(arm_field_name).test_type = 't-test';
                stats_results.session_comparisons.(arm_field_name).arm_name = arm_name;
            else
                % Perform ANOVA for multiple sessions
                [p, ~, stats] = anova1(arm_data, session_groups, 'off');
                
                stats_results.session_comparisons.(arm_field_name).p_value = p;
                stats_results.session_comparisons.(arm_field_name).significant = p < options.alpha;
                stats_results.session_comparisons.(arm_field_name).test_type = 'ANOVA';
                stats_results.session_comparisons.(arm_field_name).arm_name = arm_name;
                
                % If significant, perform post-hoc tests
                if p < options.alpha
                    [c, ~, ~, ~] = multcompare(stats, 'Alpha', options.alpha, 'Display', 'off');
                    stats_results.session_comparisons.(arm_field_name).posthoc = c;
                end
            end
        end
    end
    
    % Create summary
    stats_results.summary.significant_arm_comparisons = 0;
    stats_results.summary.significant_session_comparisons = 0;
    stats_results.summary.significant_zero_tests = 0;
    
    % Count significant comparisons
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
    
    % Count significant zero tests
    zero_fields = fieldnames(stats_results.zero_tests);
    for i = 1:length(zero_fields)
        if isfield(stats_results.zero_tests.(zero_fields{i}), 'significant')
            stats_results.summary.significant_zero_tests = stats_results.summary.significant_zero_tests + sum(stats_results.zero_tests.(zero_fields{i}).significant);
        end
    end
end

function connectCorrelationMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, valid_sessions, numArms)
    % Connect individual mouse correlation data points across sessions
    % Always uses raw correlations for display
    for arm_idx = 1:numArms
        for m = 1:length(all_mice)
            current_mouse = all_mice{m};
            
            % Collect data points for this mouse across sessions
            mouse_x = [];
            mouse_y = [];
            
            for sess_idx = 1:length(valid_sessions)
                mouse_in_session = find(strcmp(mouseIds{sess_idx}, current_mouse));
                
                if ~isempty(mouse_in_session) && arm_idx <= size(sessionData{sess_idx}, 2)
                    correlation_val = sessionData{sess_idx}(mouse_in_session, arm_idx);  % Raw correlation for display
                    
                    % Only include non-NaN correlations
                    if ~isnan(correlation_val)
                        x_pos = barPositions(arm_idx, sess_idx);
                        mouse_x = [mouse_x, x_pos];
                        mouse_y = [mouse_y, correlation_val];
                    end
                end
            end
            
            % Connect the points with a line if we have at least 2 points
            if length(mouse_x) >= 2
                plot(mouse_x, mouse_y, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5 0.3]);
            end
        end
    end
end

function addCorrelationSignificanceIndicators(stats_results, barPositions, groupedData, errorData, options, valid_sessions, arms)
    % Add significance indicators for correlation comparisons
    star_height_offset = (max(options.ylimit) - min(options.ylimit)) * 0.025;  % Slightly larger offset
    
    % 1. Skip zero test indicators (bullet points removed)
    
    % 2. Add stars for arm comparisons within sessions
    if isfield(stats_results, 'arm_comparisons')
        session_fields = fieldnames(stats_results.arm_comparisons);
        
        for s = 1:length(session_fields)
            session = session_fields{s};
            session_idx = find(strcmp(valid_sessions, session));
            
            if ~isempty(session_idx) && stats_results.arm_comparisons.(session).significant
                max_height = max(groupedData(:, session_idx) + errorData(:, session_idx));
                session_x = mean(barPositions(:, session_idx));
                star_y = max_height + star_height_offset * 6;  % Raised a bit (was * 4)
                
                p_val = stats_results.arm_comparisons.(session).p_value;
                star_symbol = getCorrelationStarSymbol(p_val);
                
                text(session_x, star_y, star_symbol, 'HorizontalAlignment', 'center', ...
                    'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
            end
        end
    end
    
    % 3. Add stars for session comparisons within arms
    if isfield(stats_results, 'session_comparisons')
        arm_fields = fieldnames(stats_results.session_comparisons);
        
        for a = 1:length(arm_fields)
            arm_field = arm_fields{a};
            if stats_results.session_comparisons.(arm_field).significant
                if isfield(stats_results.session_comparisons.(arm_field), 'arm_name')
                    arm_name = stats_results.session_comparisons.(arm_field).arm_name;
                else
                    arm_name = arm_field;
                end
                
                arm_idx = find(strcmp(arms, arm_name));
                if ~isempty(arm_idx) && arm_idx <= size(barPositions, 1)
                    p_val = stats_results.session_comparisons.(arm_field).p_value;
                    star_symbol = getCorrelationStarSymbol(p_val);
                    
                    % If we have post-hoc results, show specific pairwise comparisons
                    if isfield(stats_results.session_comparisons.(arm_field), 'posthoc') && ...
                       ~isempty(stats_results.session_comparisons.(arm_field).posthoc)
                        
                        posthoc = stats_results.session_comparisons.(arm_field).posthoc;
                        sig_comparisons = posthoc(posthoc(:, 6) < options.alpha, :);
                        
                        if ~isempty(sig_comparisons)
                            addCorrelationSessionPairwiseLines(sig_comparisons, arm_idx, barPositions, groupedData, errorData, star_height_offset, valid_sessions, options);
                        end
                    else
                        % No post-hoc data, show overall significance
                        max_height_arm = max(groupedData(arm_idx, :) + errorData(arm_idx, :));
                        arm_x = mean(barPositions(arm_idx, :));
                        star_y = max_height_arm + star_height_offset * 8;  % Raised a bit (was * 6)
                        
                        text(arm_x, star_y, star_symbol, 'HorizontalAlignment', 'center', ...
                            'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', options.stats_color);
                        
                        if length(valid_sessions) > 1
                            x_left = min(barPositions(arm_idx, :));
                            x_right = max(barPositions(arm_idx, :));
                            line_y = star_y - star_height_offset * 1.5;  % Adjusted line position
                            plot([x_left, x_right], [line_y, line_y], 'Color', options.stats_color, 'LineWidth', 2);
                        end
                    end
                end
            end
        end
    end
end

function addZeroTestIndicators(zero_tests, barPositions, groupedData, errorData, star_height_offset, valid_sessions, options)
    % Add indicators for correlations significantly different from zero
    session_fields = fieldnames(zero_tests);
    
    for s = 1:length(session_fields)
        session = session_fields{s};
        session_idx = find(strcmp(valid_sessions, session));
        
        if ~isempty(session_idx) && isfield(zero_tests.(session), 'significant')
            for a = 1:length(zero_tests.(session).significant)
                if zero_tests.(session).significant(a) && a <= size(barPositions, 1)
                    % Add small indicator near the bar
                    x_pos = barPositions(a, session_idx);
                    y_pos = groupedData(a, session_idx) + errorData(a, session_idx) + star_height_offset * 1;  % Much closer (was * 2)
                    
                    p_val = zero_tests.(session).p_values(a);
                    star_symbol = getCorrelationStarSymbol(p_val);
                    
                    % Small star to indicate significantly different from zero
                    text(x_pos, y_pos, '•', 'HorizontalAlignment', 'center', ...
                        'FontSize', options.font_size + 2, 'FontWeight', 'bold', 'Color', [0 0 0]);  % Smaller font
                end
            end
        end
    end
end

function addCorrelationSessionPairwiseLines(sig_comparisons, arm_idx, barPositions, groupedData, errorData, star_height_offset, valid_sessions, options)
    % Add pairwise comparison lines for correlation data
    base_height = max(groupedData(arm_idx, :) + errorData(arm_idx, :));
    line_spacing = star_height_offset * 2;  % Slightly larger spacing (was * 1.5)
    
    for i = 1:size(sig_comparisons, 1)
        comp = sig_comparisons(i, :);
        session1_idx = comp(1);
        session2_idx = comp(2);
        p_val = comp(6);
        
        star_symbol = getCorrelationStarSymbol(p_val);
        
        x1 = barPositions(arm_idx, session1_idx);
        x2 = barPositions(arm_idx, session2_idx);
        line_y = base_height + line_spacing * (i + 1.5);  % Slightly higher start (was i + 1)
        
        % Draw significance line
        plot([x1, x2], [line_y, line_y], 'Color', options.stats_color, 'LineWidth', 2);
        plot([x1, x1], [line_y - line_spacing*0.15, line_y], 'Color', options.stats_color, 'LineWidth', 1.5);  % Slightly longer vertical lines
        plot([x2, x2], [line_y - line_spacing*0.15, line_y], 'Color', options.stats_color, 'LineWidth', 1.5);  % Slightly longer vertical lines
        
        % Add star
        text(mean([x1, x2]), line_y + line_spacing*0.1, star_symbol, ...  % Slightly higher (was 0.05)
            'HorizontalAlignment', 'center', 'FontSize', options.font_size, 'FontWeight', 'bold', 'Color', options.stats_color);
    end
end

function star_symbol = getCorrelationStarSymbol(p_val)
    % Helper function to get star symbol based on p-value
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

function formatCorrelationPlot(arms, valid_sessions, session_labels_map, options)
    % Format the correlation plot
    ylabel('dF/F-Distance Correlation (r)', 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    title([options.title, ' - ', upper(options.group)], 'FontSize', options.font_size + 4, 'FontWeight', 'bold');
    
    % Set x-axis labels
    xticks(1:length(arms));
    xticklabels(arms);
    set(gca, 'FontSize', options.font_size);
    
    % Create session labels for legend
    sessionLabels = cell(1, length(valid_sessions));
    for i = 1:length(valid_sessions)
        if isKey(session_labels_map, valid_sessions{i})
            sessionLabels{i} = session_labels_map(valid_sessions{i});
        else
            sessionLabels{i} = strrep(valid_sessions{i}, 'sess', 'Session ');
        end
    end
    
    % Add legend
    legend(sessionLabels, 'Location', 'northeast', 'FontSize', options.font_size);
    legend('boxoff');
    
    % Set y-axis limits
    ylim(options.ylimit);
    
    % Final formatting
    box off;
    grid off;
    set(gca, 'LineWidth', 1.2);
    set(gcf, 'Color', 'white');
    set(gca, 'Color', 'white');
end

function displayCorrelationStatsSummary(stats_results)
    % Display a summary of correlation statistical results
    fprintf('\n===== Correlation Statistical Summary =====\n');
    if isfield(stats_results, 'fisher_z_used')
        if stats_results.fisher_z_used
            fprintf('Statistical tests performed using Fisher z-transformation\n');
        else
            fprintf('Statistical tests performed using raw correlations\n');
        end
    end
    fprintf('Significant correlations different from zero: %d\n', stats_results.summary.significant_zero_tests);
    fprintf('Significant arm comparisons within sessions: %d\n', stats_results.summary.significant_arm_comparisons);
    fprintf('Significant session comparisons within arms: %d\n', stats_results.summary.significant_session_comparisons);
    
    % Display zero test results (always show raw correlations for interpretation)
    if isfield(stats_results, 'zero_tests')
        fprintf('\nCorrelations significantly different from zero:\n');
        zero_fields = fieldnames(stats_results.zero_tests);
        arms = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
        
        for i = 1:length(zero_fields)
            session = zero_fields{i};
            session_display = strrep(session, 'sess', 'Session ');
            if isfield(stats_results.zero_tests.(session), 'significant')
                for a = 1:length(stats_results.zero_tests.(session).significant)
                    if stats_results.zero_tests.(session).significant(a)
                        fprintf('%s - %s: r = %.3f, p = %.4f\n', ...
                            session_display, arms{a}, ...
                            stats_results.zero_tests.(session).means(a), ...
                            stats_results.zero_tests.(session).p_values(a));
                    end
                end
            end
        end
    end
    
    % Display details for arm comparisons
    session_fields = fieldnames(stats_results.arm_comparisons);
    for i = 1:length(session_fields)
        session = session_fields{i};
        if stats_results.arm_comparisons.(session).significant
            fprintf('%s: Arms comparison p = %.4f (significant)\n', ...
                session, stats_results.arm_comparisons.(session).p_value);
        end
    end
    
    % Display details for session comparisons
    arm_fields = fieldnames(stats_results.session_comparisons);
    for i = 1:length(arm_fields)
        arm_field = arm_fields{i};
        if stats_results.session_comparisons.(arm_field).significant
            if isfield(stats_results.session_comparisons.(arm_field), 'arm_name')
                arm_display_name = stats_results.session_comparisons.(arm_field).arm_name;
            else
                arm_display_name = arm_field;
            end
            fprintf('%s: Sessions comparison p = %.4f (%s, significant)\n', ...
                arm_display_name, stats_results.session_comparisons.(arm_field).p_value, ...
                stats_results.session_comparisons.(arm_field).test_type);
        end
    end
    fprintf('===========================================\n');
end