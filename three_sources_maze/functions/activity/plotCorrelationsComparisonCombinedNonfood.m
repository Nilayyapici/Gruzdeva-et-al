function plotCorrelationsComparisonCombinedNonfood(mice_all, options)
    % Enhanced plotting function for dF/F-distance correlations with combined non-food arms
    % Creates publication-ready plots with significance indicators comparing Food vs Combined Non-Food
    % Now supports filtering by both group and memory type
    
    % Set defaults for options
    if nargin < 2
        options = struct();
    end
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'group'), options.group = 'all'; end
    if ~isfield(options, 'memory'), options.memory = 'all'; end  % NEW: Memory filtering
    if ~isfield(options, 'title'), options.title = 'dF/F-Distance Correlations: Food vs Combined Non-Food'; end
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end
    if ~isfield(options, 'dist_lim'), options.dist_lim = 200; end
    if ~isfield(options, 'dist_too_close'), options.dist_too_close = 5; end
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
    if ~isfield(options, 'use_fisher_z'), options.use_fisher_z = true; end
    if ~isfield(options, 'min_mice_per_comparison'), options.min_mice_per_comparison = 3; end
    
    % NEW: Filter mice_all by memory type before analysis
    mice_all_filtered = filterMiceByMemory(mice_all, options.memory);
    
    % Display filtering information
    if ~strcmp(options.memory, 'all')
        fprintf('Filtering for memory type: %s\n', options.memory);
        fprintf('Original data: %d entries, Filtered data: %d entries\n', ...
                size(mice_all, 1), size(mice_all_filtered, 1));
    end
    
    % Create session labels mapping
    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, ...
                                       {'Before', 'Learning', 'Test'});
    
    % First, get the standard correlation data using the existing function
    options_internal = options;
    options_internal.suppress_figures = true;  % Don't show the original plots
    options_internal.run_stats = false;       % We'll handle stats ourselves
    
    % Get the original correlation data using filtered data and existing function
    corr_data = analyzeDffDistanceCorrelations(mice_all_filtered, options_internal);
    
    if options.suppress_figures
        return; % Exit early if figures are suppressed
    end
    
    % Check if we got valid correlation data
    if isempty(corr_data) || ~isfield(corr_data, 'corr')
        error('No correlation data returned. Check if filtered data contains valid entries.');
    end
    
    % Check for valid sessions
    valid_sessions = intersect(fieldnames(corr_data.corr)', options.sessions);
    if isempty(valid_sessions)
        error('No valid sessions found for plotting correlations.');
    end
    
    % Check if we have actual data in the sessions
    has_data = false;
    for i = 1:length(valid_sessions)
        session = valid_sessions{i};
        if isfield(corr_data.corr, session) && ~isempty(corr_data.corr.(session))
            has_data = true;
            break;
        end
    end
    
    if ~has_data
        fprintf('Warning: No correlation data found for memory type "%s" and group "%s"\n', options.memory, options.group);
        fprintf('Available entries in filtered data: %d\n', size(mice_all_filtered, 1));
        return;
    end
    
    % Transform the correlation data to combine non-food arms
    [combined_corr_data, combined_z_data] = combineNonfoodArms(corr_data, valid_sessions, options);
    
    % Define arms for the combined analysis
    arms = {'Food Arm', 'Non-Food Arms'};
    numArms = length(arms);
    
    % Organize data for plotting and statistics
    sessionData = cell(1, length(valid_sessions));      % Raw correlations for display
    sessionDataStats = cell(1, length(valid_sessions)); % Data for statistics (raw or Fisher z)
    sessionMeans = cell(1, length(valid_sessions));
    sessionSEMs = cell(1, length(valid_sessions));
    mouseIds = cell(1, length(valid_sessions));
    
    for i = 1:length(valid_sessions)
        session = valid_sessions{i};
        
        % Use combined correlation data
        sessionData{i} = combined_corr_data.(session);  % Always raw for display
        
        % Choose data for statistics based on option
        if options.use_fisher_z
            sessionDataStats{i} = combined_z_data.(session);  % Use Fisher z for statistics
        else
            sessionDataStats{i} = combined_corr_data.(session);  % Use raw correlations
        end
        
        % Get mouse IDs
        mouseIds{i} = corr_data.mouse_ids.(session);
        
        % Calculate means and SEMs for display (always use raw correlations)
        sessionMeans{i} = nanmean(sessionData{i}, 1);
        sessionSEMs{i} = nanstd(sessionData{i}, 1) ./ sqrt(sum(~isnan(sessionData{i}), 1));
    end
    
    % Collect all mice across sessions
    all_mice = unique([mouseIds{:}]);
    
    % Perform statistical tests if requested
    if options.show_stats
        stats_results = performCombinedCorrelationStatisticalTests(combined_corr_data, sessionDataStats, ...
                                                                  valid_sessions, arms, options);
    else
        stats_results = [];
    end
    
    % Create the main plot with updated title including memory filter info
    createCombinedCorrelationPlot(sessionData, sessionMeans, sessionSEMs, mouseIds, all_mice, ...
                                 arms, valid_sessions, session_labels_map, options, stats_results);
    
    % Display statistical summary if requested
    if options.show_stats && ~isempty(stats_results)
        displayCombinedCorrelationStatsSummary(stats_results);
    end

    exportCombinedNonfoodCorrelationsToExcel(combined_corr_data, corr_data, valid_sessions, options);
end

% NEW FUNCTION: Filter mice_all by memory type
function mice_all_filtered = filterMiceByMemory(mice_all, memory_filter)
    % Filter the mice_all cell array by memory type
    
    if strcmp(memory_filter, 'all')
        mice_all_filtered = mice_all;
        return;
    end
    
    % Check if memory column exists (should be column 6)
    if size(mice_all, 2) < 6
        warning('Memory information not found in mice_all. Returning all data.');
        mice_all_filtered = mice_all;
        return;
    end
    
    % Find entries that match the memory filter
    memory_match = false(size(mice_all, 1), 1);
    
    for i = 1:size(mice_all, 1)
        if strcmp(mice_all{i, 6}, memory_filter)
            memory_match(i) = true;
        end
    end
    
    % Filter the data
    mice_all_filtered = mice_all(memory_match, :);
    
    if isempty(mice_all_filtered)
        warning('No entries found matching memory filter: %s', memory_filter);
    end
end

function mouse_id = extractMouseID(session_info)
    % Extract mouse ID from session info string
    % Example: 'F13_sess0' -> 'F13'
    
    if contains(session_info, '_sess')
        parts = strsplit(session_info, '_sess');
        mouse_id = parts{1};
    else
        mouse_id = session_info;
    end
end

% Include all the other existing functions from your original code

function [combined_corr_data, combined_z_data] = combineNonfoodArms(corr_data, valid_sessions, options)
    % Combine the two non-food arms into a single measure
    
    combined_corr_data = struct();
    combined_z_data = struct();
    
    for i = 1:length(valid_sessions)
        session = valid_sessions{i};
        
        if isfield(corr_data.corr, session) && ~isempty(corr_data.corr.(session))
            % Original data: [food, nonfood1, nonfood2] for each mouse (rows)
            original_corr = corr_data.corr.(session);
            original_z = corr_data.z_corr.(session);
            
            % Check if we have the expected number of columns (3 arms)
            if size(original_corr, 2) < 3
                error('Expected correlation data for 3 arms (food, nonfood1, nonfood2), but got %d columns for session %s', ...
                      size(original_corr, 2), session);
            end
            
            % Extract columns
            food_corr = original_corr(:, 1);           % Food arm correlations
            nonfood1_corr = original_corr(:, 2);       % Non-food arm 1 correlations
            nonfood2_corr = original_corr(:, 3);       % Non-food arm 2 correlations
            
            food_z = original_z(:, 1);                 % Food arm z-transformed
            nonfood1_z = original_z(:, 2);             % Non-food arm 1 z-transformed
            nonfood2_z = original_z(:, 3);             % Non-food arm 2 z-transformed
            
            % Combine non-food arms
            num_mice = size(original_corr, 1);
            combined_nonfood_corr = nan(num_mice, 1);
            combined_nonfood_z = nan(num_mice, 1);
            
            for mouse_idx = 1:num_mice
                % For raw correlations: average if both exist, or use single value
                nf1 = nonfood1_corr(mouse_idx);
                nf2 = nonfood2_corr(mouse_idx);
                
                if ~isnan(nf1) && ~isnan(nf2)
                    combined_nonfood_corr(mouse_idx) = mean([nf1, nf2]);
                elseif ~isnan(nf1)
                    combined_nonfood_corr(mouse_idx) = nf1;
                elseif ~isnan(nf2)
                    combined_nonfood_corr(mouse_idx) = nf2;
                end
                
                % For z-transformed data: average if both exist, or use single value
                nf1_z = nonfood1_z(mouse_idx);
                nf2_z = nonfood2_z(mouse_idx);
                
                if ~isnan(nf1_z) && ~isnan(nf2_z)
                    combined_nonfood_z(mouse_idx) = mean([nf1_z, nf2_z]);
                elseif ~isnan(nf1_z)
                    combined_nonfood_z(mouse_idx) = nf1_z;
                elseif ~isnan(nf2_z)
                    combined_nonfood_z(mouse_idx) = nf2_z;
                end
            end
            
            % Store combined data: [food, combined_nonfood]
            combined_corr_data.(session) = [food_corr, combined_nonfood_corr];
            combined_z_data.(session) = [food_z, combined_nonfood_z];
        else
            % No data for this session, create empty matrices
            combined_corr_data.(session) = [];
            combined_z_data.(session) = [];
        end
    end
end

function stats_results = performCombinedCorrelationStatisticalTests(combined_corr_data, sessionDataStats, ...
                                                                   valid_sessions, arms, options)
    % Perform comprehensive statistical tests for combined correlation data
    stats_results = struct();
    stats_results.arm_comparisons = struct();
    stats_results.session_comparisons = struct();
    stats_results.summary = struct();
    stats_results.zero_tests = struct();
    stats_results.fisher_z_used = options.use_fisher_z;
    
    numArms = length(arms);
    numSessions = length(valid_sessions);
    
    % Initialize summary counters
    stats_results.summary.significant_zero_tests = 0;
    stats_results.summary.significant_arm_comparisons = 0;
    stats_results.summary.significant_session_comparisons = 0;
    
    % 1. Test if correlations are significantly different from zero
    for s = 1:numSessions
        session = valid_sessions{s};
        if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s})
            session_stats_data = sessionDataStats{s};  % Use stats data (raw or Fisher z)
            
            stats_results.zero_tests.(session) = struct();
            stats_results.zero_tests.(session).means = nan(1, numArms);
            stats_results.zero_tests.(session).p_values = ones(1, numArms);
            stats_results.zero_tests.(session).significant = false(1, numArms);
            
            for a = 1:numArms
                arm_values = session_stats_data(:, a);
                valid_values = arm_values(~isnan(arm_values));
                
                if length(valid_values) >= options.min_mice_per_comparison
                    % One-sample t-test against zero
                    [~, p] = ttest(valid_values, 0);
                    stats_results.zero_tests.(session).p_values(a) = p;
                    stats_results.zero_tests.(session).significant(a) = p < options.alpha;
                    
                    % Store the display mean (always raw correlation for interpretation)
                    if isfield(combined_corr_data, session)
                        raw_values = combined_corr_data.(session)(:, a);
                        valid_raw = raw_values(~isnan(raw_values));
                        if ~isempty(valid_raw)
                            stats_results.zero_tests.(session).means(a) = mean(valid_raw);
                        end
                    end
                    
                    if stats_results.zero_tests.(session).significant(a)
                        stats_results.summary.significant_zero_tests = ...
                            stats_results.summary.significant_zero_tests + 1;
                    end
                end
            end
        end
    end
    
    % 2. Compare arms within each session (Food vs Combined Non-Food)
    for s = 1:numSessions
        session = valid_sessions{s};
        stats_results.arm_comparisons.(session) = struct();
        
        if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s}) && numArms >= 2
            session_stats_data = sessionDataStats{s};
            
            % Get data for both arms
            arm1_data = session_stats_data(:, 1);  % Food
            arm2_data = session_stats_data(:, 2);  % Combined Non-Food
            
            arm1_valid = arm1_data(~isnan(arm1_data));
            arm2_valid = arm2_data(~isnan(arm2_data));
            
            if length(arm1_valid) >= options.min_mice_per_comparison && ...
               length(arm2_valid) >= options.min_mice_per_comparison
                
                [~, p] = ttest2(arm1_valid, arm2_valid);
                stats_results.arm_comparisons.(session).p_value = p;
                stats_results.arm_comparisons.(session).significant = p < options.alpha;
                
                if stats_results.arm_comparisons.(session).significant
                    stats_results.summary.significant_arm_comparisons = ...
                        stats_results.summary.significant_arm_comparisons + 1;
                end
            else
                stats_results.arm_comparisons.(session).p_value = 1;
                stats_results.arm_comparisons.(session).significant = false;
            end
        else
            stats_results.arm_comparisons.(session).p_value = 1;
            stats_results.arm_comparisons.(session).significant = false;
        end
    end
    
    % 3. Compare sessions within each arm type
    for a = 1:numArms
        arm_name = lower(strrep(strrep(arms{a}, ' ', '_'), '-', '_'));
        field_name = sprintf('%s_comparisons', arm_name);
        stats_results.session_comparisons.(field_name) = struct();
        
        % Collect data across sessions for this arm
        arm_data_by_session = {};
        session_names = {};
        
        for s = 1:numSessions
            if s <= length(sessionDataStats) && ~isempty(sessionDataStats{s})
                arm_values = sessionDataStats{s}(:, a);
                valid_values = arm_values(~isnan(arm_values));
                
                if length(valid_values) >= options.min_mice_per_comparison
                    arm_data_by_session{end+1} = valid_values;
                    session_names{end+1} = valid_sessions{s};
                end
            end
        end
        
        % ANOVA across sessions for this arm
        if length(arm_data_by_session) >= 2
            % Prepare data for ANOVA
            all_data = [];
            group_labels = [];
            
            for sess_idx = 1:length(arm_data_by_session)
                data = arm_data_by_session{sess_idx};
                all_data = [all_data; data(:)];
                group_labels = [group_labels; repmat(sess_idx, length(data), 1)];
            end
            
            if length(unique(group_labels)) >= 2
                p = anova1(all_data, group_labels, 'off');
                stats_results.session_comparisons.(field_name).p_value = p;
                stats_results.session_comparisons.(field_name).significant = p < options.alpha;
                stats_results.session_comparisons.(field_name).sessions_compared = session_names;
                
                if stats_results.session_comparisons.(field_name).significant
                    stats_results.summary.significant_session_comparisons = ...
                        stats_results.summary.significant_session_comparisons + 1;
                end
            else
                stats_results.session_comparisons.(field_name).p_value = 1;
                stats_results.session_comparisons.(field_name).significant = false;
            end
        else
            stats_results.session_comparisons.(field_name).p_value = 1;
            stats_results.session_comparisons.(field_name).significant = false;
        end
    end
end

function createCombinedCorrelationPlot(sessionData, sessionMeans, sessionSEMs, mouseIds, all_mice, ...
                                      arms, valid_sessions, session_labels_map, options, stats_results)
    % Create the main correlation plot with combined non-food arms
    
    figure('Position', options.figure_size);
    hold on;
    
    numArms = length(arms);
    numSessions = length(valid_sessions);
    
    % Prepare data for grouped bar plot
    groupedData = zeros(numArms, numSessions);
    errorData = zeros(numArms, numSessions);
    barPositions = zeros(numArms, numSessions);
    
    for i = 1:numArms
        for j = 1:numSessions
            if i <= length(sessionMeans{j})
                groupedData(i, j) = sessionMeans{j}(i);
                errorData(i, j) = sessionSEMs{j}(i);
            end
        end
    end
    
    % Create grouped bar plot with proper styling
    h = bar(groupedData, options.bar_width);
    
    % Set colors for each session bar
    for i = 1:numSessions
        h(i).FaceColor = options.colors{mod(i-1, length(options.colors))+1};
        h(i).EdgeColor = 'none';
        h(i).FaceAlpha = 0.8;
    end
    
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
    
    % Add error bars and individual data points
    numBars = length(valid_sessions);
    groupwidth = min(options.bar_width, numBars/(numBars+1.5));
    
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
    
    % Connect points across sessions if requested
    if options.connect_points && numSessions > 1
        connectCombinedMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, ...
                                      valid_sessions, numArms);
    end
    
    % Add significance indicators if statistics were performed
    if options.show_stats && ~isempty(stats_results)
        addCombinedCorrelationSignificanceIndicators(stats_results, barPositions, groupedData, ...
                                                     errorData, options, valid_sessions, arms);
    end
    
    % Format plot
    formatCombinedCorrelationPlot(arms, valid_sessions, session_labels_map, options, h);
    
    hold off;
end

function connectCombinedMouseDataPoints(sessionData, mouseIds, barPositions, all_mice, ...
                                       valid_sessions, numArms)
    % Connect individual mouse correlation data points across sessions
    for arm_idx = 1:numArms
        for m = 1:length(all_mice)
            current_mouse = all_mice{m};
            
            % Collect data points for this mouse across sessions
            mouse_x = [];
            mouse_y = [];
            
            for sess_idx = 1:length(valid_sessions)
                mouse_in_session = find(strcmp(mouseIds{sess_idx}, current_mouse));
                
                if ~isempty(mouse_in_session) && arm_idx <= size(sessionData{sess_idx}, 2)
                    correlation_val = sessionData{sess_idx}(mouse_in_session, arm_idx);
                    
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

function addCombinedCorrelationSignificanceIndicators(stats_results, barPositions, groupedData, ...
                                                     errorData, options, valid_sessions, arms)
    % Add significance indicators for combined correlation comparisons
    star_height_offset = (max(options.ylimit) - min(options.ylimit)) * 0.025;
    
    % 1. Add stars for arm comparisons within sessions (Food vs Combined Non-Food)
    if isfield(stats_results, 'arm_comparisons')
        session_fields = fieldnames(stats_results.arm_comparisons);
        
        for s = 1:length(session_fields)
            session = session_fields{s};
            session_idx = find(strcmp(valid_sessions, session));
            
            if ~isempty(session_idx) && stats_results.arm_comparisons.(session).significant
                % Position for significance bar between arms
                x1 = barPositions(1, session_idx);  % Food arm position
                x2 = barPositions(2, session_idx);  % Combined non-food arm position
                
                % Calculate height based on data
                max_height = max([groupedData(1, session_idx) + errorData(1, session_idx), ...
                                 groupedData(2, session_idx) + errorData(2, session_idx)]);
                
                line_y = max_height + star_height_offset * 2;
                line_spacing = star_height_offset * 0.7;
                
                % Draw significance bar
                plot([x1, x1, x2, x2], [line_y, line_y + line_spacing, line_y + line_spacing, line_y], ...
                     'Color', options.stats_color, 'LineWidth', 1.5);
                
                % Add star
                p_val = stats_results.arm_comparisons.(session).p_value;
                star_symbol = getCorrelationStarSymbol(p_val);
                
                text(mean([x1, x2]), line_y + line_spacing * 1.5, star_symbol, ...
                    'HorizontalAlignment', 'center', 'FontSize', options.font_size, ...
                    'FontWeight', 'bold', 'Color', options.stats_color);
            end
        end
    end
    
    % 2. Add stars for session comparisons within arms (across sessions for same arm)
    if isfield(stats_results, 'session_comparisons')
        session_comp_fields = fieldnames(stats_results.session_comparisons);
        
        for f = 1:length(session_comp_fields)
            field_name = session_comp_fields{f};
            
            if stats_results.session_comparisons.(field_name).significant
                % Determine which arm this comparison is for
                if contains(field_name, 'food')
                    arm_idx = 1;  % Food arm
                elseif contains(field_name, 'non_food') || contains(field_name, 'nonfood')
                    arm_idx = 2;  % Combined non-food arm
                else
                    continue;  % Skip if arm can't be determined
                end
                
                % Find the highest point for this arm across all sessions
                max_height = 0;
                for sess_idx = 1:length(valid_sessions)
                    height = groupedData(arm_idx, sess_idx) + errorData(arm_idx, sess_idx);
                    if height > max_height
                        max_height = height;
                    end
                end
                
                % Position the session comparison star above the arm
                x_center = arm_idx;  % Center on the arm
                line_y = max_height + star_height_offset * 8;
                
                % Add star for session comparison
                p_val = stats_results.session_comparisons.(field_name).p_value;
                star_symbol = getCorrelationStarSymbol(p_val);
                
                text(x_center, line_y, star_symbol, ...
                    'HorizontalAlignment', 'center', 'FontSize', options.font_size, ...
                    'FontWeight', 'bold', 'Color', options.stats_color);
                
                % Add a horizontal line above the arm to indicate session comparison
                x_span = 0.2;  % Width of the line
                plot([x_center - x_span, x_center + x_span], [line_y - star_height_offset * 0.5, line_y - star_height_offset * 0.5], ...
                     'Color', options.stats_color, 'LineWidth', 2);
            end
        end
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

function displayCombinedCorrelationStatsSummary(stats_results)
    % Display a summary of correlation statistical results for combined non-food analysis
    fprintf('\n===== Combined Non-Food Correlation Statistical Summary =====\n');
    if isfield(stats_results, 'fisher_z_used')
        if stats_results.fisher_z_used
            fprintf('Statistical tests performed using Fisher z-transformation\n');
        else
            fprintf('Statistical tests performed using raw correlations\n');
        end
    end
    fprintf('Significant correlations different from zero: %d\n', stats_results.summary.significant_zero_tests);
    fprintf('Significant arm comparisons (Food vs Combined Non-Food): %d\n', stats_results.summary.significant_arm_comparisons);
    fprintf('Significant session comparisons within arms: %d\n', stats_results.summary.significant_session_comparisons);
    
    % Display zero test results
    if isfield(stats_results, 'zero_tests')
        fprintf('\nCorrelations significantly different from zero:\n');
        zero_fields = fieldnames(stats_results.zero_tests);
        arms = {'Food', 'Combined Non-Food'};
        
        for i = 1:length(zero_fields)
            session = zero_fields{i};
            session_display = strrep(session, 'sess', 'Session ');
            if isfield(stats_results.zero_tests.(session), 'significant')
                for a = 1:length(stats_results.zero_tests.(session).significant)
                    if a <= length(arms) && stats_results.zero_tests.(session).significant(a)
                        fprintf('%s - %s: r = %.3f, p = %.4f\n', ...
                            session_display, arms{a}, ...
                            stats_results.zero_tests.(session).means(a), ...
                            stats_results.zero_tests.(session).p_values(a));
                    end
                end
            end
        end
    end
    
    % Display details for arm comparisons (Food vs Combined Non-Food)
    if isfield(stats_results, 'arm_comparisons')
        session_fields = fieldnames(stats_results.arm_comparisons);
        for i = 1:length(session_fields)
            session = session_fields{i};
            if stats_results.arm_comparisons.(session).significant
                session_display = strrep(session, 'sess', 'Session ');
                fprintf('%s: Food vs Combined Non-Food comparison p = %.4f (significant)\n', ...
                    session_display, stats_results.arm_comparisons.(session).p_value);
            end
        end
    end
    
    % Display details for session comparisons
    if isfield(stats_results, 'session_comparisons')
        arm_fields = fieldnames(stats_results.session_comparisons);
        for i = 1:length(arm_fields)
            arm_field = arm_fields{i};
            if stats_results.session_comparisons.(arm_field).significant
                arm_display = strrep(strrep(arm_field, '_comparisons', ''), '_', ' ');
                arm_display = [upper(arm_display(1)), arm_display(2:end)];
                fprintf('%s across sessions: p = %.4f (significant)\n', ...
                    arm_display, stats_results.session_comparisons.(arm_field).p_value);
            end
        end
    end
    fprintf('=========================================================\n');
end

% UPDATED FUNCTION: Enhanced plot formatting with memory info
function formatCombinedCorrelationPlot(arms, valid_sessions, session_labels_map, options, h)
    % Format the correlation plot for combined non-food analysis with memory info
    ylabel('dF/F-Distance Correlation (r)', 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    % Format the group name properly - replace underscores with spaces and use title case
    formatted_group_name = formatGroupName(options.group);
    
    % Create title with both group and memory information
    title_parts = {options.title};
    if ~strcmp(options.group, 'all')
        title_parts{end+1} = formatted_group_name;
    end
    if ~strcmp(options.memory, 'all')
        formatted_memory = formatGroupName(options.memory);
        title_parts{end+1} = formatted_memory;
    end
    
    full_title = strjoin(title_parts, ' - ');
    title(full_title, 'FontSize', options.font_size + 4, 'FontWeight', 'bold');
    
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
    
    % Add legend using the bar handles only
    legend(h, sessionLabels, 'Location', 'northeast', 'FontSize', options.font_size);
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

function formatted_name = formatGroupName(group_name)
    % Helper function to format group names for display
    % Replaces underscores with spaces and applies title case
    
    if isempty(group_name)
        formatted_name = '';
        return;
    end
    
    % Replace underscores with spaces
    formatted_name = strrep(group_name, '_', ' ');
    
    % Apply title case (capitalize first letter of each word)
    words = strsplit(formatted_name, ' ');
    for i = 1:length(words)
        if ~isempty(words{i})
            words{i} = [upper(words{i}(1)), lower(words{i}(2:end))];
        end
    end
    formatted_name = strjoin(words, ' ');
end

function exportCombinedNonfoodCorrelationsToExcel(combined_corr_data, corr_data, valid_sessions, options)
    % Export per-mouse correlation values (Food Arm & Non-Food Arms) to Excel
    
    arm_labels = {'Food_Arm', 'NonFood_Arms'};
    
    % Session label map
    sess_label_map = containers.Map({'sess0','sess1','sess2'}, ...
                                    {'Session0','Session1','Session2'});
    
    % Collect all unique mouse IDs across sessions
    all_mouse_ids = {};
    for s = 1:length(valid_sessions)
        sess = valid_sessions{s};
        ids = corr_data.mouse_ids.(sess);
        all_mouse_ids = union(all_mouse_ids, ids);
    end
    all_mouse_ids = sort(all_mouse_ids);
    
    % Build column names
    col_names = {'Mouse_ID'};
    for s = 1:length(valid_sessions)
        sess = valid_sessions{s};
        if isKey(sess_label_map, sess)
            sess_str = sess_label_map(sess);
        else
            sess_str = strrep(sess, 'sess', 'Session');
        end
        for a = 1:length(arm_labels)
            col_names{end+1} = sprintf('%s_%s', sess_str, arm_labels{a});
        end
    end
    
    % Fill data matrix
    n_mice = length(all_mouse_ids);
    n_cols  = length(col_names) - 1;
    data_matrix = NaN(n_mice, n_cols);
    
    col_idx = 1;
    for s = 1:length(valid_sessions)
        sess = valid_sessions{s};
        sess_mouse_ids = corr_data.mouse_ids.(sess);
        sess_data = combined_corr_data.(sess);   % [n_mice_in_sess x 2]
        
        for a = 1:length(arm_labels)
            for m = 1:n_mice
                mouse_row = find(strcmp(sess_mouse_ids, all_mouse_ids{m}));
                if ~isempty(mouse_row) && a <= size(sess_data, 2)
                    data_matrix(m, col_idx) = sess_data(mouse_row(1), a);
                end
            end
            col_idx = col_idx + 1;
        end
    end
    
    % Build MATLAB table
    export_table = array2table(data_matrix, 'VariableNames', col_names(2:end));
    export_table.Mouse_ID = all_mouse_ids(:);
    export_table = [export_table(:,end), export_table(:,1:end-1)];
    
    % Session label for filename
    if isequal(sort(valid_sessions), sort({'sess0','sess2'}))
        sess_label = 'Test';
    elseif isequal(sort(valid_sessions), sort({'sess0','sess1'}))
        sess_label = 'Learning';
    else
        sess_label = strjoin(valid_sessions, '_');
    end
    
    excel_filename = sprintf('correlations_combined_nonfood_%s_%s_%s.xlsx', ...
        options.group, sess_label, datestr(now, 'yyyymmdd'));
    
    writetable(export_table, excel_filename);
    fprintf('Correlation data exported to: %s\n', excel_filename);
    disp(export_table);
end