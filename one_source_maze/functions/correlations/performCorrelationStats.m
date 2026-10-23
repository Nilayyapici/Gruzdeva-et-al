function performCorrelationStats(corr_data, options)
    % PERFORMCORRELATIONSTATS Perform statistical analyses on directional correlation data
    %   Compares correlations across sessions, between directions, and between groups/sources
    
    % Direction types and labels
    directions = {'towards', 'away'};
    direction_labels = {'Towards Food', 'Away From Food'};
    
    % Group and source combinations
    groups = {'fasted', 'fed'};
    sources = {'food', 'gel'};
    
    % Print header
    fprintf('\n===== Statistical Analysis of dF/F-Distance Correlations =====\n');
    
    % 1. Analyze each group-source combination
    for g = 1:length(groups)
        for src = 1:length(sources)
            group_name = groups{g};
            source_name = sources{src};
            field_name = [group_name, '_', source_name];
            
            % Skip if this combination doesn't match our filters or has no data
            if (~strcmp(options.group, 'all') && ~strcmp(group_name, options.group)) || ...
               (~strcmp(options.source, 'all') && ~strcmp(source_name, options.source)) || ...
               ~isfield(corr_data, field_name)
                continue;
            end
            
            fprintf('\n----- %s - %s -----\n', ...
                    upper(group_name(1)), group_name(2:end), ...
                    upper(source_name(1)), source_name(2:end));
            
            % Analyze each direction (towards/away)
            for d = 1:length(directions)
                direction = directions{d};
                
                % Skip if this direction doesn't exist
                if ~isfield(corr_data.(field_name), direction)
                    continue;
                end
                
                fprintf('\n* %s\n', direction_labels{d});
                fprintf('  Comparing across sessions:\n');
                
                % Extract data for all sessions
                all_session_data = [];
                session_indices = [];
                valid_sessions = {};
                
                for s = 1:length(options.sessions)
                    session = options.sessions{s};
                    
                    % Check if we have data for this session
                    has_session_data = false;
                    if isfield(corr_data.(field_name).(direction), session) && ...
                       isfield(corr_data.(field_name).(direction).(session), 'z_corr') && ...
                       ~isempty(corr_data.(field_name).(direction).(session).z_corr)
                        has_session_data = true;
                    end
                    
                    if has_session_data
                        % Get z-transformed correlation values (for better statistics)
                        z_values = corr_data.(field_name).(direction).(session).z_corr;
                        
                        % Add to combined data
                        all_session_data = [all_session_data; z_values];
                        session_indices = [session_indices; repmat(s, size(z_values))];
                        valid_sessions{end+1} = session;
                    end
                end
                
                % Skip if insufficient data
                if length(unique(session_indices)) < 2 || length(all_session_data) < 4
                    fprintf('    Insufficient data for statistical analysis.\n');
                    continue;
                end
                
                % Perform one-way ANOVA if we have enough data
                [p, tbl, stats] = anova1(all_session_data, session_indices, 'off');
                
                % Report results
                fprintf('    ANOVA p-value: %.4f\n', p);
                
                % Post-hoc test if significant
                if p < 0.05
                    fprintf('    Significant differences found between sessions. Performing Tukey HSD test.\n');
                    [c, m, h, gnames] = multcompare(stats, 'Alpha', 0.05, 'Display', 'off');
                    
                    % Print comparison results
                    fprintf('    Multiple Comparison Results:\n');
                    fprintf('    %-12s %-12s %-10s %-10s %-10s %-10s\n', ...
                            'Session 1', 'Session 2', 'Lower CI', 'Diff', 'Upper CI', 'p-value');
                    
                    for i = 1:size(c, 1)
                        sess1 = strrep(valid_sessions{c(i, 1)}, 'sess', 'Session ');
                        sess2 = strrep(valid_sessions{c(i, 2)}, 'sess', 'Session ');
                        fprintf('    %-12s %-12s %-10.4f %-10.4f %-10.4f %-10.4f\n', ...
                                sess1, sess2, c(i, 3), c(i, 4), c(i, 5), c(i, 6));
                    end
                else
                    fprintf('    No significant differences found across sessions.\n');
                end
            end
        end
    end
    
    % 2. Compare directions (towards vs away) within each group-source-session combination
    fprintf('\n\n===== Comparing Towards vs Away Movements =====\n');
    
    for g = 1:length(groups)
        for src = 1:length(sources)
            group_name = groups{g};
            source_name = sources{src};
            field_name = [group_name, '_', source_name];
            
            % Skip if this combination doesn't match our filters or has no data
            if (~strcmp(options.group, 'all') && ~strcmp(group_name, options.group)) || ...
               (~strcmp(options.source, 'all') && ~strcmp(source_name, options.source)) || ...
               ~isfield(corr_data, field_name)
                continue;
            end
            
            % Check if both directions exist
            if ~isfield(corr_data.(field_name), 'towards') || ...
               ~isfield(corr_data.(field_name), 'away')
                continue;
            end
            
            fprintf('\n----- %s - %s -----\n', ...
                    upper(group_name(1)), group_name(2:end), ...
                    upper(source_name(1)), source_name(2:end));
            
            % Process each session
            for s = 1:length(options.sessions)
                session = options.sessions{s};
                
                % Check if we have data for both directions in this session
                has_towards_data = false;
                has_away_data = false;
                
                if isfield(corr_data.(field_name).towards, session) && ...
                   isfield(corr_data.(field_name).towards.(session), 'z_corr') && ...
                   ~isempty(corr_data.(field_name).towards.(session).z_corr)
                    has_towards_data = true;
                end
                
                if isfield(corr_data.(field_name).away, session) && ...
                   isfield(corr_data.(field_name).away.(session), 'z_corr') && ...
                   ~isempty(corr_data.(field_name).away.(session).z_corr)
                    has_away_data = true;
                end
                
                if ~has_towards_data || ~has_away_data
                    continue;
                end
                
                fprintf('\n* %s\n', strrep(session, 'sess', 'Session '));
                
                % Get z-transformed correlation values
                towards_z = corr_data.(field_name).towards.(session).z_corr;
                away_z = corr_data.(field_name).away.(session).z_corr;
                
                % Get mouse IDs
                towards_ids = corr_data.(field_name).towards.(session).mouse_ids;
                away_ids = corr_data.(field_name).away.(session).mouse_ids;
                
                % Check if we can do a paired test (same mice in both directions)
                [common_ids, towards_idx, away_idx] = intersect(towards_ids, away_ids);
                
                if ~isempty(common_ids)
                    % Perform paired t-test on matched data
                    fprintf('  Found %d mice with both towards and away data. Performing paired t-test.\n', ...
                            length(common_ids));
                    
                    % Extract matched data
                    paired_towards = towards_z(towards_idx);
                    paired_away = away_z(away_idx);
                    
                    [h, p, ~, stats] = ttest(paired_towards, paired_away);
                    
                    % Report results
                    fprintf('  Paired t-test results: t(%d) = %.2f, p = %.4f\n', ...
                            stats.df, stats.tstat, p);
                    
                    if p < 0.05
                        fprintf('  Significant difference found between towards and away movements.\n');
                        
                        % Report means
                        fprintf('  Towards mean: %.4f, Away mean: %.4f\n', ...
                                mean(paired_towards), mean(paired_away));
                    else
                        fprintf('  No significant difference found between towards and away movements.\n');
                    end
                else
                    % Perform unpaired t-test
                    fprintf('  No common mice between towards and away data. Performing unpaired t-test.\n');
                    
                    [h, p, ~, stats] = ttest2(towards_z, away_z);
                    
                    % Report results
                    fprintf('  t-test results: t(%d) = %.2f, p = %.4f\n', ...
                            stats.df, stats.tstat, p);
                    if p < 0.05
                        fprintf('  Significant difference found between towards and away movements.\n');
                        
                        % Report means
                        fprintf('  Towards mean: %.4f, Away mean: %.4f\n', ...
                                mean(towards_z, 'omitnan'), mean(away_z, 'omitnan'));
                    else
                        fprintf('  No significant difference found between towards and away movements.\n');
                    end
                end
            end
        end
    end
    
    % 3. Compare across groups (fasted vs fed) for the same source and direction
    if strcmp(options.group, 'all')
        fprintf('\n\n===== Comparing Fasted vs Fed Conditions =====\n');
        
        for src = 1:length(sources)
            source_name = sources{src};
            
            % Skip if this source doesn't match our filter
            if ~strcmp(options.source, 'all') && ~strcmp(source_name, options.source)
                continue;
            end
            
            fprintf('\n----- %s -----\n', upper(source_name(1)), source_name(2:end));
            
            % Compare for each direction
            for d = 1:length(directions)
                direction = directions{d};
                
                fprintf('\n* %s\n', direction_labels{d});
                
                % Compare across each session
                for s = 1:length(options.sessions)
                    session = options.sessions{s};
                    
                    % Check if we have data for both groups in this session
                    field_name_fasted = ['fasted_', source_name];
                    field_name_fed = ['fed_', source_name];
                    
                    if isfield(corr_data, field_name_fasted) && ...
                       isfield(corr_data, field_name_fed)
                        
                        has_fasted_data = false;
                        has_fed_data = false;
                        
                        % Check fasted data
                        if isfield(corr_data.(field_name_fasted), direction) && ...
                           isfield(corr_data.(field_name_fasted).(direction), session) && ...
                           isfield(corr_data.(field_name_fasted).(direction).(session), 'z_corr') && ...
                           ~isempty(corr_data.(field_name_fasted).(direction).(session).z_corr)
                            has_fasted_data = true;
                        end
                        
                        % Check fed data
                        if isfield(corr_data.(field_name_fed), direction) && ...
                           isfield(corr_data.(field_name_fed).(direction), session) && ...
                           isfield(corr_data.(field_name_fed).(direction).(session), 'z_corr') && ...
                           ~isempty(corr_data.(field_name_fed).(direction).(session).z_corr)
                            has_fed_data = true;
                        end
                        
                        if ~has_fasted_data || ~has_fed_data
                            continue;
                        end
                        
                        fprintf('\n  %s\n', strrep(session, 'sess', 'Session '));
                        
                        % Get z-transformed correlation values
                        fasted_z = corr_data.(field_name_fasted).(direction).(session).z_corr;
                        fed_z = corr_data.(field_name_fed).(direction).(session).z_corr;
                        
                        % Perform unpaired t-test (groups are independent)
                        [h, p, ~, stats] = ttest2(fasted_z, fed_z);
                        
                        % Report results
                        fprintf('    t-test results: t(%d) = %.2f, p = %.4f\n', ...
                                stats.df, stats.tstat, p);
                        
                        if p < 0.05
                            fprintf('    Significant difference found between fasted and fed conditions.\n');
                            
                            % Report means
                            fprintf('    Fasted mean: %.4f, Fed mean: %.4f\n', ...
                                    mean(fasted_z, 'omitnan'), mean(fed_z, 'omitnan'));
                        else
                            fprintf('    No significant difference found between fasted and fed conditions.\n');
                        end
                    end
                end
            end
        end
    end
    
    % 4. Compare across sources (food vs gel) for the same group and direction
    if strcmp(options.source, 'all')
        fprintf('\n\n===== Comparing Food vs Gel Sources =====\n');
        
        for g = 1:length(groups)
            group_name = groups{g};
            
            % Skip if this group doesn't match our filter
            if ~strcmp(options.group, 'all') && ~strcmp(group_name, options.group)
                continue;
            end
            
            fprintf('\n----- %s -----\n', upper(group_name(1)), group_name(2:end));
            
            % Compare for each direction
            for d = 1:length(directions)
                direction = directions{d};
                
                fprintf('\n* %s\n', direction_labels{d});
                
                % Compare across each session
                for s = 1:length(options.sessions)
                    session = options.sessions{s};
                    
                    % Check if we have data for both sources in this session
                    field_name_food = [group_name, '_food'];
                    field_name_gel = [group_name, '_gel'];
                    
                    if isfield(corr_data, field_name_food) && ...
                       isfield(corr_data, field_name_gel)
                        
                        has_food_data = false;
                        has_gel_data = false;
                        
                        % Check food data
                        if isfield(corr_data.(field_name_food), direction) && ...
                           isfield(corr_data.(field_name_food).(direction), session) && ...
                           isfield(corr_data.(field_name_food).(direction).(session), 'z_corr') && ...
                           ~isempty(corr_data.(field_name_food).(direction).(session).z_corr)
                            has_food_data = true;
                        end
                        
                        % Check gel data
                        if isfield(corr_data.(field_name_gel), direction) && ...
                           isfield(corr_data.(field_name_gel).(direction), session) && ...
                           isfield(corr_data.(field_name_gel).(direction).(session), 'z_corr') && ...
                           ~isempty(corr_data.(field_name_gel).(direction).(session).z_corr)
                            has_gel_data = true;
                        end
                        
                        if ~has_food_data || ~has_gel_data
                            continue;
                        end
                        
                        fprintf('\n  %s\n', strrep(session, 'sess', 'Session '));
                        
                        % Get z-transformed correlation values
                        food_z = corr_data.(field_name_food).(direction).(session).z_corr;
                        gel_z = corr_data.(field_name_gel).(direction).(session).z_corr;
                        
                        % Perform unpaired t-test (sources are independent)
                        [h, p, ~, stats] = ttest2(food_z, gel_z);
                        
                        % Report results
                        fprintf('    t-test results: t(%d) = %.2f, p = %.4f\n', ...
                                stats.df, stats.tstat, p);
                        
                        if p < 0.05
                            fprintf('    Significant difference found between food and gel sources.\n');
                            
                            % Report means
                            fprintf('    Food mean: %.4f, Gel mean: %.4f\n', ...
                                    mean(food_z, 'omitnan'), mean(gel_z, 'omitnan'));
                        else
                            fprintf('    No significant difference found between food and gel sources.\n');
                        end
                    end
                end
            end
        end
    end
end