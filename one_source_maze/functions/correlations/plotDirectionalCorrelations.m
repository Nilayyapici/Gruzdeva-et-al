% Plot correlation data for towards and away movements
function plotDirectionalCorrelations(corr_data, options)
    % Direction types and labels
    directions = {'towards', 'away'};
    direction_labels = {'Towards', 'Away'};
    
    % Define session colors based on direction and source
    % TOWARDS (reddish) - sess0&2 same, sess1&3 same, food darker than gel
    towards_colors = struct();
    towards_colors.gel = {
        [1.0, 0.6, 0.6],     % sess0: Light red
        [0.8, 0.3, 0.3],     % sess1: Medium red
        [1.0, 0.6, 0.6],     % sess2: Light red (same as sess0)
        [0.8, 0.3, 0.3]      % sess3: Medium red (same as sess1)
    };
    towards_colors.food = {
        [0.8, 0.2, 0.2],     % sess0: Dark red
        [0.6, 0.0, 0.0],     % sess1: Very dark red
        [0.8, 0.2, 0.2],     % sess2: Dark red (same as sess0)
        [0.6, 0.0, 0.0]      % sess3: Very dark red (same as sess1)
    };
    
    % AWAY (blueish) - sess0&2 same, sess1&3 same, food darker than gel
    away_colors = struct();
    away_colors.gel = {
        [0.6, 0.8, 1.0],     % sess0: Light blue
        [0.3, 0.6, 0.9],     % sess1: Medium blue
        [0.6, 0.8, 1.0],     % sess2: Light blue (same as sess0)
        [0.3, 0.6, 0.9]      % sess3: Medium blue (same as sess1)
    };
    away_colors.food = {
        [0.2, 0.4, 0.8],     % sess0: Dark blue
        [0.0, 0.2, 0.6],     % sess1: Very dark blue
        [0.2, 0.4, 0.8],     % sess2: Dark blue (same as sess0)
        [0.0, 0.2, 0.6]      % sess3: Very dark blue (same as sess1)
    };
    
    % Group and source combinations
    groups = {'fasted', 'fed'};
    sources = {'food', 'gel'};
    
    % Set up figure positions
    fig_positions = [
        [100, 500, 800, 600], % fasted_food
        [950, 500, 800, 600], % fed_food
        [100, 50, 800, 600],  % fasted_gel
        [950, 50, 800, 600]   % fed_gel
    ];
    
    % Create figures for each group-source combination
    fig_idx = 1;
    
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
            
            % Check if there's any data for this combination
            has_data = false;
            for d = 1:length(directions)
                direction = directions{d};
                if isfield(corr_data.(field_name), direction)
                    for s = 1:length(options.sessions)
                        session = options.sessions{s};
                        if isfield(corr_data.(field_name).(direction), session) && ...
                           isfield(corr_data.(field_name).(direction).(session), 'corr') && ...
                           ~isempty(corr_data.(field_name).(direction).(session).corr)
                            has_data = true;
                            break;
                        end
                    end
                    if has_data
                        break;
                    end
                end
            end
            
            if ~has_data
                continue;
            end
            
            % Create figure for this group-source combination
            figure('Name', sprintf('%s-%s Correlation', group_name, source_name), ...
                   'Position', fig_positions(fig_idx, :));
            fig_idx = fig_idx + 1;
            
            % Create subplot for each direction (towards/away)
            for d = 1:length(directions)
                direction = directions{d};
                
                % Skip if this direction doesn't exist
                if ~isfield(corr_data.(field_name), direction)
                    continue;
                end
                
                % Select color scheme based on direction
                if strcmp(direction, 'towards')
                    session_colors = towards_colors.(source_name);
                else % away
                    session_colors = away_colors.(source_name);
                end
                
                % Create subplot
                subplot(1, 2, d);
                
                % Extract data for this direction across sessions
                all_sessions_data = {};
                session_labels = {};
                session_indices = [];
                
                % Process each session
                for s = 1:length(options.sessions)
                    session = options.sessions{s};
                    
                    % Check if we have data for this session
                    has_session_data = false;
                    if isfield(corr_data.(field_name).(direction), session) && ...
                       isfield(corr_data.(field_name).(direction).(session), 'corr') && ...
                       ~isempty(corr_data.(field_name).(direction).(session).corr)
                        has_session_data = true;
                    end
                    
                    if has_session_data
                        % Get correlation values
                        corr_values = corr_data.(field_name).(direction).(session).corr;
                        
                        % Store session data
                        all_sessions_data{end+1} = corr_values;
                        session_labels{end+1} = strrep(session, 'sess', 'Session ');
                        session_indices(end+1) = s;
                    end
                end
                
                % Skip if no data available for any session
                if isempty(all_sessions_data)
                    continue;
                end
                
                % Calculate means and SEMs
                num_sessions = length(all_sessions_data);
                means = zeros(num_sessions, 1);
                sems = zeros(num_sessions, 1);
                
                for s = 1:num_sessions
                    means(s) = mean(all_sessions_data{s}, 'omitnan');
                    sems(s) = std(all_sessions_data{s}, 'omitnan') / ...
                              sqrt(sum(~isnan(all_sessions_data{s})));
                end
                
                % Create the bar plot
                bar_h = bar(1:num_sessions, means, 'FaceColor', 'flat');
                
                % Set the bar colors based on direction and source
                for s = 1:num_sessions
                    % Get the session index (sess0, sess1, etc.)
                    sess_idx = session_indices(s);
                    bar_h.CData(s,:) = session_colors{sess_idx};
                end
                
                hold on;
                
                % Add error bars
                errorbar(1:num_sessions, means, sems, 'k.', 'LineWidth', 1);
                
                % Add individual data points with small jitter
                for s = 1:num_sessions
                    % Get the session index
                    sess_idx = session_indices(s);
                    
                    % Add small jitter to x-coordinates (±0.1)
                    x_jitter = s + (rand(size(all_sessions_data{s})) - 0.5) * 0.1;
                    
                    % Plot individual points
                    scatter(x_jitter, all_sessions_data{s}, 15, 'o', 'MarkerEdgeColor', 'k', ...
                            'MarkerFaceColor', session_colors{sess_idx}, ...
                            'MarkerFaceAlpha', 0.5);
                end
                
                % Add a horizontal line at y=0
                line(xlim, [0 0], 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1.5);
                
                % Format the plot
                xlabel('Session', 'FontSize', 12);
                ylabel('dF/F-Distance Correlation (r)', 'FontSize', 12);
                title(sprintf('%s', direction_labels{d}), 'FontSize', 14, 'FontWeight', 'bold');
                
                % Set y-axis limits if provided
                if ~isempty(options.ylimit)
                    ylim(options.ylimit);
                end
                
                % Set x-axis ticks and labels
                xticks(1:num_sessions);
                xticklabels({'Sess0', 'Sess1', 'Sess2', 'Sess3'});
                
                % Add grid and box
                grid off;
                box off;
                
                hold off;
            end
            
            % Add overall title for the figure
            group_title = group_name;
            group_title(1) = upper(group_title(1));
            
            source_title = source_name;
            source_title(1) = upper(source_title(1));
            
            sgtitle(sprintf('dF/F-Distance Correlations: %s - %s', group_title, source_title), ...
                   'FontSize', 16, 'FontWeight', 'bold');
        end
    end
end