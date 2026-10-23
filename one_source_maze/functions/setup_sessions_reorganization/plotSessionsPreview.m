function plotSessionsPreview(mice_all, mouse_indices)
    % PLOTSESSIONSPREVIEW Visualizes how original data will be split into sessions
    % Displays dF/F signal, door states, and discovery point for each mouse
    % to verify how reorganizeMiceData will separate the data
    %
    % Inputs:
    %   mice_all - Original mice data cell array (before reorganization)
    %   mouse_indices - Indices of mice to plot (optional, defaults to all)
    
    % Default to all mice if indices not provided
    if nargin < 2
        mouse_indices = 1:size(mice_all, 1);
    end
    
    % Column indices
    TIME_COL = 1;
    DOOR_COL = 7;
    DFF_COL = 11;
    
    % Colors for different sessions
    session_colors = {
        [0.7, 0.9, 1.0],  % sess0: Light blue (before discovery, door closed)
        [0.0, 0.3, 0.7],  % sess1: Dark blue (after discovery, door open)
        [0.3, 0.6, 0.9],  % sess2: Medium blue (second door closed)
        [0.0, 0.3, 0.7]   % sess3: Dark blue (second door open)
    };
    
    % Session labels for the legend
    session_labels = {
        'Session 0 (before discovery, door closed)',
        'Session 1 (after discovery, door open)',
        'Session 2 (second door closed)',
        'Session 3 (second door open)'
    };
    
    % Process each selected mouse
    for idx = 1:length(mouse_indices)
        i = mouse_indices(idx);
        
        % Skip if index is out of bounds
        if i > size(mice_all, 1)
            warning('Mouse index %d is out of bounds', i);
            continue;
        end
        
        % Extract data for this mouse
        mouse_id = mice_all{i, 1};
        condition = mice_all{i, 2};  % fasted/fed
        source = mice_all{i, 3};     % food/gel
        data = mice_all{i, 4};       % Actual data
        discovery = mice_all{i, 6};  % Discovery frame
        
        % Skip if discovery is invalid
        if isempty(discovery) || ~isnumeric(discovery) || discovery <= 0
            warning('Skipping mouse %s - invalid discovery frame', mouse_id);
            continue;
        end
        
        % Make sure discovery is within bounds
        discovery = min(discovery, size(data, 1));
        
        % Extract door status
        door_status = data(:, DOOR_COL);
        
        % Find transitions after discovery (where door status changes)
        post_discovery = door_status(discovery:end);
        transitions = find(diff(post_discovery) ~= 0) + discovery;
        
        % Create figure
        figure('Name', sprintf('Mouse %s (%s-%s) Session Preview', mouse_id, condition, source), ...
               'Position', [100, 100, 1200, 800]);
        
        % Create subplot for dF/F
        subplot(2, 1, 1);
        
        % Plot dF/F over time
        plot(data(:, TIME_COL), data(:, DFF_COL), 'k-', 'LineWidth', 1);
        hold on;
        
        % ---- Highlight sessions ----
        has_sessions = [false, false, false, false]; % Track which sessions exist
        
        % Session 0: Before discovery, door closed (0)
        has_session0 = false;
        session0_mask = false(size(data, 1), 1);
        session0_mask(1:discovery) = data(1:discovery, DOOR_COL) == 0;
        
        if any(session0_mask)
            has_session0 = true;
            has_sessions(1) = true;
            area(data(session0_mask, TIME_COL), data(session0_mask, DFF_COL), ...
                 'FaceColor', session_colors{1}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        
        % Calculate endpoints for session 1-3 if applicable
        if discovery < size(data, 1)
            % Session 1: After discovery, door open (1)
            has_session1 = false;
            if ~isempty(transitions)
                sess1_end = transitions(1);
            else
                sess1_end = size(data, 1);
            end
            
            session1_mask = false(size(data, 1), 1);
            session1_mask(discovery:sess1_end) = data(discovery:sess1_end, DOOR_COL) == 1;
            
            if any(session1_mask)
                has_session1 = true;
                has_sessions(2) = true;
                area(data(session1_mask, TIME_COL), data(session1_mask, DFF_COL), ...
                     'FaceColor', session_colors{2}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
            
            % Session 2: Second door closed period
            has_session2 = false;
            if length(transitions) >= 1 && transitions(1) < size(data, 1)
                if length(transitions) >= 2
                    sess2_end = transitions(2);
                else
                    sess2_end = size(data, 1);
                end
                
                session2_mask = false(size(data, 1), 1);
                session2_mask(transitions(1):sess2_end) = data(transitions(1):sess2_end, DOOR_COL) == 0;
                
                if any(session2_mask)
                    has_session2 = true;
                    has_sessions(3) = true;
                    area(data(session2_mask, TIME_COL), data(session2_mask, DFF_COL), ...
                         'FaceColor', session_colors{3}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                end
                
                % Session 3: Second door open period
                has_session3 = false;
                if length(transitions) >= 2 && transitions(2) < size(data, 1)
                    session3_mask = false(size(data, 1), 1);
                    session3_mask(transitions(2):end) = data(transitions(2):end, DOOR_COL) == 1;
                    
                    if any(session3_mask)
                        has_session3 = true;
                        has_sessions(4) = true;
                        area(data(session3_mask, TIME_COL), data(session3_mask, DFF_COL), ...
                             'FaceColor', session_colors{4}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                    end
                end
            end
        end
        
        % Add discovery line
        discovery_time = data(discovery, TIME_COL);
        y_limits = ylim;
        plot([discovery_time, discovery_time], y_limits, 'r--', 'LineWidth', 2);
        text(discovery_time, y_limits(2)*0.95, '  Discovery', 'Color', 'r', 'FontWeight', 'bold');
        
        % Add door transition lines
        for t = 1:length(transitions)
            transition_time = data(transitions(t), TIME_COL);
            
            % Check if the next index is within bounds
            if transitions(t)+1 <= size(data, 1)
                door_status_after = data(transitions(t)+1, DOOR_COL);
                
                if door_status_after == 1
                    label = '  Door Opens';
                    color = 'g';
                else
                    label = '  Door Closes';
                    color = 'b';
                end
                
                plot([transition_time, transition_time], y_limits, '--', 'Color', color, 'LineWidth', 2);
                text(transition_time, y_limits(2)*(0.9-0.05*t), label, 'Color', color, 'FontWeight', 'bold');
            end
        end
        
        % Create legend with only sessions that exist
        legend_items = {'dF/F'};
        legend_handles = [1]; % Placeholder for the dF/F line
        for s = 1:4
            if has_sessions(s)
                legend_items{end+1} = session_labels{s};
            end
        end
        legend_items{end+1} = 'Discovery';
        legend_items{end+1} = 'Door Transition';
        
        legend(legend_items, 'Location', 'best');
        title(sprintf('Mouse %s (%s - %s) dF/F with Session Coloring', ...
                     mouse_id, condition, source), 'FontSize', 14);
        ylabel('dF/F', 'FontSize', 12);
        grid on;
        
        % Create subplot for door status
        subplot(2, 1, 2);
        
        % Plot door status over time
        stairs(data(:, TIME_COL), data(:, DOOR_COL), 'k-', 'LineWidth', 2);
        hold on;
        
        % Add discovery line
        y_limits = [0, 1];
        plot([discovery_time, discovery_time], y_limits, 'r--', 'LineWidth', 2);
        
        % Add door transition lines
        for t = 1:length(transitions)
            transition_time = data(transitions(t), TIME_COL);
            plot([transition_time, transition_time], y_limits, 'g--', 'LineWidth', 2);
        end
        
        % Highlight session areas for reference
        if has_session0
            area(data(session0_mask, TIME_COL), 0.5*ones(sum(session0_mask), 1), ...
                 'FaceColor', session_colors{1}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        
        if discovery < size(data, 1)
            if has_session1
                area(data(session1_mask, TIME_COL), 0.5*ones(sum(session1_mask), 1), ...
                     'FaceColor', session_colors{2}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
            
            if length(transitions) >= 1 && has_session2
                area(data(session2_mask, TIME_COL), 0.5*ones(sum(session2_mask), 1), ...
                     'FaceColor', session_colors{3}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
            
            if length(transitions) >= 2 && has_session3
                area(data(session3_mask, TIME_COL), 0.5*ones(sum(session3_mask), 1), ...
                     'FaceColor', session_colors{4}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
        end
        
        % Format door status plot
        yticks([0, 1]);
        yticklabels({'Closed', 'Open'});
        ylim([-0.1, 1.1]);
        title('Door Status', 'FontSize', 14);
        xlabel('Time (s)', 'FontSize', 12);
        ylabel('Door State', 'FontSize', 12);
        grid on;

        % Link x-axes of both plots
        linkaxes(findobj(gcf, 'Type', 'axes'), 'x');

        % Add annotation with session counts
        annotation('textbox', [0.01, 0.01, 0.98, 0.03], ...
            'String', sprintf('Sessions: %d points in Sess0, %d in Sess1, %d in Sess2, %d in Sess3', ...
            sum(session0_mask), ...
            sum(session1_mask), ...
            sum(session2_mask), ...
            sum(session3_mask)), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end
end