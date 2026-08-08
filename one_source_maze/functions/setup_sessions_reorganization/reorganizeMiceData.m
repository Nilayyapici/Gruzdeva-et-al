function mice_all_reorganized = reorganizeMiceData(mice_all)
    % REORGANIZEMIDATA Reorganizes mice data based on door conditions
    %   Splits each mouse dataset into sessions based on door status
    %   Session 0: Door closed (0) before food discovery
    %   Session 1: Door open (1) after food discovery until second door closing
    %   Session 2: Door closed (0) for the second time
    %   Session 3: Door open (1) for the second time
    
    mice_all_reorganized = {};
    door_col = 7;
    
    % Loop through each entry in the original data
    for i = 1:size(mice_all, 1)
        % Extract relevant information
        mouse_id = mice_all{i, 1};
        side = mice_all{i, 2};     % fasted/fed
        side2 = mice_all{i, 3};    % food/gel
        data = mice_all{i, 4};     % Actual data from the file
        side3 = mice_all{i, 5};    % Signal quality indicator
        discovery = mice_all{i, 6}; % Food discovery frame
        
        % Check for valid discovery frame
        if isempty(discovery) || ~isnumeric(discovery) || discovery <= 0
            warning('Skipping mouse %s - invalid discovery frame', mouse_id);
            continue;
        end
        
        % Get door status changes after discovery
        door_status = data(:, door_col);
        
        % ---- Session 0: Before discovery, door closed (0) ----
        data_sess0 = data(1:discovery, :);
        data_sess0_closed = data_sess0(data_sess0(:, door_col) == 0, :);
        
        % Add session 0 data if not empty
        if ~isempty(data_sess0_closed)
            row_idx = size(mice_all_reorganized, 1) + 1;
            mice_all_reorganized{row_idx, 1} = [mouse_id, '_sess0']; % Add session suffix
            mice_all_reorganized{row_idx, 2} = side;
            mice_all_reorganized{row_idx, 3} = side2;
            mice_all_reorganized{row_idx, 4} = data_sess0_closed;
            mice_all_reorganized{row_idx, 5} = side3;
            mice_all_reorganized{row_idx, 6} = discovery;
        end
        
        % ---- Find transitions after discovery ----
        post_discovery = door_status(discovery:end);
        transitions = find(diff(post_discovery) ~= 0) + discovery;
        
        % ---- Session 1: After discovery, door open (1) ----
        if discovery < size(data, 1)
            % Find the end of session 1 (when door closes again)
            if ~isempty(transitions)
                sess1_end = transitions(1);
            else
                sess1_end = size(data, 1);
            end
            
            % Extract session 1 data
            data_sess1 = data(discovery:sess1_end, :);
            data_sess1_open = data_sess1(data_sess1(:, door_col) == 1, :);
            
            % Add session 1 data if not empty
            if ~isempty(data_sess1_open)
                row_idx = size(mice_all_reorganized, 1) + 1;
                mice_all_reorganized{row_idx, 1} = [mouse_id, '_sess1'];
                mice_all_reorganized{row_idx, 2} = side;
                mice_all_reorganized{row_idx, 3} = side2;
                mice_all_reorganized{row_idx, 4} = data_sess1_open;
                mice_all_reorganized{row_idx, 5} = side3;
                mice_all_reorganized{row_idx, 6} = discovery;
            end
            
            % ---- Session 2: Second door closed period ----
            if length(transitions) >= 1 && transitions(1) < size(data, 1)
                % Find the end of session 2 (when door opens again)
                if length(transitions) >= 2
                    sess2_end = transitions(2);
                else
                    sess2_end = size(data, 1);
                end
                
                % Extract session 2 data
                data_sess2 = data(transitions(1):sess2_end, :);
                data_sess2_closed = data_sess2(data_sess2(:, door_col) == 0, :);
                
                % Add session 2 data if not empty
                if ~isempty(data_sess2_closed)
                    row_idx = size(mice_all_reorganized, 1) + 1;
                    mice_all_reorganized{row_idx, 1} = [mouse_id, '_sess2'];
                    mice_all_reorganized{row_idx, 2} = side;
                    mice_all_reorganized{row_idx, 3} = side2;
                    mice_all_reorganized{row_idx, 4} = data_sess2_closed;
                    mice_all_reorganized{row_idx, 5} = side3;
                    mice_all_reorganized{row_idx, 6} = discovery;
                end
                
                % ---- Session 3: Second door open period ----
                if length(transitions) >= 2 && transitions(2) < size(data, 1)
                    % Extract session 3 data
                    data_sess3 = data(transitions(2):end, :);
                    data_sess3_open = data_sess3(data_sess3(:, door_col) == 1, :);
                    
                    % Add session 3 data if not empty
                    if ~isempty(data_sess3_open)
                        row_idx = size(mice_all_reorganized, 1) + 1;
                        mice_all_reorganized{row_idx, 1} = [mouse_id, '_sess3'];
                        mice_all_reorganized{row_idx, 2} = side;
                        mice_all_reorganized{row_idx, 3} = side2;
                        mice_all_reorganized{row_idx, 4} = data_sess3_open;
                        mice_all_reorganized{row_idx, 5} = side3;
                        mice_all_reorganized{row_idx, 6} = discovery;
                    end
                end
            end
        end
    end
    
    % Display summary of reorganization
    fprintf('Reorganization complete.\n');
    fprintf('Original data had %d entries.\n', size(mice_all, 1));
    fprintf('Reorganized data has %d entries.\n', size(mice_all_reorganized, 1));
    
    % Optional: Print session distribution
    sess_counts = zeros(1, 4);
    for i = 1:size(mice_all_reorganized, 1)
        mouse_id = mice_all_reorganized{i, 1};
        if contains(mouse_id, '_sess0')
            sess_counts(1) = sess_counts(1) + 1;
        elseif contains(mouse_id, '_sess1')
            sess_counts(2) = sess_counts(2) + 1;
        elseif contains(mouse_id, '_sess2')
            sess_counts(3) = sess_counts(3) + 1;
        elseif contains(mouse_id, '_sess3')
            sess_counts(4) = sess_counts(4) + 1;
        end
    end
    
    fprintf('Session counts:\n');
    fprintf('  Session 0 (before discovery, door closed): %d\n', sess_counts(1));
    fprintf('  Session 1 (after discovery, door open): %d\n', sess_counts(2));
    fprintf('  Session 2 (second door closed period): %d\n', sess_counts(3));
    fprintf('  Session 3 (second door open period): %d\n', sess_counts(4));
end