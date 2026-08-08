function [zone_data, mice_all_updated] = analyzeZoneTimes(mice_all, options)
    % Extract options
    group_filter = options.group;
    use_limits = options.use_limits;
    dist_lim = options.dist_lim;
    time_lim = options.time_lim * 60;
    include_horizontal = options.horizontal;
    memory_tracking = options.memory_tracking;
    
    % Initialize zone time data structure
    zone_data = struct();
    sessions = {'sess0', 'sess1', 'sess2'};
    for s = 1:length(sessions)
        session = sessions{s};
        zone_data.times.(session) = [];
        zone_data.times_percent.(session) = [];
        zone_data.groups.(session) = [];
        zone_data.mouse_ids.(session) = {};
    end
    
    mice_all_updated = mice_all;
    memory_data = struct();  % Collect from ALL mice
    
    % Loop through each entry
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};
        group = mice_all{i, 2};
        food_arm = mice_all{i, 3};
        data = mice_all{i, 4};
        
        mouse_id = extractMouseID(session_info);
        
        % Determine session type
        if contains(session_info, 'sess0')
            session_type = 'sess0';
        elseif contains(session_info, 'sess1')
            session_type = 'sess1';
        elseif contains(session_info, 'sess2')
            session_type = 'sess2';
        else
            continue;
        end
        
        % Process data based on session
        if strcmp(session_type, 'sess2')
            data = data(data(:,12) == 0, :);
        elseif strcmp(session_type, 'sess1') && ~isempty(data)
            data = data(data(:,12) == 1, :);
            if ~isempty(data)
                data(:,1) = data(:,1) - data(1,1);
            end
        elseif strcmp(session_type, 'sess0') && ~isempty(data)
            data = data(data(:,12) == 0, :);
        end
        
        if isempty(data)
            continue;
        end
        
        % Apply time limit if needed
        if use_limits && ~isempty(data(:,1))
            time_array = data(:,1);
            target_time = time_lim;
            limit_idx = find(time_array <= target_time, 1, 'last');
            
            if isempty(limit_idx)
                if target_time > max(time_array)
                    limit_idx = length(time_array);
                else
                    limit_idx = 1;
                end
            end
            data_limited = data(1:limit_idx, :);
        else
            data_limited = data;
        end
        
        if isempty(data_limited)
            continue;
        end
        
        % Calculate zone times
        time_steps = mean(diff(data_limited(:, 1)));
        total_session_time = max(data_limited(:, 1)) - min(data_limited(:, 1));
        
        if use_limits
            [zone_times, ~] = calculateZoneTimesWithLimits(data_limited, food_arm, time_steps, dist_lim, include_horizontal);
        else
            [zone_times, ~] = calculateZoneTimes(data_limited, food_arm, time_steps, include_horizontal);
        end
        
        % Store memory data for ALL mice (not just filtered group)
        if memory_tracking && (strcmp(session_type, 'sess0') || strcmp(session_type, 'sess2'))
            zone_times_percent = (zone_times / total_session_time) * 100;
            memory_data.(mouse_id).(session_type) = zone_times_percent(1);
        end
        
        % Only include in plotting data if matches group filter
        if strcmp(group_filter, 'all') || strcmp(group, group_filter)
            zone_times_percent = (zone_times / total_session_time) * 100;
            
            if isempty(zone_data.times.(session_type))
                zone_data.times.(session_type) = zone_times;
                zone_data.times_percent.(session_type) = zone_times_percent;
            else
                zone_data.times.(session_type) = [zone_data.times.(session_type); zone_times];
                zone_data.times_percent.(session_type) = [zone_data.times_percent.(session_type); zone_times_percent];
            end
            
            zone_data.groups.(session_type) = [zone_data.groups.(session_type); {group}];
            zone_data.mouse_ids.(session_type){end+1} = mouse_id;
        end
    end
    
    % Update memory strength for ALL mice
    if memory_tracking
        mice_all_updated = updateMemoryStrength(mice_all, memory_data);
    end
end
