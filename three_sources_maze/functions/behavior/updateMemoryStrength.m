%% Fixed updateMemoryStrength function - processes ALL mice regardless of group
function mice_all_updated = updateMemoryStrength(mice_all, memory_data)
    % Make a copy of mice_all
    mice_all_updated = mice_all;
    
    % Get all mouse IDs with memory data
    mouse_ids = fieldnames(memory_data);
    
    % First, ensure mice_all_updated has at least 6 columns
    if size(mice_all_updated, 2) < 6
        temp = cell(size(mice_all_updated, 1), 6);
        for i = 1:size(mice_all_updated, 1)
            for j = 1:size(mice_all_updated, 2)
                temp{i, j} = mice_all_updated{i, j};
            end
            temp{i, 6} = '';  % Initialize column 6
        end
        mice_all_updated = temp;
    end
    
    % Process each unique mouse
    processed_mice = {};
    
    for m = 1:length(mouse_ids)
        mouse_id = mouse_ids{m};
        
        % Skip if already processed
        if ismember(mouse_id, processed_mice)
            continue;
        end
        
        % Check if we have both sess0 and sess2 data
        if isfield(memory_data.(mouse_id), 'sess0') && isfield(memory_data.(mouse_id), 'sess2')
            % Get the food zone percentages
            food_zone_percent1 = memory_data.(mouse_id).sess0;
            food_zone_percent2 = memory_data.(mouse_id).sess2;
            
            % Determine memory strength
            if food_zone_percent2 > food_zone_percent1
                memory_strength = 'strong memory';
            else
                memory_strength = 'weak memory';
            end
            
            % Update ALL entries for this mouse (regardless of group)
            for i = 1:size(mice_all, 1)
                session_info = mice_all{i, 1};
                current_mouse_id = extractMouseID(session_info);
                
                if strcmp(current_mouse_id, mouse_id)
                    mice_all_updated{i, 6} = memory_strength;
                end
            end
            
            processed_mice{end+1} = mouse_id;
        end
    end
    
    fprintf('Memory tracking completed for all mice\n');
end
