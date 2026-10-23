function [zone_times, data_with_zones] = calculateZoneTimesWithLimits(data, food_arm, time_steps, dist_lim, include_horizontal)
    % Initialize times
    food_zone_time = 0;
    non_food_arm1_time = 0;
    non_food_arm2_time = 0;
    horizontal_time = 0;
    
    % Create a copy of data for zone marking
    data_with_zones = data;
    if size(data_with_zones, 2) < 17
        data_with_zones(:, 17) = zeros(size(data_with_zones, 1), 1);
    end
    
    % Process each data point
    for j = 1:size(data, 1)
        current_zone = data(j, 8);
        
        switch current_zone
            case 1  % Left arm
                if strcmp(food_arm, 'left') && data(j, 9) < dist_lim
                    food_zone_time = food_zone_time + time_steps;
                    data_with_zones(j, 17) = 1;  % Mark as food zone
                elseif strcmp(food_arm, 'center') && (data(j, 11) < dist_lim || data(j, 10) < dist_lim)
                    non_food_arm2_time = non_food_arm2_time + time_steps;
                    data_with_zones(j, 17) = 3;  % Mark as non-food arm 2
                elseif strcmp(food_arm, 'right') && (data(j, 10) < dist_lim || data(j, 11) < dist_lim)
                    non_food_arm1_time = non_food_arm1_time + time_steps;
                    data_with_zones(j, 17) = 2;  % Mark as non-food arm 1
                end
                
            case 2  % Center arm
                if strcmp(food_arm, 'center') && data(j, 9) < dist_lim
                    food_zone_time = food_zone_time + time_steps;
                    data_with_zones(j, 17) = 1;  % Mark as food zone
                elseif strcmp(food_arm, 'left') && (data(j, 10) < dist_lim || data(j, 11) < dist_lim)
                    non_food_arm2_time = non_food_arm2_time + time_steps;
                    data_with_zones(j, 17) = 3;  % Mark as non-food arm 2
                elseif strcmp(food_arm, 'right') && (data(j, 10) < dist_lim || data(j, 11) < dist_lim)
                    non_food_arm2_time = non_food_arm2_time + time_steps;
                    data_with_zones(j, 17) = 3;  % Mark as non-food arm 2
                end
                
            case 3  % Right arm
                if strcmp(food_arm, 'right') && data(j, 9) < dist_lim
                    food_zone_time = food_zone_time + time_steps;
                    data_with_zones(j, 17) = 1;  % Mark as food zone
                elseif strcmp(food_arm, 'center') && (data(j, 11) < dist_lim || data(j, 10) < dist_lim)
                    non_food_arm1_time = non_food_arm1_time + time_steps;
                    data_with_zones(j, 17) = 2;  % Mark as non-food arm 1
                elseif strcmp(food_arm, 'left') && (data(j, 11) < dist_lim || data(j, 10) < dist_lim)
                    non_food_arm1_time = non_food_arm1_time + time_steps;
                    data_with_zones(j, 17) = 2;  % Mark as non-food arm 1
                end
                
            case 4  % Horizontal arm
                horizontal_time = horizontal_time + time_steps;
                data_with_zones(j, 17) = 4;  % Mark as horizontal
                
            case 5  % Food zone
                food_zone_time = food_zone_time + time_steps;
                data_with_zones(j, 17) = 1;  % Mark as food zone
        end
    end
    
    % Return zone times
    if include_horizontal
        zone_times = [food_zone_time, non_food_arm1_time, non_food_arm2_time, horizontal_time];
    else
        zone_times = [food_zone_time, non_food_arm1_time, non_food_arm2_time];
    end
end
