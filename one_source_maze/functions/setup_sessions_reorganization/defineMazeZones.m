function mice_all = defineMazeZones(mice_all, corner_coords)
    % DEFINEMAKEZONES Define zones for all mice using polyfit boundaries
    % 
    % Inputs:
    %   mice_all - cell array with mouse data
    %   corner_coords - struct with corner coordinates (x1-x12, y1-y12, x1_door, y1_door, x2_door, y2_door)
    %
    % Output:
    %   mice_all - updated with zone column added to each mouse's data
    
    % Extract corner coordinates
    x1 = corner_coords.x1; y1 = corner_coords.y1;
    x2 = corner_coords.x2; y2 = corner_coords.y2;
    x3 = corner_coords.x3; y3 = corner_coords.y3;
    x4 = corner_coords.x4; y4 = corner_coords.y4;
    x5 = corner_coords.x5; y5 = corner_coords.y5;
    x6 = corner_coords.x6; y6 = corner_coords.y6;
    x7 = corner_coords.x7; y7 = corner_coords.y7;
    x8 = corner_coords.x8; y8 = corner_coords.y8;
    x9 = corner_coords.x9; y9 = corner_coords.y9;
    x10 = corner_coords.x10; y10 = corner_coords.y10;
    x11 = corner_coords.x11; y11 = corner_coords.y11;
    x12 = corner_coords.x12; y12 = corner_coords.y12;
    
    % Calculate boundary line coefficients using polyfit
    coefficients1 = polyfit([x1, x12], [y1, y12], 1);
    a1 = coefficients1(1); % left side out
    b1 = coefficients1(2);
    
    coefficients2 = polyfit([x2, x3], [y2, y3], 1);
    a2 = coefficients2(1); % left inner
    b2 = coefficients2(2);
    
    coefficients3 = polyfit([x4, x5], [y4, y5], 1);
    a3 = coefficients3(1); % left center
    b3 = coefficients3(2);
    
    coefficients4 = polyfit([x6, x7], [y6, y7], 1);
    a4 = coefficients4(1); % right center
    b4 = coefficients4(2);
    
    coefficients5 = polyfit([x8, x9], [y8, y9], 1);
    a5 = coefficients5(1); % right inner
    b5 = coefficients5(2);
    
    coefficients6 = polyfit([x10, x11], [y10, y11], 1);
    a6 = coefficients6(1); % right outside
    b6 = coefficients6(2);
    
    coefficients7 = polyfit([x3, x8], [y3, y8], 1);
    a7 = coefficients7(1); % horizontal inner
    b7 = coefficients7(2);
    
    coefficients8 = polyfit([x11, x12], [y11, y12], 1);
    a8 = coefficients8(1); % horizontal outside
    b8 = coefficients8(2);

    % Process each mouse
    for mouse_idx = 1:size(mice_all, 1)
        data = mice_all{mouse_idx, 4};
        mouse_id = mice_all{mouse_idx, 1};
        
        % Get x,y coordinates (columns 2,3) and distance to food (column 5)
        x_coords = data(:, 2);
        y_coords = data(:, 3);
        dist_to_food = data(:, 5);
        n_points = length(x_coords);
        
        % Find where food is located (minimum distance to food)
        [~, food_idx] = min(dist_to_food);
        food_x = x_coords(food_idx);
        food_y = y_coords(food_idx);
        
        % Determine which arm contains the food using temporary zone assignment
        if (food_y - a1*food_x) <= b1 && (food_y - a2*food_x) >= b2 && (food_y - a7*food_x) <= b7
            food_arm = 'left'; % Food is in left arm
        elseif (food_y - a5*food_x) >= b5 && (food_y - a6*food_x) < b6 && (food_y - a7*food_x) <= b7
            food_arm = 'right'; % Food is in right arm
        else
            food_arm = 'unknown';
        end
        
        % Initialize zones column
        zones = NaN(n_points, 1);
        
        % Apply zone assignment logic based on food location
        for i = 1:n_points
            x_pos = x_coords(i);
            y_pos = y_coords(i);
            
            % Determine which physical arm this point is in
            if (y_pos - a1*x_pos) <= b1 && (y_pos - a2*x_pos) >= b2 && (y_pos - a7*x_pos) <= b7
                physical_arm = 'left';
            elseif (y_pos - a3*x_pos) <= b3 && (y_pos - a4*x_pos) >= b4 && (y_pos - a7*x_pos) <= b7
                physical_arm = 'center';
            elseif (y_pos - a5*x_pos) >= b5 && (y_pos - a6*x_pos) < b6 && (y_pos - a7*x_pos) <= b7
                physical_arm = 'right';
            elseif (y_pos - a7*x_pos) > b7 && (y_pos - a8*x_pos) <= b8
                zones(i) = 4; % horizontal corridor
                continue;
            else
                zones(i) = NaN;
                continue;
            end
            
            % Assign zone numbers based on food location
            if strcmp(physical_arm, 'center')
                zones(i) = 2; % center arm always stays center
            elseif strcmp(physical_arm, food_arm)
                zones(i) = 3; % food arm
            else
                zones(i) = 1; % opposite arm (the other arm that doesn't have food)
            end
        end
        
        % Add zones as new column to data
        data_with_zones = [data, zones];
        mice_all{mouse_idx, 4} = data_with_zones;
        
        % Create scatter plot for this mouse
        figure;
        hold on;
        scatter(x_coords, y_coords, 5, zones, 'filled');
        title(['zones - ' mouse_id], 'Interpreter', 'none');
        hold off;
    end
end