function plotAllMiceInMaze(mice_all, corner_coords)
    % PLOTALLMICEINMAZE Plot all mouse trajectories overlaid on M-maze outline
    % 
    % Inputs:
    %   mice_all - cell array with mouse data
    %   corner_coords - struct with corner coordinates (x1-x12, y1-y12, etc.)
    
    % Create figure
    figure('Position', [100, 100, 1000, 800]);
    hold on;
    
    % Plot maze outline first
    plotMazeOutline(corner_coords);
    
    % Plot all mouse trajectories
    for i = 1:size(mice_all, 1)
        try
            data = mice_all{i, 4};
            
            x = data(:, 2);
            y = data(:, 3);
            
            % Sample data to reduce clutter (every 10th point)
            indices = 1:10:length(x);
            x_sample = x(indices);
            y_sample = y(indices);
            
            % Plot trajectory as scatter
            scatter(x_sample, y_sample, 3, [0.3, 0.3, 0.3], 'filled', 'MarkerFaceAlpha', 0.3);
            
        catch ME
            warning('Error plotting mouse %d (%s): %s', i, mice_all{i, 1}, ME.message);
            continue;
        end
    end
    
    % Formatting
    axis equal;
    grid on;
    xlabel('X coordinate (pixels)', 'FontSize', 12);
    ylabel('Y coordinate (pixels)', 'FontSize', 12);
    title(sprintf('All Mouse Trajectories in M-maze (N=%d mice)', size(mice_all, 1)), 'FontSize', 14, 'FontWeight', 'bold');
    
    hold off;
    
    % Print summary
    fprintf('Plotted %d mice trajectories\n', size(mice_all, 1));
end

function plotMazeOutline(corners)
    % Plot the M-maze outline using corner coordinates
    
    % Extract coordinates
    x1 = corners.x1; y1 = corners.y1;
    x2 = corners.x2; y2 = corners.y2;
    x3 = corners.x3; y3 = corners.y3;
    x4 = corners.x4; y4 = corners.y4;
    x5 = corners.x5; y5 = corners.y5;
    x6 = corners.x6; y6 = corners.y6;
    x7 = corners.x7; y7 = corners.y7;
    x8 = corners.x8; y8 = corners.y8;
    x9 = corners.x9; y9 = corners.y9;
    x10 = corners.x10; y10 = corners.y10;
    x11 = corners.x11; y11 = corners.y11;
    x12 = corners.x12; y12 = corners.y12;
    
    % Plot maze walls
    wall_color = 'k';
    wall_width = 2;
    
    plot([x1, x2], [y1, y2], wall_color, 'LineWidth', wall_width);
    plot([x2, x3], [y2, y3], wall_color, 'LineWidth', wall_width);
    plot([x3, x4], [y3, y4], wall_color, 'LineWidth', wall_width);
    plot([x4, x5], [y4, y5], wall_color, 'LineWidth', wall_width);
    plot([x5, x6], [y5, y6], wall_color, 'LineWidth', wall_width);
    plot([x6, x7], [y6, y7], wall_color, 'LineWidth', wall_width);
    plot([x7, x8], [y7, y8], wall_color, 'LineWidth', wall_width);
    plot([x8, x9], [y8, y9], wall_color, 'LineWidth', wall_width);
    plot([x9, x10], [y9, y10], wall_color, 'LineWidth', wall_width);
    plot([x10, x11], [y10, y11], wall_color, 'LineWidth', wall_width);
    plot([x11, x12], [y11, y12], wall_color, 'LineWidth', wall_width);
    plot([x12, x1], [y12, y1], wall_color, 'LineWidth', wall_width);
    
    % Plot door if coordinates provided
    if isfield(corners, 'x1_door') && isfield(corners, 'x2_door')
        plot([corners.x1_door, corners.x2_door], [corners.y1_door, corners.y2_door], wall_color, 'LineWidth', wall_width);
    end
    
    % Plot food location if provided
    if isfield(corners, 'x_food') && isfield(corners, 'y_food')
        if isfield(corners, 'r_food')
            % Draw food circle
            theta = linspace(0, 2*pi, 100);
            x_circle = corners.x_food + corners.r_food * cos(theta);
            y_circle = corners.y_food + corners.r_food * sin(theta);
            plot(x_circle, y_circle, 'r-', 'LineWidth', 2);
            fill(x_circle, y_circle, 'r', 'FaceAlpha', 0.3);
        else
            plot(corners.x_food, corners.y_food, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        end
        text(corners.x_food + 20, corners.y_food, 'Food', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    end
    
    % Mark corners with numbers
    corner_x = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12];
    corner_y = [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12];
    
    for i = 1:12
        plot(corner_x(i), corner_y(i), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'yellow');
        text(corner_x(i) + 10, corner_y(i) + 10, sprintf('%d', i), 'FontSize', 8, 'FontWeight', 'bold');
    end
end