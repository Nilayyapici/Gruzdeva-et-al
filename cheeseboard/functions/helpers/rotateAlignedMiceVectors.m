function mice_rotated = rotateAlignedMiceVectors(mice_translated, condition, target_angle)
% ROTATEALIGNEDMICEVECTORS - Rotate mice data to align center-to-food vectors
%
% Assumes centers are already aligned. Rotates each mouse's data around the
% common center to align all center-to-food vectors to the same angle.
% Vector lengths are preserved.
%
% Inputs:
%   mice_translated - Cell array with translated mouse data (centers aligned)
%   condition - String specifying condition ('saline', 'CNO', 'strong', 'weak')
%   target_angle - Target angle in degrees for all vectors (default: 0 = pointing right)
%
% Output:
%   mice_rotated - Cell array with rotated coordinates
%
% Usage:
%   mice_rotated = rotateAlignedMiceVectors(mice_translated, 'saline');
%   mice_rotated = rotateAlignedMiceVectors(mice_translated, 'saline', 0);

if nargin < 3
    target_angle = 0; % Default: point to the right
end

%% Select mice based on condition
mouse_names = mice_translated(:,1);
saline_idx = contains(mouse_names, 'saline');
cno_idx = contains(mouse_names, 'CNO');

has_memory_classification = size(mice_translated, 2) >= 7 && ~all(cellfun(@isempty, mice_translated(:,7)));
if has_memory_classification
    memory_classifications = mice_translated(:,7);
    strong_memory_idx = strcmp(memory_classifications, 'strong memory');
    weak_memory_idx = strcmp(memory_classifications, 'weak memory');
end

switch lower(condition)
    case 'saline', selected_idx = saline_idx; condition_name = 'Saline';
    case 'cno', selected_idx = cno_idx; condition_name = 'CNO';
    case 'strong', selected_idx = strong_memory_idx; condition_name = 'Strong Memory';
    case 'weak', selected_idx = weak_memory_idx; condition_name = 'Weak Memory';
    otherwise, error('Invalid condition');
end

selected_mice = find(selected_idx);
n_mice = length(selected_mice);

fprintf('\n=== ROTATING VECTORS FOR %d MICE (%s condition) ===\n', n_mice, condition_name);
fprintf('Target angle: %.1f degrees\n', target_angle);

%% Initialize rotated data
mice_rotated = mice_translated;

%% Rotate each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice_translated{mouse_idx, 1};
    food_pos = mice_translated{mouse_idx, 2};
    pre_data = mice_translated{mouse_idx, 3};
    test_data = mice_translated{mouse_idx, 4};
    center_pos = mice_translated{mouse_idx, 5};
    
    % Ensure column vectors
    food_pos = food_pos(:);
    center_pos = center_pos(:);
    
    fprintf('\nMouse %d/%d: %s\n', i, n_mice, mouse_name);
    fprintf('  Center: [%.1f, %.1f]\n', center_pos(1), center_pos(2));
    fprintf('  Original food: [%.1f, %.1f]\n', food_pos(1), food_pos(2));
    
    % Calculate current angle of center-to-food vector
    food_vector = food_pos - center_pos;
    current_angle_rad = atan2(food_vector(2), food_vector(1));
    current_angle_deg = rad2deg(current_angle_rad);
    vector_length = sqrt(food_vector(1)^2 + food_vector(2)^2);
    
    fprintf('  Current angle: %.1f degrees, Length: %.1f\n', current_angle_deg, vector_length);
    
    % Calculate rotation needed
    target_angle_rad = deg2rad(target_angle);
    rotation_angle_rad = target_angle_rad - current_angle_rad;
    rotation_angle_deg = rad2deg(rotation_angle_rad);
    
    fprintf('  Rotation to apply: %.1f degrees\n', rotation_angle_deg);
    
    % Create rotation matrix
    cos_theta = cos(rotation_angle_rad);
    sin_theta = sin(rotation_angle_rad);
    R = [cos_theta, -sin_theta;
         sin_theta, cos_theta];
    
    % Rotate food position around center
    food_rotated = R * food_vector + center_pos;
    
    % Verify
    new_angle = rad2deg(atan2(food_rotated(2) - center_pos(2), food_rotated(1) - center_pos(1)));
    new_length = sqrt((food_rotated(1) - center_pos(1))^2 + (food_rotated(2) - center_pos(2))^2);
    fprintf('  New food: [%.1f, %.1f]\n', food_rotated(1), food_rotated(2));
    fprintf('  New angle: %.1f degrees, Length: %.1f\n', new_angle, new_length);
    
    % Rotate pre data coordinates
    if ~isempty(pre_data)
        n_points = size(pre_data, 1);
        pre_data_rotated = pre_data;
        
        for j = 1:n_points
            % Get point coordinates
            x = pre_data(j, 2);
            y = pre_data(j, 3);
            
            % Vector from center to point
            point_vector = [x - center_pos(1); y - center_pos(2)];
            
            % Rotate vector
            rotated_vector = R * point_vector;
            
            % New coordinates
            pre_data_rotated(j, 2) = rotated_vector(1) + center_pos(1);
            pre_data_rotated(j, 3) = rotated_vector(2) + center_pos(2);
        end
    else
        pre_data_rotated = pre_data;
    end
    
    % Rotate test data coordinates
    if ~isempty(test_data)
        n_points = size(test_data, 1);
        test_data_rotated = test_data;
        
        for j = 1:n_points
            % Get point coordinates
            x = test_data(j, 2);
            y = test_data(j, 3);
            
            % Vector from center to point
            point_vector = [x - center_pos(1); y - center_pos(2)];
            
            % Rotate vector
            rotated_vector = R * point_vector;
            
            % New coordinates
            test_data_rotated(j, 2) = rotated_vector(1) + center_pos(1);
            test_data_rotated(j, 3) = rotated_vector(2) + center_pos(2);
        end
    else
        test_data_rotated = test_data;
    end
    
    % Store rotated data
    mice_rotated{mouse_idx, 2} = food_rotated;
    mice_rotated{mouse_idx, 3} = pre_data_rotated;
    mice_rotated{mouse_idx, 4} = test_data_rotated;
    % Center stays the same
end

fprintf('\n=== ROTATION COMPLETE ===\n');

%% Visualize results
figure('Name', sprintf('Rotated Data - %s', condition_name), ...
       'Position', [100, 100, 1000, 800], 'Color', 'white');

hold on;

% Generate colors
colors = lines(n_mice);

% Collect food positions for statistics
all_food_angles = zeros(n_mice, 1);
all_food_distances = zeros(n_mice, 1);

% Plot each mouse after rotation
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice_rotated{mouse_idx, 1};
    food_pos = mice_rotated{mouse_idx, 2};
    pre_data = mice_rotated{mouse_idx, 3};
    center_pos = mice_rotated{mouse_idx, 5};
    
    % Calculate angle and distance
    dx = food_pos(1) - center_pos(1);
    dy = food_pos(2) - center_pos(2);
    all_food_angles(i) = rad2deg(atan2(dy, dx));
    all_food_distances(i) = sqrt(dx^2 + dy^2);
    
    % Plot trajectory (subsampled)
    if size(pre_data, 1) > 100
        sample_idx = 1:50:size(pre_data, 1);
        x_traj = pre_data(sample_idx, 2);
        y_traj = pre_data(sample_idx, 3);
    else
        x_traj = pre_data(:, 2);
        y_traj = pre_data(:, 3);
    end
    
    plot(x_traj, y_traj, '-', 'Color', [colors(i,:), 0.3], 'LineWidth', 0.5);
    
    % Plot center
    scatter(center_pos(1), center_pos(2), 150, colors(i,:), 'x', 'LineWidth', 3);
    
    % Plot food
    scatter(food_pos(1), food_pos(2), 150, colors(i,:), 'o', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Plot vector from center to food
    quiver(center_pos(1), center_pos(2), dx, dy, 0, ...
           'Color', colors(i,:), 'LineWidth', 2, 'MaxHeadSize', 1);
    
    % Add label
    text(food_pos(1) + 10, food_pos(2) + 10, mouse_name, ...
         'FontSize', 9, 'Color', colors(i,:), 'FontWeight', 'bold');
end

% Draw target angle reference line
common_center = mice_rotated{selected_mice(1), 5};
target_line_length = mean(all_food_distances) * 1.2;
target_x = [common_center(1), common_center(1) + target_line_length * cos(deg2rad(target_angle))];
target_y = [common_center(2), common_center(2) + target_line_length * sin(deg2rad(target_angle))];
plot(target_x, target_y, 'r--', 'LineWidth', 3, 'DisplayName', sprintf('Target: %.0f°', target_angle));

axis equal;
grid on;
xlabel('X Position', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position', 'FontSize', 12, 'FontWeight', 'bold');

% Calculate statistics
mean_angle = mean(all_food_angles);
std_angle = std(all_food_angles);
mean_distance = mean(all_food_distances);
std_distance = std(all_food_distances);

title(sprintf('Rotated Data - All Vectors Aligned - %s (n=%d)\nAngle: %.1f° (±%.2f°), Distance: %.1f (±%.1f)', ...
      condition_name, n_mice, mean_angle, std_angle, mean_distance, std_distance), ...
      'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

hold off;

fprintf('\n=== ALIGNMENT STATISTICS ===\n');
fprintf('Mean angle: %.2f degrees (target: %.1f)\n', mean_angle, target_angle);
fprintf('Angle std dev: %.2f degrees\n', std_angle);
fprintf('Mean distance: %.2f\n', mean_distance);
fprintf('Distance std dev: %.2f\n', std_distance);

end