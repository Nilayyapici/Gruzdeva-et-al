function visualizeOriginalMiceData(mice, condition)
% VISUALIZEORIGINALMICEDATA - Plot all mice before rotation to verify data
%
% Shows pre-session trajectories, centers, food positions, and center-to-food vectors
% for all mice in the specified condition
%
% Usage:
%   visualizeOriginalMiceData(mice, 'saline')
%   visualizeOriginalMiceData(mice, 'CNO')

%% Select mice based on condition
mouse_names = mice(:,1);
saline_idx = contains(mouse_names, 'saline');
cno_idx = contains(mouse_names, 'CNO');

has_memory_classification = size(mice, 2) >= 7 && ~all(cellfun(@isempty, mice(:,7)));
if has_memory_classification
    memory_classifications = mice(:,7);
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

fprintf('\n=== VISUALIZING ORIGINAL DATA ===\n');
fprintf('Condition: %s (%d mice)\n', condition_name, n_mice);

%% Create figure
figure('Name', sprintf('Original Data - %s', condition_name), ...
       'Position', [100, 100, 1000, 800], 'Color', 'white');

hold on;

% Generate colors for each mouse
colors = lines(n_mice);

%% Plot each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice{mouse_idx, 1};
    food_pos = mice{mouse_idx, 2};
    pre_data = mice{mouse_idx, 3};
    center_pos = mice{mouse_idx, 5};
    
    fprintf('\nMouse %d: %s\n', i, mouse_name);
    fprintf('  Food: [%.1f, %.1f]\n', food_pos(1), food_pos(2));
    fprintf('  Center: [%.1f, %.1f]\n', center_pos(1), center_pos(2));
    
    % Calculate angle
    dx = food_pos(1) - center_pos(1);
    dy = food_pos(2) - center_pos(2);
    angle = rad2deg(atan2(dy, dx));
    distance = sqrt(dx^2 + dy^2);
    fprintf('  Angle: %.1f degrees, Distance: %.1f\n', angle, distance);
    
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

%% Format plot
axis equal;
grid on;
xlabel('X Position', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Original Data Before Rotation - %s (n=%d)', condition_name, n_mice), ...
      'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

hold off;

fprintf('\n=== VISUALIZATION COMPLETE ===\n');

end