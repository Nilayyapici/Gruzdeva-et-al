function mice_translated = translateMiceCenters(mice, condition, target_center)
% TRANSLATEMICECENTERS - Translate all mice data to align centers at target position
%
% Simple translation (no rotation) to move all centers to the same location
%
% Inputs:
%   mice - Cell array with mouse data
%   condition - String specifying condition ('saline', 'CNO', 'strong', 'weak')
%   target_center - Target center position [x, y] (default: [350, 250])
%
% Output:
%   mice_translated - Cell array with translated coordinates
%
% Usage:
%   mice_translated = translateMiceCenters(mice, 'saline');
%   mice_translated = translateMiceCenters(mice, 'CNO', [350, 250]);

if nargin < 3
    target_center = [350; 250]; % Default target
end

% Ensure column vector
target_center = target_center(:);

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

fprintf('\n=== TRANSLATING DATA FOR %d MICE (%s condition) ===\n', n_mice, condition_name);
fprintf('Target center: [%.1f, %.1f]\n', target_center(1), target_center(2));

%% Initialize translated data
mice_translated = mice;

%% Translate each mouse
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice{mouse_idx, 1};
    food_pos = mice{mouse_idx, 2};
    pre_data = mice{mouse_idx, 3};
    test_data = mice{mouse_idx, 4};
    center_pos = mice{mouse_idx, 5};
    
    % Ensure column vectors
    food_pos = food_pos(:);
    center_pos = center_pos(:);
    
    fprintf('\nMouse %d/%d: %s\n', i, n_mice, mouse_name);
    fprintf('  Original center: [%.1f, %.1f]\n', center_pos(1), center_pos(2));
    
    % Calculate translation offset
    offset = target_center - center_pos;
    
    fprintf('  Translation offset: [%.1f, %.1f]\n', offset(1), offset(2));
    
    % Translate food position
    food_translated = food_pos + offset;
    
    % Translate pre data coordinates (columns 2 and 3 are x, y)
    if ~isempty(pre_data)
        pre_data_translated = pre_data;
        pre_data_translated(:, 2) = pre_data(:, 2) + offset(1);
        pre_data_translated(:, 3) = pre_data(:, 3) + offset(2);
    else
        pre_data_translated = pre_data;
    end
    
    % Translate test data coordinates
    if ~isempty(test_data)
        test_data_translated = test_data;
        test_data_translated(:, 2) = test_data(:, 2) + offset(1);
        test_data_translated(:, 3) = test_data(:, 3) + offset(2);
    else
        test_data_translated = test_data;
    end
    
    % The center itself becomes the target
    center_translated = target_center;
    
    fprintf('  New center: [%.1f, %.1f]\n', center_translated(1), center_translated(2));
    fprintf('  New food: [%.1f, %.1f]\n', food_translated(1), food_translated(2));
    
    % Store translated data
    mice_translated{mouse_idx, 2} = food_translated;
    mice_translated{mouse_idx, 3} = pre_data_translated;
    mice_translated{mouse_idx, 4} = test_data_translated;
    mice_translated{mouse_idx, 5} = center_translated;
end

fprintf('\n=== TRANSLATION COMPLETE ===\n');

%% Visualize results
figure('Name', sprintf('Translated Data - %s', condition_name), ...
       'Position', [100, 100, 1000, 800], 'Color', 'white');

hold on;

% Generate colors
colors = lines(n_mice);

% Plot each mouse after translation
for i = 1:n_mice
    mouse_idx = selected_mice(i);
    mouse_name = mice_translated{mouse_idx, 1};
    food_pos = mice_translated{mouse_idx, 2};
    pre_data = mice_translated{mouse_idx, 3};
    center_pos = mice_translated{mouse_idx, 5};
    
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
    
    % Plot center (should all be at target now)
    scatter(center_pos(1), center_pos(2), 150, colors(i,:), 'x', 'LineWidth', 3);
    
    % Plot food
    scatter(food_pos(1), food_pos(2), 150, colors(i,:), 'o', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
    
    % Plot vector from center to food
    dx = food_pos(1) - center_pos(1);
    dy = food_pos(2) - center_pos(2);
    quiver(center_pos(1), center_pos(2), dx, dy, 0, ...
           'Color', colors(i,:), 'LineWidth', 2, 'MaxHeadSize', 1);
    
    % Add label
    text(food_pos(1) + 10, food_pos(2) + 10, mouse_name, ...
         'FontSize', 9, 'Color', colors(i,:), 'FontWeight', 'bold');
end

% Mark the target center with a large marker
scatter(target_center(1), target_center(2), 300, 'k', 'x', 'LineWidth', 4);
text(target_center(1) + 15, target_center(2) + 15, 'Target Center', ...
     'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

axis equal;
grid on;
xlabel('X Position', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Translated Data - All Centers at [%.0f, %.0f] - %s (n=%d)', ...
      target_center(1), target_center(2), condition_name, n_mice), ...
      'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

hold off;

end