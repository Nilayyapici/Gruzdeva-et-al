function plotMouseDffWithInteractions(mice_all, options)
% PLOTMOUSEDFFWITHINTERACTIONS Creates plots of raw DFF signals with highlighted
% interaction periods for each mouse across before, water, and object sessions
%
% Parameters:
%   mice_all - Cell array with mouse data
%   options - Struct with visualization parameters:
%     - lim_dist: distance threshold for interaction in cm (default: 3)

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_dist_water = 9;  % Distance to water
col_dist_object = 10; % Distance to object
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)
col_water_int = 15;  % Water interaction indicator (added column)
col_object_int = 16; % Object interaction indicator (added column)

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'lim_dist')
    options.lim_dist = 3; % Default distance threshold in cm
end

% Extract all unique mouse IDs from the dataset
mouse_ids = {};
for i = 1:size(mice_all, 1)
    current_id = mice_all{i, 1};
    if ischar(current_id)
        % Extract the mouse ID (e.g., 'F29' from 'F29_water_object_data_water')
        parts = strsplit(current_id, '_');
        mouse_id = parts{1};

        % Add to list if not already there
        if ~ismember(mouse_id, mouse_ids)
            mouse_ids{end+1} = mouse_id;
        end
    end
end

fprintf('Found %d unique mice in the dataset\n', length(mouse_ids));

% Process each mouse individually
for m = 1:length(mouse_ids)
    current_mouse = mouse_ids{m};
    fprintf('Processing mouse %s (%d of %d)\n', current_mouse, m, length(mouse_ids));

    % Find sessions for this mouse: before, water, and object
    before_idx = NaN;
    water_idx = NaN;
    object_idx = NaN;

    for i = 1:size(mice_all, 1)
        session_id = mice_all{i, 1};
        if ischar(session_id) && startsWith(session_id, current_mouse)
            if contains(session_id, 'before')
                before_idx = i;
            elseif contains(session_id, 'water')
                water_idx = i;
            elseif contains(session_id, 'object')
                object_idx = i;
            end
        end
    end

    % Create figure for this mouse
    figure('Position', [100, 100, 1200, 900]);
    sgtitle(['Mouse ' current_mouse ' - DFF Signals with Interactions (dist < ' num2str(options.lim_dist) ' cm)']);

    % Process before session (subplot 1)
    if ~isnan(before_idx)
        subplot(3, 1, 1);
        plotSessionWithInteractions(mice_all, before_idx, 'before', options);
    else
        subplot(3, 1, 1);
        title('Before Session - No Data Available');
        axis off;
    end

    % Process water session (subplot 2)
    if ~isnan(water_idx)
        subplot(3, 1, 2);
        plotSessionWithInteractions(mice_all, water_idx, 'water', options);
    else
        subplot(3, 1, 2);
        title('Water Session - No Data Available');
        axis off;
    end

    % Process object session (subplot 3)
    if ~isnan(object_idx)
        subplot(3, 1, 3);
        plotSessionWithInteractions(mice_all, object_idx, 'object', options);
    else
        subplot(3, 1, 3);
        title('Object Session - No Data Available');
        axis off;
    end

    % Add a better legend with colored boxes
    legend_pos = [0.72, 0.93, 0.25, 0.06];
    annotation('textbox', legend_pos, 'String', '', 'EdgeColor', 'k', 'BackgroundColor', 'w');

    % Add colored boxes with text for each item
    box_width = 0.01;
    box_height = 0.01;
    text_offset = 0.025;
    start_y = legend_pos(2) - 0.01;
    start_x = legend_pos(1) + 0.01;

    % DFF Signal
    annotation('line', [start_x, start_x+box_width], [start_y, start_y], 'Color', 'k', 'LineWidth', 2);
    annotation('textbox', [start_x+text_offset, start_y-0.005, 0.1, 0.01], 'String', 'DFF Signal', ...
        'EdgeColor', 'none', 'FontSize', 9);

    % Grooming
    start_y = start_y - 0.012;
    annotation('rectangle', [start_x, start_y, box_width, box_height], 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    annotation('textbox', [start_x+text_offset, start_y-0.005, 0.1, 0.01], 'String', 'Grooming', ...
        'EdgeColor', 'none', 'FontSize', 9);

    % Climbing
    start_y = start_y - 0.012;
    annotation('rectangle', [start_x, start_y, box_width, box_height], 'FaceColor', 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    annotation('textbox', [start_x+text_offset, start_y-0.005, 0.1, 0.01], 'String', 'Climbing', ...
        'EdgeColor', 'none', 'FontSize', 9);

    % Water Interaction
    start_y = start_y - 0.012;
    annotation('rectangle', [start_x, start_y, box_width, box_height], 'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    annotation('textbox', [start_x+text_offset, start_y-0.005, 0.1, 0.01], 'String', 'Water Interaction', ...
        'EdgeColor', 'none', 'FontSize', 9);

    % Object Interaction
    start_y = start_y - 0.012;
    annotation('rectangle', [start_x, start_y, box_width, box_height], 'FaceColor', [1 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    annotation('textbox', [start_x+text_offset, start_y-0.005, 0.1, 0.01], 'String', 'Object Interaction', ...
        'EdgeColor', 'none', 'FontSize', 9);

    % Print session info
    fprintf('Mouse %s: before=%d, water=%d, object=%d\n', current_mouse, ...
        ~isnan(before_idx), ~isnan(water_idx), ~isnan(object_idx));
end
end

function plotSessionWithInteractions(mice_all, session_idx, session_type, options)
% Helper function to plot DFF for a single session with interactions highlighted

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)
col_water_int = 15;  % Water interaction indicator (added column)
col_object_int = 16; % Object interaction indicator (added column)

% Get session data
data = mice_all{session_idx, 4};
session_id = mice_all{session_idx, 1};

% Get the raw time and DFF data
time = data(:, col_time);
dff = data(:, col_dff);

% Get interaction data
if size(data, 2) >= col_water_int
    water_interaction = data(:, col_water_int);
else
    water_interaction = zeros(size(time));
end

if size(data, 2) >= col_object_int
    object_interaction = data(:, col_object_int);
else
    object_interaction = zeros(size(time));
end

% Create the base DFF plot
plot(time, dff, 'k-', 'LineWidth', 1);
hold on;

% Calculate y-axis limits for visualization
y_limits = [min(dff) - 0.1*abs(min(dff)), max(dff) + 0.1*abs(max(dff))];
if isnan(y_limits(1)) || isnan(y_limits(2))
    y_limits = [-1, 1]; % Default if NaN
end

% Highlight grooming periods with pink
if size(data, 2) >= col_grooming
    grooming = data(:, col_grooming);
    if any(grooming > 0)
        % Find contiguous grooming periods
        grooming_starts = find(diff([0; grooming]) > 0);
        grooming_ends = find(diff([grooming; 0]) < 0);

        for i = 1:length(grooming_starts)
            t_start = time(grooming_starts(i));
            t_end = time(grooming_ends(i));

            % Create pink background for grooming
            patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    end
end

% Highlight climbing periods with green
if size(data, 2) >= col_climbing
    climbing = data(:, col_climbing);
    if any(climbing > 0)
        % Find contiguous climbing periods
        climbing_starts = find(diff([0; climbing]) > 0);
        climbing_ends = find(diff([climbing; 0]) < 0);

        for i = 1:length(climbing_starts)
            t_start = time(climbing_starts(i));
            t_end = time(climbing_ends(i));

            % Create green background for climbing
            patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    end
end

% Highlight water interaction periods with blue
if any(water_interaction > 0)
    % Find contiguous interaction periods
    water_starts = find(diff([0; water_interaction]) > 0);
    water_ends = find(diff([water_interaction; 0]) < 0);

    for i = 1:length(water_starts)
        t_start = time(water_starts(i));
        t_end = time(water_ends(i));

        % Create blue background for water interaction
        patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Print number of water interactions
    fprintf('Session %s: %d water interaction periods found\n', session_id, length(water_starts));
else
    fprintf('Session %s: No water interactions found\n', session_id);
end

% Highlight object interaction periods with orange
if any(object_interaction > 0)
    % Find contiguous interaction periods
    object_starts = find(diff([0; object_interaction]) > 0);
    object_ends = find(diff([object_interaction; 0]) < 0);

    for i = 1:length(object_starts)
        t_start = time(object_starts(i));
        t_end = time(object_ends(i));

        % Create orange background for object interaction
        patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            [1 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Print number of object interactions
    fprintf('Session %s: %d object interaction periods found\n', session_id, length(object_starts));
else
    fprintf('Session %s: No object interactions found\n', session_id);
end

% Format plot
xlabel('Time (s)');
ylabel('ΔF/F');
title([upper(session_type(1)) session_type(2:end) ' Session']);
grid on;

% Add session ID as subtitle
title_parts = strsplit(session_id, '_');
if length(title_parts) >= 2
    subtitle_text = [title_parts{1}];
    if length(title_parts) >= 4
        subtitle_text = [subtitle_text ' - ' title_parts{end}];
    end
    subtitle(subtitle_text);
end

% Set y-axis limits consistently
ylim(y_limits);

% Add text showing interaction counts
water_count = sum(water_interaction);
object_count = sum(object_interaction);
text_str = sprintf('Total frames: Water interactions: %d, Object interactions: %d', ...
    water_count, object_count);
text(0.5, 0.95, text_str, 'Units', 'normalized', 'HorizontalAlignment', 'center');
end