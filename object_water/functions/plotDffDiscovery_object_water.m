function plotDffDiscovery_object_water(mice_all)
% PLOTMOUSEDFFWITHDISCOVERIES Creates plots of raw DFF signals over time
% for each mouse with three subplots (before, water, object) and discovery markers
%
% Parameters:
%   mice_all - Cell array with mouse data

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)

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
    figure('Position', [100, 100, 1200, 600]);
    sgtitle(['Mouse ' current_mouse ' - DFF Signals with Discovery, Grooming and Climbing']);
    
    % Process before session (subplot 1)
    if ~isnan(before_idx)
        subplot(3, 1, 1);
        plotSessionDff(mice_all, before_idx, 'before');
    end
    
    % Process water session (subplot 2)
    if ~isnan(water_idx)
        subplot(3, 1, 2);
        plotSessionDff(mice_all, water_idx, 'water');
    end
    
    % Process object session (subplot 3)
    if ~isnan(object_idx)
        subplot(3, 1, 3);
        plotSessionDff(mice_all, object_idx, 'object');
    end
    
    % Add a legend as a textbox
    legnd = {'DFF Signal', 'Grooming', 'Climbing', 'Discovery Event'};
    legndh = nan(1, length(legnd));
    
    % Add legend manually since overlapping patches make automatic legend difficult
    annotation('textbox', [0.75, 0.85, 0.2, 0.1], 'String', legnd, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
        'EdgeColor', 'black', 'FontSize', 9);
    
    % Print session info
    fprintf('Mouse %s: before=%d, water=%d, object=%d\n', current_mouse, ...
        ~isnan(before_idx), ~isnan(water_idx), ~isnan(object_idx));
end
end

function plotSessionDff(mice_all, session_idx, session_type)
% Helper function to plot DFF for a single session

% Constants
col_time = 1;     % Time column
col_dff = 6;      % DFF data column
col_grooming = 13;   % Grooming indicator
col_climbing = 14;   % Climbing indicator (if available)

% Get session data
data = mice_all{session_idx, 4};
discovery = mice_all{session_idx, 6};
session_id = mice_all{session_idx, 1};

% Get the raw time and DFF data
time = data(:, col_time);
dff = data(:, col_dff);

% Create the base DFF plot
plot(time, dff, 'k-', 'LineWidth', 1);
hold on;

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
            y_limits = ylim();
            
            % Create pink background for grooming
            patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
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
            y_limits = ylim();
            
            % Create green background for climbing
            patch([t_start t_end t_end t_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
end

% Add discovery line if applicable (not for 'before' sessions)
if ~strcmp(session_type, 'before')
    if strcmp(session_type, 'water')
        disc_idx = discovery(1); % Water discovery
        disc_color = 'b';        % Blue for water
    else % object
        disc_idx = discovery(2); % Object discovery
        disc_color = 'r';        % Red for object
    end
    
    if ~isnan(disc_idx)
        discovery_time = time(disc_idx);
        xline(discovery_time, '--', 'Color', disc_color, 'LineWidth', 2);
    end
end

% Format plot
xlabel('Time (s)');
ylabel('ΔF/F');
title([upper(session_type(1)) session_type(2:end) ' Session']);
grid off;

% Add session ID as subtitle
title_parts = strsplit(session_id, '_');
if length(title_parts) >= 4
    subtitle([title_parts{1} ' - ' title_parts{end}]);
end

% Adjust y-axis to be symmetric if DFF is centered around zero
y_range = max(abs(dff)) * 1.1;
if mean(dff) < y_range/5
    ylim([-y_range y_range]);
end
end