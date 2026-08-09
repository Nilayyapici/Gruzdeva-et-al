function plotDoorEncounterAligned(mice_all, options)
% PLOTDOORENCOUNTERALIGNED Creates door encounter-aligned visualization with raster and average plots
%
% Door encounter is defined as the first time the mouse crosses below a distance threshold
% to the food location during closed door sessions (sess0 and sess2 only).
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - sessions: 'sess0', 'sess2', or 'both' (default: 'both')
%     - distance_threshold: distance threshold for door encounter in cm (default: 4)
%     - time_window: two-element vector specifying time range around encounter [before after] in seconds (default: [-20 60])
%     - ylim_avg: y-axis limits for average plot [min max] (default: auto)
%     - caxis_range: two-element vector specifying color axis range for DFF (default: [-2 2]) - DEPRECATED, use clim
%     - clim: two-element vector specifying color limits for raster plot [min max] (default: auto)
%     - z_score: whether to z-score the DFF signals for each encounter individually (default: false)
%     - smooth_window: window size for smoothing in frames (default: 15)
%     - min_time_between_encounters: minimum time between encounters to consider separate events in seconds (default: 30)
%     - show_validation: whether to show validation plots with trajectories (default: false)
%     - colormap: colormap for raster plots (default: orangewhiteblue(256))

% Constants
COL_TIME = 1;     % Time column
COL_X = 2;        % X coordinate
COL_Y = 3;        % Y coordinate
COL_DIST = 5;     % Distance to food column
COL_DFF = 11;     % DFF data column

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'state')
    options.state = 'all';
end

if ~isfield(options, 'source')
    options.source = 'all';
end

if ~isfield(options, 'sessions')
    options.sessions = 'both'; % Default to both sess0 and sess2
end

if ~isfield(options, 'distance_threshold')
    options.distance_threshold = 4; % Default 4 cm
end

if ~isfield(options, 'time_window')
    options.time_window = [-20, 60]; % Default in seconds
end

if ~isfield(options, 'caxis_range')
    options.caxis_range = [-2, 2];
end

if ~isfield(options, 'z_score')
    options.z_score = false;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_time_between_encounters')
    options.min_time_between_encounters = 30; % Default 30 seconds
end

if ~isfield(options, 'show_validation')
    options.show_validation = false;
end

if ~isfield(options, 'colormap')
    options.colormap = orangewhiteblue(256); % Default colormap
end

if ~isfield(options, 'clim')
    options.clim = []; % Empty means auto-adjust
end

% Filter mice based on state, source, and session (only sess0 and/or sess2)
n = size(mice_all, 1);
mask = true(n, 1);

% Apply state filter
if ~strcmpi(options.state, 'all')
    mask = mask & strcmp(mice_all(:, 2), options.state);
end

% Apply source filter
if ~strcmpi(options.source, 'all')
    mask = mask & strcmp(mice_all(:, 3), options.source);
end

% Apply session filter (only sess0 and/or sess2 - closed door sessions)
session_mask = false(n, 1);
for i = 1:n
    session_info = mice_all{i, 1};
    if strcmp(options.sessions, 'sess0') && contains(session_info, 'sess0')
        session_mask(i) = true;
    elseif strcmp(options.sessions, 'sess2') && contains(session_info, 'sess2')
        session_mask(i) = true;
    elseif strcmp(options.sessions, 'both') && (contains(session_info, 'sess0') || contains(session_info, 'sess2'))
        session_mask(i) = true;
    end
end
mask = mask & session_mask;

% Extract filtered mice
filtered_mice = mice_all(mask, :);
num_mice = size(filtered_mice, 1);

if num_mice == 0
    error('No mice found with the specified criteria for sess0/sess2');
end

fprintf('Processing %d mice for door encounter analysis (%s)\n', num_mice, options.sessions);
fprintf('Distance threshold: %.1f cm\n', options.distance_threshold);

% Initialize cell arrays to store aligned data
all_mouse_aligned_dff = cell(num_mice, 1);
all_mouse_encounter_counts = zeros(num_mice, 1);
all_validation_data = cell(num_mice, 1); % For validation plots

% Process each mouse to find door encounters and align data
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    distance = data(:, COL_DIST);
    
    % Quick validation of distance data
    min_dist = min(distance);
    max_dist = max(distance);
    if i <= 3  % Show first 3 mice for debugging
        fprintf('Mouse %d (%s): distance range %.1f - %.1f cm\n', ...
                i, filtered_mice{i, 1}, min_dist, max_dist);
    end
    
    % Detect door encounters (threshold crossings)
    encounter_indices = detectDoorEncounters(distance, time, options.distance_threshold, options.min_time_between_encounters);
    
    num_encounters = length(encounter_indices);
    all_mouse_encounter_counts(i) = num_encounters;
    
    % Store validation data
    all_validation_data{i} = struct('data', data, 'encounters', encounter_indices, 'mouse_id', filtered_mice{i, 1});
    
    if num_encounters == 0
        fprintf('No door encounters found for mouse %d (%s)\n', i, filtered_mice{i, 1});
        continue;
    end
    
    fprintf('Found %d door encounters for mouse %d (%s)\n', num_encounters, i, filtered_mice{i, 1});
    
    % Initialize cell array to store aligned data for each encounter
    all_encounters_aligned_dff = cell(num_encounters, 1);
    all_encounters_aligned_time = cell(num_encounters, 1);
    
    % Process each encounter
    for j = 1:num_encounters
        encounter_idx = encounter_indices(j);
        
        % Calculate time range to extract (in indices)
        encounter_time = time(encounter_idx);
        
        % Find indices within the time window
        time_mask = time >= (encounter_time + options.time_window(1)) & ...
                    time <= (encounter_time + options.time_window(2));
        
        % Extract aligned data
        aligned_time = time(time_mask) - encounter_time;
        aligned_dff = dff(time_mask);
        
        % Apply z-scoring if requested (same as plotFoodEventsAligned)
        if isfield(options, 'z_score') && options.z_score
            aligned_dff = (aligned_dff - nanmean(aligned_dff)) / nanstd(aligned_dff);
        end
        
        % Apply smoothing if requested
        if options.smooth_window > 0
            aligned_dff = movmean(aligned_dff, options.smooth_window, 'omitnan');
        end
        
        % Store aligned data for this encounter
        all_encounters_aligned_dff{j} = aligned_dff;
        all_encounters_aligned_time{j} = aligned_time;
    end
    
    % Store all encounters for this mouse
    all_mouse_aligned_dff{i} = struct('encounters_dff', {all_encounters_aligned_dff}, ...
                                      'encounters_time', {all_encounters_aligned_time});
end

% Show validation plots if requested
if options.show_validation
    plotDoorEncounterValidation(all_validation_data, options);
end

% Determine the common time grid for plotting
time_min = options.time_window(1);
time_max = options.time_window(2);
time_grid = linspace(time_min, time_max, 1000); % 1000 time points for smooth visualization

% Interpolate all encounters onto common time grid for each mouse
mouse_avg_dff = nan(num_mice, length(time_grid));

for i = 1:num_mice
    if all_mouse_encounter_counts(i) == 0
        continue;
    end
    
    % Initialize matrix to store interpolated data for all encounters
    encounters_interp = nan(all_mouse_encounter_counts(i), length(time_grid));
    
    % Interpolate each encounter onto common time grid
    for j = 1:all_mouse_encounter_counts(i)
        encounter_time = all_mouse_aligned_dff{i}.encounters_time{j};
        encounter_dff = all_mouse_aligned_dff{i}.encounters_dff{j};
        
        if isempty(encounter_time) || isempty(encounter_dff)
            continue;
        end
        
        % Interpolate this encounter
        encounters_interp(j, :) = interp1(encounter_time, encounter_dff, time_grid, 'linear', NaN);
    end
    
    % Calculate mouse average across all encounters
    mouse_avg_dff(i, :) = nanmean(encounters_interp, 1);
    
    % Validate mouse data quality
    if all_mouse_encounter_counts(i) > 0
        % Check for constant values (suspicious)
        dff_range = nanmax(mouse_avg_dff(i, :)) - nanmin(mouse_avg_dff(i, :));
        if dff_range < 0.001  % Essentially constant
            fprintf('  WARNING: Mouse %d (%s) has nearly constant dF/F (range=%.6f)\n', ...
                    i, filtered_mice{i, 1}, dff_range);
        end
        
        % Check data coverage
        valid_points = sum(~isnan(mouse_avg_dff(i, :)));
        coverage = valid_points / length(time_grid);
        if coverage < 0.5
            fprintf('  WARNING: Mouse %d (%s) has poor data coverage (%.1f%%)\n', ...
                    i, filtered_mice{i, 1}, coverage * 100);
        end
    end
end

% No additional z-scoring needed - already applied to individual encounters if requested

% Create figure with two subplots
figure('Position', [100, 100, 1000, 600]);

% 1. Average DFF plot (left subplot)
subplot(1, 2, 1);
hold on;

% Calculate overall average DFF across mice
avg_dff = nanmean(mouse_avg_dff, 1);
sem_dff = nanstd(mouse_avg_dff, 0, 1) ./ sqrt(sum(~isnan(mouse_avg_dff), 1));

% Add encounter onset line
xline(0, '--k', 'LineWidth', 1.5);

% Use orange/yellow colors for door encounters
encounter_color = [0.8, 0.5, 0.0]; % Orange
sem_color = [1.0, 0.8, 0.6]; % Light orange

% Add SEM shading
ciplot(avg_dff - sem_dff, avg_dff + sem_dff, time_grid, sem_color, 0.5);

% Plot the mean line on top
plot(time_grid, avg_dff, 'Color', encounter_color, 'LineWidth', 2);

% Format plot
xlabel('Time from Door Encounter (seconds)', 'FontSize', 12);
if isfield(options, 'z_score') && options.z_score
    ylabel('ΔF/F (z-scored)', 'FontSize', 12);
else
    ylabel('ΔF/F', 'FontSize', 12);
end
title({'ΔF/F aligned to Door Encounter', ...
       sprintf('Distance < %.1f cm', options.distance_threshold)});
grid off;
box off;

% Set x-axis limits
xlim(options.time_window);

% Set y-axis limits if provided
if isfield(options, 'ylim_avg')
    ylim(options.ylim_avg);
end

% 2. Raster plot (right subplot) - Filter for mice with >1 encounters
subplot(1, 2, 2);

% Filter for mice with more than 1 encounter
mice_with_sufficient_encounters = all_mouse_encounter_counts > 1;
filtered_mouse_avg_dff = mouse_avg_dff(mice_with_sufficient_encounters, :);
filtered_mice_info = filtered_mice(mice_with_sufficient_encounters, :);
filtered_encounter_counts = all_mouse_encounter_counts(mice_with_sufficient_encounters);
num_mice_raster = sum(mice_with_sufficient_encounters);

if num_mice_raster == 0
    % If no mice have >1 encounters, show a message
    text(0.5, 0.5, {'No mice with >1 encounters', 'Cannot display raster plot'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'center', ...
         'FontSize', 14, 'Units', 'normalized');
    xlabel('Time from Door Encounter (seconds)', 'FontSize', 12);
    ylabel('Mouse #', 'FontSize', 12);
    title({'Activity across mice', ...
           'No mice with >1 encounters'});
else
    % Create raster image with filtered data
    imagesc(time_grid, 1:num_mice_raster, filtered_mouse_avg_dff);
    
    % Add encounter onset line
    hold on;
    xline(0, '--k', 'LineWidth', 1.5);
    
    % Format plot
    xlabel('Time from Door Encounter (seconds)', 'FontSize', 12);
    ylabel('Mouse #', 'FontSize', 12);
    title({'Activity across mice', ...
           ['n = ' num2str(num_mice_raster) ' mice, ' ...
            num2str(sum(filtered_encounter_counts)) ' encounters']});
    
    % Use specified colormap
    colormap(options.colormap);
    h = colorbar;
    if isfield(options, 'z_score') && options.z_score
        ylabel(h, 'ΔF/F (z-scored)', 'FontSize', 12);
    else
        ylabel(h, 'ΔF/F', 'FontSize', 12);
    end
    
    % Set color limits - use clim if provided, otherwise use caxis_range
    if ~isempty(options.clim)
        caxis(options.clim);
    else
        caxis(options.caxis_range);
    end
    
    % Add mouse IDs and count of encounters for filtered mice
    if num_mice_raster <= 30
        mouse_labels = cell(num_mice_raster, 1);
        for i = 1:num_mice_raster
            if ischar(filtered_mice_info{i, 1})
                id_parts = strsplit(filtered_mice_info{i, 1}, '_');
                if length(id_parts) > 1
                    mouse_id = id_parts{1};
                else
                    mouse_id = filtered_mice_info{i, 1};
                end
            else
                mouse_id = num2str(i);
            end
            mouse_labels{i} = [mouse_id, ' (', num2str(filtered_encounter_counts(i)), ')'];
        end
        set(gca, 'YTick', 1:num_mice_raster, 'YTickLabel', mouse_labels);
    else
        % Just show mouse numbers if there are too many
        set(gca, 'YTick', round(linspace(1, num_mice_raster, 10)));
    end
    
    % Set x-axis limits
    xlim(options.time_window);
end

% Add overall title with threshold information
if isfield(options, 'z_score') && options.z_score
    z_score_text = 'Z-scored';
else
    z_score_text = 'Raw';
end
sgtitle({[z_score_text, ' ΔF/F values - Door Encounters (< ', ...
         num2str(options.distance_threshold), ' cm) - ', ...
         options.state, ' mice, ', options.source], ...
         ['Sessions: ', options.sessions]}, 'FontWeight', 'bold');

% Output statistics
total_encounters = sum(all_mouse_encounter_counts);
mice_with_encounters = sum(all_mouse_encounter_counts > 0);

fprintf('\n=== ANALYSIS SUMMARY ===\n');
fprintf('Distance threshold: %.1f cm\n', options.distance_threshold);
fprintf('Mice processed: %d\n', num_mice);
fprintf('Mice with encounters: %d (%.1f%%)\n', mice_with_encounters, 100 * mice_with_encounters / num_mice);
fprintf('Total encounters: %d\n', total_encounters);
fprintf('Mice in raster plot: %d (>1 encounter)\n', num_mice_raster);

% Validation warnings
if options.distance_threshold > 50
    fprintf('WARNING: Distance threshold (%.1f cm) seems unusually high!\n', options.distance_threshold);
end

if mice_with_encounters == 0
    fprintf('WARNING: No encounters detected - check distance threshold\n');
elseif mice_with_encounters / num_mice < 0.1
    fprintf('WARNING: Very few mice (%.1f%%) have encounters - threshold may be too strict\n', ...
            100 * mice_with_encounters / num_mice);
end

if total_encounters > 0
    fprintf('Average pre-encounter DFF: %.3f\n', nanmean(avg_dff(time_grid < 0)));
    fprintf('Average during-encounter DFF: %.3f\n', nanmean(avg_dff(time_grid >= 0)));
end
fprintf('========================\n');

% =========================================================
% SAVE FIGURE DATA TO EXCEL
% =========================================================
% Build a filename based on filter options
excel_filename = sprintf('DoorEncounter_%s_%s_%s_thr%.0fcm.xlsx', ...
    options.state, options.source, options.sessions, options.distance_threshold);

% --- Sheet 1: Average trace with SEM ---
% avg_dff and sem_dff are already computed above
trace_table = table(time_grid(:), avg_dff(:), sem_dff(:), ...
    avg_dff(:) - sem_dff(:), avg_dff(:) + sem_dff(:), ...
    'VariableNames', {'Time_s', 'Mean_dFF', 'SEM_dFF', 'Lower_CI', 'Upper_CI'});

writetable(trace_table, excel_filename, 'Sheet', 'Average_Trace');

% --- Sheet 2: Per-mouse raster (one row per mouse) ---
% Use the full mouse_avg_dff matrix (all mice, before the >1 encounter filter)
% Rows = mice, columns = time points
% Build a header row: mouse ID + session info, then one column per time point
raster_mouse_labels = cell(num_mice, 1);
for i = 1:num_mice
    raster_mouse_labels{i} = sprintf('%s_%s_%s', ...
        filtered_mice{i,1}, filtered_mice{i,2}, filtered_mice{i,3});
end

% Combine mouse labels with DFF data
time_header = arrayfun(@(t) sprintf('t=%.3f', t), time_grid, 'UniformOutput', false);
col_names   = [{'Mouse_ID'}, time_header];

raster_data = [raster_mouse_labels, num2cell(mouse_avg_dff)];
raster_table = cell2table(raster_data, 'VariableNames', col_names);

writetable(raster_table, excel_filename, 'Sheet', 'Raster_PerMouse');

fprintf('Figure data saved to: %s\n', excel_filename);
% =========================================================

% Optional: return the aligned data if requested
if nargout > 0
    aligned_data = struct();
    aligned_data.mouse_avg_dff = mouse_avg_dff;
    aligned_data.filtered_mouse_avg_dff = filtered_mouse_avg_dff;
    aligned_data.mice_with_sufficient_encounters = mice_with_sufficient_encounters;
    aligned_data.time_grid = time_grid;
    aligned_data.encounter_counts = all_mouse_encounter_counts;
    aligned_data.mouse_info = filtered_mice;
    aligned_data.all_encounters = all_mouse_aligned_dff;
    aligned_data.validation_data = all_validation_data;
    varargout{1} = aligned_data;
end
end

function encounter_indices = detectDoorEncounters(distance, time, threshold, min_time_between)
% DETECTDOORENCOUNTERS Detect when mouse first crosses below distance threshold
%
% Parameters:
%   distance - Vector of distances to food
%   time - Vector of time points
%   threshold - Distance threshold for encounter
%   min_time_between - Minimum time between encounters (seconds)
%
% Returns:
%   encounter_indices - Indices where encounters occurred

encounter_indices = [];

% Check if any points are actually above threshold first
if all(distance <= threshold)
    fprintf('  No encounters detected: all distances <= %.1f cm threshold\n', threshold);
    return;
end

if all(distance > threshold)
    fprintf('  No encounters detected: all distances > %.1f cm threshold\n', threshold);
    return;
end

% Find all points where distance crosses below threshold
% (was above threshold, now below threshold)
above_threshold = distance > threshold;

% Find crossing points (from above to below) - fix the artificial crossing bug
crossings = [];
for i = 2:length(above_threshold)
    if above_threshold(i-1) && ~above_threshold(i)  % Crossing from above to below
        crossings = [crossings; i];
    end
end

if isempty(crossings)
    fprintf('  No threshold crossings detected\n');
    return;
end

% Filter crossings to ensure minimum time between encounters
if length(crossings) == 1
    encounter_indices = crossings;
else
    encounter_indices = crossings(1); % Always include first crossing
    
    for i = 2:length(crossings)
        % Check if enough time has passed since last encounter
        time_since_last = time(crossings(i)) - time(encounter_indices(end));
        
        if time_since_last >= min_time_between
            encounter_indices = [encounter_indices; crossings(i)];
        end
    end
end

fprintf('  Detected %d threshold crossings, filtered to %d encounters (min %.1fs apart)\n', ...
        length(crossings), length(encounter_indices), min_time_between);
end

function plotDoorEncounterValidation(validation_data, options)
% PLOTDOORENCOUNTERVALIDATION Create validation plots showing trajectories with encounter markers
%
% Parameters:
%   validation_data - Cell array with validation data for each mouse
%   options - Options struct with distance_threshold

% Constants
COL_TIME = 1;
COL_X = 2;
COL_Y = 3;
COL_DIST = 5;

% Create figure for validation plots
figure('Position', [100, 100, 1400, 1000]);

% Calculate number of subplots needed
num_mice = length(validation_data);
num_cols = min(4, num_mice);
num_rows = ceil(num_mice / num_cols);

for i = 1:num_mice
    if isempty(validation_data{i})
        continue;
    end
    
    data = validation_data{i}.data;
    encounters = validation_data{i}.encounters;
    mouse_id = validation_data{i}.mouse_id;
    
    % Extract trajectory data
    x = data(:, COL_X);
    y = data(:, COL_Y);
    time = data(:, COL_TIME);
    distance = data(:, COL_DIST);
    
    % Create subplot
    subplot(num_rows, num_cols, i);
    
    % Plot trajectory colored by distance to food
    scatter(x, y, 10, distance, 'filled');
    hold on;
    
    % Mark encounter points
    if ~isempty(encounters)
        encounter_x = x(encounters);
        encounter_y = y(encounters);
        plot(encounter_x, encounter_y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        
        % Add numbers to encounters
        for j = 1:length(encounters)
            text(encounter_x(j), encounter_y(j), num2str(j), ...
                 'Color', 'white', 'FontSize', 8, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end
    end
    
    % Format plot
    colorbar;
    caxis([0, max(distance)]);
    colormap(jet);
    title(sprintf('%s\n%d encounters', mouse_id, length(encounters)));
    xlabel('X position');
    ylabel('Y position');
    axis equal;
    grid on;
    
    % Find actual food location (where distance is minimum)
    [min_dist, min_idx] = min(distance);
    food_x = x(min_idx);
    food_y = y(min_idx);
    
    % Mark the actual food location
    plot(food_x, food_y, 'ks', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'yellow');
    
    % Add threshold circle around the ACTUAL food location
    theta = linspace(0, 2*pi, 100);
    threshold_x = food_x + options.distance_threshold * cos(theta);
    threshold_y = food_y + options.distance_threshold * sin(theta);
    plot(threshold_x, threshold_y, 'k--', 'LineWidth', 2);
    
    % Add legend for clarity
    if i == 1  % Only add legend to first subplot
        legend('Trajectory', 'Encounters', 'Food Location', 'Threshold Circle', 'Location', 'best');
    end
end

sgtitle(sprintf('Door Encounter Validation - Threshold: %.1f cm', options.distance_threshold));

% Create a second figure showing distance vs time for validation
figure('Position', [100, 100, 1400, 800]);

for i = 1:min(6, num_mice) % Show only first 6 mice to avoid crowding
    if isempty(validation_data{i})
        continue;
    end
    
    data = validation_data{i}.data;
    encounters = validation_data{i}.encounters;
    mouse_id = validation_data{i}.mouse_id;
    
    time = data(:, COL_TIME);
    distance = data(:, COL_DIST);
    
    % Create subplot
    subplot(2, 3, i);
    
    % Plot distance vs time
    plot(time, distance, 'b-', 'LineWidth', 1);
    hold on;
    
    % Add threshold line
    yline(options.distance_threshold, 'k--', 'Threshold', 'LineWidth', 2);
    
    % Mark encounters
    if ~isempty(encounters)
        encounter_times = time(encounters);
        encounter_distances = distance(encounters);
        plot(encounter_times, encounter_distances, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        
        % Add encounter numbers
        for j = 1:length(encounters)
            text(encounter_times(j), encounter_distances(j) + 1, num2str(j), ...
                 'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end
    end
    
    % Format plot
    title(sprintf('%s - %d encounters', mouse_id, length(encounters)));
    xlabel('Time (s)');
    ylabel('Distance to Food (cm)');
    ylim([0, max(distance) * 1.1]);
    grid on;
end

sgtitle('Distance vs Time - Door Encounter Detection');
end

function ciplot(lower, upper, x, color, alpha)
% CIPLOT creates a shaded area to visualize confidence intervals
x_polygon = [x, fliplr(x)];
y_polygon = [upper, fliplr(lower)];

% Remove any NaN values
nan_indices = isnan(x_polygon) | isnan(y_polygon);
x_polygon(nan_indices) = [];
y_polygon(nan_indices) = [];

% Create the polygon
h = fill(x_polygon, y_polygon, color);
set(h, 'EdgeColor', 'none');
set(h, 'FaceAlpha', alpha);
end

function cmap = orangewhiteblue(m)
% ORANGEWHITEBLUE creates an orange-to-white-to-blue colormap for door encounters
if nargin < 1
   m = 256;
end

% Define key colors
colors = [
    0.0, 0.4, 0.8;  % Blue
    0.4, 0.7, 1.0;  % Light blue
    1.0, 1.0, 1.0;  % White (center)
    1.0, 0.7, 0.3;  % Light orange
    0.8, 0.4, 0.0   % Dark orange
];

% Positions for the colors
positions = [0, 0.25, 0.5, 0.75, 1.0];

% Interpolate
n = size(colors, 1);
cmap = zeros(m, 3);
idx = round(positions * (m-1)) + 1;

for i = 1:n-1
    i1 = idx(i);
    i2 = idx(i+1);
    ni = i2 - i1 + 1;
    
    r = linspace(colors(i,1), colors(i+1,1), ni);
    g = linspace(colors(i,2), colors(i+1,2), ni);
    b = linspace(colors(i,3), colors(i+1,3), ni);
    
    cmap(i1:i2, :) = [r', g', b'];
end
end