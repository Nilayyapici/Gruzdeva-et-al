
clear all
folderPath = '';
cd(folderPath); % Change current directory to the folder path

% Get a list of all .mat files in the directory
files = dir('*.mat');

% Initialize a cell array to hold the data
mice_all = {};

% Define mappings for side2 based on varName substrings
% Order matters now - more specific (longer) names should come first
keySet = {'FDRE14CNO','FDRE15CNO','FDRE16CNO2','FDRE4CNO2','MDRE6CNO','MDRE8CNO','MDRI17CNO',...
    'FCRI15opto2','FCRI16opto2',...
    'MCRI10optotest','MCRI12optotest','FCRI15optotest','FCRI16optotest','MCRI18optotest','MCRI21optotest','MCRI22optotest',...
    'MCRI10control2','MCRI12control','FCRI15control','FCRI16control','MCRI18control','FCRI24control','FCRI26control',...
    'MCRI10control3','MCRI12control2','FCRI15control2','FCRI16control2','MCRI18control2','MCRI21control2','MCRI22control2','FCRI24control2','FCRI26control2',...
    'MCRI10opto','FCRI15opto','FCRI16opto','MCRI12opto','MCRI18opto','MCRI21opto','MCRI22opto','FCRI24opto','FCRI26opto',...
    'FDRE16_2','M21_2',...
    'F11','F12','F13','F14','F15','F29','F30','FDRE14','FDRE15','FDRE4','FDRI31','FDRI32','M21','M22','M24','M28','MDRE8','MDRI17','MCRI8control','MCRI10control','MCRI21','MCRI22','FDRE34',...
    'MCRI12darkcontrol','FCRI15darkcontrol','FCRI16darkcontrol','MCRI21darkcontrol','MCRI22darkcontrol', 'FCRI24darkcontrol','FCRI26darkcontrol',...
    'FCRI15darkcontrol2','FCRI16darkcontrol2','MCRI21darkcontrol2','MCRI22darkcontrol2', 'FCRI24darkcontrol2','FCRI26darkcontrol2',...
    'FCRI15darkopto','FCRI16darkopto','MCRI21darkopto','MCRI22darkopto','FCRI24darkopto','FCRI26darkopto'};
% Where food was availiable for the specific session in keySet
valueSet = {'left','left','center','left','right','left','center',...
    'left','right',...
    'center','left','right','left','left','center','center',...
    'left','center','center','center','left','center','left',...
    'right','right','left','right','center','left','left','left','right',...
    'right','right','left','right','center','left','left','right','left',...
    'left','left',...
    'center','left','right','left','right','left','left','center','center','right','right','right','left','right','right','left','center','left','right','right','right','right','right',...
    'left','center','left','right','left','center','center',...
    'center','left','left','right','right','center',...
    'right','right','center','center','left','right'};

% Verify that keySet and valueSet have the same length
if length(keySet) ~= length(valueSet)
    error('keySet and valueSet must have the same length');
end

sideMap = containers.Map(keySet, valueSet);

% Loop through each file in the directory
for i = 1:length(files)
    fileName = files(i).name;

    % Load the file
    data = load(fileName);

    % Assume the variable inside the .mat file has the same name as the file
    varName = erase(fileName, '.mat'); % Remove the .mat extension to get the variable name

    % Determine the 'side' based on the name convention
    % Check for most specific patterns first
    if contains(varName, 'opto') && contains(varName, 'CRI') && ~contains(varName, 'dark')
        if contains(varName, 'optotest')
            side = 'chrimson_optotest';
        elseif contains(varName, 'opto2')
            side = 'light_opto';
        else
            side = 'light_opto';
        end
    elseif contains(varName, 'opto') && ~contains(varName, 'dark')
        side = 'opto'; % For other types of opto mice
    elseif contains(varName, 'CNO') && ~contains(varName, 'dark')
        side = 'CNO';
    elseif contains(varName, 'CRI') && contains(varName, 'dark')
        if contains (varName, 'control2')
            side = 'dark_control';
        elseif contains (varName, 'control')
            side = 'dark_control';
        elseif contains(varName, 'opto')
            side = 'dark_opto';
        end
    elseif contains(varName, 'CRI') && contains(varName, 'control2') && ~contains(varName, 'dark') && ~contains(varName, 'MCRI10control2')
        side = 'light_control';
    elseif contains(varName, 'CRI') && contains(varName, 'control') && ~contains(varName, 'dark') && ~contains(varName, 'MCRI10control') && ~contains(varName, 'MCRI8control')
        side = 'light_control';
    elseif contains(varName, 'CRI') && ~contains(varName, 'dark') && ~contains(varName, 'light') && ~contains(varName, 'opto')
        side = 'unknown';
    elseif contains(varName, 'MCRI10control') && contains(varName, 'MCRI10control2')&& contains(varName, 'MCRI8control')
        side = 'unknown';
    else
        side = 'control';
    end

    % Determine 'side2' using the map
    % Sort keys by length in descending order to check longest matches first
    sortedKeys = keySet;
    [~, idx] = sort(cellfun(@length, sortedKeys), 'descend');
    sortedKeys = sortedKeys(idx);

    side2 = 'unknown'; % Default value if no match is found
    for k = 1:length(sortedKeys)
        if contains(varName, sortedKeys{k})
            side2 = sideMap(sortedKeys{k});
            break; % Stop searching once a match is found
        end
    end

    % Determine the signal quality based on the name convention (check the signal from all the mice in the preprocessing)
    if contains(varName, 'F12') || contains(varName, 'F11') || contains(varName, 'F15') || contains(varName, 'FDRE4')
        side3 = 'bad_signal';
    else
        side3 = 'ok_signal';
    end

    % Append the loaded data and associated info into the cell array
    % Ensure the variable name is valid and exists in 'data'
    if isfield(data, varName)
        mice_all{end+1, 1} = varName; % Name or identifier of the mouse session
        mice_all{end, 2} = side; % Treatment group (CNO, control, chrimson, chrimson opto, etc.)
        mice_all{end, 3} = side2; % Food arm location (center, left, right)
        mice_all{end, 4} = data.(varName); % Actual data from the file
        mice_all{end, 5} = side3; % Signal quality (ok_signal, bad_signal)
        % mice_all{end, 6} = side4; % Behavior quality (ok_signal, bad_signal)

        % Debug output to verify correct assignment
        fprintf('Loaded: %s -> Group: %s, Food arm: %s, Signal: %s\n', ...
            varName, side, side2, side3);
    end
end

% % save the assembled cell array for later use
% save('mice_all_data.mat', 'mice_all');

%data (mice_all{end, 4}):
% 1-time
% 2-x
% 3-y
% 4-465
% 5-405
% 6-dF/F
% 7-speed
% 8-zones (1 - )
% 9-distance to food
% 10-distance to food diff2 (left for center, center for side)
% 11-distance to food diff1 (right for center, opposite for side)
% 12-door 0 closed/1 open/2 closed for the second time
% 13-grooming

% Only mice with good signal
% Find indices of mice with 'ok_signal'
ok_signal_indices = strcmpi(mice_all(:, 5), 'ok_signal'); % Case-insensitive comparison

% Extract only the rows corresponding to 'ok_signal'
mice_good = mice_all(ok_signal_indices, :);
mice_all = mice_good;

%% Reorganize sessions - split sess1 into sess0 and sess1 based on door condition
% Create a new cell array to hold the reorganized data
mice_all_reorganized = {};
% Loop through each entry in the original data
for i = 1:size(mice_all, 1)
    session_info = mice_all{i, 1}; % Mouse session identifier
    side = mice_all{i, 2}; % condition
    side2 = mice_all{i, 3}; % Specific side (center, left, right)
    data = mice_all{i, 4}; % Actual data from the file
    side3 = mice_all{i, 5}; % Signal quality indicator
    % Determine the session number
    if contains(session_info, 'sess1')
        session_number = 1;
    elseif contains(session_info, 'sess2')
        session_number = 2;
    else
        session_number = NaN;
    end
    % If this is session 1, we need to split it into sess0 and sess1
    if session_number == 1
        % Extract data where door=0 (closed) for sess0
        data_sess0 = data(data(:,12) == 0, :);
        % Extract data where door=1 (open) for sess1
        data_sess1 = data(data(:,12) == 1, :);
        % Create a new entry for sess0 (modify the original session name)
        if ~isempty(data_sess0)
            sess0_name = strrep(session_info, 'sess1', 'sess0');
            row_idx = size(mice_all_reorganized, 1) + 1;
            mice_all_reorganized{row_idx, 1} = sess0_name;
            mice_all_reorganized{row_idx, 2} = side;
            mice_all_reorganized{row_idx, 3} = side2;
            mice_all_reorganized{row_idx, 4} = data_sess0;
            mice_all_reorganized{row_idx, 5} = side3;
        end
        % Keep the original entry for sess1 but with filtered data
        if ~isempty(data_sess1)
            row_idx = size(mice_all_reorganized, 1) + 1;
            mice_all_reorganized{row_idx, 1} = session_info;
            mice_all_reorganized{row_idx, 2} = side;
            mice_all_reorganized{row_idx, 3} = side2;
            mice_all_reorganized{row_idx, 4} = data_sess1;
            mice_all_reorganized{row_idx, 5} = side3;
        end
    else
        % For session 2 or any other, just add it as is
        row_idx = size(mice_all_reorganized, 1) + 1;
        mice_all_reorganized{row_idx, 1} = session_info;
        mice_all_reorganized{row_idx, 2} = side;
        mice_all_reorganized{row_idx, 3} = side2;
        mice_all_reorganized{row_idx, 4} = data;
        mice_all_reorganized{row_idx, 5} = side3;
    end
end

% Replace the original array with the reorganized one
mice_all = mice_all_reorganized;

% %Clean up temporary variables
% clear mice_all_reorganized session_info side side2 data side3 session_number;
% clear data_sess0 data_sess1 sess0_name row_idx i;

% % Optional: Save the reorganized cell array
% save('mice_all_reorganized.mat', 'mice_all');

%% Behavior
%% Time only in the zone limited by dist_lim and time_lim

options = struct();
options.plot_type = 'time';  % or 'raw' for absolute time
options.sessions = {'sess0','sess2'};  % Will show as "Before", "Learning", "Test"
options.group = 'light_control';  % or 'control' or 'CNO'
options.title = 'Time Spent in Maze Arms Across Sessions';
options.use_limits = true;  % Set to true if you want distance/time limits
options.memory_tracking = true;
options.dist_lim = 90;
options.time_lim = 10;  % minutes
options.horizontal = false;  % Set to true to include horizontal arm
options.connect_points = true;  % Connect individual mouse data points
options.show_stats = true;  % Show significance indicators
options.alpha = 0.05;  % Significance level
options.colors = {[0.6 0.6 0.6],[0.2 0.4 0.8]};  % Colors for sessions 0 and 2 Control
options.figure_size = [100, 100,400, 600];
options.font_size = 14;
options.title = 'Time in Maze Arms';
options.ylim = [0 500];
options.combine_nonfood_arms = true;

% Create the enhanced session comparison plot
[data_table, stats_results] = plotSessionComparisonWithStats(mice_all, options);

% Export individual mouse data to Excel
% Re-run analyzeZoneTimes to get the raw per-mouse data
[zone_data_export, ~] = analyzeZoneTimes(mice_all, options);

if strcmp(options.plot_type, 'time')
    data_to_export = zone_data_export.times_percent;
else
    data_to_export = zone_data_export.times;
end

% Arm labels
if isfield(options, 'combine_nonfood_arms') && options.combine_nonfood_arms
    arm_labels = {'Food_Arm', 'NonFood_Arms'};
elseif isfield(options, 'horizontal') && options.horizontal
    arm_labels = {'Food_Arm', 'NonFood_Arm1', 'NonFood_Arm2', 'Horizontal'};
else
    arm_labels = {'Food_Arm', 'NonFood_Arm1', 'NonFood_Arm2'};
end

% Session label map
sess_label_map = containers.Map({'sess0','sess1','sess2'}, {'Session0','Session1','Session2'});

% Build a combined table: rows = mice, columns = Session_Arm combinations
valid_sessions = intersect(fieldnames(data_to_export)', options.sessions, 'stable');

% Collect all unique mouse IDs across sessions
all_mouse_ids = {};
for s = 1:length(valid_sessions)
    sess = valid_sessions{s};
    ids = zone_data_export.mouse_ids.(sess);
    all_mouse_ids = union(all_mouse_ids, ids);
end
all_mouse_ids = sort(all_mouse_ids);

% Build column names and data matrix
col_names = {'Mouse_ID'};
for s = 1:length(valid_sessions)
    sess = valid_sessions{s};
    if isKey(sess_label_map, sess)
        sess_str = sess_label_map(sess);
    else
        sess_str = strrep(sess, 'sess', 'Session');
    end
    for a = 1:length(arm_labels)
        col_names{end+1} = sprintf('%s_%s', sess_str, arm_labels{a});
    end
end

% Fill data matrix
n_mice = length(all_mouse_ids);
n_cols = length(col_names) - 1;
data_matrix = NaN(n_mice, n_cols);

col_idx = 1;
for s = 1:length(valid_sessions)
    sess = valid_sessions{s};
    sess_data = data_to_export.(sess);          % [n_mice_in_sess x n_arms]
    sess_mouse_ids = zone_data_export.mouse_ids.(sess);

    for a = 1:length(arm_labels)
        for m = 1:n_mice
            mouse_row = find(strcmp(sess_mouse_ids, all_mouse_ids{m}));
            if ~isempty(mouse_row) && a <= size(sess_data, 2)
                data_matrix(m, col_idx) = sess_data(mouse_row(1), a);
            end
        end
        col_idx = col_idx + 1;
    end
end

% Create MATLAB table and write to Excel
export_table = array2table(data_matrix, 'VariableNames', col_names(2:end));
export_table.Mouse_ID = all_mouse_ids(:);
export_table = [export_table(:,end), export_table(:,1:end-1)];  % Mouse_ID first

% Save
% Determine session label for filename
if isequal(sort(options.sessions), sort({'sess0','sess2'}))
    sess_label = 'Test';
elseif isequal(sort(options.sessions), sort({'sess0','sess1'}))
    sess_label = 'Learning';
else
    sess_label = strjoin(options.sessions, '_');
end

excel_filename = sprintf('individual_mouse_data_%s_%s_%s.xlsx', options.group, sess_label, datestr(now,'yyyymmdd'));writetable(export_table, excel_filename);
fprintf('Individual mouse data exported to: %s\n', excel_filename);
disp(export_table);

%% Time between conditions

options = struct();
options.use_limits = true;

options.dist_lim = 50;
options.time_lim =10;  % minutes
options.ylim = [0 100];
options.plot_type = 'percent';
options.sessions = {'sess0', 'sess1'};
options.conditions = {'light_control','dark_control'};
options.show_stats = true;
options.figure_size = [100, 100, 400, 600];
plotConditionComparisonWithStats(mice_all, options);

%% Heatmaps differences

options = struct();
options.diff_clim = [-0.025 0.025];        % Custom color limits ±5%
options.exclude_grooming = true;   % Remove grooming timepoints
options.bin_size = 3;              % Higher resolution
options.smooth_factor = 3;
plot_spatial_difference_heatmaps(mice_all, 'control', slanCM('viola', 64), options);

%% Activity analysis
%% Pearson corr with distance

options = struct();
options.sessions = {'sess0','sess1'};  % Will show as "Before", "Learning", "Test"
options.group = 'control';  % or 'dark_control'
options.memory = 'all';
options.title = 'dF/F-Distance Correlations';
options.time_lim =10;  % minutes
options.dist_lim =155;  % for positions as far as this limit
options.dist_too_close =5;  % limit for how close mouse can be to food
options.remove_grooming = true;  % Remove grooming periods
options.speed_threshold = 0;  % Only include speeds ≥ 0
options.connect_points = true;  % Connect individual mouse data points
options.show_stats = true;  % Show significance indicators
options.ylimit = [-0.55, 0.55];  % Custom y-axis limits
options.colors = {[0.6 0.6 0.6],[0.2 0.4 0.8]};  % Before, Learning, Test
options.font_size = 14;
options.stats_color = [0 0 0];  % Green for significance indicators
options.show_zero_line = false;  % Show horizontal line at r=0
options.use_fisher_z = false;  % Use raw correlations for statistics
options.figure_size = [100, 100, 450, 600];

% % Create the enhanced correlation plot
% plotCorrelationsComparisonWithStats(mice_all, options);

% options.figure_size = [100, 100, 450, 600];
plotCorrelationsComparisonCombinedNonfood(mice_all, options);

%% Compare conditions

options = struct();
options.conditions = {'light_control', 'dark_control'};  % Conditions to compare
options.sessions = {'sess0', 'sess2'};
options.show_stats = true;
options.connect_points = true;
options.speed_threshold = 5;
options.dist_too_close = 0;
options.dist_lim = 150;
options.use_fisher_z = false;  % Use Fisher z for statistics
% options.group_spacing = 0.8;  % Control bar spacing
options.figure_size = [100, 100, 400, 600];
options.ylimit = [-0.3, 0.85];  % Custom y-axis limits

% Run the function
plotCorrelationConditionComparison(mice_all, options);

%% Heamaps of activity difference

options.diff_clim = [-0.7,0.7];
options.bin_size = 7;
options.smooth_factor = 2;  % imgaussfilt smoothing
options.interpolate_display = true;  % Interpolate for display
options.interp_factor =5;  % 4x finer grid
options.speed_threshold = 5;

plot_activity_difference_heatmaps_smoothM3(mice_all, 'light_control', blueWhiteRed, options);

%% Plot an example

mouse_index = 'FDRI32';  % The index of your mouse in mice_all
% session_type = 'sess1';  % 'sess0', 'sess1', or 'sess2'
arm_type = 'food';  % 'food', 'nonfood1', or 'nonfood2'
options.group = 'control';  % Only show control group mice

% Optional custom settings
options = struct();
options.dist_limit =0;  % Minimum distance thresholdval
options.remove_grooming = true;  % Remove grooming periods
options.smooth_window = 1;  % Smoothing window size
options.colormap = blueWhiteRed;  % Colormap for heatmap
options.time_xlim = [100, 250];  % Set x-axis limits for time plots (in seconds)

options.sigma = 1;
options.resolution = 30;

plotMouseMultiSessionData(mice_all, mouse_index, arm_type, options);

%% Unified food run detection (same for 1 source and 3 sources)

food_area = 20;
cfg.validate_plots = 0;
cfg.min_run_extent = food_area;
cfg.paradigm = 'three';
cfg.group_filter = 'control';
cfg.time_lim = 10;  % minutes
run_data_three = detect_food_runs(mice_all, food_area, cfg);
run_data =run_data_three;

%% Slopes for combined nonfood arms

options = struct();
options.sessions = {'sess0', 'sess1', 'sess2'};
options.slope_range = 100; % Test -100 to 0 (towards) and 0 to +100 (away)
options.max_distance = 100;
options.slope_alpha = 0.05;
options.plot_results = true;
options.group_by = 'arm'; % session or arm

slope_results = analyze_combined_nonfood_slopes(run_data, options);

% Visualize results
options.plot_types = {'barplot'}; %{'heatmap', 'scatter'};
visualize_combined_nonfood_slopes(slope_results, options);

%% Separate for sessions Slopes for combined nonfood arms

opts.food_color    = [0.1 0.4 0.8];
opts.nonfood_color = [0.8 0.3 0.1];
opts.plot_results  = false;
opts.slope_range   = 100;
opts.max_distance  = 100;
opts.excel_file    = 'Slopes_by_session_maze3.xlsx';
opts.save_excel    = true;

slope_results = analyze_combined_nonfood_slopes(run_data, opts);
plot_slope_bars_by_session(slope_results, opts);

%% dFF of runs to distance (Gradient)

options = struct();
options.sessions = {'sess0', 'sess1', 'sess2'};  % Analyze all three sessions
options.smoothing = 40;                          % Use 15-point smoothing window
options.y_limits = [-0.7, 0.7];
options.signal_type = 'dff';  % '465' or '405' or 'dff'
options.max_distance = 100;
options.apply_zscore = true;

plot_normalized_dff_across_mice(run_data, options);

%% Nonfood combined gradient

options = struct();
options.sessions = {'sess0', 'sess1', 'sess2'};
options.smoothing = 30;
options.y_limits = [-0.85, 0.7];
options.signal_type = 'dff';
options.max_distance = 100;
options.apply_zscore = true;
options.run_stats = true;
options.figure_width = 470;   % Make it wider
options.figure_height = 800;  % Make it taller
options.title = 'dF/F to distance'; % Your custom title
options.axis_label_font_size = 12;    % Large font for axis labels
options.tick_font_size = 8;          % Small font for tick numbers
options.subplot_title_font_size = 18; % Large subplot titles
options.x_tick_spacing = 25;
options.y_tick_spacing = 0.25;         % Every 0.2 units on y-axis (e.g., -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)

plot_combined_nonfood_dff_across_mice(run_data, options);

%% Compare sessions for each arm Nonfood combined

options = struct();
options.sessions = {'sess0', 'sess1', 'sess2'};
options.smoothing = 35;
options.ylim = [-0.85, 0.7];
options.max_distance = 100;
options.run_stats = true;
options.colors = {[0.6 0.6 0.6], [0.1 0.7 0.9],[0.2 0.4 0.8]}; % control
options.axis_label_font_size = 12;
options.tick_font_size = 8;
options.subplot_title_font_size = 14;
options.main_title_font_size = 14;

% Tick spacing
options.x_tick_spacing = 25;
options.y_tick_spacing = 0.25;

% Figure size
options.figure_width = 470;
options.figure_height = 500;

% Custom title
% options.title = 'Session Comparison: Food vs Non-Food Arms';
options.export_vector = true;
options.export_filename = 'my_figure';
plot_combined_nonfood_compare_sessions(run_data, options);
fig = gcf;

% Wait for complete rendering
pause(0.5);

% Export as vector PDF
export_vector_svg(fig, 'session_comparison_vectorized');


%% Traces by sessions

opts.food_color    = [0.1 0.4 0.8];
opts.nonfood_color = [0.8 0.3 0.1];
opts.smoothing = 30;
opts.max_distance  = 100;
opts.ylim          = [-0.9 0.7];
opts.excel_file    = 'traces_maze3.xlsx';

plot_traces_by_session(run_data, opts);