%% Data input and Reorginize data if needed
clear all
% cd('')
% data from all the recorded mice (with both types of food (pellet=food and gel), the format for each mouse from the preprocessing_synch_photometry_individual)

% Get a list of all .mat files in the directory
files = dir('*.mat');

% Initialize a cell array to hold the data
mice_all = {};

% Define mappings for side2 based on varName substrings
keySet = {'F2','F9I_GLM_fasted_food','F9I_GLM_fed_food','F10_GLM_fed_food','F10I_GLM_fasted_food','F11_GLM_fasted_food','F11_GLM_fasted_gel','F11tg_GLM_fed_gel','F12_GLM_fasted_food','F12_GLM_fasted_gel','F12_GLM_fed_gel','F13_GLM_fasted_food','F13_GLM_fasted_gel','F13tg_GLM_fed_gel',...
    'F14_GLM_fasted_food','F14_GLM_fasted_gel','F15_GLM_fasted_food','F15_GLM_fasted_gel','F24_GLM_fasted_food','F24_GLM_fed_food','FDRE14_GLM_fasted_gel','FDRE16_GLM_fasted_gel','M15_GLM_fasted_food','M26_GLM_fasted_food','M26_GLM_fed_food','MDRE27_GLM_fasted_gel','MPhEl1_GLM_fasted_food'};

% Loop through each file in the directory
for i = 1:length(files)
    fileName = files(i).name;
    % Load the file
    data = load(fileName);
    % Assume the variable inside the .mat file has the same name as the file
    varName = erase(fileName, '.mat');  % Remove the .mat extension to get the variable name

    % Determine the 'side' based on the name convention (you might need to adjust this logic)
    if contains(varName, 'gel')
        side2 = 'gel';
    elseif contains(varName, 'food')
        side2 = 'food';
    end

    if contains(varName, 'fasted')
        side = 'fasted';
    elseif contains(varName, 'fed')
        side = 'fed';
    end

        % Determine the 'side' based on the name convention (you might need to adjust this logic)
    if contains(varName, 'MDRE27') || contains(varName, 'F11')
        side3 = 'bad_signal';
    else
        side3 = 'ok_signal';
    end


    % Append the loaded data and associated info into the cell array
    % Ensure the variable name is valid and exists in 'data'
    if isfield(data, varName)
        mice_all{end+1, 1} = varName; % Name or identifier of the mouse session
        mice_all{end, 2} = side;      % 
        mice_all{end, 3} = side2;     % Specific side (center or left)
        mice_all{end, 4} = data.(varName); % Actual data from the file
        mice_all{end, 5} = side3;  % Actual data from the file
    end
end

pix_to_cm = 0.17; 

% Loop through each entry in mice_all to recalculate speed
for i = 1:size(mice_all, 1)
    data = mice_all{i, 4};

    time = data(:, 1);
    x = data(:, 2);
    y = data(:, 3);

    dt = diff(time);          % Time differences
    dx = diff(x);             % X position differences  
    dy = diff(y);             % Y position differences

    % Calculate distance traveled between consecutive points
    distance = sqrt(dx.^2 + dy.^2)*pix_to_cm;

    % Calculate speed (distance/time)
    % Handle division by zero in case of identical timestamps
    speed_calc = zeros(size(distance));
    valid_dt = dt > 0;
    speed_calc(valid_dt) = distance(valid_dt) ./ dt(valid_dt);
    % Smooth the calculated speed (using a simple moving average)
    window_size = min(5, length(speed_calc));
    if window_size > 1
        smoothed_speed = movmean(speed_calc, window_size);
    else
        smoothed_speed = instant_speed;
    end
    
    % Match lengths (append first value at the beginning)
    speed_new = [smoothed_speed(1); smoothed_speed];

    % Replace the speed column (column 7) with recalculated values
    data(:, 4) = speed_new;

    % Update the data back into mice_all
    mice_all{i, 4} = data;
end

speed_threshold = 70; % cm/s

% Loop through each entry in mice_all to filter high speeds
for i = 1:size(mice_all, 1)
    data = mice_all{i, 4};    
    % Extract speed column (column 7)
    speed = data(:, 4);

    % Find indices where speed is <= threshold
    valid_speed_indices = speed <= speed_threshold;

    % Filter data to keep only valid speed points
    data_filtered = data(valid_speed_indices, :);

    % Update the data back into mice_all
    mice_all{i, 4} = data_filtered;

end

% save the assembled cell array for later use
% save('mice_all_data.mat', 'mice_all');

%data (mice_all{end, 4}):
% 1-time
% 2-x
% 3-y
% 4-speed
% 5-dist to food (shortest)
% 6-path from the last visit
% 7-door closed-0/open-1
% 8-food visit
% 9-eating 0/1
% 10-grooming 0/1
% 11-dff

%% Find coordinats that will work for all mice

corners.x1 = 45; corners.y1 = 10;
corners.x2 = 155; corners.y2 = 10;
corners.x3 = 170; corners.y3 = 420;
corners.x4 = 210; corners.y4 = 420;
corners.x5 = 200; corners.y5 = 10;
corners.x6 = 320; corners.y6 = 7;
corners.x7 = 330; corners.y7 = 420;
corners.x8 = 365; corners.y8 = 420;
corners.x9 = 370; corners.y9 = 10;
corners.x10 = 490; corners.y10 = 2;
corners.x11 = 490; corners.y11 = 510;
corners.x12 = 50; corners.y12 = 508;

% Plot all mice
% plotAllMiceInMaze(mice_all, corners);

mice_all = defineMazeZones(mice_all, corners);

close all

%% Clean the tracking artefacts 

mice_all = remove_maze_tracking_artifacts(mice_all, corners, struct('example', true));

%% Only mice with good signal

% Find indices of mice with 'ok_signal'
ok_signal_indices = strcmpi(mice_all(:, 5), 'ok_signal'); % Case-insensitive comparison

% Extract only the rows corresponding to 'ok_signal'
mice_good = mice_all(ok_signal_indices, :);

% Save the filtered data
% save('mice_good_data.mat', 'mice_good');

mice_all = mice_good;

%% Find food discovery

n =size(mice_all,1);

for ii = 1:n
    data = mice_all{ii,4};

    %find food discovery
    a1 = data(:,8)>0;
    b1 =find(a1>0);
    a2 = data(:,9)>0;
    b2 =find(a2>0);
    if ii == 26 || b2(1)<=b1(1) % for fasted smelly ADD if ii == 30
        discovery = b2(1);
    elseif b1(1)<b2(1)
        discovery = b1(1);
    end
    mice_all{ii,6} = discovery;
end

%% Add acceleration - column 15 to the data

opts.smooth_win_sec = 2;
opts.example = 0;
opts.n_example = 3;
opts.example_rows = [1,10,13];

mice_all = add_acceleration_column(mice_all, opts);

%% Add head direction - an angle beween head direction vector and the vector to food - column 16 to the data

opts = struct;
opts.example = true;
opts.n_example = 0;
opts.n_vectors =50;

mice_all = add_head_direction_column(mice_all, opts);

%% Add angle beween head direction vector and the vector to the closest path - column 17 to the data

opts = struct;
opts.example = true;
opts.n_example = 0;
opts.n_vectors =50;

mice_all = add_path_direction_column(mice_all, opts);

%% Add elapsed time from food discovery - column 18 to the data

opts = struct;
opts.n_examples = 0;
[mice_all, info] = add_elapsed_distance_columns(mice_all, opts);

%% Raster plot food discovery

options = struct();
% options.sort_by = 'response_magnitude'; % Sort mice by peak response
options.state = 'fasted';
options.source = 'food';
options.smooth_window = 1;            % Apply more smoothing
options.z_score = true;               % Use raw DFF values
options.time_window = [-20, 20]; % 50 seconds before to 200 seconds after discovery
options.caxis_range = [-3, 3];          % Adjust color range for raw values
options.ylim_avg = [-1.5, 2]; 
options.plot_distance = false;
options.bar_window = 15;

plotDiscoveryAlignedData(mice_all, options);

%% Comparing fed and fasted - discovery

options.source = 'food';
options.window = 15;
compareDiscoveryResponse(mice_all,options);

%% plot eating/food visits

options = struct();
options.mouse_index = 1; %6, 19,24
options.xlim = [800 1500];  % Optional: zoom to first 500 seconds
options.ylim = [-10 15];

plotActivityWithBehavior(mice_all, options);

%% Raster plot food visits or eating without food discovery

options = struct();
options.state = 'fasted';  % or 'fed' or 'all'
options.source = 'food';   % or 'gel' or 'all'
options.event_type = 'visit'; %or 'eating'
options.smooth_window = 1;
options.ylim_avg = [-1.5, 2]; 
options.caxis_range = [-1.5, 2];          % Adjust color range for raw values
options.time_window = [-20, 20];
options.z_score = true;
options.min_events = 1;
options.bar_window = 15;

plotFoodEventsAligned(mice_all, options);

%% Comparing fed and fasted - events

options.source = 'food';
options.event_type = 'visit';
options.window = 15;
compareEventResponse(mice_all,options);  

%% Compare correlations before and after food discovery (with no direction)
% If want to check correlation to the opposite arm - change in the function
% column distance for 14

options = struct();
options.state = {'fasted'};  % Include both states for analysis
options.source = {'food'};   % Include both sources for analysis
options.dist_limit =5;             % Minimum distance threshold
options.speed_threshold = 5;
options.remove_grooming = true;     % Remove grooming periods
options.use_fisher_z = false;        % Use Fisher z-transformation
options.plot_combined = true;       % Enable combined plotting
options.combined_state = 'fasted';  % Choose which state to plot (fed or fasted)

compareBeforeAfterCorrelationsEnhanced(mice_all, options);

%% Validate distance limit

validation_options.dist_limit = 5;
validation_options.state = 'fasted';
validation_options.source = 'food';
validateDistanceLimit(mice_all, validation_options);

%% Plot an example

options = struct();
options.dist_limit = 0;             % Minimum distance threshold
options.remove_grooming = false;     % Remove grooming periods
options.remove_eating = false;
options.time_before_limits = [398.6, 768.6];  
options.time_after_limits = [398.6, 768.6];
options.colormap = blueWhiteRed;
% options.colormap_limits = [-1, 1];
options.sigma =4;
options.resolution = 25;
options.difference_limits = [-0.8,0.8];

options.bin_size = 4;  % 4 cm bins
options.smooth_factor = 2;  % imgaussfilt smoothing
options.interpolate_display = 1;
options.interp_factor = 2;  % 4x finer display

fasted_food_indices = [];
for i = 1:size(mice_all, 1)
    if strcmp(mice_all{i, 2}, 'fasted') && strcmp(mice_all{i, 3}, 'gel')
        fasted_food_indices = [fasted_food_indices, i];
        fprintf('Mouse %d: ID=%s, Condition=%s, Stimulus=%s\n', ...
                i, mice_all{i, 1}, mice_all{i, 2}, mice_all{i, 3});
    end
end
% options.dist_limit = 5;
% options.remove_grooming = true;

%Adjust smoothness of color transitions
options.contour_levels = 30; % More levels = very smooth

% Standard resolution works fine with contourf
% options.resolution = 35;  % Default, works great
% options.sigma = 1;         % Standard smoothing
% Select a specific mouse
if ~isempty(fasted_food_indices)
    mouse_index = fasted_food_indices(4); % SELECT HERE
    % Plot the df/f vs distance for this mouse
    plotMouseDffVsDistance(mice_all, mouse_index, options);
    plotMouseDffDifferenceHeatmap(mice_all, mouse_index, options); %% Figure 1 heatmaps

else
    disp('No mice found with fasted condition and food stimulus.');
end

%% Reorganize sessions - split sess1 into sess0 and sess1 based on door condition

%view specific mice
plotSessionsPreview(mice_all, [1, 3, 5]);
%
mice_all_reorganized = reorganizeMiceData(mice_all);
% save('mice_all_reorganized.mat', 'mice_all_reorganized');
mice_all = mice_all_reorganized;

%% Unified run detection (same for 1 source and 3 sources)

food_area = 20;
cfg1.validate_plots = 0;
cfg1.paradigm      = 'single';
cfg1.group_filter  = 'fasted';
cfg1.source_filter = 'food';
cfg2.validate_plots = 0;
cfg2.paradigm      = 'single';
cfg2.group_filter  = 'fasted';
cfg2.source_filter = 'gel';
runs_fasted_food = detect_food_runs(mice_all, food_area, cfg1);
runs_fasted_gel = detect_food_runs(mice_all, food_area, cfg2);

%% Correlations bar plots

options = struct();
options.sessions = {'sess0', 'sess1'};
options.run_category = 'food'; %  'not_food' or 'food'
options.separate_plots = true;  % or omit this line (default is false)
options.use_fisher_z = false;  % Use Fisher z-transformation
plot_correlation_comparison(runs_fasted_food, options);

%% dFF of runs to distance (Gradient)
% Create options struct
options = struct();
options.sessions = {'sess0', 'sess1'};  % Select sessions 0 and 1
% options.conditions = {'fasted_food','fasted_gel' };  % Compare fasted vs fed for food
options.smoothing = 20;  % Smoothing window size
options.xlim = [0,85];
% options.xlim_away = [0,100];
options.ylim = [-0.65 0.65];
options.figure_size = [500 600];
options.run_category = 'food';

% Generate the plot
% plot_dff_distance_single_source(runs_fed_food, runs_fasted_food, runs_fed_gel, runs_fasted_gel, options);
plot_dff_distance_fasted_only(runs_fasted_food, runs_fasted_gel, options)
print(gcf, '-dpdf', 'traces80cm.pdf', '-painters');

%% Analysis of slopes

options = struct();
options.slope_range = 210;        % Distance range for testing (±50 units)
options.max_distance = 210;     % Maximum distance to analyze
options.slope_alpha = 0.05;     % Significance level
options.bin_width =0.05;
options.plot_results = true;    % Generate plots
slope_results = analyze_single_arm_slopes(runs_fasted_food, options);
options.separate_plots = true;
options.plot_types = {'barplot'};
% options.session_colors = {[0.2, 0.4, 0.8], [0.8, 0.2, 0.2]}; % Blue, Red
options.save_plots = 1;
visualize_single_arm_slopes(slope_results, options);

%% Runs after food visit vs after eating

food_area = 20;
validate_plots = 0;
max_time_gap = 1;
max_spatial_dist = 200;
group_filter = 'fasted';
source_filter = 'food';
lookback_window = 15;

runs_fasted_food_beh = analyze_runs_with_behavior_classification(...
    mice_all, food_area, group_filter, source_filter, validate_plots, max_time_gap, max_spatial_dist, lookback_window);

% Plotting Runs after food visit vs after eating

options = struct();
options.sessions = {'sess1'};
options.lookback_window = 5;  % seconds
options.smoothing = 25;
options.xlim_towards = [0 200];
options.xlim_away = [0 200];
options.ylim = [-1 1];
% options.save_dir = '';

plot_runs_by_behavior(runs_fasted_food_beh, options);

%% Slopes (for the behavior - visit vs eating)

slope_opts = struct();
opts.sessions = {'sess1'};
opts.directions = {'towards', 'away'};
opts.behaviors = {'eating', 'visit'};
opts.min_runs_per_mouse = 2;
opts.slope_range = 150;
% opts.save_dir = '';

behavior_slopes = analyze_behavior_slopes(runs_fasted_food_beh, opts);

opts.plot_types = {'barplot'};
opts.separate_directions = true;

visualize_behavior_slopes(behavior_slopes, opts);

%% Averaged both types of runs dynamic

options = struct();
options.max_distance = 200;    % Set maximum distance to plot
options.min_mice = 3;         
options.max_runs = 20;
% options.ylim_dff = [-0.5, ]; % Set custom y-axis limits
options.save_data   = true;
options.save_path   = '';   % base name, no extension (optional)
options.save_format = 'xlsx';                          % 'mat' (default), 'csv', or 'both'
res = plot_dff_vs_distance_session_comparison(runs_fasted_food, options);

%% Plot activity vs head angle

opts = struct;
opts.exclude_eating = 1;
opts.exclude_grooming = 1;
opts.n_bins = 36;
opts.sessions = 1;
opts.state = 'fasted';
opts.source = 'food';
opts.dist_threshold = 80;
opts.speed_threshold = 0;
opts.dist_remove = 0;
opts.weight_by_occupancy = 1;
opts.sector_ylim = [-1.3,2];
opts.angle_col = 16; % 17 - to shortest path, 16 - absolute to food

plot_activity_polar(mice_all, opts);

%% GLM per mouse models comparision

opts.state         = 'fasted';
opts.source        = 'food';
opts.baseline_correct = 'sliding_pct';
opts.baseline_win_sec  = 60;    % 60s window (default)
opts.baseline_pct      = 8;     % 8th percentile (default)
opts.cv_type   = 'chunked'; 
opts.shuf_type = 'circular';   % roll dF/F by random offset
opts.n_cv_iter     = 100;      % CV iterations (default 100)
opts.cv_train_frac = 0.80;     % 80/20 split (default)
opts.excl_abs_time = true;
opts.tsf_type = 'temporal_distance'; % 'time_since_food'/'time_to_food'/'temporal_distance'
opts.headdirection    = 0;
glm_results = glm_dff_model_comparison(mice_all, opts);

%% Plotting behavioral traces 

opts = struct();
opts.state = 'fasted';
opts.source = 'food';
opts.session = 1;
opts.save_dir = 'behavior_trace_examples';

plot_behavior_traces_per_mouse(mice_all, opts);