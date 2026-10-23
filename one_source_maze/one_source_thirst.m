%% Data input and Reorginize data if needed
clear all
% cd('')
% data from all the recorded mice (with both types of food (pellet=food and gel), the format for each mouse from the preprocessing_synch_photometry_individual)

% Get a list of all .mat files in the directory
files = dir('*.mat');

% Initialize a cell array to hold the data
mice_all = {};

% Define mappings for side2 based on varName substrings
% keySet = {'F9I_GLM_thirsty_1','F9I_GLM_thirsty_2','F9I_GLM_thirsty_3','F10I_GLM_thirsty_1','F10I_GLM_thirsty_2','FPhEl6_GLM_thirsty_1'};

keySet = {'F9I_GLM_thirsty_2','F10I_GLM_thirsty_2','FPhEl6_GLM_thirsty_1'}; % 3 unique mice

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
    else
        side2 = 'water';
    end

    if contains(varName, 'fasted')
        side = 'fasted';
    elseif contains(varName, 'fed')
        side = 'fed';
    elseif contains(varName, 'thirsty')
        side = 'thirsty';
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

%% All the data before door closed

mice_all_sess = mice_all;
for i = 1:size(mice_all, 1)
    d = mice_all{i, 4};
    first_open = find(d(:, 7) >= 1, 1);                   % first 1
    if isempty(first_open), continue; end
    reclose = find(d(first_open:end, 7) < 1, 1);           % first 0 after the 1s
    if ~isempty(reclose)
        mice_all_sess{i, 4} = d(1:first_open + reclose - 2, :);   % keep up to last 1
    end
end

figure;plot(mice_all_sess{3, 4}(:,7))

mice_all = mice_all_sess;

%% Find water discovery

n =size(mice_all,1);

for ii = 1:n
    data = mice_all{ii,4};

    %find food discovery
    a1 = data(:,8)>0;
    b1 =find(a1>0);
    discovery = b1(1);
    mice_all{ii,6} = discovery;
end

%% Raster plot water discovery

options = struct();
% options.sort_by = 'response_magnitude'; % Sort mice by peak response
options.state = 'thirsty';
options.source = 'water';
options.smooth_window = 1;            % Apply more smoothing
options.z_score = true;               % Use raw DFF values
options.time_window = [-20, 20]; % 50 seconds before to 200 seconds after discovery
options.caxis_range = [-3, 3];          % Adjust color range for raw values
options.ylim_avg = [-1.5, 2]; 
options.plot_distance = false;
options.bar_window = 15;

plotDiscoveryAlignedData(mice_all, options);

%% Raster plot water visits or drinking without food discovery

options = struct();
options.state = 'thirsty';  % or 'fed' or 'all'
options.source = 'water';   % or 'gel' or 'all'
options.event_type = 'visit'; %or 'eating'
options.smooth_window = 1;
options.ylim_avg = [-1.5, 2]; 
options.caxis_range = [-1.5, 2];          % Adjust color range for raw values
options.time_window = [-20, 20];
options.z_score = true;
options.min_events = 1;
options.bar_window = 15;

plotFoodEventsAligned(mice_all, options);

%% Compare correlations before and after water discovery (with no direction)
% If want to check correlation to the opposite arm - change in the function
% column distance for 14

options = struct();
options.state = {'thirsty'};  % Include both states for analysis
options.source = {'water'};   % Include both sources for analysis
options.dist_limit =5;             % Minimum distance threshold
options.speed_threshold = 5;
options.remove_grooming = true;     % Remove grooming periods
options.use_fisher_z = false;        % Use Fisher z-transformation
options.plot_combined =0;       % Enable combined plotting
options.combined_state = 'thirsty';  % Choose which state to plot (fed or fasted)

compareBeforeAfterCorrelationsEnhanced(mice_all, options);

%% Plot an example

options = struct();
options.dist_limit = 0;             % Minimum distance threshold
options.remove_grooming = false;     % Remove grooming periods
options.remove_eating = false;
options.time_before_limits = [50, 420];  
options.time_after_limits = [0, 370];
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
    if strcmp(mice_all{i, 2}, 'thirsty') && strcmp(mice_all{i, 3}, 'water')
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
    mouse_index = fasted_food_indices(3); % SELECT HERE
    % Plot the df/f vs distance for this mouse
    plotMouseDffVsDistance(mice_all, mouse_index, options);
    plotMouseDffDifferenceHeatmap(mice_all, mouse_index, options); %% Figure 1 heatmaps

else
    disp('No mice found with fasted condition and food stimulus.');
end
