% Define the folder containing the matrix files
clear all
folderPath = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\M_maze\water_object';
cd(folderPath);  % Change current directory to the folder path

% Get a list of all .mat files in the directory
files = dir('*.mat');

% Initialize a cell array to hold the data
mice_all = {};
% 
% Define mappings for side2 based on varName substrings
keySet = {'MDRI29','MDRE29','FDRE34', 'F30','F29','M28'};
valueSet = {'right','right','right','right','right','right'};
sideMap = containers.Map(keySet, valueSet);

% Loop through each file in the directory
for i = 1:length(files)
    fileName = files(i).name;
    % Load the file
    data = load(fileName);
    % Assume the variable inside the .mat file has the same name as the file
    varName = erase(fileName, '.mat');  % Remove the .mat extension to get the variable name

    % Determine the 'side' based on the name convention (you might need to adjust this logic)
    if contains(varName, 'CNO')
        side = 'CNO';
    elseif contains(varName, 'MCRI')
        side = 'chrimson';
    else
        side = 'control';
    end

    % Determine 'side2' using the map
    side2 = 'unknown';  % Default value if no match is found
    for k = keySet
        if contains(varName, k{1})  % Access k as a string using k{1}
            side2 = sideMap(k{1});  % Use the string directly
            break;  % Stop searching once a match is found
        end
    end

        % Determine the 'side' based on the name convention (you might need to adjust this logic)
    if contains(varName, 'MDRI29')
        side3 = 'bad_signal';
    else
        side3 = 'ok_signal';
    end


    % Append the loaded data and associated info into the cell array
    % Ensure the variable name is valid and exists in 'data'
    if isfield(data, varName)
        mice_all{end+1, 1} = varName; % Name or identifier of the mouse session
        mice_all{end, 2} = side;      % CNO or control side
        mice_all{end, 3} = side2;     % Specific side (center or left)
        mice_all{end, 4} = data.(varName); % Actual data from the file
        mice_all{end, 5} = side3;  % Actual data from the file
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
% 8-zones (1 - 4)
% 9-distance to water
% 10-distance to object
% 11- nothing
% 12-door 0 closed/1 open/2 closed for the second time
% 13-grooming
% 14-climbing or nothing if no climbing

%% Correction for distance - coefficient pix_to_cm

old_pix_to_cm = 75/400;
pix_to_cm = 0.17; % mean for 63.5/385 (horizontal) and 82.55/470 (vertical)
% Loop through each entry in mice_all
for i = 1:size(mice_all, 1)
    data = mice_all{i, 4};
    data(:,9:10) =  data(:,9:10)./old_pix_to_cm.*pix_to_cm;
    mice_all{i, 4} = data;
end

%% Reorganize sessions - split sess1 into sess0 and sess1 based on door condition
% Create a new cell array to hold the reorganized data
mice_all_reorganized = {};

% Loop through each entry in the original data
for i = 1:size(mice_all, 1)
    mouse_id = mice_all{i, 1}(1,1:6); % Mouse session identifier

    side = mice_all{i, 2};         % CNO or control side
    side2 = mice_all{i, 3};        % Specific side (center, left, right)
    data = mice_all{i, 4};         % Actual data from the file
    side3 = mice_all{i, 5};        % Signal quality indicator

    door_col = 12;

    % Get door status changes after discovery
    door_status = data(:, door_col);

    % ---- Find transitions after discovery ----
    transitions_idx = find(diff(door_status) ~= 0);
    
    % Adjust for diff index shift
    transitions = transitions_idx + 1;

    % Extract data where door=0 (closed) for sess0
    before = data(1:transitions(1),:);
    water = data(transitions(1):transitions(2),:);
    object = data(transitions(3):end,:);

    % Create a new entry for sess0 (modify the original session name)
    if ~isempty(before)
        row_idx = size(mice_all_reorganized, 1) + 1;
        mice_all_reorganized{row_idx, 1} = [mouse_id, '_before'];
        mice_all_reorganized{row_idx, 2} = side;
        mice_all_reorganized{row_idx, 3} = side2;
        mice_all_reorganized{row_idx, 4} = before;
        mice_all_reorganized{row_idx, 5} = side3;
    end

    % Keep the original entry for sess1 but with filtered data
    if ~isempty(water)
        row_idx = size(mice_all_reorganized, 1) + 1;
        mice_all_reorganized{row_idx, 1} = [mouse_id, '_water'];
        mice_all_reorganized{row_idx, 2} = side;
        mice_all_reorganized{row_idx, 3} = side2;
        mice_all_reorganized{row_idx, 4} = water;
        mice_all_reorganized{row_idx, 5} = side3;
    end
    if ~isempty(object)
        % For session 2 or any other, just add it as is
        row_idx = size(mice_all_reorganized, 1) + 1;
        mice_all_reorganized{row_idx, 1} = [mouse_id, '_object'];
        mice_all_reorganized{row_idx, 2} = side;
        mice_all_reorganized{row_idx, 3} = side2;
        mice_all_reorganized{row_idx, 4} = object;
        mice_all_reorganized{row_idx, 5} = side3;
    end
end

%
% figure
% hold on
% plot(mice_all_reorganized{1,4}(:,1),mice_all_reorganized{1,4}(:,6))
% plot(mice_all_reorganized{2,4}(:,1),mice_all_reorganized{2,4}(:,6))
% plot(mice_all_reorganized{3,4}(:,1), mice_all_reorganized{3,4}(:,6))

% Replace the original array with the reorganized one
mice_all = mice_all_reorganized;

%% Only mice with good signal

% Find indices of mice with 'ok_signal'
ok_signal_indices = strcmpi(mice_all(:, 5), 'ok_signal'); % Case-insensitive comparison

% Extract only the rows corresponding to 'ok_signal'
mice_good = mice_all(ok_signal_indices, :);

% Save the filtered data
save('mice_good_data.mat', 'mice_good');

% mice_all = mice_good;

%% find interaction with object/water
% when the x,y is closer than 2 cm to the object/water --> interaction
% data(:,15) - water 0(no)/1(yes), data(:,16) - object 0(no)/1(yes)
lim_dist = 3.5;

for i = 1:size(mice_all, 1)
    % Get the data matrix for this mouse/session
    data = mice_all{i, 4};
    mouse_id = mice_all{i, 1};
    
    % Initialize new columns with zeros
    water_interaction = zeros(size(data, 1), 1);
    object_interaction = zeros(size(data, 1), 1);
    
    % Apply interaction rules based on mouse_id
    if contains(mouse_id, 'before')
        % No interactions in "before" sessions (doors are closed)
        water_interaction(:) = 0;
        object_interaction(:) = 0;
    elseif contains(mouse_id, 'water')
        % Only water interactions possible in "water" sessions
        water_interaction(data(:, 9) < lim_dist) = 1; % Column 9 is distance to water
        object_interaction(:) = 0; % No object interactions
    elseif contains(mouse_id, 'object')
        % Only object interactions possible in "object" sessions
        water_interaction(:) = 0; % No water interactions
        object_interaction(data(:, 10) < lim_dist) = 1; % Column 10 is distance to object
    else
        % Default case - shouldn't occur with your naming convention
        water_interaction(data(:, 9) < lim_dist) = 1;
        object_interaction(data(:, 10) < lim_dist) = 1;
    end

    % Add new columns to the data matrix

    data(:, 15) = water_interaction;
    data(:, 16) = object_interaction;

    % Update the data in mice_all
    mice_all{i, 4} = data;

    % % validation - plotting
    % plot_mouse_interactions(data, mouse_id, 'water');
    % plot_mouse_interactions(data, mouse_id, 'object');

    %find "discovery" first transition for interaction
    indx_water=find(data(:,15)>0);
    indx_object=find(data(:,16)>0);
    if ~isempty(indx_water)
        discovery_water = indx_water(1);
    else
        discovery_water = NaN;
    end
    if ~isempty(indx_object)
        discovery_object = indx_object(1);
    else
        discovery_object = NaN;
    end

    mice_all{i, 6} = [discovery_water,discovery_object];

end

%% Raster plot discovery

% Define options (optional)
options = struct();
options.time_window = [-20, 20];
options.caxis_range = [-3, 3];
options.z_score = true;
options.smooth_window = 1;

% Plot discovery-aligned activity for a specific mouse
% plotDiscovery_object_water(mice_all, options);

plotDiscovery_object_water_with_previous(mice_all, options);

% % Plot individual DFF traces for all mice
% plotDffDiscovery_object_water(mice_all);

%% find Interactions 

options = struct();
options.lim_dist = 3.5; % 3cm threshold for interactions
options.time_window = [-20, 20]; % 5 seconds before, 10 seconds after
options.z_score = true;
options.smooth_window = 10;
options.min_duration =5; % Minimum 5 frames to count as an interaction

% Run analysis
plotInteractionAligned(mice_all, options);

% % Validation
% options = struct();
% options.lim_dist = 3.5; % Use the same distance threshold as in your interaction detection
% 
% plotMouseDffWithInteractions(mice_all, options);

%%  Correlations with distance

options = struct();
options.sessions = {'before', 'water'};           % Only compare these sessions
options.group = 'control';                       % Only analyze this group
options.time_lim = 10;                            % time limit
options.dist_lim = 250;                          % dfor positions as far as this limit
options.dist_too_close = 5;                      % limit for how close mouse can be to food
options.filter_order_dist = 1;                   % order filter for distance
options.cutoff_freq_dist = 0.08;                 % Different cutoff for dF/F
options.filter_order_dff = 1;                    % order filter for distance
options.cutoff_freq_dff = 0.05;                  % Different cutoff for dF/F
options.suppress_figures = true;  % Add this line to turn off figures
options.remove_grooming = true;                  % Remove grooming periods
options.speed_threshold = 0;                     % Only include speeds ≥ 0
options.run_stats = true;

corr_data = analyzeDffDistanceCorr_water_object(mice_all, options);
plotBeforeAfterCorr_water_obj(corr_data, options)

%% Heatmaps&examples

options.dist_limit = 0;  % Your interaction threshold from the main script
options.remove_grooming = false;
options.time_before_limits = [150, 500];  
options.time_water_limits = [150, 500];
options.time_object_limits = [150, 500];
options.colormap = blueWhiteRed;
options.sigma = 4;
options.resolution = 35;

% Get unique mouse prefixes from your data
mouse_prefixes = {'F29', 'F30', 'FDRE34', 'M28', 'MDRE29', 'MDRI29'};

% Run analysis for each mouse
for i =1%1:length(mouse_prefixes) 2,5 - ok
    fprintf('Analyzing mouse %s...\n', mouse_prefixes{i});
    try
        plotMouseWaterObjectAnalysis(mice_all, mouse_prefixes{i}, options);
    catch ME
        fprintf('Error analyzing mouse %s: %s\n', mouse_prefixes{i}, ME.message);
    end
end