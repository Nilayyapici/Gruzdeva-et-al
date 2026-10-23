function validate_path_lengths(baseDir, mouseNames, options)
% Validate path length calculations by plotting each mouse individually
%
% Usage:
% mouseNames = {'MDRE8', 'F13', 'M21', 'FDRE14'};
% baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';
% options.num_trials = 10;
% options.food_discovery_threshold = 1.2;
% validate_path_lengths(baseDir, mouseNames, options);

if nargin < 3
    options = struct();
end

if ~isfield(options, 'num_trials'), options.num_trials = 10; end
if ~isfield(options, 'food_discovery_threshold'), options.food_discovery_threshold = 1.2; end

% Load all mouse data
fprintf('Loading data for validation...\n');
allMiceData = {};
for mouseIdx = 1:length(mouseNames)
    mouseName = mouseNames{mouseIdx};
    mouseData = loadMouseData(baseDir, mouseName);
    if ~isempty(mouseData)
        allMiceData{end+1} = mouseData;
        fprintf('  Loaded %s\n', mouseName);
    end
end

if isempty(allMiceData)
    error('No mouse data found!');
end

% Create figure with subplots for each mouse
nMice = length(allMiceData);
figure('Name', 'Path Length Validation - Individual Mice', ...
    'Position', [50, 50, 400*min(nMice, 2), 400*ceil(nMice/2)]);

% Process each mouse
for mouseIdx = 1:nMice
    subplot(ceil(nMice/2), min(nMice, 2), mouseIdx);
    
    mouseData = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    
    % Storage for this mouse - initialize with NaN for all trials
    pathLengths = NaN(1, options.num_trials);
    timesToFood = NaN(1, options.num_trials);
    numPoints = NaN(1, options.num_trials);
    
    % Process each trial
    for trialIdx = 1:size(learningData, 1)
        trialName = learningData{trialIdx, 1};
        trialData = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};
        
        if ~isempty(trialData) && ~isempty(foodCoords)
            % Extract trial number
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                
                if trialNum >= 1 && trialNum <= options.num_trials
                    % Calculate path length
                    [timeToFood, pathLength, ~, nPoints] = calculatePathToFood(trialData, foodCoords, options.food_discovery_threshold);
                    
                    if ~isnan(pathLength)
                        pathLengths(trialNum) = pathLength;
                        timesToFood(trialNum) = timeToFood;
                        numPoints(trialNum) = nPoints;
                    end
                end
            end
        end
    end
    
    % Count valid trials
    validTrials = ~isnan(pathLengths);
    nValidTrials = sum(validTrials);
    
    % Plot path lengths - this will now connect trials 1-2-3-4 in order
    hold on;
    plot(1:options.num_trials, pathLengths, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.3 0.6 0.9], 'Color', [0.3 0.6 0.9]);
    
    % Add labels and formatting
    xlabel('Trial Number');
    ylabel('Path Length to Food (cm)');
    title(sprintf('%s (n=%d trials)', mouseData.name, nValidTrials));
    grid on;
    xlim([0.5, options.num_trials + 0.5]);
    
    % Print statistics for this mouse
    if nValidTrials > 0
        validPathLengths = pathLengths(validTrials);
        validTimes = timesToFood(validTrials);
        validPoints = numPoints(validTrials);
        
        fprintf('\n%s:\n', mouseData.name);
        fprintf('  Mean path length: %.1f ± %.1f cm\n', mean(validPathLengths), std(validPathLengths));
        fprintf('  Range: %.1f - %.1f cm\n', min(validPathLengths), max(validPathLengths));
        fprintf('  Mean time to food: %.1f ± %.1f s\n', mean(validTimes), std(validTimes));
        fprintf('  Mean points per trial: %.1f\n', mean(validPoints));
        
        % Print individual trial details
        fprintf('  Trial details:\n');
        for trialNum = 1:options.num_trials
            if validTrials(trialNum)
                fprintf('    Trial %d: %.1f cm, %.1f s, %d points\n', ...
                    trialNum, pathLengths(trialNum), timesToFood(trialNum), numPoints(trialNum));
            else
                fprintf('    Trial %d: No data\n', trialNum);
            end
        end
    end
    
    hold off;
end

sgtitle('Path Length to Food Discovery - Individual Mouse Validation', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Overall statistics
fprintf('\n=== OVERALL STATISTICS ===\n');
allPathLengths = [];
allTimes = [];
for mouseIdx = 1:nMice
    mouseData = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    
    for trialIdx = 1:size(learningData, 1)
        trialName = learningData{trialIdx, 1};
        trialData = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};
        
        if ~isempty(trialData) && ~isempty(foodCoords)
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                if trialNum >= 1 && trialNum <= options.num_trials
                    [timeToFood, pathLength, ~, ~] = calculatePathToFood(trialData, foodCoords, options.food_discovery_threshold);
                    if ~isnan(pathLength)
                        allPathLengths(end+1) = pathLength;
                        allTimes(end+1) = timeToFood;
                    end
                end
            end
        end
    end
end

if ~isempty(allPathLengths)
    fprintf('All mice combined:\n');
    fprintf('  Mean path length: %.1f ± %.1f cm\n', mean(allPathLengths), std(allPathLengths));
    fprintf('  Median path length: %.1f cm\n', median(allPathLengths));
    fprintf('  Range: %.1f - %.1f cm\n', min(allPathLengths), max(allPathLengths));
    fprintf('  Mean time to food: %.1f ± %.1f s\n', mean(allTimes), std(allTimes));
    fprintf('  Total valid trials: %d\n', length(allPathLengths));
end
fprintf('==========================\n\n');

end

function mouseData = loadMouseData(baseDir, mouseName)
% Load learning data for a single mouse from saline condition only
mouseData = [];
mouseFolders = dir(fullfile(baseDir, sprintf('%s*', mouseName)));

for i = 1:length(mouseFolders)
    if mouseFolders(i).isdir
        mouseFolder = fullfile(baseDir, mouseFolders(i).name);
        salineFolders = dir(fullfile(mouseFolder, '*saline*'));
        salineFolders = salineFolders([salineFolders.isdir]);

        for j = 1:length(salineFolders)
            salineFolder = fullfile(mouseFolder, salineFolders(j).name);
            learningFile = fullfile(salineFolder, sprintf('%s_learning.mat', mouseName));

            if exist(learningFile, 'file')
                try
                    data = load(learningFile);
                    learningData = eval(sprintf('data.%s_learning', mouseName));
                    mouseData.name = mouseName;
                    mouseData.folder = salineFolder;
                    mouseData.learningData = learningData;
                    mouseData.condition = salineFolders(j).name;
                    return;
                catch ME
                    fprintf('    Error loading %s: %s\n', learningFile, ME.message);
                end
            end
        end
    end
end
end

function [timeToFood, pathLength, meanSpeed, nPoints] = calculatePathToFood(trialData, foodCoords, foodDiscoveryThreshold)
% Calculate path metrics until food discovery with detailed output

timeToFood = NaN;
pathLength = NaN;
meanSpeed = NaN;
nPoints = 0;

if isempty(trialData) || isempty(foodCoords)
    return;
end

% Extract data
x_coords = trialData(:, 2);
y_coords = trialData(:, 3);
time_coords = trialData(:, 1);
speed_values = trialData(:, 7);

% Remove NaN values
valid_idx = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(time_coords);
x_clean = x_coords(valid_idx);
y_clean = y_coords(valid_idx);
time_clean = time_coords(valid_idx);
speed_clean = speed_values(valid_idx);

if length(x_clean) < 2 || length(foodCoords) < 2
    return;
end

% Calculate distances to food and find discovery point
pix_to_cm = 75/400;
distances_to_food = sqrt((x_clean - foodCoords(1)).^2 + (y_clean - foodCoords(2)).^2) * pix_to_cm;
foodDiscoveryIdx = find(distances_to_food < foodDiscoveryThreshold, 1, 'first');

% If food was discovered, calculate metrics up to that point
if ~isempty(foodDiscoveryIdx)
    % Time to food discovery
    timeToFood = time_clean(foodDiscoveryIdx);
    
    % Number of points until discovery
    nPoints = foodDiscoveryIdx;
    
    % Path length until food discovery
    x_to_food = x_clean(1:foodDiscoveryIdx);
    y_to_food = y_clean(1:foodDiscoveryIdx);
    dx = diff(x_to_food);
    dy = diff(y_to_food);
    distances = sqrt(dx.^2 + dy.^2);
    pathLength = sum(distances) * pix_to_cm;
    
    % Mean speed until food discovery
    speed_to_food = speed_clean(1:foodDiscoveryIdx);
    valid_speed = speed_to_food(~isnan(speed_to_food));
    if ~isempty(valid_speed)
        meanSpeed = mean(valid_speed);
    end
end
end