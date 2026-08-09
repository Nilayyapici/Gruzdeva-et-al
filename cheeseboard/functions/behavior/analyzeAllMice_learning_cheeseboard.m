function allMiceData = analyzeAllMice_learning_cheeseboard(baseDir, mouseNames, options)
% Comprehensive analysis of all mice learning data with customizable trial number
%
% Usage:
% mouseNames = {'MDRE8', 'F13', 'M21', 'FDRE14'};
% baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';
%
% options.use_zscore = true;           % Use z-scored dF/F (default: false)
% options.distance_limit = 15;        % Maximum distance to analyze in cm (default: 20)
% options.bin_size = 0.5;             % Bin size for distance binning in cm (default: 0.5)
% options.dff_ylim = [-2, 2];         % Y-axis limits for dF/F plots (default: auto)
% options.food_discovery_threshold = 1.2; % Distance threshold for food discovery in cm (default: 1.2)
% options.num_trials = 10;            % Number of trials to analyze (default: 10, max: 12)
% options.smoothing_method = 'movmean'; % Smoothing method: 'movmean', 'sgolay', or 'none' (default: 'movmean')
% options.smoothing_window = 3;       % Smoothing window size (default: 3)
% options.sgolay_order = 2;           % Savitzky-Golay polynomial order (default: 2, only used if smoothing_method = 'sgolay')
% options.early_trials = [1, 2, 3];   % Trials to include in early phase (default: [1, 2, 3])
% options.middle_trials = [4, 5, 6, 7]; % Trials to include in middle phase (default: [4, 5, 6, 7])
% options.late_trials = [8, 9, 10];   % Trials to include in late phase (default: [8, 9, 10])
% options.plot_phases = {'early', 'middle', 'late'}; % Which phases to plot (default: all three)
% options.save_xlsx = true;           % Save analysis data to Excel file (default: true)
% options.xlsx_path = '';             % Path for xlsx output (default: baseDir/analysis_results.xlsx)
% options.speed_threshold = 50;       % Exclude points where speed > threshold in cm/s (default: 50)
% options.min_speed_threshold = 0;    % Exclude points where speed < threshold in cm/s (default: 0, disabled)
% options.exclude_first_n_frames = 0; % Exclude first N frames from each trial (default: 0)
%
% analyzeAllMice_learning_cheeseboard(baseDir, mouseNames, options);

if nargin < 2
    mouseNames = {'MDRE8', 'F13', 'M21', 'FDRE14'};
    baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';
end

if nargin < 3
    options = struct();
end

% Set default options
if ~isfield(options, 'use_zscore'), options.use_zscore = false; end
if ~isfield(options, 'distance_limit'), options.distance_limit = 20; end
if ~isfield(options, 'bin_size'), options.bin_size = 0.5; end
if ~isfield(options, 'dff_ylim'), options.dff_ylim = []; end
if ~isfield(options, 'food_discovery_threshold'), options.food_discovery_threshold = 1.2; end
if ~isfield(options, 'num_trials'), options.num_trials = 10; end
if ~isfield(options, 'smoothing_method'), options.smoothing_method = 'movmean'; end
if ~isfield(options, 'smoothing_window'), options.smoothing_window = 3; end
if ~isfield(options, 'sgolay_order'), options.sgolay_order = 2; end
if ~isfield(options, 'early_trials'), options.early_trials = [1, 2, 3]; end
if ~isfield(options, 'middle_trials'), options.middle_trials = [4, 5, 6, 7]; end
if ~isfield(options, 'late_trials'), options.late_trials = [8, 9, 10]; end
if ~isfield(options, 'plot_phases'), options.plot_phases = {'early', 'middle', 'late'}; end
if ~isfield(options, 'save_xlsx'), options.save_xlsx = true; end
if ~isfield(options, 'xlsx_path'), options.xlsx_path = ''; end
if ~isfield(options, 'speed_threshold'), options.speed_threshold = 50; end
if ~isfield(options, 'min_speed_threshold'), options.min_speed_threshold = 0; end
if ~isfield(options, 'exclude_first_n_frames'), options.exclude_first_n_frames = 0; end

% Validate options
if options.num_trials < 1 || options.num_trials > 12
    error('Number of trials must be between 1 and 12');
end
if ~ismember(options.smoothing_method, {'movmean', 'sgolay', 'none'})
    error('Smoothing method must be "movmean", "sgolay", or "none"');
end
if options.exclude_first_n_frames < 0
    error('exclude_first_n_frames must be >= 0');
end
valid_phases = {'early', 'middle', 'late'};
if ~iscell(options.plot_phases)
    options.plot_phases = {options.plot_phases};
end
for i = 1:length(options.plot_phases)
    if ~ismember(options.plot_phases{i}, valid_phases)
        error('plot_phases must contain only "early", "middle", or "late"');
    end
end

fprintf('Analysis options:\n');
fprintf('  Z-scored dF/F: %s\n', mat2str(options.use_zscore));
fprintf('  Distance limit: %.1f cm\n', options.distance_limit);
fprintf('  Bin size: %.1f cm\n', options.bin_size);
fprintf('  Food discovery threshold: %.1f cm\n', options.food_discovery_threshold);
fprintf('  Number of trials to analyze: %d\n', options.num_trials);
fprintf('  Speed threshold: %.1f cm/s\n', options.speed_threshold);
fprintf('  Min speed threshold: %.1f cm/s\n', options.min_speed_threshold);
fprintf('  Exclude first N frames: %d\n', options.exclude_first_n_frames);
fprintf('  Smoothing method (dF/F vs distance only): %s\n', options.smoothing_method);
if ~strcmp(options.smoothing_method, 'none')
    fprintf('  Smoothing window: %d\n', options.smoothing_window);
    if strcmp(options.smoothing_method, 'sgolay')
        fprintf('  Savitzky-Golay order: %d\n', options.sgolay_order);
    end
end
fprintf('  Early trials: [%s]\n', num2str(options.early_trials));
fprintf('  Middle trials: [%s]\n', num2str(options.middle_trials));
fprintf('  Late trials: [%s]\n', num2str(options.late_trials));
fprintf('  Phases to plot: %s\n', strjoin(options.plot_phases, ', '));
if ~isempty(options.dff_ylim)
    fprintf('  dF/F Y-limits: [%.1f, %.1f]\n', options.dff_ylim(1), options.dff_ylim(2));
else
    fprintf('  dF/F Y-limits: Auto\n');
end

% Initialize data storage
allMiceData = {};

% Load all mouse data
fprintf('Loading data for %d mice (saline condition only)...\n', length(mouseNames));
for mouseIdx = 1:length(mouseNames)
    mouseName = mouseNames{mouseIdx};
    fprintf('  Loading %s... \n', mouseName);
    try
        mouseData = loadMouseData(baseDir, mouseName);
        if ~isempty(mouseData)
            allMiceData{end+1} = mouseData;
            fprintf('    ✓ Success: %d trials from %s\n', size(mouseData.learningData, 1), mouseData.condition);
        else
            fprintf('    ✗ No saline learning data found\n');
        end
    catch ME
        fprintf('    ✗ Error: %s\n', ME.message);
    end
end

if isempty(allMiceData)
    error('No mouse data found!');
end

fprintf('\nSuccessfully loaded %d mice\n\n', length(allMiceData));

% 1. Learning Curves Analysis
lcData = plotLearningCurves(allMiceData, options);

% 2. dFF vs Distance Analysis - Individual Trials
[trialDFFData, distanceBinCenters] = analyzeDFFvsDistance_IndividualTrials(allMiceData, options);

% 3. dFF vs Distance Analysis - Trial Phases (Early, Middle, Late)
phaseDFFData = analyzeDFFvsDistance_TrialPhases(allMiceData, options);

% 4. Summary Statistics
summaryData = printSummaryStats(allMiceData, options);

% 5. Save to xlsx if requested
if options.save_xlsx
    if isempty(options.xlsx_path)
        xlsxFile = fullfile(baseDir, sprintf('analysis_results_%s.xlsx', datestr(now, 'yyyymmdd_HHMMSS')));
    else
        xlsxFile = options.xlsx_path;
    end
    saveAnalysisToXLSX(xlsxFile, allMiceData, lcData, trialDFFData, phaseDFFData, ...
        summaryData, distanceBinCenters, options);
end

fprintf('\nReturning allMiceData for further analysis...\n');
end

% =========================================================================
%  SHARED FRAME/SPEED FILTER HELPER
% =========================================================================

function data = applyFrameAndSpeedFilter(data, options)
% Exclude first N frames, then remove rows where speed > max threshold
% or speed < min threshold. Rows with NaN speed are kept.
n = options.exclude_first_n_frames;
if n > 0 && size(data, 1) > n
    data = data(n+1:end, :);
end
speed_vals = data(:, 7);
keep = isnan(speed_vals) ...
     | (speed_vals <= options.speed_threshold & speed_vals >= options.min_speed_threshold);
data = data(keep, :);
end

% =========================================================================
%  XLSX EXPORT
% =========================================================================

function saveAnalysisToXLSX(xlsxFile, allMiceData, lcData, trialDFFData, phaseDFFData, ...
    summaryData, binCenters, options)

fprintf('\nSaving results to: %s\n', xlsxFile);

maxTrials = options.num_trials;
nBins     = length(binCenters);

% --- Sheet 1: Learning Curves -----------------------------------------
mouseNames  = cellfun(@(m) m.name, allMiceData, 'UniformOutput', false);
trialLabels = arrayfun(@(t) sprintf('Trial_%d', t), 1:maxTrials, 'UniformOutput', false);

ttf_header = [{'Mouse'}, strcat('TimeToFood_', trialLabels)];
ttf_data   = [mouseNames', num2cell(lcData.timeToFood)];
pl_header  = [{'Mouse'}, strcat('PathLength_cm_', trialLabels)];
pl_data    = [mouseNames', num2cell(lcData.pathLength)];
sp_header  = [{'Mouse'}, strcat('MeanSpeed_cmps_', trialLabels)];
sp_data    = [mouseNames', num2cell(lcData.speed)];

mean_ttf = nanmean(lcData.timeToFood, 1);
sem_ttf  = nanstd(lcData.timeToFood, 0, 1) ./ sqrt(sum(~isnan(lcData.timeToFood), 1));
mean_pl  = nanmean(lcData.pathLength, 1);
sem_pl   = nanstd(lcData.pathLength, 0, 1) ./ sqrt(sum(~isnan(lcData.pathLength), 1));
mean_sp  = nanmean(lcData.speed, 1);
sem_sp   = nanstd(lcData.speed, 0, 1) ./ sqrt(sum(~isnan(lcData.speed), 1));

ttf_data = [ttf_data; [{'Mean'}, num2cell(mean_ttf)]; [{'SEM'}, num2cell(sem_ttf)]];
pl_data  = [pl_data;  [{'Mean'}, num2cell(mean_pl)];  [{'SEM'}, num2cell(sem_pl)]];
sp_data  = [sp_data;  [{'Mean'}, num2cell(mean_sp)];  [{'SEM'}, num2cell(sem_sp)]];

blank = cell(1, maxTrials + 1);
lc_table = [ttf_header; ttf_data; blank; ...
            {'--- Path Length (cm) ---'}, blank(2:end); ...
            pl_header; pl_data; blank; ...
            {'--- Mean Speed (cm/s) ---'}, blank(2:end); ...
            sp_header; sp_data];

writecell(lc_table, xlsxFile, 'Sheet', 'Learning_Curves');
fprintf('  Sheet written: Learning_Curves\n');

% --- Sheet 2: DFF vs Distance - Individual Trials ---------------------
trial_headers = {};
for t = 1:maxTrials
    trial_headers{end+1} = sprintf('Trial%d_Mean', t);
    trial_headers{end+1} = sprintf('Trial%d_SEM',  t);
    trial_headers{end+1} = sprintf('Trial%d_nDataPts', t);
end
dff_trial_header = [{'Distance_cm'}, trial_headers];
dff_trial_data   = num2cell(binCenters');
for t = 1:maxTrials
    if ~isempty(trialDFFData{t}.meanDFF)
        dff_trial_data = [dff_trial_data, ...
            num2cell(trialDFFData{t}.meanDFF'), ...
            num2cell(trialDFFData{t}.semDFF'), ...
            num2cell(trialDFFData{t}.counts')];
    else
        nanCol  = num2cell(NaN(nBins, 1));
        zeroCol = num2cell(zeros(nBins, 1));
        dff_trial_data = [dff_trial_data, nanCol, nanCol, zeroCol];
    end
end
writecell([dff_trial_header; dff_trial_data], xlsxFile, 'Sheet', 'DFF_vs_Distance_Trials');
fprintf('  Sheet written: DFF_vs_Distance_Trials\n');

% --- Sheet 3: DFF vs Distance - Phases --------------------------------
phaseNames = fieldnames(phaseDFFData);
phase_headers = {};
for i = 1:length(phaseNames)
    pn = phaseNames{i};
    phase_headers{end+1} = sprintf('%s_Mean', pn);
    phase_headers{end+1} = sprintf('%s_SEM',  pn);
    phase_headers{end+1} = sprintf('%s_nDataPts', pn);
end
dff_phase_header = [{'Distance_cm'}, phase_headers];
dff_phase_data   = num2cell(binCenters');
for i = 1:length(phaseNames)
    pn = phaseNames{i};
    if ~isempty(phaseDFFData.(pn).meanDFF)
        dff_phase_data = [dff_phase_data, ...
            num2cell(phaseDFFData.(pn).meanDFF'), ...
            num2cell(phaseDFFData.(pn).semDFF'), ...
            num2cell(phaseDFFData.(pn).counts')];
    else
        nanCol  = num2cell(NaN(nBins, 1));
        zeroCol = num2cell(zeros(nBins, 1));
        dff_phase_data = [dff_phase_data, nanCol, nanCol, zeroCol];
    end
end
writecell([dff_phase_header; dff_phase_data], xlsxFile, 'Sheet', 'DFF_vs_Distance_Phases');
fprintf('  Sheet written: DFF_vs_Distance_Phases\n');

% --- Sheet 4: Summary Stats -------------------------------------------
sum_header = {'Mouse', 'Condition', 'TrialsAttempted', 'TrialsSuccessful', 'SuccessRate_pct'};
writecell([sum_header; summaryData.rows], xlsxFile, 'Sheet', 'Summary_Stats');
fprintf('  Sheet written: Summary_Stats\n');

fprintf('Excel file saved successfully.\n');
end

% =========================================================================
%  DATA LOADING
% =========================================================================

function mouseData = loadMouseData(baseDir, mouseName)
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
                fprintf('    Found: %s\n', salineFolder);
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
if isempty(mouseData)
    fprintf('    No saline folder with learning file found for %s\n', mouseName);
end
end

% =========================================================================
%  LEARNING CURVES
% =========================================================================

function lcData = plotLearningCurves(allMiceData, options)
maxTrials        = options.num_trials;
timeToFoodMatrix = NaN(length(allMiceData), maxTrials);
pathLengthMatrix = NaN(length(allMiceData), maxTrials);
speedMatrix      = NaN(length(allMiceData), maxTrials);

for mouseIdx = 1:length(allMiceData)
    mouseData    = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    for trialIdx = 1:size(learningData, 1)
        trialName  = learningData{trialIdx, 1};
        trialData  = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};
        if ~isempty(trialData) && ~isempty(foodCoords)
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                if trialNum >= 1 && trialNum <= maxTrials
                    [timeToFood, pathLength, meanSpeed] = calculateTrialMetrics(trialData, foodCoords, options);
                    timeToFoodMatrix(mouseIdx, trialNum) = timeToFood;
                    pathLengthMatrix(mouseIdx, trialNum) = pathLength;
                    speedMatrix(mouseIdx, trialNum)      = meanSpeed;
                end
            end
        end
    end
end

lcData.timeToFood = timeToFoodMatrix;
lcData.pathLength = pathLengthMatrix;
lcData.speed      = speedMatrix;

figure('Name', sprintf('Group Learning Curves (%d Trials)', maxTrials), 'Position', [100, 100, 1200, 800]);

subplot(3, 1, 1);
plotMeanWithSEM(1:maxTrials, timeToFoodMatrix, [0.8 0.2 0.8]);
xlabel('Trial Number'); ylabel('Time to Food Discovery (s)');
title('Learning Curve: Time to Food Discovery');
xlim([0.5, maxTrials + 0.5]); grid off;

subplot(3, 1, 2);
plotMeanWithSEM(1:maxTrials, pathLengthMatrix, [0.3 0.6 0.9]);
xlabel('Trial Number'); ylabel('Path Length (cm)');
title('Learning Curve: Path Length');
xlim([0.5, maxTrials + 0.5]); grid off;

subplot(3, 1, 3);
plotMeanWithSEM(1:maxTrials, speedMatrix, [0.9 0.4 0.3]);
xlabel('Trial Number'); ylabel('Mean Speed (cm/s)');
title('Learning Curve: Mean Speed');
xlim([0.5, maxTrials + 0.5]); grid off;

sgtitle(sprintf('Group Learning Curves - %d Trials', maxTrials), 'FontSize', 16, 'FontWeight', 'bold');
end

% =========================================================================
%  DFF vs DISTANCE - INDIVIDUAL TRIALS
% =========================================================================

function [trialDFFData, binCenters] = analyzeDFFvsDistance_IndividualTrials(allMiceData, options)
maxTrials    = options.num_trials;
distanceBins = 0:options.bin_size:options.distance_limit;
binCenters   = distanceBins(1:end-1) + options.bin_size/2;

if options.use_zscore
    fprintf('\nZ-scoring dF/F data across all mice and trials...\n');
    allMiceData = zscore_dff_data(allMiceData, options);
end

if maxTrials <= 4,      rows = 2; cols = 2;
elseif maxTrials <= 6,  rows = 2; cols = 3;
elseif maxTrials <= 9,  rows = 3; cols = 3;
elseif maxTrials <= 12, rows = 3; cols = 4;
else,                   rows = 4; cols = ceil(maxTrials/4);
end

figure('Name', sprintf('dF/F vs Distance: Individual Trials (1-%d)', maxTrials), ...
    'Position', [50, 50, 300*cols, 300*rows]);

rawData = cell(maxTrials, 1);
for t = 1:maxTrials
    rawData{t}.distances   = [];
    rawData{t}.dff         = [];
    rawData{t}.mouseNames  = {};
    rawData{t}.mouseLabels = []; % per-point mouse index for bin-level mouse counting
end

fprintf('\nProcessing dF/F vs distance for individual trials (1-%d)...\n', maxTrials);
for mouseIdx = 1:length(allMiceData)
    mouseData    = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    fprintf('  Processing %s...\n', mouseData.name);
    for trialIdx = 1:size(learningData, 1)
        trialName  = learningData{trialIdx, 1};
        data       = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};
        if ~isempty(data) && ~isempty(foodCoords)
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                if trialNum >= 1 && trialNum <= maxTrials
                    [dffValues, distances] = extractTrajectoryToFood(data, foodCoords, options);
                    if ~isempty(dffValues) && ~isempty(distances)
                        rawData{trialNum}.distances   = [rawData{trialNum}.distances;   distances];
                        rawData{trialNum}.dff         = [rawData{trialNum}.dff;          dffValues];
                        rawData{trialNum}.mouseNames{end+1} = mouseData.name;
                        rawData{trialNum}.mouseLabels = [rawData{trialNum}.mouseLabels; ...
                                                         repmat(mouseIdx, length(distances), 1)];
                        fprintf('    Trial %d: %d data points\n', trialNum, length(dffValues));
                    end
                end
            end
        end
    end
end

all_dff_values = [];
for t = 1:maxTrials
    if ~isempty(rawData{t}.dff), all_dff_values = [all_dff_values; rawData{t}.dff]; end
end
if ~isempty(options.dff_ylim)
    dff_ylim = options.dff_ylim;
elseif ~isempty(all_dff_values)
    dff_ylim = [prctile(all_dff_values, 5), prctile(all_dff_values, 95)];
    rp = (dff_ylim(2) - dff_ylim(1)) * 0.1;
    dff_ylim = [dff_ylim(1) - rp, dff_ylim(2) + rp];
else
    dff_ylim = [-2, 2];
end
fprintf('Using dF/F Y-axis limits: [%.2f, %.2f]\n', dff_ylim(1), dff_ylim(2));

trialDFFData = cell(maxTrials, 1);
for trialNum = 1:maxTrials
    trialDFFData{trialNum}.meanDFF = NaN(1, length(binCenters));
    trialDFFData{trialNum}.semDFF  = NaN(1, length(binCenters));
    trialDFFData{trialNum}.counts  = zeros(1, length(binCenters));
end

for trialNum = 1:maxTrials
    subplot(rows, cols, trialNum);
    if ~isempty(rawData{trialNum}.distances)
        [meanDFF, semDFF, counts] = binDFFByDistance(rawData{trialNum}.distances, ...
            rawData{trialNum}.dff, distanceBins);
        trialDFFData{trialNum}.meanDFF = meanDFF;
        trialDFFData{trialNum}.semDFF  = semDFF;
        trialDFFData{trialNum}.counts  = counts;
        if ~strcmp(options.smoothing_method, 'none')
            meanDFF = applySmoothingToArray(meanDFF, options);
        end
        % Require at least 2 mice per bin (not just 2 data points)
        mouseCounts = countMicePerBin(rawData{trialNum}.distances, rawData{trialNum}.mouseLabels, distanceBins);
        validBins = mouseCounts >= 2 & ~isnan(meanDFF);
        if sum(validBins) > 0
            colors     = flip(slanCM('gem', maxTrials));
            trialColor = colors(trialNum, :);
            plotMeanWithSEMShade(binCenters(validBins), meanDFF(validBins), semDFF(validBins), trialColor);
            xlabel('Distance to Food (cm)');
            if options.use_zscore, ylabel('Z-scored dF/F'); else, ylabel('dF/F (%)'); end
            title(sprintf('Trial %d (n=%d mice)', trialNum, length(unique(rawData{trialNum}.mouseNames))));
            grid off; xlim([0, options.distance_limit]); ylim(dff_ylim);
        else
            text(0.5, 0.5, sprintf('Trial %d\nInsufficient data\n(n=%d mice)', ...
                trialNum, length(unique(rawData{trialNum}.mouseNames))), ...
                'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
            axis off;
        end
    else
        text(0.5, 0.5, sprintf('Trial %d\nNo data', trialNum), ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10, 'Color', [0.5 0.5 0.5]);
        axis off;
    end
end
sgtitle('dF/F vs Distance to Food: Individual Trials (Until Food Discovery)', 'FontSize', 16, 'FontWeight', 'bold');

% Summary overlay figure
figure('Name', sprintf('dF/F vs Distance: Trial Comparison (1-%d)', maxTrials), 'Position', [100, 100, 1000, 600]);
colors        = flip(slanCM('gem', maxTrials));
plotHandles   = [];
legendEntries = {};
hold on;
for trialNum = 1:maxTrials
    if ~isempty(rawData{trialNum}.distances)
        [meanDFF, semDFF, counts] = binDFFByDistance(rawData{trialNum}.distances, rawData{trialNum}.dff, distanceBins);
        if ~strcmp(options.smoothing_method, 'none')
            meanDFF = applySmoothingToArray(meanDFF, options);
            semDFF  = applySmoothingToArray(semDFF, options);
        end
        validBins = countMicePerBin(rawData{trialNum}.distances, rawData{trialNum}.mouseLabels, distanceBins) >= 2 ...
                    & ~isnan(meanDFF);
        if sum(validBins) > 0
            x_valid    = binCenters(validBins);
            mean_valid = meanDFF(validBins);
            sem_valid  = semDFF(validBins);
            if length(x_valid) > 1
                x_fill = [x_valid, fliplr(x_valid)];
                y_fill = [mean_valid + sem_valid, fliplr(mean_valid - sem_valid)];
                vfi = ~isnan(y_fill);
                if sum(vfi) > 0
                    fill(x_fill(vfi), y_fill(vfi), colors(trialNum, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
            end
            h = plot(x_valid, mean_valid, '-', 'LineWidth', 2, 'Color', colors(trialNum, :), ...
                'DisplayName', sprintf('Trial %d (n=%d)', trialNum, length(unique(rawData{trialNum}.mouseNames))));
            plotHandles(end+1)   = h;
            legendEntries{end+1} = sprintf('Trial %d (n=%d)', trialNum, length(unique(rawData{trialNum}.mouseNames)));
            fprintf('Summary plot: Trial %d plotted with %d valid bins\n', trialNum, sum(validBins));
        else
            fprintf('Summary plot: Trial %d skipped - insufficient data\n', trialNum);
        end
    else
        fprintf('Summary plot: Trial %d skipped - no data\n', trialNum);
    end
end
plot([0, options.distance_limit], [0, 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Distance to Food (cm)');
if options.use_zscore, ylabel('Z-scored dF/F'); else, ylabel('dF/F (%)'); end
title('dF/F vs Distance: All Trials Comparison (Until Food Discovery)');
if ~isempty(plotHandles)
    legend(plotHandles, legendEntries, 'Location', 'eastoutside');
else
    text(0.5, 0.5, 'No data to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
end
legend('box', 'off'); grid off;
xlim([0, options.distance_limit]); ylim(dff_ylim);
hold off;

fprintf('\n=== dF/F vs DISTANCE ANALYSIS SUMMARY ===\n');
for trialNum = 1:maxTrials
    if ~isempty(rawData{trialNum}.distances)
        fprintf('Trial %d: %d mice, %d total data points\n', trialNum, ...
            length(unique(rawData{trialNum}.mouseNames)), length(rawData{trialNum}.distances));
    else
        fprintf('Trial %d: No data\n', trialNum);
    end
end
fprintf('==========================================\n');
end

% =========================================================================
%  DFF vs DISTANCE - PHASES
% =========================================================================

function phaseDFFData = analyzeDFFvsDistance_TrialPhases(allMiceData, options)
fprintf('\n=== ANALYZING TRIAL PHASES ===\n');

% Apply z-scoring if requested (must match what IndividualTrials does)
if options.use_zscore
    fprintf('Z-scoring dF/F data across all mice and trials...\n');
    allMiceData = zscore_dff_data(allMiceData, options);
end

distanceBins = 0:options.bin_size:options.distance_limit;
binCenters   = distanceBins(1:end-1) + options.bin_size/2;
nBins        = length(binCenters);

phases = struct();
if ismember('early',  options.plot_phases)
    phases.early.trials  = options.early_trials;
    phases.early.name    = 'Early';
    phases.early.color   = [0.2220 0.8008 0.9653];
end
if ismember('middle', options.plot_phases)
    phases.middle.trials = options.middle_trials;
    phases.middle.name   = 'Middle';
    phases.middle.color  = [0.5057 0.3757 0.9646];
end
if ismember('late',   options.plot_phases)
    phases.late.trials   = options.late_trials;
    phases.late.name     = 'Late';
    phases.late.color    = [0.4562 0.0588 0.3805];
end

phaseNames = fieldnames(phases);

rawData = struct();
for i = 1:length(phaseNames)
    pn = phaseNames{i};
    rawData.(pn).distances   = [];
    rawData.(pn).dff         = [];
    rawData.(pn).mouseNames  = {};
    rawData.(pn).mouseLabels = [];
end

fprintf('Processing dF/F vs distance for trial phases...\n');
for mouseIdx = 1:length(allMiceData)
    mouseData    = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    fprintf('  Processing %s...\n', mouseData.name);
    for trialIdx = 1:size(learningData, 1)
        trialName  = learningData{trialIdx, 1};
        data       = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};
        if ~isempty(data) && ~isempty(foodCoords)
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                for i = 1:length(phaseNames)
                    pn = phaseNames{i};
                    if ismember(trialNum, phases.(pn).trials)
                        [dffValues, distances] = extractTrajectoryToFood(data, foodCoords, options);
                        if ~isempty(dffValues) && ~isempty(distances)
                            rawData.(pn).distances   = [rawData.(pn).distances;   distances];
                            rawData.(pn).dff         = [rawData.(pn).dff;          dffValues];
                            rawData.(pn).mouseNames{end+1} = mouseData.name;
                            rawData.(pn).mouseLabels = [rawData.(pn).mouseLabels; ...
                                                        repmat(mouseIdx, length(distances), 1)];
                            fprintf('    Trial %d (%s phase): %d data points\n', trialNum, phases.(pn).name, length(dffValues));
                        end
                        break;
                    end
                end
            end
        end
    end
end

all_dff_values = [];
for i = 1:length(phaseNames)
    pn = phaseNames{i};
    if ~isempty(rawData.(pn).dff), all_dff_values = [all_dff_values; rawData.(pn).dff]; end
end
if ~isempty(options.dff_ylim)
    dff_ylim = options.dff_ylim;
elseif ~isempty(all_dff_values)
    dff_ylim = [prctile(all_dff_values, 5), prctile(all_dff_values, 95)];
    rp = (dff_ylim(2) - dff_ylim(1)) * 0.1;
    dff_ylim = [dff_ylim(1) - rp, dff_ylim(2) + rp];
else
    dff_ylim = [-2, 2];
end

phaseDFFData = struct();
for i = 1:length(phaseNames)
    pn = phaseNames{i};
    phaseDFFData.(pn).meanDFF = NaN(1, nBins);
    phaseDFFData.(pn).semDFF  = NaN(1, nBins);
    phaseDFFData.(pn).counts  = zeros(1, nBins);
end

figure('Name', 'dF/F vs Distance: Trial Phases Comparison', 'Position', [100, 100, 1000, 600]);
hold on;
plotHandles   = [];
legendEntries = {};

for i = 1:length(phaseNames)
    pn = phaseNames{i};
    if ~isempty(rawData.(pn).distances)
        [meanDFF, semDFF, counts] = binDFFByDistance(rawData.(pn).distances, rawData.(pn).dff, distanceBins);
        phaseDFFData.(pn).meanDFF = meanDFF;
        phaseDFFData.(pn).semDFF  = semDFF;
        phaseDFFData.(pn).counts  = counts;
        if ~strcmp(options.smoothing_method, 'none')
            meanDFF = applySmoothingToArray(meanDFF, options);
            semDFF  = applySmoothingToArray(semDFF, options);
        end
        validBins = countMicePerBin(rawData.(pn).distances, rawData.(pn).mouseLabels, distanceBins) >= 2 ...
                    & ~isnan(meanDFF);
        if sum(validBins) > 0
            x_valid    = binCenters(validBins);
            mean_valid = meanDFF(validBins);
            sem_valid  = semDFF(validBins);
            if length(x_valid) > 1
                x_fill = [x_valid, fliplr(x_valid)];
                y_fill = [mean_valid + sem_valid, fliplr(mean_valid - sem_valid)];
                vfi = ~isnan(y_fill);
                if sum(vfi) > 0
                    fill(x_fill(vfi), y_fill(vfi), phases.(pn).color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
            end
            nMice      = length(unique(rawData.(pn).mouseNames));
            trialsList = sprintf('%d', phases.(pn).trials(1));
            if length(phases.(pn).trials) > 1
                trialsList = sprintf('%d-%d', phases.(pn).trials(1), phases.(pn).trials(end));
            end
            h = plot(x_valid, mean_valid, '-', 'LineWidth', 3, 'Color', phases.(pn).color, ...
                'DisplayName', sprintf('%s (Trials %s, n=%d)', phases.(pn).name, trialsList, nMice));
            plotHandles(end+1)   = h;
            legendEntries{end+1} = sprintf('%s (Trials %s, n=%d)', phases.(pn).name, trialsList, nMice);
            fprintf('%s phase: %d mice, %d data points, %d valid bins\n', ...
                phases.(pn).name, nMice, length(rawData.(pn).distances), sum(validBins));
        end
    end
end

plot([0, options.distance_limit], [0, 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Distance to Food (cm)', 'FontSize', 12);
if options.use_zscore, ylabel('Z-scored dF/F', 'FontSize', 12); else, ylabel('dF/F (%)', 'FontSize', 12); end
title('dF/F vs Distance: Learning Phase Comparison', 'FontSize', 14, 'FontWeight', 'bold');
if ~isempty(plotHandles)
    legend(plotHandles, legendEntries, 'Location', 'northeast', 'Box', 'off');
else
    text(0.5, 0.5, 'No data to plot', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 14);
end
grid off; xlim([0, options.distance_limit]); ylim(dff_ylim);
hold off;
fprintf('=================================\n\n');
end

% =========================================================================
%  SUMMARY STATS
% =========================================================================

function summaryData = printSummaryStats(allMiceData, options)
maxTrials = options.num_trials;
fprintf('\n=== GROUP ANALYSIS SUMMARY (SALINE CONDITION - %d TRIALS) ===\n', maxTrials);
fprintf('Total mice analyzed: %d\n', length(allMiceData));

totalTrials      = 0;
successfulTrials = 0;
rows             = {};

for mouseIdx = 1:length(allMiceData)
    mouseData    = allMiceData{mouseIdx};
    learningData = mouseData.learningData;
    mouseTrials  = 0;
    mouseSucc    = 0;
    for trialIdx = 1:size(learningData, 1)
        trialName = learningData{trialIdx, 1};
        numbers   = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
        if ~isempty(numbers)
            trialNum = str2double(numbers{end}{2});
            if trialNum >= 1 && trialNum <= maxTrials
                mouseTrials = mouseTrials + 1;
                if ~isempty(learningData{trialIdx, 3})
                    mouseSucc = mouseSucc + 1;
                end
            end
        end
    end
    totalTrials      = totalTrials      + mouseTrials;
    successfulTrials = successfulTrials + mouseSucc;
    successRate      = 0;
    if mouseTrials > 0, successRate = 100 * mouseSucc / mouseTrials; end
    rows{end+1} = {mouseData.name, mouseData.condition, mouseTrials, mouseSucc, successRate};
    fprintf('  %s (%s): %d/%d trials successful (trials 1-%d)\n', ...
        mouseData.name, mouseData.condition, mouseSucc, mouseTrials, maxTrials);
end

fprintf('\nOverall: %d/%d trials successful (%.1f%%) for trials 1-%d\n', ...
    successfulTrials, totalTrials, 100*successfulTrials/totalTrials, maxTrials);
fprintf('Smoothing method used: %s\n', options.smoothing_method);
if ~strcmp(options.smoothing_method, 'none')
    fprintf('Smoothing window: %d\n', options.smoothing_window);
    if strcmp(options.smoothing_method, 'sgolay')
        fprintf('Savitzky-Golay order: %d\n', options.sgolay_order);
    end
end
fprintf('Speed threshold: %.1f cm/s\n', options.speed_threshold);
fprintf('Min speed threshold: %.1f cm/s\n', options.min_speed_threshold);
fprintf('Excluded first frames: %d\n', options.exclude_first_n_frames);
fprintf('==========================================\n\n');

totalRate = 0;
if totalTrials > 0, totalRate = 100 * successfulTrials / totalTrials; end
rows{end+1} = {'TOTAL', 'all', totalTrials, successfulTrials, totalRate};
summaryData.rows = vertcat(rows{:});
end

% =========================================================================
%  HELPER FUNCTIONS
% =========================================================================

function plotMeanWithSEM(x, dataMatrix, color)
meanData = nanmean(dataMatrix, 1);
semData  = nanstd(dataMatrix, 0, 1) ./ sqrt(sum(~isnan(dataMatrix), 1));
hold on;
plot(x, meanData, '-o', 'MarkerSize', 2, 'Color', color, 'MarkerFaceColor', color);
validIdx = ~isnan(meanData) & ~isnan(semData);
if sum(validIdx) > 0
    x_valid    = x(validIdx);
    mean_valid = meanData(validIdx);
    sem_valid  = semData(validIdx);
    patch([x_valid, fliplr(x_valid)], [mean_valid + sem_valid, fliplr(mean_valid - sem_valid)], ...
        color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
hold off;
end

function [timeToFood, pathLength, meanSpeed] = calculateTrialMetrics(trialData, foodCoords, options)
% For path length: only exclude first N frames (no speed filter) so trajectory
% gaps from removed high-speed frames don't artificially shorten the path.
% For timeToFood and meanSpeed: apply full filter.
n = options.exclude_first_n_frames;
if n > 0 && size(trialData, 1) > n
    trialDataForPath = trialData(n+1:end, :);
else
    trialDataForPath = trialData;
end
trialDataFiltered = applyFrameAndSpeedFilter(trialData, options);

timeToFood = NaN; pathLength = NaN; meanSpeed = NaN;
if isempty(trialDataFiltered) || isempty(foodCoords), return; end

pix_to_cm = 0.16;

% --- timeToFood and meanSpeed from fully filtered data ---
x_f      = trialDataFiltered(:, 2);
y_f      = trialDataFiltered(:, 3);
t_f      = trialDataFiltered(:, 1);
spd_f    = trialDataFiltered(:, 7);
vi       = ~isnan(x_f) & ~isnan(y_f) & ~isnan(t_f);
x_fc     = x_f(vi); y_fc = y_f(vi); t_fc = t_f(vi); spd_fc = spd_f(vi);
if length(x_fc) < 2 || length(foodCoords) < 2, return; end
dist_f   = sqrt((x_fc - foodCoords(1)).^2 + (y_fc - foodCoords(2)).^2) * pix_to_cm;
fdIdx    = find(dist_f < options.food_discovery_threshold, 1, 'first');
if ~isempty(fdIdx)
    timeToFood  = t_fc(fdIdx);
    valid_speed = spd_fc(1:fdIdx);
    valid_speed = valid_speed(~isnan(valid_speed));
    if ~isempty(valid_speed), meanSpeed = mean(valid_speed); end
end

% --- pathLength from frame-excluded but NOT speed-filtered data ---
x_p   = trialDataForPath(:, 2);
y_p   = trialDataForPath(:, 3);
vip   = ~isnan(x_p) & ~isnan(y_p);
x_pc  = x_p(vip); y_pc = y_p(vip);
if length(x_pc) < 2, return; end
dist_p = sqrt((x_pc - foodCoords(1)).^2 + (y_pc - foodCoords(2)).^2) * pix_to_cm;
fdIdxP = find(dist_p < options.food_discovery_threshold, 1, 'first');
if ~isempty(fdIdxP)
    dx = diff(x_pc(1:fdIdxP));
    dy = diff(y_pc(1:fdIdxP));
    pathLength = sum(sqrt(dx.^2 + dy.^2)) * pix_to_cm;
end
end

function smoothedArray = applySmoothingToArray(dataArray, options)
smoothedArray = dataArray;
validIdx = ~isnan(dataArray);
if sum(validIdx) >= options.smoothing_window
    validData      = dataArray(validIdx);
    validPositions = find(validIdx);
    switch options.smoothing_method
        case 'movmean'
            if length(validData) >= options.smoothing_window
                smoothedArray(validPositions) = movmean(validData, options.smoothing_window);
            end
        case 'sgolay'
            if length(validData) >= options.smoothing_window && options.smoothing_window > options.sgolay_order
                try
                    smoothedArray(validPositions) = sgolayfilt(validData, options.sgolay_order, options.smoothing_window);
                catch
                    smoothedArray(validPositions) = movmean(validData, options.smoothing_window);
                end
            end
    end
end
end

function [dffValues, distances] = extractTrajectoryToFood(trialData, foodCoords, options)
% Apply frame and speed filters before extracting trajectory
trialData  = applyFrameAndSpeedFilter(trialData, options);
dffValues  = []; distances = [];
x_coords   = trialData(:, 2);
y_coords   = trialData(:, 3);
dff_vals   = trialData(:, 6);
valid_idx  = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(dff_vals);
x_clean    = x_coords(valid_idx);
y_clean    = y_coords(valid_idx);
dff_clean  = dff_vals(valid_idx);
if length(x_clean) < 2 || length(foodCoords) < 2, return; end
pix_to_cm        = 0.16;
dist_to_food     = sqrt((x_clean - foodCoords(1)).^2 + (y_clean - foodCoords(2)).^2) * pix_to_cm;
foodDiscoveryIdx = find(dist_to_food < options.food_discovery_threshold, 1, 'first');
if ~isempty(foodDiscoveryIdx)
    distances = dist_to_food(1:foodDiscoveryIdx);
    dffValues = dff_clean(1:foodDiscoveryIdx);
else
    withinRange = dist_to_food <= options.distance_limit;
    distances   = dist_to_food(withinRange);
    dffValues   = dff_clean(withinRange);
end
end

function [meanDFF, semDFF, counts] = binDFFByDistance(distances, dffValues, distanceBins)
meanDFF = NaN(1, length(distanceBins)-1);
semDFF  = NaN(1, length(distanceBins)-1);
counts  = zeros(1, length(distanceBins)-1);
for i = 1:length(distanceBins)-1
    binIdx   = distances >= distanceBins(i) & distances < distanceBins(i+1);
    validDFF = dffValues(binIdx & ~isnan(dffValues));
    if ~isempty(validDFF)
        meanDFF(i) = mean(validDFF);
        semDFF(i)  = 0;
        if length(validDFF) > 1, semDFF(i) = std(validDFF) / sqrt(length(validDFF)); end
        counts(i)  = length(validDFF);
    end
end
end

function zscoredMiceData = zscore_dff_data(allMiceData, options)
zscoredMiceData = allMiceData;
for mouseIdx = 1:length(allMiceData)
    mouseData     = allMiceData{mouseIdx};
    learningData  = mouseData.learningData;
    all_mouse_dff = [];
    for trialIdx = 1:size(learningData, 1)
        trialData = learningData{trialIdx, 3};
        if ~isempty(trialData)
            % Apply filters before computing z-score parameters so the
            % baseline is computed on the same data that will be analysed
            trialData = applyFrameAndSpeedFilter(trialData, options);
            dff_vals  = trialData(:, 6);
            all_mouse_dff = [all_mouse_dff; dff_vals(~isnan(dff_vals))];
        end
    end
    if ~isempty(all_mouse_dff)
        dff_mean = mean(all_mouse_dff);
        dff_std  = std(all_mouse_dff);
        fprintf('  %s: Mean=%.3f, STD=%.3f (n=%d points)\n', mouseData.name, dff_mean, dff_std, length(all_mouse_dff));
        for trialIdx = 1:size(learningData, 1)
            if ~isempty(learningData{trialIdx, 3})
                orig = zscoredMiceData{mouseIdx}.learningData{trialIdx, 3}(:, 6);
                zscoredMiceData{mouseIdx}.learningData{trialIdx, 3}(:, 6) = (orig - dff_mean) / dff_std;
            end
        end
    else
        fprintf('  %s: No valid dF/F data found\n', mouseData.name);
    end
end
end

function plotMeanWithSEMShade(x, meanData, semData, color)
hold on;
if length(x) > 1
    x_fill = [x, fliplr(x)];
    y_fill = [meanData + semData, fliplr(meanData - semData)];
    vfi = ~isnan(y_fill);
    if sum(vfi) > 0
        fill(x_fill(vfi), y_fill(vfi), color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end
plot(x, meanData, '-', 'LineWidth', 2, 'Color', color);
plot([min(x), max(x)], [0 0], 'k--', 'LineWidth', 1);
hold off;
end

function mouseCounts = countMicePerBin(distances, mouseLabels, distanceBins)
% Count the number of unique mice contributing data to each distance bin.
mouseCounts = zeros(1, length(distanceBins)-1);
for i = 1:length(distanceBins)-1
    binIdx = distances >= distanceBins(i) & distances < distanceBins(i+1);
    if any(binIdx)
        mouseCounts(i) = length(unique(mouseLabels(binIdx)));
    end
end
end