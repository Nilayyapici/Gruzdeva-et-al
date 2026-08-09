function processBatchTrials(mouseName, mouseDir, x_food, y_food)
    % Batch process all trials for a single mouse
    % This script processes all trial folders (A1, A2, B4, etc.) for any mouse
    %
    % Usage: processBatchTrials('MDRE8', 'C:\Users\Anna\...\MDRE8\MDRE8_saline2', 301, 184)
    
    % Initialize the learning cell array with dynamic name
    learningData = {};
    
    % Get all trial folders (A1, A2, B4, etc.)
    trialFolders = dir(mouseDir);
    trialFolders = trialFolders([trialFolders.isdir] & ~startsWith({trialFolders.name}, '.'));
    
    % Filter to get only trial folders (those that start with A, B, C, etc.)
    trialNames = {};
    for i = 1:length(trialFolders)
        if ~isempty(regexp(trialFolders(i).name, '^[A-Z]\d+$', 'once'))
            trialNames{end+1} = trialFolders(i).name;
        end
    end
    
    % Sort trial names by letter then number (A1, A2, A10, B1, B2, etc.)
    if ~isempty(trialNames)
        % Extract letters and numbers for proper sorting
        trialParts = cell(length(trialNames), 2);
        for i = 1:length(trialNames)
            % Use regex to separate letter(s) and number(s)
            tokens = regexp(trialNames{i}, '^([A-Z]+)(\d+)$', 'tokens');
            if ~isempty(tokens)
                trialParts{i, 1} = tokens{1}{1};  % Letter part
                trialParts{i, 2} = str2double(tokens{1}{2});  % Number part
            else
                trialParts{i, 1} = trialNames{i};  % Fallback
                trialParts{i, 2} = 0;
            end
        end
        
        % Sort by letter first, then by number
        [~, sortIdx] = sortrows(trialParts, [1, 2]);
        trialNames = trialNames(sortIdx);
    end
    
    fprintf('Found %d trials to process: %s\n', length(trialNames), strjoin(trialNames, ', '));
    fprintf('Food location: x=%.1f, y=%.1f\n', x_food, y_food);
    
    % Process each trial
    for trialIdx = 1:length(trialNames)
        trialName = trialNames{trialIdx};
        trialPath = fullfile(mouseDir, trialName);
        
        fprintf('\nProcessing trial %d/%d: %s\n', trialIdx, length(trialNames), trialName);
        
        try
            % Change to trial directory
            cd(trialPath);
            
            % Process this trial
            [trialData, foodCoords] = processSingleTrial(mouseName, trialName, x_food, y_food);
            
            % Store in learning array
            learningData{trialIdx, 1} = sprintf('%s_%s', mouseName, trialName);
            learningData{trialIdx, 2} = foodCoords;
            learningData{trialIdx, 3} = trialData;
            
            fprintf('  - Successfully processed %s\n', trialName);
            
        catch ME
            fprintf('  - Error processing %s: %s\n', trialName, ME.message);
            % Store empty data for failed trials
            learningData{trialIdx, 1} = sprintf('%s_%s', mouseName, trialName);
            learningData{trialIdx, 2} = [];
            learningData{trialIdx, 3} = [];
        end
    end
    
    % Save the complete learning array with dynamic variable name
    cd(mouseDir);
    
    % Create variable with dynamic name (mouseName_learning)
    eval(sprintf('%s_learning = learningData;', mouseName));
    
    % Save with dynamic filename
    save(sprintf('%s_learning.mat', mouseName), sprintf('%s_learning', mouseName));
    
    fprintf('\nBatch processing complete! Saved %s_learning.mat with %d trials\n', mouseName, length(trialNames));
    
    % Display summary
    fprintf('\nSummary:\n');
    for i = 1:size(learningData, 1)
        if ~isempty(learningData{i, 3})
            fprintf('  %s: %d data points\n', learningData{i, 1}, size(learningData{i, 3}, 1));
        else
            fprintf('  %s: FAILED\n', learningData{i, 1});
        end
    end
end

function [txy_zones, foodCoords] = processSingleTrial(mouseName, trialName, x_food, y_food)
    % Process a single trial - adapted from your original script
    % Now accepts food coordinates as parameters
    
    % Set smoothing parameters
    order = 3;
    framelen = 501;
    
    % Find files automatically - flexible patterns for complex naming
    % Try multiple patterns to handle different naming conventions
    
    % Pattern 1: Standard naming (MDRE8_A1_*)
    trackingFiles = dir(sprintf('%s_%s20*.csv', mouseName, trialName));
    analogFiles = dir(sprintf('%s_%s_AI20*.csv', mouseName, trialName));
    dlcFiles = dir(sprintf('%s_%s20*DLC*.csv', mouseName, trialName));
    
    % Pattern 2: With "saline" or other text (F13_saline_A1_*)
    if isempty(trackingFiles)
        trackingFiles = dir(sprintf('%s*%s20*.csv', mouseName, trialName));
    end
    if isempty(analogFiles)
        analogFiles = dir(sprintf('%s*%s_AI20*.csv', mouseName, trialName));
    end
    if isempty(dlcFiles)
        dlcFiles = dir(sprintf('%s*%s20*DLC*.csv', mouseName, trialName));
    end
    
    % Pattern 3: Most flexible - search for any file containing mouse name and trial
    if isempty(trackingFiles)
        allCSV = dir('*.csv');
        for i = 1:length(allCSV)
            fname = allCSV(i).name;
            if contains(fname, mouseName) && contains(fname, trialName) && ~contains(fname, 'AI') && ~contains(fname, 'DLC')
                trackingFiles = allCSV(i);
                break;
            end
        end
    end
    if isempty(analogFiles)
        allCSV = dir('*AI*.csv');
        for i = 1:length(allCSV)
            fname = allCSV(i).name;
            if contains(fname, mouseName) && contains(fname, trialName)
                analogFiles = allCSV(i);
                break;
            end
        end
    end
    if isempty(dlcFiles)
        allCSV = dir('*DLC*.csv');
        for i = 1:length(allCSV)
            fname = allCSV(i).name;
            if contains(fname, mouseName) && contains(fname, trialName)
                dlcFiles = allCSV(i);
                break;
            end
        end
    end
    
    % Photometry and AIN files - flexible patterns to handle "saline" etc.
    photFiles = dir(sprintf('photometry_%s_%s_*.mat', mouseName, trialName));
    ainFiles = dir(sprintf('ain_%s_%s_*.mat', mouseName, trialName));
    
    % If not found, try with wildcard for middle text (like "saline")
    if isempty(photFiles)
        photFiles = dir(sprintf('photometry_%s*%s_*.mat', mouseName, trialName));
    end
    if isempty(ainFiles)
        ainFiles = dir(sprintf('ain_%s*%s_*.mat', mouseName, trialName));
    end
    
    % If still not found, search for any photometry/ain files containing mouse and trial
    if isempty(photFiles)
        allMat = dir('photometry*.mat');
        for i = 1:length(allMat)
            fname = allMat(i).name;
            if contains(fname, mouseName) && contains(fname, trialName)
                photFiles = allMat(i);
                break;
            end
        end
    end
    if isempty(ainFiles)
        allMat = dir('ain*.mat');
        for i = 1:length(allMat)
            fname = allMat(i).name;
            if contains(fname, mouseName) && contains(fname, trialName)
                ainFiles = allMat(i);
                break;
            end
        end
    end
    
    % Debug: Show what files were found
    fprintf('    Found files:\n');
    if ~isempty(trackingFiles), fprintf('      Tracking: %s\n', trackingFiles(1).name); end
    if ~isempty(analogFiles), fprintf('      Analog: %s\n', analogFiles(1).name); end
    if ~isempty(dlcFiles), fprintf('      DLC: %s\n', dlcFiles(1).name); end
    if ~isempty(photFiles), fprintf('      Photometry: %s\n', photFiles(1).name); end
    if ~isempty(ainFiles), fprintf('      AIN: %s\n', ainFiles(1).name); end
    
    % Check if all required files exist
    if isempty(trackingFiles) || isempty(analogFiles) || isempty(dlcFiles) || ...
       isempty(photFiles) || isempty(ainFiles)
        missingFiles = {};
        if isempty(trackingFiles), missingFiles{end+1} = 'tracking'; end
        if isempty(analogFiles), missingFiles{end+1} = 'analog'; end
        if isempty(dlcFiles), missingFiles{end+1} = 'DLC'; end
        if isempty(photFiles), missingFiles{end+1} = 'photometry'; end
        if isempty(ainFiles), missingFiles{end+1} = 'ain'; end
        error('Missing required files for %s_%s: %s', mouseName, trialName, strjoin(missingFiles, ', '));
    end
    
    % Load data files
    tracking = readtable(trackingFiles(1).name);
    Analog_table = readtable(analogFiles(1).name);
    DLC = readtable(dlcFiles(1).name);
    
    load(photFiles(1).name); % loads 'photometry'
    load(ainFiles(1).name);   % loads 'ain'
    
    %% Handle DLC/tracking frame mismatch
    % Check if DLC has exactly 1 more frame than tracking
    trackingRows = size(tracking, 1);
    dlcRows = size(DLC, 1);
    
    if dlcRows == trackingRows + 1
        fprintf('    Warning: DLC file has 1 extra frame. Removing first frame from DLC.\n');
        DLC(1, :) = [];  % Remove first row from DLC
    elseif dlcRows ~= trackingRows
        error('Frame count mismatch: Tracking has %d frames, DLC has %d frames. Difference is %d (only 1-frame difference can be auto-corrected).', ...
            trackingRows, dlcRows, abs(dlcRows - trackingRows));
    end
    
    %% Process tracking data
    txy = [];
    txy(:,1) = table2array(tracking(:,1));
    txy(:,2:5) = table2array(DLC(:,[2,3,5,6]));
    
    Analog_copy = table2array(Analog_table);
    Analog_time = Analog_copy(:,1);
    Analog_data = Analog_copy(:,2);
    
    %% Process analog data
    Copy_photometry = find(Analog_data <= 10);
    Analog_time(Copy_photometry) = [];
    Analog_data(Copy_photometry) = [];
    Analog_bonsai = [Analog_time, Analog_data];
    
    %% Synchronization
    start_phot = Analog_bonsai(2,1);
    frame_beh = findnearest(start_phot, txy(:,1));
    dt_synch = (start_phot - txy(1, 1))/1000;
    
    txy(:, 1) = (txy(:, 1) - txy(1, 1))/1000;
    Analog_bonsai(:, 1) = (Analog_bonsai(:, 1) - Analog_bonsai(1, 1))/1000;
    
    time_tracking = txy(:,1);
    txy(1:frame_beh-1, :) = [];
    txy(:, 1) = (txy(:, 1) - txy(1, 1));
    
    %% Process photometry
    photometry(:,1) = photometry(:,1)*1.000138;
    ain(:,1) = ain(:,1)*1.000138;
    
    %% Validation figure - check synchronization
    figure('Name', sprintf('Validation: %s_%s', mouseName, trialName), 'Position', [100, 100, 800, 400]);
    hold all
    plot(Analog_bonsai(:,1), Analog_bonsai(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analog TTL');
    plot(photometry(:,1), photometry(:,3)*100, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Isosbestic x100');
    plot(ain(:,1), ain(:,2)*100, 'g-', 'LineWidth', 1.5, 'DisplayName', 'AIN x100');
    xlabel('Time (s)');
    ylabel('Signal');
    title(sprintf('Synchronization Check: %s_%s', mouseName, trialName));
    legend('Location', 'best');
    grid on;
    hold off;
    
    % Pause to allow inspection
    fprintf('    Press any key to continue to next trial...\n');
    pause;
    
    %% Arena coordinates
    r_food = 20;
    pix_to_cm = 75/400;  % Pixel to cm conversion factor
    
    % Food coordinates are now passed as parameters
    % x_food and y_food are function parameters
    
    %% Calculate speed
    % Calculate displacement in pixels
    dx = diff(txy(:, 4));  % centroid x displacement
    dy = diff(txy(:, 5));  % centroid y displacement
    dt = diff(txy(:, 1));  % time difference
    
    % Calculate speed in pixels/second
    displacement_pixels = sqrt(dx.^2 + dy.^2);
    speed_pixels_per_sec = displacement_pixels ./ dt;
    
    % Convert to cm/sec
    speed_cm_per_sec = speed_pixels_per_sec * pix_to_cm;
    
    % Create speed array - speed has one less point than position
    xy_speed = txy(1:end-1, :);  % Remove last point to match speed array length
    xy_speed(:, 6) = speed_cm_per_sec;
    
    % Remove first point
    xy_speed(1, :) = [];
    
    %% Remove tracking artifacts
    xy_corrected = xy_speed;
    for i = 1:length(xy_speed)
        if xy_speed(i, 6) >= 60 
           xy_corrected(i, :) = NaN;
        end
    end
    xy_corrected = rmmissing(xy_corrected);
    
    %% Process photometry data
    auto = photometry(:, 3);
    gcamp = photometry(:, 2);
    time_phot = photometry(:, 1);
    
    % Fitting and dFF calculation
    reg = polyfit(auto, gcamp, 1);
    a = reg(1);
    b = reg(2);
    controlFit = a.*auto + b;
    dff = (gcamp - controlFit)./controlFit;
    dff = dff * 100;
    
    tautodff = [time_phot, auto, gcamp, dff];
    
    % Apply smoothing
    tautodff(:,3) = sgolayfilt(tautodff(:, 3), order, framelen);
    tautodff(:,2) = sgolayfilt(tautodff(:, 2), order, framelen);
    
    % Recalculate dFF after smoothing
    reg = polyfit(tautodff(:,2), tautodff(:,3), 1);
    a = reg(1);
    b = reg(2);
    controlFit = a.*tautodff(:,2) + b;
    dff = (tautodff(:,3) - controlFit)./controlFit;
    dff = dff * 100;
    tautodff(:,4) = dff;
    
    %% Interpolation
    txy_phot_speed = [];
    interp_coord = interp1(xy_corrected(:, 1), xy_corrected(:, [2,3,6]), tautodff(:,1));
    txy_phot_speed = time_phot;
    txy_phot_speed(:, [2,3,7]) = interp_coord;
    txy_phot_speed(:, 4) = tautodff(:,3);
    txy_phot_speed(:, 5) = tautodff(:,2);
    txy_phot_speed(:, 6) = tautodff(:,4);
    interp_coord_centroid = interp1(xy_corrected(:, 1), xy_corrected(:, [4,5]), tautodff(:,1));
    txy_phot_speed(:, [8 9]) = interp_coord_centroid;
    
    % Clean up timing if needed
    start_index = findnearest(5, txy_phot_speed(:,1));
    txy_phot_speed(1:start_index, :) = [];
    
    %% Binning
    n = 0.05;
    time_step = mean(diff(txy_phot_speed(:,1)));
    
    txy_phot_sp_bin = zeros([round(length(txy_phot_speed)/round(n/time_step))-1, size(txy_phot_speed,2)]);
    for i = round(n/time_step):round(n/time_step):length(txy_phot_speed)
        txy_phot_sp_bin(i/round(n/time_step),1) = txy_phot_speed(i,1);
        txy_phot_sp_bin(i/round(n/time_step),[2,3]) = txy_phot_speed(i,[2,3]);
        txy_phot_sp_bin(i/round(n/time_step),4:end) = mean(txy_phot_speed(i-round(n/time_step)+1:i,4:end));
    end
    
    %% Final processing
    txy_zones = txy_phot_sp_bin;
    
    % Calculate distance to food (9th column) using passed coordinates
    for i = 1:length(txy_zones)
        txy_zones(i,9) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2);
    end
    
    % Return food coordinates
    foodCoords = [x_food; y_food];
    
    % Save individual trial data
    save('txy_phot_speed.mat', 'txy_phot_speed');
end