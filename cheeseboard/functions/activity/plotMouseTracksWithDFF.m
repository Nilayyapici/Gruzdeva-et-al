function plotMouseTracksWithDFF(mouseName, mouseDir)
    % Plot mouse tracks with dFF values overlayed to check data quality
    % 
    % Usage: plotMouseTracksWithDFF('MDRE8', 'C:\Users\Anna\...\MDRE8\MDRE8_saline2')
    
    if nargin < 2
        % Default values if not provided
        mouseName = 'MDRE8';
        mouseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard\MDRE8\MDRE8_saline2';
    end
    
    % Load the learning data
    cd(mouseDir);
    dataFile = sprintf('%s_learning.mat', mouseName);
    
    if ~exist(dataFile, 'file')
        error('File %s not found in %s', dataFile, mouseDir);
    end
    
    load(dataFile);
    
    % Get the learning data variable dynamically
    learningData = eval(sprintf('%s_learning', mouseName));
    
    % Get number of trials
    numTrials = size(learningData, 1);
    
    % Set up subplot grid - adjust based on number of trials
    if numTrials <= 4
        rows = 2; cols = 2;
    elseif numTrials <= 6
        rows = 2; cols = 3;
    elseif numTrials <= 9
        rows = 3; cols = 3;
    elseif numTrials <= 12
        rows = 3; cols = 4;
    else
        rows = 4; cols = 4;
    end
    
    % Create figure
    figure('Name', sprintf('%s Learning Tracks with dFF', mouseName), ...
           'Position', [100, 100, 1200, 800]);
    
    % Arena parameters (adjust if needed)
    x0 = 336; y0 = 262;
    arena_radius = 227;  % External radius
    
    % Create array to store trial data by position number
    trialsByPosition = cell(12, 3); % Max 12 positions, 3 columns (name, food, data)
    
    % Extract trial numbers and place trials in correct positions
    fprintf('Organizing trials by position:\n');
    for i = 1:numTrials
        trialName = learningData{i, 1};
        foodCoords = learningData{i, 2};
        trialData = learningData{i, 3};
        
        % Extract number from trial name - find the last number in the string
        numbers = regexp(trialName, '\d+', 'match');
        if ~isempty(numbers)
            % Take the last number found (in case there are multiple numbers)
            trialNumber = str2double(numbers{end});
            fprintf('  %s -> Position %d\n', trialName, trialNumber);
            if trialNumber >= 1 && trialNumber <= 12
                trialsByPosition{trialNumber, 1} = trialName;
                trialsByPosition{trialNumber, 2} = foodCoords;
                trialsByPosition{trialNumber, 3} = trialData;
            end
        else
            fprintf('  %s -> Could not extract number\n', trialName);
        end
    end
    
    % Check which positions have data
    fprintf('Positions with data: ');
    for pos = 1:12
        if ~isempty(trialsByPosition{pos, 1})
            fprintf('%d ', pos);
        end
    end
    fprintf('\n');
    
    % Calculate global dFF range for consistent color scaling
    allDFFValues = [];
    for pos = 1:12
        if ~isempty(trialsByPosition{pos, 3})
            trialData = trialsByPosition{pos, 3};
            dff_values = trialData(:, 6);
            valid_dff = dff_values(~isnan(dff_values));
            allDFFValues = [allDFFValues; valid_dff];
        end
    end
    
    if ~isempty(allDFFValues)
        globalMin = prctile(allDFFValues, 5);  % Use 5th percentile to avoid outliers
        globalMax = prctile(allDFFValues, 95); % Use 95th percentile to avoid outliers
        fprintf('Global dFF range: %.2f to %.2f%%\n', globalMin, globalMax);
    else
        globalMin = -5;
        globalMax = 5;
        fprintf('No dFF data found, using default range: %.1f to %.1f%%\n', globalMin, globalMax);
    end
    
    % Plot trials in their numbered positions
    for position = 1:min(12, rows*cols)
        subplot(rows, cols, position);
        
        trialName = trialsByPosition{position, 1};
        foodCoords = trialsByPosition{position, 2};
        trialData = trialsByPosition{position, 3};
        
        % Handle empty positions
        if isempty(trialName)
            text(0.5, 0.5, sprintf('Trial %d\nNot Found', position), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
            axis([0 1 0 1]);
            title(sprintf('Trial %d - Missing', position), 'Color', [0.5 0.5 0.5]);
            continue;
        end
        
        % Fix title - replace underscores with spaces
        titleText = strrep(trialName, '_', ' ');
        
        if isempty(trialData)
            axis off
            continue;
        end
        
        % Extract coordinates and dFF
        x_coords = trialData(:, 2);  % x coordinates (nose)
        y_coords = trialData(:, 3);  % y coordinates (nose)
        dff_values = trialData(:, 6); % dFF values
        
        % Remove NaN values
        valid_idx = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(dff_values);
        x_coords = x_coords(valid_idx);
        y_coords = y_coords(valid_idx);
        dff_values = dff_values(valid_idx);
        
        if isempty(x_coords)
            text(200, 225, sprintf('%s\nNo valid data', titleText), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'red');
            axis([0 400 0 450]);
            title(titleText, 'Color', 'red');
            continue;
        end
        
        % Find food discovery point (nose within 0.5 cm of food)
        foodDiscoveryIdx = [];
        if ~isempty(foodCoords) && length(foodCoords) >= 2
            pix_to_cm = 75/400;
            distances_to_food = sqrt((x_coords - foodCoords(1)).^2 + (y_coords - foodCoords(2)).^2) * pix_to_cm;
            foodDiscoveryIdx = find(distances_to_food < 1.2, 1, 'first');
        end
        
        % Truncate trajectory at food discovery (if found)
        if ~isempty(foodDiscoveryIdx)
            x_coords = x_coords(1:foodDiscoveryIdx);
            y_coords = y_coords(1:foodDiscoveryIdx);
            dff_values = dff_values(1:foodDiscoveryIdx);
            trajectoryType = 'to food';
            fprintf('    %s: Trajectory until food discovery (%d points)\n', titleText, length(x_coords));
        else
            trajectoryType = 'full';
            fprintf('    %s: Full trajectory - no food discovery (%d points)\n', titleText, length(x_coords));
        end
        
        % Plot arena boundary
        hold on;
        theta = linspace(0, 2*pi, 100);
        arena_x = x0 + arena_radius * cos(theta);
        arena_y = y0 + arena_radius * sin(theta);
        plot(arena_x, arena_y, 'k-', 'LineWidth', 2);
        box off
        axis off
        
        % Create scatter plot with dFF color coding
        h = scatter(x_coords, y_coords, 8, dff_values, 'filled');
        
        % Set alpha (transparency) if supported by MATLAB version
        try
            h.MarkerFaceAlpha = 0.7;
        catch
            % Alpha not supported in older MATLAB versions - skip transparency
        end
        
        % Set colormap and consistent color axis
        colormap(blueWhiteRed);
        caxis([globalMin, globalMax]);  % Use global range for all plots
        
        % Set axis properties
        axis equal;
        xlim([0, 650]);
        ylim([0, 500]);
        xlabel('X coordinate (pixels)');
        ylabel('Y coordinate (pixels)');

        % Plot food location if available
        if ~isempty(foodCoords) && length(foodCoords) >= 2
            scatter(foodCoords(1), foodCoords(2), 60, 'y', 'filled', 's', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 2);
        end
        
        % Add start point (green circle)
        if ~isempty(x_coords)
            scatter(x_coords(1), y_coords(1), 30, 'g', 'filled', 'o', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1);
        end
        
        % % Add end point (red circle if food found, blue if not)
        % if ~isempty(x_coords)
        %     if strcmp(trajectoryType, 'to food')
        %         scatter(x_coords(end), y_coords(end), 30, 'r', 'filled', 'o', ...
        %             'MarkerEdgeColor', 'k', 'LineWidth', 1);
        %     else
        %         scatter(x_coords(end), y_coords(end), 30, 'b', 'filled', 'o', ...
        %             'MarkerEdgeColor', 'k', 'LineWidth', 1);
        %     end
        % end
        
        title(sprintf('%s (%s)', titleText, trajectoryType));
        
        % % Add statistics text
        % if ~isempty(foodDiscoveryIdx)
        %     validIdx = valid_idx(1:foodDiscoveryIdx);
        %     timeToFood = trialData(foodDiscoveryIdx, 1);
        %     text(10, 480, sprintf('Food found: %.1fs', timeToFood), ...
        %         'FontSize', 8, 'BackgroundColor', 'white', 'Color', 'red');
        % else
        %     text(10, 480, sprintf('No food discovery'), ...
        %         'FontSize', 8, 'BackgroundColor', 'white', 'Color', 'blue');
        % end
        
        hold off;
    end
    
    % Add single colorbar for the entire figure
    c = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
    c.Label.String = 'dF/F (%)';
    
    % Add overall title
    sgtitle(sprintf('%s Learning Sessions', mouseName), ...
        'FontSize', 16, 'FontWeight', 'bold');
    
    % Create summary figure
    createSummaryPlots(learningData, mouseName);
end


function createSummaryPlots(learningData, mouseName)
    % Create summary plots for time to food discovery, path length and mean speed across trials 1-12
    
    figure('Name', sprintf('%s Learning Progression', mouseName), ...
           'Position', [150, 150, 350, 800]);
    
    % Initialize arrays for 12 trials
    trialNumbers = 1:12;
    timeToFood = NaN(12, 1);
    pathLengths = NaN(12, 1);
    meanSpeeds = NaN(12, 1);
    trialLabels = cell(12, 1);
    
    % Extract time to food discovery, path length and mean speed for each trial
    numTrials = size(learningData, 1);
    
    for i = 1:numTrials
        trialName = learningData{i, 1};
        trialData = learningData{i, 3};
        foodCoords = learningData{i, 2};
        
        % Extract trial number
        numbers = regexp(trialName, '\d+', 'match');
        if ~isempty(numbers)
            trialNum = str2double(numbers{end});
            
            if trialNum >= 1 && trialNum <= 12 && ~isempty(trialData) && ~isempty(foodCoords)
                % Calculate time to food discovery
                x_coords = trialData(:, 2);
                y_coords = trialData(:, 3);
                time_coords = trialData(:, 1);
                
                % Remove NaN values
                valid_idx = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(time_coords);
                x_clean = x_coords(valid_idx);
                y_clean = y_coords(valid_idx);
                time_clean = time_coords(valid_idx);
                
                if ~isempty(x_clean) && length(foodCoords) >= 2
                    % Calculate distance to food for each timepoint
                    pix_to_cm = 75/400;
                    distances_to_food = sqrt((x_clean - foodCoords(1)).^2 + (y_clean - foodCoords(2)).^2) * pix_to_cm;
                    
                    % Find first time when distance < 0.5 cm
                    foodDiscoveryIdx = find(distances_to_food < 1.2, 1, 'first');
                    
                    if ~isempty(foodDiscoveryIdx)
                        timeToFood(trialNum) = time_clean(foodDiscoveryIdx);
                    end
                end
                
                % Calculate path length
                if length(x_clean) > 1
                    % Calculate total path length in cm
                    dx = diff(x_clean);
                    dy = diff(y_clean);
                    distances = sqrt(dx.^2 + dy.^2);
                    pix_to_cm = 75/400;
                    pathLengths(trialNum) = sum(distances) * pix_to_cm;
                end
                
                % Calculate mean speed
                speed_values = trialData(:, 7);
                valid_speed = speed_values(~isnan(speed_values));
                
                if ~isempty(valid_speed)
                    meanSpeeds(trialNum) = mean(valid_speed);
                end
                
                % Create clean trial label (extract just the trial part)
                trialPart = regexp(trialName, '([A-Z]+\d+)$', 'tokens');
                if ~isempty(trialPart)
                    trialLabels{trialNum} = trialPart{1}{1};
                else
                    trialLabels{trialNum} = sprintf('Trial %d', trialNum);
                end
            end
        end
    end
    
    % Fill missing trial labels
    for i = 1:12
        if isempty(trialLabels{i})
            trialLabels{i} = sprintf('%d', i);
        end
    end
    
    % Plot 1: Time to Food Discovery across trials
    subplot(3, 1, 1);
    % Find valid data points
    validTimeIdx = ~isnan(timeToFood);
    
    % Plot all points (including NaN as gaps)
    plot(trialNumbers, timeToFood, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0.8 0.2 0.8], 'MarkerFaceColor', [0.8 0.2 0.8], 'MarkerEdgeColor', 'k');
    
    xlabel('Trial Number');
    ylabel('Time to Food Discovery (s)');
    title('Time to Food Discovery Across Learning Trials');
    set(gca, 'XTick', 1:12, 'XTickLabel', trialLabels, 'XTickLabelRotation', 45);
    grid off;
    box off
    xlim([0.5, 12.5]);
    
    % Plot 2: Path Length across trials
    subplot(3, 1, 2);
    % Find valid data points
    validPathIdx = ~isnan(pathLengths);
    
    % Plot all points (including NaN as gaps)
    plot(trialNumbers, pathLengths, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0.3 0.6 0.9], 'MarkerFaceColor', [0.3 0.6 0.9], 'MarkerEdgeColor', 'k');
    
    xlabel('Trial Number');
    ylabel('Path Length (cm)');
    title('Path Length Across Learning Trials');
    set(gca, 'XTick', 1:12, 'XTickLabel', trialLabels, 'XTickLabelRotation', 45);
    grid off;
    box off
    xlim([0.5, 12.5]);
    
    % Plot 3: Mean Speed across trials
    subplot(3, 1, 3);
    % Find valid data points
    validSpeedIdx = ~isnan(meanSpeeds);
    
    % Plot all points (including NaN as gaps)
    plot(trialNumbers, meanSpeeds, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0.9 0.4 0.3], 'MarkerFaceColor', [0.9 0.4 0.3], 'MarkerEdgeColor', 'k');
    
    xlabel('Trial Number');
    ylabel('Mean Speed (cm/s)');
    title('Mean Speed Across Learning Trials');
    set(gca, 'XTick', 1:12, 'XTickLabel', trialLabels, 'XTickLabelRotation', 45);
    grid off;
    box off
    xlim([0.5, 12.5]);
    
    sgtitle(sprintf('%s Learning Curves', mouseName), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Print summary to command window
    fprintf('\n=== %s LEARNING PROGRESSION ===\n', mouseName);
    fprintf('Trial\tTime to Food (s)\tPath Length (cm)\tMean Speed (cm/s)\n');
    fprintf('-----\t---------------\t---------------\t-----------------\n');
    
    for i = 1:12
        if ~isnan(timeToFood(i)) && ~isnan(pathLengths(i)) && ~isnan(meanSpeeds(i))
            fprintf('%d\t%.1f\t\t%.1f\t\t%.2f\n', i, timeToFood(i), pathLengths(i), meanSpeeds(i));
        elseif ~isnan(timeToFood(i)) && ~isnan(pathLengths(i))
            fprintf('%d\t%.1f\t\t%.1f\t\tN/A\n', i, timeToFood(i), pathLengths(i));
        elseif ~isnan(timeToFood(i)) && ~isnan(meanSpeeds(i))
            fprintf('%d\t%.1f\t\tN/A\t\t%.2f\n', i, timeToFood(i), meanSpeeds(i));
        elseif ~isnan(pathLengths(i)) && ~isnan(meanSpeeds(i))
            fprintf('%d\tN/A\t\t%.1f\t\t%.2f\n', i, pathLengths(i), meanSpeeds(i));
        elseif ~isnan(timeToFood(i))
            fprintf('%d\t%.1f\t\tN/A\t\tN/A\n', i, timeToFood(i));
        elseif ~isnan(pathLengths(i))
            fprintf('%d\tN/A\t\t%.1f\t\tN/A\n', i, pathLengths(i));
        elseif ~isnan(meanSpeeds(i))
            fprintf('%d\tN/A\t\tN/A\t\t%.2f\n', i, meanSpeeds(i));
        else
            fprintf('%d\tN/A\t\tN/A\t\tN/A\n', i);
        end
    end
    
    % Calculate learning metrics
    validTimes = timeToFood(~isnan(timeToFood));
    validPaths = pathLengths(~isnan(pathLengths));
    validSpeeds = meanSpeeds(~isnan(meanSpeeds));
    
    if length(validTimes) > 1
        timeCorr = corrcoef(find(~isnan(timeToFood)), validTimes);
        fprintf('\nTime to food trend: r = %.3f (negative = learning)\n', timeCorr(1,2));
    end
    
    if length(validPaths) > 1
        pathCorr = corrcoef(find(~isnan(pathLengths)), validPaths);
        fprintf('Path length trend: r = %.3f (negative = learning)\n', pathCorr(1,2));
    end
    
    if length(validSpeeds) > 1
        speedCorr = corrcoef(find(~isnan(meanSpeeds)), validSpeeds);
        fprintf('Mean speed trend: r = %.3f\n', speedCorr(1,2));
    end
    
    fprintf('==============================\n\n');
end