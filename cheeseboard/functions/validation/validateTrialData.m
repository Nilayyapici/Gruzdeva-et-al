function validateTrialData(allMiceData, options)
    % Validate specific trial data by plotting individual trajectories and averaged dF/F
    %
    % Usage: 
    % options.trial_number = 1;               % Which trial to analyze (default: 1)
    % options.use_zscore = false;             % Use z-scored dF/F (default: false)
    % options.distance_limit = 20;            % Maximum distance to analyze in cm (default: 20)
    % options.bin_size = 0.5;                 % Bin size for distance binning in cm (default: 0.5)
    % options.dff_ylim = [];                  % Y-axis limits for dF/F plots (default: auto)
    % options.food_discovery_threshold = 1.2; % Distance threshold for food discovery in cm (default: 1.2)
    % options.speed_threshold = 50;           % Exclude points where speed > threshold in cm/s (default: 50)
    % options.exclude_first_n_frames = 0;     % Exclude first N frames from each trial trajectory (default: 0)
    %
    % validateTrialData(allMiceData, options)
    %
    % This function creates detailed plots for the specified trial to check data quality:
    % 1. All mouse trajectories overlaid on the same axes (color-coded by mouse)
    % 2. Individual mouse trajectories with dF/F color-coding
    % 3. Individual dF/F vs distance plots for each mouse
    % 4. Averaged dF/F vs distance plot with all mice combined

    if nargin < 2
        options = struct();
    end

    % Set default options
    if ~isfield(options, 'trial_number'), options.trial_number = 1; end
    if ~isfield(options, 'use_zscore'), options.use_zscore = false; end
    if ~isfield(options, 'distance_limit'), options.distance_limit = 20; end
    if ~isfield(options, 'bin_size'), options.bin_size = 0.5; end
    if ~isfield(options, 'dff_ylim'), options.dff_ylim = []; end
    if ~isfield(options, 'food_discovery_threshold'), options.food_discovery_threshold = 1.2; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 50; end
    if ~isfield(options, 'exclude_first_n_frames'), options.exclude_first_n_frames = 0; end

    % Validate options
    if options.trial_number < 1 || options.trial_number > 12
        error('Trial number must be between 1 and 12');
    end
    if options.exclude_first_n_frames < 0
        error('exclude_first_n_frames must be >= 0');
    end

    fprintf('\n=== VALIDATING TRIAL %d DATA ===\n', options.trial_number);
    fprintf('Options:\n');
    fprintf('  Trial number: %d\n', options.trial_number);
    fprintf('  Z-scored dF/F: %s\n', mat2str(options.use_zscore));
    fprintf('  Distance limit: %.1f cm\n', options.distance_limit);
    fprintf('  Bin size: %.1f cm\n', options.bin_size);
    fprintf('  Food discovery threshold: %.1f cm\n', options.food_discovery_threshold);
    fprintf('  Speed threshold: %.1f cm/s\n', options.speed_threshold);
    fprintf('  Exclude first N frames: %d\n', options.exclude_first_n_frames);
    fprintf('  pix_to_cm: 0.16\n');
    if ~isempty(options.dff_ylim)
        fprintf('  dF/F Y-limits: [%.1f, %.1f]\n', options.dff_ylim(1), options.dff_ylim(2));
    else
        fprintf('  dF/F Y-limits: Auto\n');
    end

    % Extract specified trial data for all mice
    trial_data = [];
    trial_mice_names = {};

    for mouseIdx = 1:length(allMiceData)
        mouseData    = allMiceData{mouseIdx};
        learningData = mouseData.learningData;

        for trialIdx = 1:size(learningData, 1)
            trialName  = learningData{trialIdx, 1};
            data       = learningData{trialIdx, 3};
            foodCoords = learningData{trialIdx, 2};

            if ~isempty(data) && ~isempty(foodCoords)
                numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
                if ~isempty(numbers)
                    trialNum = str2double(numbers{end}{2});

                    if trialNum == options.trial_number
                        trial_entry = struct();
                        trial_entry.mouseName  = mouseData.name;
                        trial_entry.data       = data;
                        trial_entry.foodCoords = foodCoords;
                        trial_entry.trialName  = trialName;

                        trial_data = [trial_data; trial_entry];
                        trial_mice_names{end+1} = mouseData.name;

                        fprintf('  Found Trial %d for %s: %d data points\n', ...
                                options.trial_number, mouseData.name, size(data, 1));
                        break;
                    end
                end
            end
        end
    end

    if isempty(trial_data)
        error('No Trial %d data found! Available trials may be different.', options.trial_number);
    end

    fprintf('Total mice with Trial %d data: %d\n', options.trial_number, length(trial_data));

    % Z-score data if requested
    if options.use_zscore
        fprintf('Z-scoring dF/F data...\n');
        trial_data = zscore_trial_data(trial_data);
    end

    % Figure 0: All trajectories overlaid on same axes
    createOverlaidTrajectoryPlot(trial_data, options);

    % Figures 1 & 2: Individual trajectory and dF/F plots
    createValidationPlots(trial_data, options);

    % Figure 3: Averaged analysis
    createAveragedAnalysis(trial_data, options);
end

% =========================================================================
%  SHARED FRAME/SPEED FILTER HELPER
% =========================================================================

function data = applyFrameAndSpeedFilter(data, options)
    % Exclude first N frames
    n = options.exclude_first_n_frames;
    if n > 0 && size(data, 1) > n
        data = data(n+1:end, :);
    end
    % Exclude rows where speed > threshold (column 7); keep NaN speed rows
    speed_vals = data(:, 7);
    keep = isnan(speed_vals) | speed_vals <= options.speed_threshold;
    data = data(keep, :);
end

% =========================================================================
%  ALL TRAJECTORIES OVERLAID
% =========================================================================

function createOverlaidTrajectoryPlot(trial_data, options)
    pix_to_cm = 0.16;
    nMice     = length(trial_data);
    trialNum  = options.trial_number;
    colors    = lines(nMice);

    figure('Name', sprintf('Trial %d: All Trajectories Overlaid', trialNum), ...
           'Position', [50, 50, 800, 700]);
    hold on;

    legendHandles = gobjects(nMice, 1);

    for i = 1:nMice
        mouse_data = trial_data(i);

        % Apply frame and speed filters
        data_filt  = applyFrameAndSpeedFilter(mouse_data.data, options);

        x_coords  = data_filt(:, 2);
        y_coords  = data_filt(:, 3);
        valid_idx = ~isnan(x_coords) & ~isnan(y_coords);
        x_clean   = x_coords(valid_idx) * pix_to_cm;
        y_clean   = y_coords(valid_idx) * pix_to_cm;

        if length(x_clean) < 2, continue; end

        foodX_cm     = mouse_data.foodCoords(1) * pix_to_cm;
        foodY_cm     = mouse_data.foodCoords(2) * pix_to_cm;
        dist_to_food = sqrt((x_clean - foodX_cm).^2 + (y_clean - foodY_cm).^2);
        discoveryIdx = find(dist_to_food < options.food_discovery_threshold, 1, 'first');

        if ~isempty(discoveryIdx)
            x_plot = x_clean(1:discoveryIdx);
            y_plot = y_clean(1:discoveryIdx);
        else
            x_plot = x_clean;
            y_plot = y_clean;
        end

        h = plot(x_plot, y_plot, '-', 'Color', colors(i, :), 'LineWidth', 1.5, ...
                 'DisplayName', mouse_data.mouseName);
        legendHandles(i) = h;

        plot(x_plot(1),   y_plot(1),   'o', 'Color', colors(i, :), ...
             'MarkerSize', 8, 'MarkerFaceColor', colors(i, :), 'HandleVisibility', 'off');
        plot(x_plot(end), y_plot(end), 's', 'Color', colors(i, :), ...
             'MarkerSize', 8, 'MarkerFaceColor', colors(i, :), 'HandleVisibility', 'off');
        plot(foodX_cm, foodY_cm, 'p', 'Color', colors(i, :), ...
             'MarkerSize', 14, 'MarkerFaceColor', 'yellow', 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end

    hold off;
    axis equal;
    xlabel('X (cm)');
    ylabel('Y (cm)');
    legend(legendHandles, 'Location', 'eastoutside', 'Box', 'off');
    title(sprintf('Trial %d: All Mouse Trajectories Overlaid', trialNum), ...
          'FontSize', 14, 'FontWeight', 'bold');
    annotation('textbox', [0.01 0.01 0.35 0.07], 'String', ...
               'Circle = start  |  Square = end  |  Star = food location', ...
               'EdgeColor', 'none', 'FontSize', 8, 'FitBoxToText', 'on');
    grid off;
    fprintf('  Overlaid trajectory plot created for %d mice\n', nMice);
end

% =========================================================================
%  INDIVIDUAL VALIDATION PLOTS
% =========================================================================

function zscoredData = zscore_trial_data(trial_data)
    zscoredData = trial_data;
    for i = 1:length(trial_data)
        data      = trial_data(i).data;
        dff_vals  = data(:, 6);
        valid_dff = dff_vals(~isnan(dff_vals));
        if ~isempty(valid_dff)
            dff_mean = mean(valid_dff);
            dff_std  = std(valid_dff);
            zscoredData(i).data(:, 6) = (data(:, 6) - dff_mean) / dff_std;
            fprintf('    %s: Mean=%.3f, STD=%.3f\n', trial_data(i).mouseName, dff_mean, dff_std);
        end
    end
end

function createValidationPlots(trial_data, options)
    nMice    = length(trial_data);
    trialNum = options.trial_number;
    pix_to_cm = 0.16;

    if nMice <= 6,       rows = 2; cols = 3;
    elseif nMice <= 9,   rows = 3; cols = 3;
    elseif nMice <= 12,  rows = 3; cols = 4;
    elseif nMice <= 16,  rows = 4; cols = 4;
    else,                rows = 4; cols = 5;
    end

    % Global dF/F limits
    all_dff = [];
    for i = 1:nMice
        [dffValues, ~] = extractTrajectoryToFood_validation(trial_data(i), options);
        if ~isempty(dffValues), all_dff = [all_dff; dffValues]; end
    end

    if ~isempty(options.dff_ylim)
        dff_clim = options.dff_ylim;
    elseif ~isempty(all_dff)
        dff_clim = [prctile(all_dff, 5), prctile(all_dff, 95)];
    else
        dff_clim = [-2, 2];
    end
    fprintf('Using dF/F color limits: [%.2f, %.2f]\n', dff_clim(1), dff_clim(2));

    % --- Figure 1: Trajectory plots with dF/F color-coding ---------------
    figure('Name', sprintf('Trial %d: Mouse Trajectories with dF/F', trialNum), ...
           'Position', [50, 50, 1400, 1000]);

    for i = 1:min(nMice, rows*cols)
        subplot(rows, cols, i);

        mouse_data = trial_data(i);
        [dffValues, ~] = extractTrajectoryToFood_validation(mouse_data, options);

        if ~isempty(dffValues)
            % Apply frame and speed filters
            data_filt  = applyFrameAndSpeedFilter(mouse_data.data, options);

            x_coords   = data_filt(:, 2);
            y_coords   = data_filt(:, 3);
            valid_idx  = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(data_filt(:, 6));
            x_clean    = x_coords(valid_idx) * pix_to_cm;
            y_clean    = y_coords(valid_idx) * pix_to_cm;
            dff_clean  = data_filt(valid_idx, 6);

            foodX_cm     = mouse_data.foodCoords(1) * pix_to_cm;
            foodY_cm     = mouse_data.foodCoords(2) * pix_to_cm;
            dist_to_food = sqrt((x_clean - foodX_cm).^2 + (y_clean - foodY_cm).^2);
            foodDiscoveryIdx = find(dist_to_food < options.food_discovery_threshold, 1, 'first');

            if ~isempty(foodDiscoveryIdx)
                x_plot   = x_clean(1:foodDiscoveryIdx);
                y_plot   = y_clean(1:foodDiscoveryIdx);
                dff_plot = dff_clean(1:foodDiscoveryIdx);
                title_suffix = ' (until food discovery)';
            else
                withinLimit = dist_to_food <= options.distance_limit;
                x_plot   = x_clean(withinLimit);
                y_plot   = y_clean(withinLimit);
                dff_plot = dff_clean(withinLimit);
                title_suffix = ' (no food discovery)';
            end

            scatter(x_plot, y_plot, 20, dff_plot, 'filled');
            colormap(blueWhiteRed(256));
            caxis(dff_clim);

            hold on;
            plot(foodX_cm, foodY_cm, 'ks', 'MarkerSize', 12, 'LineWidth', 3, 'MarkerFaceColor', 'yellow');
            if ~isempty(x_plot)
                plot(x_plot(1), y_plot(1), 'go', 'MarkerSize', 8, 'LineWidth', 2);
            end
            plot(x_plot, y_plot, 'k-', 'LineWidth', 0.5);
            hold off;

            title(sprintf('%s%s', mouse_data.mouseName, title_suffix));
            xlabel('X (cm)'); ylabel('Y (cm)');
            xlim([10, 90]);
            ylim([0, 80]);
            fprintf('  %s: %d points plotted\n', mouse_data.mouseName, length(x_plot));
        else
            text(0.5, 0.5, sprintf('%s\nNo valid data', mouse_data.mouseName), ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized');
            axis off;
        end
    end

    c = colorbar('Position', [0.92, 0.15, 0.02, 0.7]);
    if options.use_zscore
        c.Label.String = 'Z-scored dF/F';
    else
        c.Label.String = 'dF/F (%)';
    end
    sgtitle(sprintf('Trial %d: Individual Mouse Trajectories with dF/F Color-Coding', trialNum), ...
            'FontSize', 16, 'FontWeight', 'bold');

    % --- Figure 2: Individual dF/F vs distance plots ----------------------
    figure('Name', sprintf('Trial %d: Individual dF/F vs Distance', trialNum), ...
           'Position', [100, 100, 1400, 1000]);

    for i = 1:min(nMice, rows*cols)
        subplot(rows, cols, i);

        mouse_data = trial_data(i);
        [dffValues, distances] = extractTrajectoryToFood_validation(mouse_data, options);

        if ~isempty(dffValues) && ~isempty(distances)
            scatter(distances, dffValues, 15, 'b', 'filled');
            hold on;

            if length(distances) > 10
                try
                    [dist_sorted, sort_idx] = sort(distances);
                    dff_sorted  = dffValues(sort_idx);
                    window_size = min(20, floor(length(dist_sorted)/5));
                    if window_size >= 3
                        dff_smooth = movmean(dff_sorted, window_size);
                        plot(dist_sorted, dff_smooth, 'r-', 'LineWidth', 2);
                    end
                catch
                end
            end

            plot([0, max(distances)], [0, 0], 'k--', 'LineWidth', 1);
            hold off;

            xlabel('Distance to Food (cm)');
            if options.use_zscore, ylabel('Z-scored dF/F'); else, ylabel('dF/F (%)'); end
            title(sprintf('%s (n=%d)', mouse_data.mouseName, length(dffValues)));
            grid on;
            xlim([0, options.distance_limit]);
            ylim(dff_clim);
        else
            text(0.5, 0.5, sprintf('%s\nNo valid data', mouse_data.mouseName), ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized');
            axis off;
        end
    end

    sgtitle(sprintf('Trial %d: Individual Mouse dF/F vs Distance to Food', trialNum), ...
            'FontSize', 16, 'FontWeight', 'bold');
end

% =========================================================================
%  AVERAGED ANALYSIS
% =========================================================================

function createAveragedAnalysis(trial_data, options)
    trialNum      = options.trial_number;
    all_distances = [];
    all_dff       = [];
    mouse_labels  = [];

    for i = 1:length(trial_data)
        [dffValues, distances] = extractTrajectoryToFood_validation(trial_data(i), options);
        if ~isempty(dffValues) && ~isempty(distances)
            all_distances = [all_distances; distances];
            all_dff       = [all_dff;       dffValues];
            mouse_labels  = [mouse_labels;  repmat(i, length(distances), 1)];
        end
    end

    if isempty(all_distances)
        fprintf('No valid data for averaged analysis of Trial %d\n', trialNum);
        return;
    end

    distanceBins = 0:options.bin_size:options.distance_limit;
    binCenters   = distanceBins(1:end-1) + options.bin_size/2;
    [meanDFF, semDFF, counts] = binDFFByDistance_validation(all_distances, all_dff, distanceBins);

    figure('Name', sprintf('Trial %d: Averaged dF/F vs Distance Analysis', trialNum), ...
           'Position', [150, 150, 1200, 800]);

    % Plot 1: Raw data scatter
    subplot(2, 2, 1);
    colors = lines(length(trial_data));
    for i = 1:length(trial_data)
        mouse_mask = mouse_labels == i;
        if sum(mouse_mask) > 0
            scatter(all_distances(mouse_mask), all_dff(mouse_mask), 10, colors(i, :), 'filled', ...
                    'DisplayName', trial_data(i).mouseName);
            hold on;
        end
    end
    xlabel('Distance to Food (cm)');
    if options.use_zscore, ylabel('Z-scored dF/F'); else, ylabel('dF/F (%)'); end
    title('Raw Data: All Mice Combined');
    legend('Location', 'eastoutside');
    grid on; xlim([0, options.distance_limit]);
    hold off;

    % Plot 2: Binned averages with SEM
    subplot(2, 2, 2);
    validBins = counts >= 2 & ~isnan(meanDFF);
    if sum(validBins) > 0
        x_valid    = binCenters(validBins);
        mean_valid = meanDFF(validBins);
        sem_valid  = semDFF(validBins);
        fill([x_valid, fliplr(x_valid)], [mean_valid + sem_valid, fliplr(mean_valid - sem_valid)], ...
             'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;
        plot(x_valid, mean_valid, 'b-', 'LineWidth', 3);
        plot([0, options.distance_limit], [0, 0], 'k--', 'LineWidth', 1);
        hold off;
    end
    xlabel('Distance to Food (cm)');
    if options.use_zscore, ylabel('Z-scored dF/F'); else, ylabel('dF/F (%)'); end
    title(sprintf('Averaged dF/F (%.1f cm bins)', options.bin_size));
    grid on; xlim([0, options.distance_limit]);

    % Plot 3: Data point counts per bin
    subplot(2, 2, 3);
    bar(binCenters, counts, 'FaceColor', [0.7 0.7 0.7]);
    xlabel('Distance to Food (cm)');
    ylabel('Number of Data Points');
    title('Data Points per Distance Bin');
    grid on; xlim([0, options.distance_limit]);

    % Plot 4: Statistics summary
    subplot(2, 2, 4);
    axis off;
    totalPoints       = length(all_distances);
    nMice             = length(trial_data);
    avgPointsPerMouse = totalPoints / nMice;
    validBinCount     = sum(validBins);

    stats_text = {
        sprintf('TRIAL %d STATISTICS:', trialNum);
        '';
        sprintf('Total mice: %d', nMice);
        sprintf('Total data points: %d', totalPoints);
        sprintf('Avg points per mouse: %.1f', avgPointsPerMouse);
        sprintf('Distance range: 0-%.1f cm', options.distance_limit);
        sprintf('Bin size: %.1f cm', options.bin_size);
        sprintf('Valid bins (>=2 points): %d/%d', validBinCount, length(binCenters));
        '';
        'dF/F Statistics:';
        sprintf('Mean: %.3f', nanmean(all_dff));
        sprintf('STD: %.3f', nanstd(all_dff));
        sprintf('Range: %.3f to %.3f', min(all_dff), max(all_dff));
        '';
        'Distance Statistics:';
        sprintf('Mean: %.2f cm', nanmean(all_distances));
        sprintf('STD: %.2f cm', nanstd(all_distances));
        sprintf('Range: %.2f to %.2f cm', min(all_distances), max(all_distances));
        '';
        'Analysis Options:';
        sprintf('Z-scored: %s', mat2str(options.use_zscore));
        sprintf('Food threshold: %.1f cm', options.food_discovery_threshold);
        sprintf('Speed threshold: %.1f cm/s', options.speed_threshold);
        sprintf('Excluded first frames: %d', options.exclude_first_n_frames);
        sprintf('pix_to_cm: 0.16');
    };

    text(0.1, 0.9, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top', 'Units', 'normalized');
    sgtitle(sprintf('Trial %d: Comprehensive Analysis Summary', trialNum), 'FontSize', 16, 'FontWeight', 'bold');

    fprintf('\n=== TRIAL %d VALIDATION SUMMARY ===\n', trialNum);
    fprintf('Mice analyzed: %d\n', nMice);
    fprintf('Total data points: %d\n', totalPoints);
    fprintf('Average points per mouse: %.1f\n', avgPointsPerMouse);
    fprintf('Valid distance bins: %d/%d\n', validBinCount, length(binCenters));
    fprintf('dF/F range: %.3f to %.3f\n', min(all_dff), max(all_dff));
    fprintf('Distance range: %.2f to %.2f cm\n', min(all_distances), max(all_distances));
    fprintf('===================================\n');
end

% =========================================================================
%  TRAJECTORY EXTRACTION (uses frame + speed filters)
% =========================================================================

function [dffValues, distances] = extractTrajectoryToFood_validation(mouse_data, options)
    dffValues  = [];
    distances  = [];
    pix_to_cm  = 0.16;

    % Apply frame and speed filters first
    data = applyFrameAndSpeedFilter(mouse_data.data, options);

    foodCoords = mouse_data.foodCoords;
    x_coords  = data(:, 2);
    y_coords  = data(:, 3);
    dff_vals  = data(:, 6);
    valid_idx = ~isnan(x_coords) & ~isnan(y_coords) & ~isnan(dff_vals);
    x_clean   = x_coords(valid_idx);
    y_clean   = y_coords(valid_idx);
    dff_clean = dff_vals(valid_idx);

    if length(x_clean) < 2 || length(foodCoords) < 2, return; end

    dist_to_food = sqrt((x_clean - foodCoords(1)).^2 + (y_clean - foodCoords(2)).^2) * pix_to_cm;
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

% =========================================================================
%  BINNING
% =========================================================================

function [meanDFF, semDFF, counts] = binDFFByDistance_validation(distances, dffValues, distanceBins)
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