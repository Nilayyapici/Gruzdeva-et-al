function plotMouseTrajectories_allTrials(allMiceData, options)
% Plot a figure per mouse showing trajectories for trials 1-N in a grid layout.
%
% Usage:
% options.num_trials = 10;               % Number of trials to plot (default: 10)
% options.speed_threshold = 50;          % Exclude points where speed > threshold cm/s (default: 50)
% options.min_speed_threshold = 5;       % Exclude points where speed < threshold cm/s (default: 5)
% options.exclude_first_n_frames = 20;    % Exclude first N frames from each trial (default: 20)
% options.food_discovery_threshold = 1.2;% Distance threshold for food discovery in cm (default: 1.2)
% options.xlim = [10, 90];              % X-axis limits in cm (default: [10, 90])
% options.ylim = [0, 80];               % Y-axis limits in cm (default: [0, 80])
% options.color_by = 'trial';           % Color trajectories by: 'trial' or 'dff' (default: 'trial')
% options.dff_ylim = [];                % Color limits for dF/F coloring (default: auto)
% options.rows = 2;                     % Number of subplot rows (default: 2)
%
% plotMouseTrajectories_allTrials(allMiceData, options)

if nargin < 2
    options = struct();
end

% Set default options
if ~isfield(options, 'num_trials'),              options.num_trials = 10; end
if ~isfield(options, 'speed_threshold'),         options.speed_threshold = 50; end
if ~isfield(options, 'min_speed_threshold'),     options.min_speed_threshold = 5; end
if ~isfield(options, 'exclude_first_n_frames'),  options.exclude_first_n_frames = 5; end
if ~isfield(options, 'food_discovery_threshold'),options.food_discovery_threshold = 1; end
if ~isfield(options, 'xlim'),                    options.xlim = [15, 95]; end
if ~isfield(options, 'ylim'),                    options.ylim = [0, 80]; end
if ~isfield(options, 'color_by'),                options.color_by = 'dff'; end
if ~isfield(options, 'dff_ylim'),                options.dff_ylim = []; end
if ~isfield(options, 'rows'),                    options.rows = 2; end

% Validate
if ~ismember(options.color_by, {'trial', 'dff'})
    error('color_by must be "trial" or "dff"');
end

pix_to_cm = 0.16;
maxTrials = options.num_trials;
cols      = ceil(maxTrials / options.rows);

% Colors for 'trial' mode - one per trial
trialColors = lines(maxTrials);

fprintf('Plotting trajectories for %d mice, %d trials each...\n', length(allMiceData), maxTrials);

for mouseIdx = 1:length(allMiceData)
    mouseData    = allMiceData{mouseIdx};
    learningData = mouseData.learningData;

    fprintf('  Mouse %s...\n', mouseData.name);

    % Collect all trial data indexed by trial number
    trialSlots = cell(maxTrials, 1);  % trialSlots{t} = struct with data and foodCoords

    for trialIdx = 1:size(learningData, 1)
        trialName  = learningData{trialIdx, 1};
        data       = learningData{trialIdx, 3};
        foodCoords = learningData{trialIdx, 2};

        if ~isempty(data) && ~isempty(foodCoords)
            numbers = regexp(trialName, '([A-Z])(\d+)(?!.*\d)', 'tokens');
            if ~isempty(numbers)
                trialNum = str2double(numbers{end}{2});
                if trialNum >= 1 && trialNum <= maxTrials
                    trialSlots{trialNum}.data       = data;
                    trialSlots{trialNum}.foodCoords = foodCoords;
                    trialSlots{trialNum}.name       = trialName;
                end
            end
        end
    end

    % Compute global dF/F color limits for this mouse (if color_by = 'dff')
    if strcmp(options.color_by, 'dff')
        if ~isempty(options.dff_ylim)
            dff_clim = options.dff_ylim;
        else
            all_dff = [];
            for t = 1:maxTrials
                if ~isempty(trialSlots{t})
                    d = applyFrameAndSpeedFilter(trialSlots{t}.data, options);
                    dff_vals = d(:, 6);
                    all_dff  = [all_dff; dff_vals(~isnan(dff_vals))];
                end
            end
            if ~isempty(all_dff)
                dff_clim = [prctile(all_dff, 5), prctile(all_dff, 95)];
            else
                dff_clim = [-2, 2];
            end
        end
    end

    % Create figure for this mouse
    fig = figure('Name', sprintf('%s: Trajectories Trials 1-%d', mouseData.name, maxTrials), ...
                 'Position', [50, 50, 220*cols, 220*options.rows]);

    for t = 1:maxTrials
        subplot(options.rows, cols, t);

        if isempty(trialSlots{t})
            % No data for this trial
            text(0.5, 0.5, sprintf('Trial %d\nNo data', t), ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
                 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
            axis off;
            continue;
        end

        % Apply frame and speed filter
        data_filt  = applyFrameAndSpeedFilter(trialSlots{t}.data, options);
        foodCoords = trialSlots{t}.foodCoords;

        x_coords  = data_filt(:, 2);
        y_coords  = data_filt(:, 3);
        dff_vals  = data_filt(:, 6);
        valid_idx = ~isnan(x_coords) & ~isnan(y_coords);
        x_clean   = x_coords(valid_idx) * pix_to_cm;
        y_clean   = y_coords(valid_idx) * pix_to_cm;
        dff_clean = dff_vals(valid_idx);

        if length(x_clean) < 2
            text(0.5, 0.5, sprintf('Trial %d\nInsufficient data', t), ...
                 'HorizontalAlignment', 'center', 'Units', 'normalized', ...
                 'FontSize', 9, 'Color', [0.5 0.5 0.5]);
            axis off;
            continue;
        end

        foodX_cm = foodCoords(1) * pix_to_cm;
        foodY_cm = foodCoords(2) * pix_to_cm;

        % Find food discovery point
        dist_to_food     = sqrt((x_clean - foodX_cm).^2 + (y_clean - foodY_cm).^2);
        foodDiscoveryIdx = find(dist_to_food < options.food_discovery_threshold, 1, 'first');

        if ~isempty(foodDiscoveryIdx)
            x_plot   = x_clean(1:foodDiscoveryIdx);
            y_plot   = y_clean(1:foodDiscoveryIdx);
            dff_plot = dff_clean(1:foodDiscoveryIdx);
            discovered = true;
        else
            x_plot   = x_clean;
            y_plot   = y_clean;
            dff_plot = dff_clean;
            discovered = false;
        end

        hold on;

        % Draw arena boundary
        x0_cm          = 336 * pix_to_cm;
        y0_cm          = 262 * pix_to_cm;
        arena_r_cm     = 227 * pix_to_cm;
        theta          = linspace(0, 2*pi, 360);
        plot(x0_cm + arena_r_cm * cos(theta), y0_cm + arena_r_cm * sin(theta), ...
             'k-', 'LineWidth', 1.2);

        if strcmp(options.color_by, 'dff') && ~all(isnan(dff_plot))
            % Color by dF/F value
            scatter(x_plot, y_plot, 8, dff_plot, 'filled');
            colormap(blueWhiteRed(256));
            clim(dff_clim);
        else
            % Color by trial number (solid line)
            plot(x_plot, y_plot, '-', 'Color', trialColors(t, :), 'LineWidth', 1.2);
        end

        % Start marker
        plot(x_plot(1), y_plot(1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.2 0.7 0.2], ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

        % Food location marker
        plot(foodX_cm, foodY_cm, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow', ...
             'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

        % End marker - different if food found vs not
        if discovered
            plot(x_plot(end), y_plot(end), 's', 'MarkerSize', 6, 'MarkerFaceColor', 'yellow', ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        else
            plot(x_plot(end), y_plot(end), 's', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.2 0.2], ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        end

        hold off;

        % Axes formatting
        xlim(options.xlim);
        ylim(options.ylim);
        xlabel('X (cm)', 'FontSize', 7);
        ylabel('Y (cm)', 'FontSize', 7);
        set(gca, 'FontSize', 7);

        if discovered
            title(sprintf('Trial %d ✓', t), 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1]);
        else
            title(sprintf('Trial %d ✗', t), 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.6 0.1 0.1]);
        end

        grid off;
        box on;
        axis off;
    end

    % Add colorbar for dF/F mode
    if strcmp(options.color_by, 'dff')
        cb = colorbar('Position', [0.93, 0.1, 0.015, 0.8]);
        cb.Label.String = 'dF/F (%)';
        cb.FontSize = 8;
    end

    sgtitle(sprintf('%s  —  Trials 1–%d', mouseData.name, maxTrials), ...
            'FontSize', 13, 'FontWeight', 'bold');

    % Small legend annotation
    annotation(fig, 'textbox', [0.01, 0.01, 0.4, 0.04], ...
               'String', 'Circle = start  |  Star = food  |  Square (green) = found  |  Square (red) = not found', ...
               'EdgeColor', 'none', 'FontSize', 6, 'FitBoxToText', 'on');

    fprintf('    Done (%d/%d trials with data)\n', ...
            sum(~cellfun(@isempty, trialSlots)), maxTrials);
end

fprintf('Done. One figure created per mouse.\n');
end

% =========================================================================
%  SHARED FRAME/SPEED FILTER  (same as in main analysis function)
% =========================================================================

function data = applyFrameAndSpeedFilter(data, options)
n = options.exclude_first_n_frames;
if n > 0 && size(data, 1) > n
    data = data(n+1:end, :);
end
speed_vals = data(:, 7);
keep = isnan(speed_vals) ...
     | (speed_vals <= options.speed_threshold & speed_vals >= options.min_speed_threshold);
data = data(keep, :);
end