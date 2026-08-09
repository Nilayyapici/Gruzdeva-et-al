%% Load data
clear all
cd C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Data_tables_all_tests\Cheeseboard
load('Cheeseboard.mat')

%% Reorganize the data

mice = reorganize_cheeseboard_data(Cheeseboard);

%data:
% data(:,1) - time
% data(:,[2,3]) - x,y
% data(:,4) - 465
% data(:,5) - 405
% data(:,6) - dff
% data(:,7) - speed
% data(:,8) - zones - 0-outside, 1-food, 2- area2, 3 - area3
% data(:,9) - distance
% data(:,10) - grooming

%% fixing the scale

old_cm_to_pix = 75/400;
pix_to_cm = 0.16; 
% Loop through each entry in mice_all
for i = 1:size(mice, 1)
    data_pre = mice{i, 3};
    data_test = mice{i, 4};
    data_pre(:,9) =  data_pre(:,9)./old_cm_to_pix.*pix_to_cm;
    data_test(:,9) =  data_test(:,9)./old_cm_to_pix.*pix_to_cm;
    mice{i, 3} = data_pre;
    mice{i, 4} = data_test;
end

%% fixing the scale

old_cm_to_pix = 75/400;
pix_to_cm = 0.16; % mean for 63.5/385 (horizontal) and 82.55/470 (vertical)
% Loop through each entry in mice_all
for i = 1:size(mice, 1)
     data_pre = mice{i, 3};
     data_test = mice{i, 4};
     data_pre(:,9) =  data_pre(:,9)./old_cm_to_pix.*pix_to_cm;
     data_test(:,9) =  data_test(:,9)./old_cm_to_pix.*pix_to_cm;
     mice{i, 3} = data_pre;
     mice{i, 4} = data_test;
end

% Loop through each entry in mice_all to recalculate speed
for i = 1:size(mice, 1)
    for j = 3:4
        data = mice{i, j};

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
            smoothed_speed = speed_calc;
        end

        % Match lengths (append first value at the beginning)
        speed_new = [smoothed_speed(1); smoothed_speed];

        % Replace the speed column (column 7) with recalculated values
        data(:, 7) = speed_new;

        % Update the data back into mice_all
        mice{i, j} = data;
    end
end

speed_threshold = 80; % cm/s

% Loop through each entry in mice_all to filter high speeds
for i = 1:size(mice, 1)
    for j = 3:4
        data = mice{i, j};
        % Extract speed column (column 7)
        speed = data(:, 7);

        % Find indices where speed is <= threshold
        valid_speed_indices = speed <= speed_threshold;

        % Filter data to keep only valid speed points
        data_filtered = data(valid_speed_indices, :);

        % Update the data back into mice_all
        mice{i, j} = data_filtered;
    end
end

%% Time in the zone by food
options = struct();
options.group = {'saline','CNO'};
options.time_limit = 10;   % first 5 minutes
options.plot_trajectories = true;

results = cheeseboard_behavioral_zones(mice, options);
mice = results.mice;

%% Time in zones

options = struct();
options.group = 'saline';
options.time_limit = 10;
options.figure_size = [10, 3];  % Width x Height in inches
options.font_size = 12;
options.save_figure = false;
results_saline = cheeseboard_zones_flexible(mice, options);

%% Latencies plot

options = struct();
options.group = 'saline';
options.time_limit = 10;
options.figure_size = [3,4];
options.discovery_distance = 2;
results_saline = cheeseboard_latency_flexible(mice, options);

%% Remove mice without signal

% mice MDRI17, MDRE6 from the all observation and looking at the raw data
% Initialize an array to store indices of rows to keep
keepIndices = true(size(mice, 1), 1);

for i = 1:size(mice, 1)
    if contains(mice{i, 1}, 'MDRE6', 'IgnoreCase', true)|| contains(mice{i, 1}, 'MDRI17', 'IgnoreCase', true)
        keepIndices(i) = false;  % Mark this row not to be kept
    end
end

% Create a new cell array excluding the specified rows
Mice_good = mice(keepIndices, :);
mice = Mice_good;

%% overall Correlation of dF/F with distance

options.group = {'saline'};
options.time_limit = 10;
options.speed_threshold =15;
plot_dff_distance_correlation_clean(mice, 'saline',options);

%% Trajectories of behavior

options = struct();
options.time_limit = 10;
options.speed_threshold = 0;
options.exclude_grooming = false;
trajectoryCheeseboard(mice, 'saline', options);

%% example of heatmaps and traces
% Plot an example mouse
mouse_name = 'M22_saline'; % Full mouse name from mice{i,1}

% Set options
options = struct();
options.dist_limit = 0;           % Minimum distance threshold
options.remove_grooming = true;   % Remove grooming periods
options.smooth_window = 30;       % Smoothing window size
options.sigma = 1;                % Gaussian smoothing for heatmap
options.resolution = 30;          % Heatmap resolution
options.time_xlim = [0, 525];     % Time axis limits (optional)
options.colormap = blueWhiteRed;  % Colormap for heatmap
options.dff_ylim = [-2, 4];  % dF/F range
options.dist_ylim = [0, 60];     % Distance range (cm)

options.group_filter = 'all';     % 'all', 'saline', or 'CNO'

plotMouseCheeseboardData(mice, mouse_name, options);

%% Rotation of the data and averaging across mice

visualizeOriginalMiceData(mice, 'saline');
mice_translated = translateMiceCenters(mice, 'saline', [350, 250]);
mice_rotated = rotateAlignedMiceVectors(mice_translated, 'saline', 0);

%% Averaged Occupancy

options = struct();
options.bin_size = 5;
options.exclude_grooming = true;
options.speed_threshold = 5;
options.colormap_name = slanCM('viola', 64);
plotAveragedOccupancyHeatmaps(mice_rotated, 'saline', options);

%% Averaged Activity

options = struct();
options.exclude_grooming = true;
options.resolution = 40;
options.speed_threshold = 5;
options.sigma = 2;
options.clim_activity = [-0.4 0.4];
options.clim_diff = [-0.4 0.4];
options.colormap_name = blueWhiteRed;

plotAveragedActivityHeatmaps(mice_rotated, 'saline', options);

% Smoothed heatmaps of activity

options.time_limit = 10;            % Analyze only first 10 minutes
options.exclude_grooming = false;   % Exclude grooming periods
options.speed_threshold = 0;        % Minimum speed threshold
options.bin_size = 5;               % Spatial bin size in cm (default: 4)
options.smooth_factor = 2.0;        % Gaussian smoothing factor (default: 2.0)
options.interpolate_display = true; % Interpolate to finer grid (default: true)
options.interp_factor = 4;          % Interpolation factor (default: 4)
options.colormap_name = blueWhiteRed;      % Colormap
options.clim_activity = [-0.4 0.4];    % Color limits for activity
options.clim_diff = [-0.35 0.35];        % Color limits for difference

plotAveragedActivityHeatmaps_smooth(mice_rotated, 'saline', options);

%% Runs towards and away

options = struct();
options.food_area = 6;              % 2 cm from food hole
options.group_filter = 'saline';       % Analyze all mice
options.session_filter = 'both';    % Both pre and test
options.exclude_grooming = true;    % Remove grooming
options.validate_plots = false;      % Show validation plots
options.threshold = 0;            % Distance velocity threshold
options.time_limit = 10;
options.memory = 'all';
options.colormap = blueWhiteRed;
options.caxis = [-0.8, 0.8];
options.plot_specific_mouse = 'FDRE14';  % Show plots only for this mouse

% Run analysis
run_data = analyze_cheeseboard_runs(Mice_good, options);

%% Correlations bar plots

options = struct();
results_saline = plot_runs_dff_distance_correlation(run_data,'saline');

%% Gradients

options = struct();
options.max_distance = 40; % Plot up to 35 cm from food
options.ylim = [-0.25 0.45]; % Y-axis limits for z-scored dF/F
options.smoothing = 15; % Smoothing window size
options.run_stats = true; % Include statistical tests
plot_cheeseboard_dff_distance_sessions(run_data, options);

%% Slopes
options = struct();
options.slope_range = 1000;        % Distance range for testing (±50 units)
options.max_distance = 35;     % Maximum distance to analyze
options.bin_width =0.01;
options.plot_results = true;    % Generate plots
options.plot_types = {'barplot'};

[slope_results] = analyze_cheeseboard_slopes(run_data, options);
visualize_cheeseboard_slopes(slope_results, options);

%% Compare speed between towards and away runs
% Works on run_data output from analyze_cheeseboard_runs

sessions   = {'pre', 'test'};
sess_colors = {[0.2 0.5 0.9], [0.9 0.3 0.2]};  % blue=pre, red=test

figure('Position', [100 100 900 500]);

for si = 1:2
    sess = sessions{si};

    towards_speed = [];
    away_speed    = [];
    towards_mouse = {};
    away_mouse    = {};

    for i = 1:length(run_data)
        if ~strcmp(run_data(i).session, sess), continue; end

        for j = 1:length(run_data(i).runs)
            run = run_data(i).runs(j);
            mean_spd = mean(run.speed, 'omitnan');
            if strcmp(run.type, 'towards')
                towards_speed(end+1) = mean_spd;
                towards_mouse{end+1} = run_data(i).mouse_id;
            else
                away_speed(end+1) = mean_spd;
                away_mouse{end+1}  = run_data(i).mouse_id;
            end
        end
    end

    % Per-animal means
    unique_ids   = unique([towards_mouse, away_mouse]);
    n_animals    = length(unique_ids);
    mean_towards = NaN(n_animals, 1);
    mean_away    = NaN(n_animals, 1);

    for a = 1:n_animals
        aid = unique_ids{a};
        t_idx  = strcmp(towards_mouse, aid);
        aw_idx = strcmp(away_mouse,    aid);
        if any(t_idx),  mean_towards(a) = mean(towards_speed(t_idx),  'omitnan'); end
        if any(aw_idx), mean_away(a)    = mean(away_speed(aw_idx),    'omitnan'); end
    end

    valid = ~isnan(mean_towards) & ~isnan(mean_away);

    % Stats
    [~, p, ~, stats] = ttest(mean_towards(valid), mean_away(valid));
    fprintf('\n--- %s session ---\n', sess);
    fprintf('Towards:  mean=%.3f  SEM=%.3f\n', ...
        mean(mean_towards(valid)), std(mean_towards(valid))/sqrt(sum(valid)));
    fprintf('Away:     mean=%.3f  SEM=%.3f\n', ...
        mean(mean_away(valid)),    std(mean_away(valid))/sqrt(sum(valid)));
    fprintf('Paired t-test: t(%d)=%.3f, p=%.4f\n', stats.df, stats.tstat, p);

    % Plot
    subplot(1, 2, si); hold on;
    col = sess_colors{si};

    % Individual lines
    for a = 1:n_animals
        if valid(a)
            plot([1 2], [mean_towards(a), mean_away(a)], ...
                'Color', [0.75 0.75 0.75], 'LineWidth', 1);
        end
    end

    % Bars
    bar_data = [mean(mean_towards(valid)), mean(mean_away(valid))];
    bar_sem  = [std(mean_towards(valid))/sqrt(sum(valid)), ...
        std(mean_away(valid))/sqrt(sum(valid))];
    b = bar(1:2, bar_data, 0.5);
    b.FaceColor = 'flat';
    b.CData     = [col * 0.7; col];   % slightly different shade per direction
    b.FaceAlpha = 0.6;
    errorbar(1:2, bar_data, bar_sem, 'k.', 'LineWidth', 1.5, 'CapSize', 8);

    % Individual dots
    scatter(1 + (rand(sum(valid),1)-0.5)*0.08, mean_towards(valid), ...
        40, col*0.7, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter(2 + (rand(sum(valid),1)-0.5)*0.08, mean_away(valid), ...
        40, col,     'filled', 'MarkerFaceAlpha', 0.8);

    % Significance bracket
    y_max = max([mean_towards(valid); mean_away(valid)]) * 1.18;
    if     p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = 'n.s.'; end
    plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1);
    text(1.5, y_max * 1.03, sprintf('%s  (p=%.3f)', sig, p), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');

    set(gca, 'XTick', 1:2, 'XTickLabel', {'Towards', 'Away'}, 'FontSize', 11);
    ylabel('Mean speed (cm/s)');
    title(sprintf('%s session  (n=%d animals)', sess, sum(valid)));
    box off; xlim([0.4 2.6]);
end

sgtitle('Run speed: towards vs away', 'FontSize', 13);