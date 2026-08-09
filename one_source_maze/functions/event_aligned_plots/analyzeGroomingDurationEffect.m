function analyzeGroomingDurationEffect(mice_all, options)
% ANALYZEGROOMINQDURATIONEFFECT Analyzes the correlation between grooming duration
% and the increase in neural activity from pre-grooming baseline
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters (see plotGroomingAlignedData)
%     Additional options:
%     - baseline_window: time window for baseline calculation in seconds [start end] (default: [-15 -5])
%     - response_window: time window for response calculation in seconds [start end] (default: [2 12])
%     - max_duration: maximum grooming duration to include in correlation (default: 30 seconds)

% Constants
COL_TIME = 1;     % Time column
COL_DFF = 11;     % DFF data column
COL_GROOM = 10;   % Grooming column (0/1)

% Set default options if not provided
if ~exist('options', 'var')
    options = struct();
end

if ~isfield(options, 'state')
    options.state = 'all';
end

if ~isfield(options, 'source')
    options.source = 'all';
end

if ~isfield(options, 'z_score')
    options.z_score = true;
end

if ~isfield(options, 'smooth_window')
    options.smooth_window = 15;
end

if ~isfield(options, 'min_grooming_duration')
    options.min_grooming_duration = 2; % Default 2 seconds
end

if ~isfield(options, 'max_duration')
    options.max_duration = 30; % Default 30 seconds (exclude longer events as outliers)
end

% Time windows for calculating baseline and response
if ~isfield(options, 'baseline_window')
    options.baseline_window = [-15, -5]; % 15-5 seconds before grooming onset
end

if ~isfield(options, 'response_window')
    options.response_window = [2, 12]; % 2-12 seconds after grooming onset
end

% Filter mice based on state and source
n = size(mice_all, 1);
mask = true(n, 1);

if ~strcmpi(options.state, 'all')
    mask = mask & strcmp(mice_all(:, 2), options.state);
end

if ~strcmpi(options.source, 'all')
    mask = mask & strcmp(mice_all(:, 3), options.source);
end

% Extract filtered mice
filtered_mice = mice_all(mask, :);
num_mice = size(filtered_mice, 1);

if num_mice == 0
    error('No mice found with the specified criteria');
end

% Initialize arrays to store results
all_grooming_durations = [];
all_dff_increases = [];
all_baseline_dffs = [];
all_response_dffs = [];
all_mouse_ids = {};

% Track outliers for reporting
outlier_durations = [];
outlier_increases = [];
outlier_mouse_ids = {};

% Process each mouse to find grooming events and calculate DFF changes
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    grooming = data(:, COL_GROOM);
    
    % Apply z-scoring if requested
    if options.z_score
        dff = (dff - nanmean(dff)) / nanstd(dff);
    end
    
    % Apply smoothing if requested
    if options.smooth_window > 0
        dff = movmean(dff, options.smooth_window, 'omitnan');
    end
    
    % Find onsets and offsets of grooming events
    groom_changes = diff([0; grooming; 0]); % Add padding and compute differences
    groom_onsets = find(groom_changes == 1);
    groom_offsets = find(groom_changes == -1) - 1;
    
    % Get mouse ID for labeling
    if ischar(filtered_mice{i, 1})
        id_parts = strsplit(filtered_mice{i, 1}, '_');
        if length(id_parts) > 1
            mouse_id = id_parts{1};
        else
            mouse_id = filtered_mice{i, 1};
        end
    else
        mouse_id = num2str(i);
    end
    
    % Calculate durations in seconds
    groom_durations = time(groom_offsets) - time(groom_onsets);
    
    % Filter for minimum duration
    valid_events = groom_durations >= options.min_grooming_duration;
    groom_onsets = groom_onsets(valid_events);
    groom_offsets = groom_offsets(valid_events);
    groom_durations = groom_durations(valid_events);
    
    num_events = length(groom_onsets);
    
    if num_events == 0
        fprintf('No valid grooming events found for mouse %s\n', mouse_id);
        continue;
    end
    
    % Process each grooming event
    for j = 1:num_events
        onset = groom_onsets(j);
        onset_time = time(onset);
        duration = groom_durations(j);
        
        % Calculate baseline DFF (before grooming)
        baseline_mask = time >= (onset_time + options.baseline_window(1)) & ...
                        time <= (onset_time + options.baseline_window(2));
        if sum(baseline_mask) > 0
            baseline_dff = nanmean(dff(baseline_mask));
        else
            continue; % Skip if no baseline period
        end
        
        % Calculate response DFF (during grooming)
        response_mask = time >= (onset_time + options.response_window(1)) & ...
                        time <= (onset_time + options.response_window(2));
        if sum(response_mask) > 0
            response_dff = nanmean(dff(response_mask));
        else
            continue; % Skip if no response period
        end
        
        % Calculate DFF increase
        dff_increase = response_dff - baseline_dff;
        
        % Store results, separating outliers from main dataset
        if duration > options.max_duration
            % Store as outlier
            outlier_durations = [outlier_durations; duration];
            outlier_increases = [outlier_increases; dff_increase];
            outlier_mouse_ids = [outlier_mouse_ids; {mouse_id}];
        else
            % Store in main dataset for correlation
            all_grooming_durations = [all_grooming_durations; duration];
            all_dff_increases = [all_dff_increases; dff_increase];
            all_baseline_dffs = [all_baseline_dffs; baseline_dff];
            all_response_dffs = [all_response_dffs; response_dff];
            all_mouse_ids = [all_mouse_ids; {mouse_id}];
        end
    end
end

% Check if we have enough data points
if isempty(all_grooming_durations)
    error('No valid grooming events found with the specified criteria');
end

% Report on outliers
num_outliers = length(outlier_durations);
if num_outliers > 0
    fprintf('Found %d grooming events longer than %d seconds (excluded from correlation)\n', ...
            num_outliers, options.max_duration);
end

% Compute correlation (using only non-outlier events)
[r, p] = corrcoef(all_grooming_durations, all_dff_increases);
correlation = r(1,2);
p_value = p(1,2);

% Create figure
figure('Position', [100, 100, 900, 700]);

% 1. Scatter plot with linear fit
subplot(2, 2, [1, 3]);
scatter(all_grooming_durations, all_dff_increases, 50, 'filled', 'MarkerFaceAlpha', 0.7);
hold on;

% Plot outliers with different styling if they exist
if ~isempty(outlier_durations)
    scatter(outlier_durations, outlier_increases, 50, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
end

% Add linear fit (using only non-outlier events)
p = polyfit(all_grooming_durations, all_dff_increases, 1);
x_fit = linspace(min(all_grooming_durations), max(all_grooming_durations), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Add zero reference line
yline(0, '--k', 'LineWidth', 1);

% Format plot
xlabel('Grooming Duration (seconds)', 'FontSize', 12);
ylabel('ΔF/F Increase', 'FontSize', 12);
if options.z_score
    ylabel('Z-scored ΔF/F Increase', 'FontSize', 12);
end
title({
    'Correlation between Grooming Duration and Neural Activity Increase',
    sprintf('r = %.3f, p = %.4f, n = %d events (outliers excluded)', correlation, p_value, length(all_grooming_durations))
}, 'FontSize', 14);
grid on;

% Add legend if outliers exist
if ~isempty(outlier_durations)
    legend({'Included events', sprintf('Outliers (>%ds)', options.max_duration), 'Linear fit'}, 'Location', 'best');
end

% Add text annotation explaining the calculation
text_x = min(all_grooming_durations) + 0.05 * (max(all_grooming_durations) - min(all_grooming_durations));
text_y = max(all_dff_increases) - 0.2 * (max(all_dff_increases) - min(all_dff_increases));
text_str = {
    sprintf('Baseline window: [%d, %d] sec', options.baseline_window(1), options.baseline_window(2)),
    sprintf('Response window: [%d, %d] sec', options.response_window(1), options.response_window(2)),
    'Increase = mean(response) - mean(baseline)',
};
text(text_x, text_y, text_str, 'FontSize', 10, 'BackgroundColor', [0.9, 0.9, 0.9, 0.5]);

% 2. Histogram of grooming durations
subplot(2, 2, 2);
histogram(all_grooming_durations, 'FaceColor', [0.4, 0.4, 0.8], 'EdgeColor', 'none');
hold on;
if ~isempty(outlier_durations)
    histogram(outlier_durations, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'none');
end
xlabel('Grooming Duration (seconds)', 'FontSize', 12);
ylabel('Count', 'FontSize', 12);
title('Distribution of Grooming Durations', 'FontSize', 12);
if ~isempty(outlier_durations)
    legend({'Included', 'Outliers'});
end
grid on;

% 3. Histogram of DFF increases
subplot(2, 2, 4);
histogram(all_dff_increases, 'FaceColor', [0.4, 0.6, 0.4], 'EdgeColor', 'none');
xlabel('ΔF/F Increase', 'FontSize', 12);
if options.z_score
    xlabel('Z-scored ΔF/F Increase', 'FontSize', 12);
end
ylabel('Count', 'FontSize', 12);
title('Distribution of Neural Activity Increases', 'FontSize', 12);
grid on;

% Add overall title
sgtitle({
    ['Relationship Between Grooming Duration and Neural Activity: ', ...
     options.state, ' mice, ', options.source, ' stimulus'],
}, 'FontWeight', 'bold');

% Print some statistics
fprintf('Analysis complete: %d grooming events analyzed (excluding %d outliers)\n', ...
        length(all_grooming_durations), num_outliers);
fprintf('Mean grooming duration: %.2f seconds (SD: %.2f)\n', mean(all_grooming_durations), std(all_grooming_durations));
fprintf('Mean ΔF/F increase: %.3f (SD: %.3f)\n', mean(all_dff_increases), std(all_dff_increases));
fprintf('Correlation: r = %.3f, p = %.4f\n', correlation, p_value);

% Calculate quartiles of grooming duration
duration_quartiles = quantile(all_grooming_durations, [0.25, 0.5, 0.75]);
fprintf('Grooming duration quartiles: Q1=%.1f, Median=%.1f, Q3=%.1f seconds\n', ...
        duration_quartiles(1), duration_quartiles(2), duration_quartiles(3));

% Optional: return the analysis results if requested
if nargout > 0
    results = struct();
    results.grooming_durations = all_grooming_durations;
    results.dff_increases = all_dff_increases;
    results.baseline_dffs = all_baseline_dffs;
    results.response_dffs = all_response_dffs;
    results.mouse_ids = all_mouse_ids;
    results.correlation = correlation;
    results.p_value = p_value;
    results.fit_params = p;
    results.outliers.durations = outlier_durations;
    results.outliers.increases = outlier_increases;
    results.outliers.mouse_ids = outlier_mouse_ids;
    varargout{1} = results;
end
end