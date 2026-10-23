function [varargout] = plotEventDurationCorrelation(mice_all, options)
% PLOTEVENTDURATIONCORRELATION Plots correlation between event duration and maximum DFF change
%
% Parameters:
%   mice_all - Cell array with mouse data as described in the spec
%   options - Struct with visualization parameters:
%     - state: 'fed', 'fasted', or 'all' (default: 'all')
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - event_type: 'visit' or 'eating' (default: 'visit')
%     - z_score: whether to z-score the DFF signals (default: false)
%     - plot_individual_mice: whether to plot individual mice (default: true)
%     - ylim: y-axis limits [min max] (default: auto)
%     - baseline_window: seconds before event onset to use for baseline [-X, 0] (default: [-5, 0])
%     - dff_measure: how to measure DFF change (default: 'max')
%       - 'max': maximum DFF value during event - baseline
%       - 'mean': mean DFF during event - baseline
%       - 'abs_max': largest absolute change (positive or negative)
%       - 'min': minimum DFF value during event - baseline
%       - 'min_change': minimum DFF when DFF decreases (omits events with increases)
%     - max_duration: maximum event duration in seconds to include (default: Inf)

% Constants
COL_TIME = 1;     % Time column
COL_VISIT = 8;    % Food visit column
COL_EATING = 9;   % Eating column
COL_DFF = 11;     % DFF data column

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

if ~isfield(options, 'event_type')
    options.event_type = 'visit';
end

if ~isfield(options, 'z_score')
    options.z_score = false;
end

if ~isfield(options, 'plot_individual_mice')
    options.plot_individual_mice = true;
end

if ~isfield(options, 'baseline_window')
    options.baseline_window = [-5, 0]; % Default: 5 seconds before onset
end

if ~isfield(options, 'dff_measure')
    options.dff_measure = 'max'; % Default: maximum DFF during event
    % Other options: 'mean', 'abs_max', 'min', 'min_change'
end

if ~isfield(options, 'max_duration')
    options.max_duration = Inf; % Default: no maximum duration limit
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

% Determine which column to use based on event_type
if strcmpi(options.event_type, 'visit')
    event_col = COL_VISIT;
    event_name = 'Food Visit';
    color_map = winter(num_mice); % Blue colormap for visits
elseif strcmpi(options.event_type, 'eating')
    event_col = COL_EATING;
    event_name = 'Eating';
    color_map = autumn(num_mice); % Red/orange colormap for eating
else
    error('Invalid event_type. Must be either "visit" or "eating"');
end

% Initialize arrays to store event data
all_event_durations = [];
all_dff_changes = [];
mouse_event_durations = cell(num_mice, 1);
mouse_dff_changes = cell(num_mice, 1);
mouse_ids = cell(num_mice, 1);
all_mouse_event_counts = zeros(num_mice, 1);

% Process each mouse to calculate event durations and DFF changes
for i = 1:num_mice
    data = filtered_mice{i, 4};
    time = data(:, COL_TIME);
    dff = data(:, COL_DFF);
    event_signal = data(:, event_col);
    
    % Get the food discovery frame
    discovery_frame = filtered_mice{i, 6};
    discovery_time_sec = data(discovery_frame, COL_TIME);
    
    % Store mouse ID
    if ischar(filtered_mice{i, 1})
        id_parts = strsplit(filtered_mice{i, 1}, '_');
        if length(id_parts) > 1
            mouse_ids{i} = id_parts{1};
        else
            mouse_ids{i} = filtered_mice{i, 1};
        end
    else
        mouse_ids{i} = num2str(i);
    end
    
    % Find onsets and offsets of events
    event_changes = diff([0; event_signal; 0]); % Add padding and compute differences
    event_onsets = find(event_changes == 1);
    event_offsets = find(event_changes == -1) - 1;
    
    % Calculate durations in seconds
    event_durations = time(event_offsets) - time(event_onsets);
    
    % Filter for events after discovery and within maximum duration
    valid_events = time(event_onsets) > discovery_time_sec & ... % Only include events after discovery
                   event_durations <= options.max_duration;      % Only include events shorter than max_duration
    
    event_onsets = event_onsets(valid_events);
    event_offsets = event_offsets(valid_events);
    event_durations = event_durations(valid_events);
    
    num_events = length(event_onsets);
    all_mouse_event_counts(i) = num_events;
    
    if num_events == 0
        fprintf('No valid %s events found for mouse %d after food discovery\n', lower(event_name), i);
        continue;
    end
    
    % Calculate DFF changes for each event
    event_dff_changes = zeros(num_events, 1);
    
    for j = 1:num_events
        onset_frame = event_onsets(j);
        offset_frame = event_offsets(j);
        onset_time = time(onset_frame);
        offset_time = time(offset_frame);
        
        % Calculate baseline DFF from window before event
        pre_mask = time >= (onset_time + options.baseline_window(1)) & ...
                  time <= (onset_time + options.baseline_window(2));
        pre_dff = dff(pre_mask);
        
        if isempty(pre_dff)
            % If no baseline window available, use average of available data before event
            pre_mask = time < onset_time;
            pre_dff = dff(pre_mask);
            if isempty(pre_dff)
                % If still no data, skip this event
                event_dff_changes(j) = NaN;
                continue;
            end
        end
        
        % Get DFF during the event (from onset to offset)
        during_mask = (time >= onset_time) & (time <= offset_time);
        during_dff = dff(during_mask);
        
        if isempty(during_dff)
            event_dff_changes(j) = NaN;
            continue;
        end
        
        % Apply z-scoring if requested (to each segment separately)
        if isfield(options, 'z_score') && options.z_score
            if ~isempty(pre_dff)
                pre_mean = mean(pre_dff, 'omitnan');
                pre_std = std(pre_dff, 'omitnan');
                if pre_std > 0
                    pre_dff = (pre_dff - pre_mean) / pre_std;
                    during_dff = (during_dff - pre_mean) / pre_std;
                end
            end
        end
        
        % Calculate baseline
        baseline_dff = mean(pre_dff, 'omitnan');
        
        % Measure DFF based on selected method
        if strcmpi(options.dff_measure, 'max')
            % Maximum DFF during event
            [max_dff, max_idx] = max(during_dff);
            event_dff_changes(j) = max_dff - baseline_dff;
        elseif strcmpi(options.dff_measure, 'mean')
            % Mean DFF during event
            mean_dff = mean(during_dff, 'omitnan');
            event_dff_changes(j) = mean_dff - baseline_dff;
        elseif strcmpi(options.dff_measure, 'abs_max')
            % Maximum absolute change (could be positive or negative)
            [pos_max, pos_idx] = max(during_dff);
            [neg_max, neg_idx] = min(during_dff);
            pos_change = pos_max - baseline_dff;
            neg_change = neg_max - baseline_dff;
            
            if abs(pos_change) > abs(neg_change)
                event_dff_changes(j) = pos_change;
            else
                event_dff_changes(j) = neg_change;
            end
        elseif strcmpi(options.dff_measure, 'min')
            % Minimum DFF during event
            [min_dff, min_idx] = min(during_dff);
            event_dff_changes(j) = min_dff - baseline_dff;
        elseif strcmpi(options.dff_measure, 'min_change')
            % Minimum DFF value during event (when DFF decreases)
            [min_dff, min_idx] = min(during_dff);
            event_dff_changes(j) = min_dff - baseline_dff;
            
            % If the change is positive (no decrease), set to NaN
            if event_dff_changes(j) > 0
                event_dff_changes(j) = NaN; % Skip events with no decrease
            end
        else
            % Default to max if invalid option
            [max_dff, max_idx] = max(during_dff);
            event_dff_changes(j) = max_dff - baseline_dff;
        end
    end
    
    % Store event data for this mouse
    mouse_event_durations{i} = event_durations;
    mouse_dff_changes{i} = event_dff_changes;
    
    % Add to the overall arrays
    all_event_durations = [all_event_durations; event_durations];
    all_dff_changes = [all_dff_changes; event_dff_changes];
end

% Create figure
figure('Position', [100, 100, 800, 600]);
hold on;

% Plot individual mice if requested
if options.plot_individual_mice
    for i = 1:num_mice
        if ~isempty(mouse_event_durations{i})
            % Remove NaN values
            valid_idx = ~isnan(mouse_dff_changes{i});
            durations = mouse_event_durations{i}(valid_idx);
            changes = mouse_dff_changes{i}(valid_idx);
            
            if ~isempty(durations)
                scatter(durations, changes, 50, ...
                    color_map(i,:), 'filled', 'MarkerFaceAlpha', 0.7, ...
                    'DisplayName', sprintf('%s (n=%d)', mouse_ids{i}, sum(valid_idx)));
                
                % Linear fit for each mouse if enough data points
                if length(durations) > 2
                    p = polyfit(durations, changes, 1);
                    x_range = linspace(min(durations), max(durations), 100);
                    y_fit = polyval(p, x_range);
                    plot(x_range, y_fit, '-', 'Color', color_map(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
        end
    end
end

% For correlation calculation, we need valid values
valid_idx = ~isnan(all_dff_changes);
valid_durations = all_event_durations(valid_idx);
valid_changes = all_dff_changes(valid_idx);

% Calculate and plot the overall correlation
if ~isempty(valid_durations) && length(valid_durations) > 2
    % Linear fit for all data
    p_all = polyfit(valid_durations, valid_changes, 1);
    x_fit = linspace(min(valid_durations), max(valid_durations), 100);
    y_fit = polyval(p_all, x_fit);
    
    % Plot the overall fit line
    plot(x_fit, y_fit, 'k-', 'LineWidth', 3, 'DisplayName', 'Overall fit');
    
    % Calculate correlation coefficient and p-value
    [r, p_val] = corr(valid_durations, valid_changes, 'rows', 'complete');
    
    % Add correlation stats to the plot
    if p_val < 0.001
        p_text = 'p < 0.001***';
    elseif p_val < 0.01
        p_text = 'p < 0.01**';
    elseif p_val < 0.05
        p_text = 'p < 0.05*';
    else
        p_text = sprintf('p = %.3f', p_val);
    end
    
    % Add annotation
    text(0.05, 0.95, sprintf('r = %.3f, %s', r, p_text), ...
        'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.7]);
    
    % Add equation for the fit line
    if p_all(2) >= 0
        eq_text = sprintf('y = %.3fx + %.3f', p_all(1), p_all(2));
    else
        eq_text = sprintf('y = %.3fx - %.3f', p_all(1), abs(p_all(2)));
    end
    
    text(0.05, 0.89, eq_text, 'Units', 'normalized', 'FontSize', 12, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', [1 1 1 0.7]);
end

% Plot formatting
xlabel(sprintf('%s Duration (seconds)', event_name), 'FontSize', 14);

% Y-axis label based on DFF measure type
if strcmpi(options.dff_measure, 'max')
    measure_text = 'Maximum';
elseif strcmpi(options.dff_measure, 'mean')
    measure_text = 'Mean';
elseif strcmpi(options.dff_measure, 'abs_max')
    measure_text = 'Maximum Absolute';
elseif strcmpi(options.dff_measure, 'min')
    measure_text = 'Minimum';
elseif strcmpi(options.dff_measure, 'min_change')
    measure_text = 'Minimum (Decrease Only)';
else
    measure_text = 'Maximum';
end

if isfield(options, 'z_score') && options.z_score
    ylabel(sprintf('%s ΔF/F Change (z-scored)', measure_text), 'FontSize', 14);
else
    ylabel(sprintf('%s ΔF/F Change', measure_text), 'FontSize', 14);
end

title({
    sprintf('Correlation between %s Duration and %s Neural Activity Change', event_name, measure_text),
    sprintf('State: %s, Source: %s (n = %d mice, %d events)', ...
        upper(options.state(1)), options.state(2:end), ...
        sum(all_mouse_event_counts > 0), ...
        sum(~isnan(all_dff_changes)))}, 'FontSize', 16);

% Set y-axis limits if provided
if isfield(options, 'ylim')
    ylim(options.ylim);
end

% Add grid and legend
grid on;
if options.plot_individual_mice && sum(all_mouse_event_counts > 0) > 1
    legend('show', 'Location', 'best');
end

% Optional: return the correlation data if requested
if nargout > 0
    corr_data = struct();
    corr_data.all_durations = all_event_durations;
    corr_data.all_dff_changes = all_dff_changes;
    corr_data.mouse_durations = mouse_event_durations;
    corr_data.mouse_dff_changes = mouse_dff_changes;
    corr_data.mouse_ids = mouse_ids;
    corr_data.valid_durations = valid_durations;
    corr_data.valid_changes = valid_changes;
    if exist('r', 'var')
        corr_data.r = r;
        corr_data.p = p_val;
        corr_data.fit = p_all;
    end
    varargout{1} = corr_data;
end
end