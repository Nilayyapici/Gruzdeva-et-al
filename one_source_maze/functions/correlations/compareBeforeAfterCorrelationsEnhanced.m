function compareBeforeAfterCorrelationsEnhanced(mice, options)
% This enhanced function includes both individual state/source plots
% and combined plots with 4 bars for a selected state
% NOW WITH SPEED THRESHOLD FILTERING
%
% Usage:
%   compareBeforeAfterCorrelationsEnhanced(mice, options)
%
% Parameters:
%   mice: cell array with mouse data
%   options: struct with fields:
%     - state: cell array of states to compare {'fed', 'fasted'}
%     - source: cell array of sources {'food', 'gel'}
%     - dist_limit: minimum distance threshold (default: 5)
%     - speed_threshold: minimum speed threshold in cm/s (default: 0.5)
%     - plot_combined: boolean, whether to plot combined 4-bar charts (default: false)
%     - combined_state: string, which state to use for combined plots ('fed' or 'fasted')

% Set default speed threshold if not provided
if ~isfield(options, 'speed_threshold')
    options.speed_threshold = 0; % Default minimum speed threshold in cm/s
end

% Use existing function for basic analysis (but we need to modify it to include speed)
% For now, we'll do the complete analysis here with speed filtering

% Define constants
COL_SPEED = 4;    % Speed
COL_DIST = 5;     % Distance to food
COL_DOOR = 7;     % Door status
COL_FOOD_INT = 8; % Food interaction
COL_EATING = 9;   % Eating
COL_GROOM = 10;   % Grooming
COL_DFF = 11;     % DFF data

% Set default options if not provided
if ~isfield(options, 'dist_limit')
    options.dist_limit = 5;
end

% Option to remove grooming periods (default: true)
if ~isfield(options, 'remove_grooming')
    options.remove_grooming = true;
end

% Option to use Fisher z-transformation (default: true)
if ~isfield(options, 'use_fisher_z')
    options.use_fisher_z = true;
end

% Convert string source to cell array for consistent handling
if ischar(options.source)
    options.source = {options.source};
end

% Initialize results structure
n = size(mice, 1);
results = struct(...
    'mouse_id', cell(n, 1), ...
    'condition', cell(n, 1), ...
    'stimulus', cell(n, 1), ...
    'discovery_frame', zeros(n, 1), ...
    'corr_before', zeros(n, 2), ... % [rho, pval]
    'corr_after', zeros(n, 2), ...  % [rho, pval]
    'z_before', zeros(n, 1), ...    % Fisher z-transformed correlation
    'z_after', zeros(n, 1), ...     % Fisher z-transformed correlation
    'delta_z', zeros(n, 1) ...      % Change in z (after - before)
    );

% Calculate correlations for each mouse with speed filtering
for i = 1:n
    % Store metadata
    results(i).mouse_id = mice{i, 1};
    results(i).condition = mice{i, 2};
    results(i).stimulus = mice{i, 3};

    % Get data and discovery frame
    data = mice{i, 4};
    discovery = mice{i, 6};
    results(i).discovery_frame = discovery;

    % Find end frame (second closed door or end of data)
    if length(data) > 11000
        closed_indices = find(data(discovery:end, COL_DOOR) < 1);
        if ~isempty(closed_indices)
            end_frame = discovery + closed_indices(1) - 1;
        else
            end_frame = length(data);
        end
    else
        end_frame = length(data);
    end

    % Calculate correlation before discovery WITH SPEED FILTERING
    if options.remove_grooming
        % Remove rows with grooming, eating, food interactions, and apply speed/distance thresholds
        valid_before = data(1:discovery, COL_SPEED) > options.speed_threshold & ...
            data(1:discovery, COL_DIST) > options.dist_limit & ...
            data(1:discovery, COL_GROOM) == 0 & ...
            data(1:discovery, COL_EATING) == 0 & ...
            data(1:discovery, COL_FOOD_INT) == 0;
    else
        % Only filter by speed, distance thresholds, and eating (keep grooming)
        valid_before = data(1:discovery, COL_SPEED) > options.speed_threshold & ...
            data(1:discovery, COL_DIST) > options.dist_limit & ...
            data(1:discovery, COL_EATING) == 0 & ...
            data(1:discovery, COL_FOOD_INT) == 0;
    end

    if any(valid_before)
        [rho_before, pval_before] = corr(...
            data(valid_before, COL_DFF), ...
            data(valid_before, COL_DIST), ...
            'Type', 'Pearson');
    else
        rho_before = NaN;
        pval_before = NaN;
    end

    results(i).corr_before = [rho_before, pval_before];

    % Apply Fisher z-transformation if enabled
    if options.use_fisher_z
        results(i).z_before = fisher_z(rho_before);
    else
        results(i).z_before = rho_before; % Store the raw correlation
    end

    % Calculate correlation after discovery WITH SPEED FILTERING
    if options.remove_grooming
        % Remove rows with grooming, eating, food interactions, and apply speed/distance thresholds
        valid_after = data(discovery:end_frame, COL_SPEED) > options.speed_threshold & ...
            data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
            data(discovery:end_frame, COL_GROOM) == 0 & ...
            data(discovery:end_frame, COL_EATING) == 0 & ...
            data(discovery:end_frame, COL_FOOD_INT) == 0;
    else
        % Only filter by speed, distance thresholds, and eating (keep grooming)
        valid_after = data(discovery:end_frame, COL_SPEED) > options.speed_threshold & ...
            data(discovery:end_frame, COL_DIST) > options.dist_limit & ...
            data(discovery:end_frame, COL_EATING) == 0 & ...
            data(discovery:end_frame, COL_FOOD_INT) == 0;
    end

    if any(valid_after)
        dff_values = data(discovery:end_frame, COL_DFF);
        dist_values = data(discovery:end_frame, COL_DIST);

        [rho_after, pval_after] = corr(...
            dff_values(valid_after), ...
            dist_values(valid_after), ...
            'Type', 'Pearson');
    else
        rho_after = NaN;
        pval_after = NaN;
    end

    results(i).corr_after = [rho_after, pval_after];

    % Apply Fisher z-transformation if enabled
    if options.use_fisher_z
        results(i).z_after = fisher_z(rho_after);
    else
        results(i).z_after = rho_after; % Store the raw correlation
    end

    % Calculate the change in correlation (after - before)
    results(i).delta_z = results(i).z_after - results(i).z_before;
end

% Group data by condition and stimulus
grouped = groupResultsByCondition(results);

% Plot before vs after comparisons (individual plots)
plotBeforeAfterComparisons(grouped, options);

% If plot_combined is set to true, generate the combined plots
if isfield(options, 'plot_combined') && options.plot_combined
    % Determine which state to use for combined plot
    if isfield(options, 'combined_state')
        state_to_plot = options.combined_state;
    else
        % Default to the first state in the options
        state_to_plot = options.state{1};
    end

    % Create the combined plot for the selected state
    plotBeforeAfterBothSources(grouped, state_to_plot);
end

% Print summary with speed threshold information
fprintf('Analysis complete with speed threshold: %.2f cm/s\n', options.speed_threshold);

% Count valid data points before and after speed filtering
total_before = 0;
total_after = 0;
valid_before = 0;
valid_after = 0;

for i = 1:length(results)
    if ~isnan(results(i).z_before)
        valid_before = valid_before + 1;
    end
    if ~isnan(results(i).z_after)
        valid_after = valid_after + 1;
    end
    total_before = total_before + 1;
    total_after = total_after + 1;
end

% ── Export individual mouse before/after correlations to Excel ──────────
% Determine state and source labels for filename
if iscell(options.state)
    state_label = strjoin(options.state, '_');
else
    state_label = options.state;
end
if iscell(options.source)
    source_label = strjoin(options.source, '_');
else
    source_label = options.source;
end

% Filter results to selected state(s) and source(s)
state_list  = cellstr(options.state);
source_list = cellstr(options.source);

keep = false(length(results), 1);
for i = 1:length(results)
    keep(i) = ismember(results(i).condition, state_list) && ...
        ismember(results(i).stimulus,  source_list);
end
filtered_results = results(keep);

if ~isempty(filtered_results)
    n_exp = length(filtered_results);

    mouse_ids = reshape({filtered_results.mouse_id}',   n_exp, 1);
    conditions = reshape({filtered_results.condition}', n_exp, 1);
    stimuli   = reshape({filtered_results.stimulus}',   n_exp, 1);
    r_before  = reshape(arrayfun(@(x) x.corr_before(1), filtered_results), n_exp, 1);
    r_after   = reshape(arrayfun(@(x) x.corr_after(1),  filtered_results), n_exp, 1);
    z_before  = reshape(arrayfun(@(x) x.z_before, filtered_results),       n_exp, 1);
    z_after   = reshape(arrayfun(@(x) x.z_after,  filtered_results),       n_exp, 1);
    delta_z   = reshape(arrayfun(@(x) x.delta_z,  filtered_results),       n_exp, 1);

    if options.use_fisher_z
        z_col_before = 'FisherZ_Before';
        z_col_after  = 'FisherZ_After';
        z_col_delta  = 'Delta_FisherZ';
    else
        z_col_before = 'r_stat_Before';
        z_col_after  = 'r_stat_After';
        z_col_delta  = 'Delta_r';
    end

    if options.use_fisher_z
        export_table = table(mouse_ids, conditions, stimuli, ...
            r_before, r_after, z_before, z_after, delta_z, ...
            'VariableNames', {'Mouse_ID', 'State', 'Source', ...
            'r_Before', 'r_After', ...
            'FisherZ_Before', 'FisherZ_After', 'Delta_FisherZ'});
    else
        export_table = table(mouse_ids, conditions, stimuli, ...
            r_before, r_after, delta_z, ...
            'VariableNames', {'Mouse_ID', 'State', 'Source', ...
            'r_Before', 'r_After', 'Delta_r'});
    end

    excel_filename = sprintf('before_after_%s_%s_%s.xlsx', ...
        state_label, source_label, datestr(now, 'yyyymmdd'));

    writetable(export_table, excel_filename);
    fprintf('Data exported to: %s\n', excel_filename);
    disp(export_table);
else
    fprintf('No data matched state/source selection for Excel export.\n');
end

fprintf('Valid correlations: Before discovery: %d/%d (%.1f%%), After discovery: %d/%d (%.1f%%)\n', ...
    valid_before, total_before, 100*valid_before/total_before, ...
    valid_after, total_after, 100*valid_after/total_after);
end