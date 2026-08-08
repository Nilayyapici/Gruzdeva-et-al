function compareBeforeAfterCorrelations_Speed(mice, options)
    % COMPAREBEFOREAFTERCORRELATIONS_SPEED Compares correlations of DFF with speed before vs after food discovery
    % Now calculates speed from x, y coordinates and time instead of using pre-existing speed column
    %   mice: cell array with mouse data
    %   options: struct with fields:
    %     - state: cell array of states to compare {'fed', 'fasted'}
    %     - source: cell array of sources {'food', 'gel'}
    %     - speed_threshold: minimum speed threshold (default: 0.5)
    %     - x_col: column index for x coordinate (default: 2)
    %     - y_col: column index for y coordinate (default: 3)
    %     - time_col: column index for time (default: 1, if empty uses frame numbers)
    %     - frame_rate: sampling rate in Hz (default: 30, used if time_col is empty)
    %     - alpha: significance threshold for zero-comparison t-tests (default: 0.05)

    % Define constants
    COL_DIST = 14;     % Distance to food
    COL_DOOR = 7;      % Door status
    COL_FOOD_INT = 8;  % Food interaction
    COL_EATING = 9;    % Eating
    COL_GROOM = 10;    % Grooming
    COL_DFF = 11;      % DFF data

    % Set default options if not provided
    if ~isfield(options, 'speed_threshold')
        options.speed_threshold = 0.5;
    end
    if ~isfield(options, 'x_col')
        options.x_col = 2;
    end
    if ~isfield(options, 'y_col')
        options.y_col = 3;
    end
    if ~isfield(options, 'time_col')
        options.time_col = [];
    end
    if ~isfield(options, 'frame_rate')
        options.frame_rate = 30;
    end
    if ~isfield(options, 'alpha')
        options.alpha = 0.05;
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

    % Calculate correlations for each mouse
    for i = 1:n
        results(i).mouse_id = mice{i, 1};
        results(i).condition = mice{i, 2};
        results(i).stimulus = mice{i, 3};

        data = mice{i, 4};
        discovery = mice{i, 6};
        results(i).discovery_frame = discovery;

        % Calculate speed from x, y coordinates and time
        speed = calculateSpeed(data, options);

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

        % Calculate correlation before discovery
        valid_before = speed(1:discovery) > options.speed_threshold & ...
                      data(1:discovery, COL_GROOM) == 0 & ...
                      data(1:discovery, COL_EATING) == 0 & ...
                      data(1:discovery, COL_FOOD_INT) == 0 & ...
                      ~isnan(speed(1:discovery));

        if any(valid_before)
            [rho_before, pval_before] = corr(...
                data(valid_before, COL_DFF), ...
                speed(valid_before), ...
                'Type', 'Pearson');
        else
            rho_before = NaN;
            pval_before = NaN;
        end
        results(i).corr_before = [rho_before, pval_before];
        results(i).z_before = fisher_z(rho_before);

        % Calculate correlation after discovery
        valid_after = speed(discovery:end_frame) > options.speed_threshold & ...
                     data(discovery:end_frame, COL_GROOM) == 0 & ...
                     data(discovery:end_frame, COL_EATING) == 0 & ...
                     data(discovery:end_frame, COL_FOOD_INT) == 0 & ...
                     ~isnan(speed(discovery:end_frame));

        if any(valid_after)
            dff_values = data(discovery:end_frame, COL_DFF);
            speed_values = speed(discovery:end_frame);
            [rho_after, pval_after] = corr(...
                dff_values(valid_after), ...
                speed_values(valid_after), ...
                'Type', 'Pearson');
        else
            rho_after = NaN;
            pval_after = NaN;
        end
        results(i).corr_after = [rho_after, pval_after];
        results(i).z_after = fisher_z(rho_after);

        % Calculate the change in correlation (after - before)
        results(i).delta_z = results(i).z_after - results(i).z_before;
    end

    % Group data by condition and stimulus
    grouped = groupResultsByCondition(results);

    % Test whether before/after correlations differ from zero (stored separately
    % to avoid struct-array assignment errors)
    zero_tests = testCorrelationsAgainstZero(grouped, options.alpha);

    % Plot before vs after comparisons
    plotBeforeAfterComparisonSpeed(grouped, options);
end


% =========================================================
% testCorrelationsAgainstZero
% One-sample t-tests vs 0 for each group, before and after.
% Results returned as a plain struct (separate from grouped,
% which is a struct array and cannot be written into directly).
% =========================================================
function zero_tests = testCorrelationsAgainstZero(grouped, alpha)
    group_names = fieldnames(grouped);
    zero_tests  = struct();

    fprintf('\n========================================================\n');
    fprintf('  One-sample t-tests: correlations vs. zero (alpha=%.2f)\n', alpha);
    fprintf('========================================================\n');

    for g = 1:length(group_names)
        gname = group_names{g};
        grp   = grouped.(gname);

        z_before = [grp.z_before];
        z_after  = [grp.z_after];

        z_before_valid = z_before(~isnan(z_before));
        z_after_valid  = z_after(~isnan(z_after));

        fprintf('\nGroup: %s  (n_before=%d, n_after=%d)\n', ...
            gname, numel(z_before_valid), numel(z_after_valid));

        % --- BEFORE ---
        if numel(z_before_valid) > 1
            [~, p_before, ci_before, stats_before] = ttest(z_before_valid, 0);
            sig_before    = p_before < alpha;
            mean_r_before = tanh(mean(z_before_valid));

            fprintf('  BEFORE  r=%.3f (Fisher z=%.3f), t(%d)=%.3f, p=%.4f%s\n', ...
                mean_r_before, mean(z_before_valid), stats_before.df, ...
                stats_before.tstat, p_before, sigLabel(sig_before));

            zero_tests.(gname).before.p      = p_before;
            zero_tests.(gname).before.tstat  = stats_before.tstat;
            zero_tests.(gname).before.df     = stats_before.df;
            zero_tests.(gname).before.ci     = ci_before;
            zero_tests.(gname).before.sig    = sig_before;
            zero_tests.(gname).before.mean_r = mean_r_before;
        else
            fprintf('  BEFORE  insufficient data (n=%d)\n', numel(z_before_valid));
            zero_tests.(gname).before.p      = NaN;
            zero_tests.(gname).before.tstat  = NaN;
            zero_tests.(gname).before.df     = NaN;
            zero_tests.(gname).before.ci     = [NaN NaN];
            zero_tests.(gname).before.sig    = false;
            zero_tests.(gname).before.mean_r = NaN;
        end

        % --- AFTER ---
        if numel(z_after_valid) > 1
            [~, p_after, ci_after, stats_after] = ttest(z_after_valid, 0);
            sig_after    = p_after < alpha;
            mean_r_after = tanh(mean(z_after_valid));

            fprintf('  AFTER   r=%.3f (Fisher z=%.3f), t(%d)=%.3f, p=%.4f%s\n', ...
                mean_r_after, mean(z_after_valid), stats_after.df, ...
                stats_after.tstat, p_after, sigLabel(sig_after));

            zero_tests.(gname).after.p      = p_after;
            zero_tests.(gname).after.tstat  = stats_after.tstat;
            zero_tests.(gname).after.df     = stats_after.df;
            zero_tests.(gname).after.ci     = ci_after;
            zero_tests.(gname).after.sig    = sig_after;
            zero_tests.(gname).after.mean_r = mean_r_after;
        else
            fprintf('  AFTER   insufficient data (n=%d)\n', numel(z_after_valid));
            zero_tests.(gname).after.p      = NaN;
            zero_tests.(gname).after.tstat  = NaN;
            zero_tests.(gname).after.df     = NaN;
            zero_tests.(gname).after.ci     = [NaN NaN];
            zero_tests.(gname).after.sig    = false;
            zero_tests.(gname).after.mean_r = NaN;
        end
    end
    fprintf('========================================================\n\n');
end


function label = sigLabel(is_sig)
    if is_sig
        label = '  *';
    else
        label = '  ns';
    end
end


% =========================================================
% calculateSpeed — unchanged
% =========================================================
function speed = calculateSpeed(data, options)
    % Calculate speed from x, y coordinates and time
    % Returns speed in units/second

    x = data(:, options.x_col);
    y = data(:, options.y_col);

    if isempty(options.time_col)
        dt = 1 / options.frame_rate;
        time_diff = ones(size(x)) * dt;
    else
        time = data(:, options.time_col);
        time_diff = [NaN; diff(time)];
    end

    dx = [NaN; diff(x)];
    dy = [NaN; diff(y)];
    distance = sqrt(dx.^2 + dy.^2);
    speed = distance ./ time_diff;

    speed(1) = NaN;
    speed(time_diff <= 0) = NaN;
    speed(isinf(speed)) = NaN;
end