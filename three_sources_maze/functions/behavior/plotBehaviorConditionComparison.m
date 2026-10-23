function results = plotBehaviorConditionComparison(mice_all, options)
% PLOTBEHAVIORCONDITIONCOMPARISON  Compare average speed and/or total distance
% travelled between conditions across sessions.
%
% Layout (same as plotCorrelationConditionComparison):
%   x-axis        = conditions
%   coloured bars = sessions
%   dots          = individual mice, grey lines connect the same mouse across sessions
%
% Usage:
%   opts.conditions = {'light_control', 'dark_control', 'light_opto', 'dark_opto'};
%   opts.metric     = {'speed', 'speed_xy', 'distance'};
%                     % 'speed'    = mean of speed column 7 (cm/s)
%                     % 'speed_xy' = speed recalculated from x/y/time, px -> cm with pixel_to_cm
%                     % 'distance' = total path length from x/y: sum of step lengths, px -> cm
%   (distance = sum(speed column * dt), so it always matches average speed)
%   opts.time_lim_s  = 600; % recommended for distance: sessions differ in length
%   opts.excel_file = 'C:\data\behavior_results.xlsx';   % default: behavior_results.xlsx in current folder
%   res = plotBehaviorConditionComparison(mice_all, opts);
%
% ---------------------------------------------------------------------------
% DATA FORMAT
%   mice_all is an N x 5 cell array (a table also works), one row per
%   mouse-session:
%     {1} name        e.g. 'FCRI15control_sess1'  (mouse ID = name without '_sessX')
%     {2} group       e.g. 'light_control'        (matched against options.conditions)
%     {3} food arm    'left' / 'right' / 'center' (not used here)
%     {4} data        [nFrames x 13] double
%     {5} signal flag e.g. 'ok_signal'            (only used if require_ok_signal = true)
%
%   Data matrix columns (defaults):
%     1 time | 2 x | 3 y | 4 465 | 5 405 | 6 dF/F | 7 speed | 8 zones
%     9-11 distances to food | 13 grooming
%
% A fully custom extractor can still be passed as options.data_fn:
%   [values, ids] = data_fn(mice_all, condition, session, metric, options)
% ---------------------------------------------------------------------------

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'conditions'), error('Must specify conditions to compare in options.conditions'); end

    % ----- Metric / extraction options -----
    options = setDefault(options, 'metric', 'speed');            % 'speed', 'distance' or both in a cell
    options = setDefault(options, 'sessions', {'sess0', 'sess1', 'sess2'});
    options = setDefault(options, 'data_fn', []);                % custom extractor (see above)
    options = setDefault(options, 'col_time', 1);                % column with time (or use fps)
    options = setDefault(options, 'col_x', 2);                   % column with x position
    options = setDefault(options, 'col_y', 3);                   % column with y position
    options = setDefault(options, 'col_speed', 7);               % column with precomputed speed
    options = setDefault(options, 'use_speed_column', true);     % speed metric: use col 7 (true) or compute from x/y (false)
    options = setDefault(options, 'time_to_s', 1);               % e.g. 0.001 if the time column is in ms
    options = setDefault(options, 'fps', []);                    % used if col_time is empty
    options = setDefault(options, 'col_grooming', 13);           % grooming column, nonzero = grooming
    options = setDefault(options, 'require_ok_signal', false);   % only keep rows flagged 'ok_signal'
    options = setDefault(options, 'pixel_to_cm', 0.173);         % cm per pixel for x/y (columns 2-3)
    options = setDefault(options, 'max_gap_factor', 3);          % dt > factor*median(dt) = missing frames, not counted
    options = setDefault(options, 'time_lim_s', Inf);            % only use first N seconds of each session
    options = setDefault(options, 'remove_grooming', true);      % exclude grooming frames
    options = setDefault(options, 'speed_threshold', 0);         % cm/s; slower steps treated as stationary (removes jitter)
    options = setDefault(options, 'max_speed', Inf);             % cm/s; faster steps treated as tracking artefacts
    options = setDefault(options, 'smooth_window', 0);           % samples; moving-average on x/y (0 = off)
    options = setDefault(options, 'distance_method', 'xy');      % 'xy' = sum of sqrt((x1-x0)^2+(y1-y0)^2) * pixel_to_cm; 'speed_column' = sum(speed*dt)
    options = setDefault(options, 'speed_col_scale', 1);         % multiply column 7 by this to get cm/s (1 if already cm/s)
    options = setDefault(options, 'check_consistency', true);    % warn if x/y path length disagrees with speed column
    options = setDefault(options, 'filter_order_dist', 1);       % Butterworth low-pass on x/y (same as correlation pipeline)
    options = setDefault(options, 'cutoff_freq_dist', 0);        % low-pass on x/y: 0 = OFF (default); e.g. 0.08 to filter
    options = setDefault(options, 'distance_units', 'm');        % 'm' or 'cm'

    % ----- Statistics options -----
    options = setDefault(options, 'show_stats', true);
    options = setDefault(options, 'alpha', 0.05);
    options = setDefault(options, 'paired_sessions', true);      % same mice across sessions -> paired tests

    % ----- Plot options -----
    options = setDefault(options, 'connect_points', true);
    options = setDefault(options, 'ylimit', []);                 % [] = automatic
    options = setDefault(options, 'ylabel', []);                 % [] = automatic
    options = setDefault(options, 'title', []);                  % [] = automatic
    options = setDefault(options, 'session_colors', {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.8 0.2], [0.8 0.6 0.2]});
    options = setDefault(options, 'figure_size', [100, 100, 800, 600]);
    options = setDefault(options, 'font_size', 12);
    options = setDefault(options, 'bar_width', 0.8);
    options = setDefault(options, 'legend_location', 'northeast');

    % ----- Export -----
    options = setDefault(options, 'excel_file', 'behavior_results.xlsx'); % saved in current folder unless a full path is given; '' = no export

    session_labels_map = containers.Map({'sess0', 'sess1', 'sess2'}, {'Before', 'Learning', 'Test'});

    metrics = cellstr(options.metric);
    nC = numel(options.conditions);
    nS = numel(options.sessions);
    results = struct();

    for k = 1:numel(metrics)
        metric = lower(metrics{k});
        if ~ismember(metric, {'speed', 'speed_xy', 'distance'})
            error('Unknown metric "%s". Use ''speed'', ''speed_xy'' or ''distance''.', metric);
        end

        % Collect per-mouse values: rows = conditions, columns = sessions
        values = cell(nC, nS);
        ids = cell(nC, nS);
        for c = 1:nC
            for s = 1:nS
                if isempty(options.data_fn)
                    [v, id] = extractBehaviorMetric(mice_all, options.conditions{c}, options.sessions{s}, metric, options);
                else
                    [v, id] = options.data_fn(mice_all, options.conditions{c}, options.sessions{s}, metric, options);
                end
                values{c, s} = double(v(:));
                ids{c, s} = reshape(cellstr(id), [], 1);
            end
        end

        [stats, groupedData, errorData] = createBehaviorFigure(values, ids, metric, options, session_labels_map);

        results.(metric).values = values;
        results.(metric).mouse_ids = ids;
        results.(metric).mean = groupedData;
        results.(metric).sem = errorData;
        results.(metric).stats = stats;
    end

    if ~isempty(options.excel_file)
        exportResultsToExcel(results, options, getSessionLabels(options.sessions, session_labels_map));
    end
end

%% ========================================================================
%  Data extraction
%  ========================================================================
function [vals, ids] = extractBehaviorMetric(mice_all, condition, session, metric, options)
    % Rows of mice_all: {name, group, food_arm, data, signal_flag}
    if istable(mice_all), mice_all = table2cell(mice_all); end

    vals = zeros(0, 1);
    ids = cell(0, 1);
    for r = 1:size(mice_all, 1)
        if ~strcmpi(char(mice_all{r, 2}), condition), continue; end

        % Split 'FCRI15control_sess1' -> mouse 'FCRI15control', session 'sess1'
        tok = regexp(char(mice_all{r, 1}), '^(.*)_(sess\d+)$', 'tokens', 'once');
        if isempty(tok) || ~strcmp(tok{2}, session), continue; end

        if options.require_ok_signal && size(mice_all, 2) >= 5 && ...
           ~strcmp(char(mice_all{r, 5}), 'ok_signal')
            continue;
        end

        data = mice_all{r, 4};
        if isempty(data), continue; end

        vals(end+1, 1) = computeMetricFromTrack(data, metric, options); %#ok<AGROW>
        ids{end+1, 1} = tok{1};                                         %#ok<AGROW>
    end
end

function v = computeMetricFromTrack(data, metric, options)
    x = double(data(:, options.col_x));
    y = double(data(:, options.col_y));
    n = numel(x);

    % Time (seconds)
    if ~isempty(options.col_time)
        t = double(data(:, options.col_time)) * options.time_to_s;
    elseif ~isempty(options.fps)
        t = (0:n-1)' / options.fps;
    else
        error('Set options.col_time (time column) or options.fps (frame rate).');
    end
    t = t - t(1);

    % Frame-wise inclusion mask
    frameOK = t <= options.time_lim_s;                         % time window
    if options.remove_grooming && ~isempty(options.col_grooming)
        g = data(:, options.col_grooming);
        frameOK = frameOK & ~(~isnan(g) & g ~= 0);             % drop grooming
    end

    % ---- Speed column: average speed, or distance = sum(speed * dt) ----
    use_col = (strcmp(metric, 'speed') && options.use_speed_column) || ...
              (strcmp(metric, 'distance') && strcmpi(options.distance_method, 'speed_column'));
    if use_col
        sp = double(data(:, options.col_speed)) * options.speed_col_scale;   % cm/s
        dt = [diff(t); 0];                                     % duration of each frame (s)
        dt_med = median(dt(dt > 0));
        valid = frameOK & ~isnan(sp) & dt > 0 & dt <= options.max_gap_factor * dt_med & ...
                sp >= options.speed_threshold & sp <= options.max_speed;
        if ~any(valid), v = NaN; return; end
        if strcmp(metric, 'speed')
            v = mean(sp(valid));
        else
            v = sum(sp(valid) .* dt(valid));
            if strcmpi(options.distance_units, 'm'), v = v / 100; end
        end
        return;
    end

    % ---- Otherwise compute from the x/y trajectory ----
    if options.cutoff_freq_dist > 0
        x = lowpassTrack(x, options.filter_order_dist, options.cutoff_freq_dist);
        y = lowpassTrack(y, options.filter_order_dist, options.cutoff_freq_dist);
    end
    if options.smooth_window > 1
        x = movmean(x, options.smooth_window, 'omitnan');
        y = movmean(y, options.smooth_window, 'omitnan');
    end
    step = hypot(diff(x), diff(y)) * options.pixel_to_cm;
    dt = diff(t);
    inst_speed = step ./ dt;

    dt_med = median(dt(dt > 0));
    valid = frameOK(1:end-1) & frameOK(2:end) & ~isnan(step) & dt > 0 & ...
            dt <= options.max_gap_factor * dt_med & ...           % skip jumps across missing frames
            inst_speed >= options.speed_threshold & inst_speed <= options.max_speed;
    if ~any(valid), v = NaN; return; end

    if options.check_consistency && strcmp(metric, 'distance') && ~isempty(options.col_speed)
        sp = double(data(1:end-1, options.col_speed)) * options.speed_col_scale;
        ok = valid & ~isnan(sp);
        d_speed = sum(sp(ok) .* dt(ok));
        d_xy = sum(step(ok));
        if d_speed > 0 && (d_xy / d_speed > 1.5 || d_xy / d_speed < 0.67)
            warning(['x/y path length (%.0f) is %.1fx the distance implied by the speed column (%.0f cm). ', ...
                     'Check pixel_to_cm, tracking jumps (max_speed) or use distance_method = ''speed_column''.'], ...
                     d_xy, d_xy / d_speed, d_speed);
        end
    end

    switch metric
        case {'speed', 'speed_xy'}
            % Speed recalculated from x/y (px -> cm via pixel_to_cm) and time:
            % total path length / total time over valid steps (cm/s)
            v = sum(step(valid)) / sum(dt(valid));
        case 'distance'
            v = sum(step(valid));
            if strcmpi(options.distance_units, 'm'), v = v / 100; end
    end
end

function xf = lowpassTrack(x, order, cutoff)
    % Zero-phase Butterworth low-pass; NaN gaps are interpolated for filtering
    % and put back afterwards.
    nanIdx = isnan(x);
    if all(nanIdx) || numel(x) < 3 * (order + 1), xf = x; return; end
    xf = fillmissing(x, 'linear', 'EndValues', 'nearest');
    [b, a] = butter(order, cutoff, 'low');
    xf = filtfilt(b, a, xf);
    xf(nanIdx) = NaN;
end

function [stats, groupedData, errorData] = createBehaviorFigure(values, ids, metric, options, session_labels_map)
    [nC, nS] = size(values);
    sessionLabels = getSessionLabels(options.sessions, session_labels_map);
    conditionLabels = cellfun(@formatGroupName, options.conditions, 'UniformOutput', false);

    % Mean and SEM
    groupedData = nan(nC, nS);
    errorData = zeros(nC, nS);
    nData = zeros(nC, nS);
    for c = 1:nC
        for s = 1:nS
            v = values{c, s};
            v = v(~isnan(v));
            nData(c, s) = numel(v);
            if ~isempty(v)
                groupedData(c, s) = mean(v);
                if numel(v) > 1, errorData(c, s) = std(v) / sqrt(numel(v)); end
            end
        end
    end

    % Statistics
    stats = struct('session_comparisons', {cell(nC, 1)}, 'condition_comparisons', {cell(nS, 1)});
    if options.show_stats
        stats = performBehaviorStatisticalTests(values, ids, options, sessionLabels, conditionLabels);
    end

    % Bars
    figure('Position', options.figure_size, 'Color', 'white');
    plotData = groupedData;
    plotData(isnan(plotData)) = 0;
    if nC == 1
        h = bar(1, plotData, options.bar_width);   % single group of session bars
    else
        h = bar(plotData, options.bar_width);
    end
    hold on;
    for s = 1:nS
        h(s).FaceColor = options.session_colors{mod(s-1, numel(options.session_colors)) + 1};
        h(s).EdgeColor = 'none';
        h(s).FaceAlpha = 0.8;
    end

    barPositions = getBarPositions(h, nC, nS, options);

    % Error bars and individual points
    for s = 1:nS
        errorbar(barPositions(:, s), groupedData(:, s), errorData(:, s), 'k', ...
                 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');
        for c = 1:nC
            v = values{c, s};
            v = v(~isnan(v));
            if ~isempty(v)
                scatter(repmat(barPositions(c, s), numel(v), 1), v, 25, 'k', 'filled', ...
                        'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
            end
        end
    end

    if options.connect_points && nS > 1
        connectMouseDataPointsAcrossSessions(values, ids, barPositions);
    end

    % Vertical extent of the data, used to place significance brackets
    allPts = vertcat(values{:});
    allPts = [allPts(~isnan(allPts)); groupedData(~isnan(groupedData)) + errorData(~isnan(groupedData))];
    if isempty(allPts), allPts = 0; end
    dataMin = min([0; allPts]);
    dataMax = max(allPts);
    yspan = dataMax - dataMin;
    if yspan <= 0, yspan = 1; end
    offset = 0.06 * yspan;

    y_top = dataMax;
    if options.show_stats
        y_top = addBehaviorSignificance(stats, barPositions, groupedData, errorData, values, options, offset);
    end

    % Formatting
    [defaultYLabel, defaultTitle] = getMetricLabels(metric, options);
    if isempty(options.ylabel), yl = defaultYLabel; else, yl = options.ylabel; end
    if isempty(options.title), tl = defaultTitle; else, tl = options.title; end
    ylabel(yl, 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    title(tl, 'FontSize', options.font_size + 4, 'FontWeight', 'bold');

    xticks(1:nC);
    xticklabels(conditionLabels);
    xlim([0.5, nC + 0.5]);
    if isempty(options.ylimit)
        ylim([dataMin, y_top + 2 * offset]);
    else
        ylim(options.ylimit);
    end

    legend(h, sessionLabels, 'Location', options.legend_location, 'FontSize', options.font_size);
    legend('boxoff');
    set(gca, 'FontSize', options.font_size, 'LineWidth', 1.2, 'Color', 'white');
    box off;
    grid off;
    hold off;

    displayBehaviorStatsSummary(stats, groupedData, errorData, nData, metric, options, ...
                                sessionLabels, conditionLabels);
end

function pos = getBarPositions(h, nC, nS, options)
    pos = zeros(nC, nS);
    try
        for s = 1:nS
            pos(:, s) = h(s).XEndPoints(:);   % R2019b+
        end
    catch
        groupwidth = min(options.bar_width, nS / (nS + 1.5));
        for c = 1:nC
            for s = 1:nS
                if nS == 1
                    pos(c, s) = c;
                else
                    pos(c, s) = c + (s - (nS + 1) / 2) * groupwidth / nS;
                end
            end
        end
    end
end

function connectMouseDataPointsAcrossSessions(values, ids, barPositions)
    [nC, nS] = size(values);
    for c = 1:nC
        all_mice = unique(vertcat(ids{c, :}));
        for m = 1:numel(all_mice)
            mx = [];
            my = [];
            for s = 1:nS
                idx = find(strcmp(ids{c, s}, all_mice{m}), 1);
                if ~isempty(idx) && ~isnan(values{c, s}(idx))
                    mx(end+1) = barPositions(c, s);    %#ok<AGROW>
                    my(end+1) = values{c, s}(idx);     %#ok<AGROW>
                end
            end
            if numel(mx) >= 2
                plot(mx, my, '-', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5 0.4], 'HandleVisibility', 'off');
            end
        end
    end
end

%% ========================================================================
%  Statistics
%  ========================================================================
function stats = performBehaviorStatisticalTests(values, ids, options, sessionLabels, conditionLabels)
    [nC, nS] = size(values);
    stats = struct('session_comparisons', {cell(nC, 1)}, 'condition_comparisons', {cell(nS, 1)});
    stats.paired_sessions = options.paired_sessions;

    % 1. Sessions within each condition
    for c = 1:nC
        sIdx = find(cellfun(@(v) any(~isnan(v)), values(c, :)));
        if numel(sIdx) < 2, continue; end

        if options.paired_sessions
            M = buildPairedMatrix(values(c, sIdx), ids(c, sIdx));
            if size(M, 1) < 2
                warning('%s: fewer than 2 mice with data in all sessions - session test skipped.', conditionLabels{c});
                continue;
            end
            if numel(sIdx) == 2
                [~, p] = ttest(M(:, 1), M(:, 2));
                stats.session_comparisons{c} = makeTest(conditionLabels{c}, 'paired t-test', p, ...
                                                        size(M, 1), [sIdx(1), sIdx(2), p], options.alpha);
            else
                p = rmAnova1(M);
                pw = zeros(0, 3);
                if p < options.alpha
                    pairs = nchoosek(1:numel(sIdx), 2);
                    for i = 1:size(pairs, 1)
                        [~, pp] = ttest(M(:, pairs(i, 1)), M(:, pairs(i, 2)));
                        pw(end+1, :) = [sIdx(pairs(i, 1)), sIdx(pairs(i, 2)), min(1, pp * size(pairs, 1))]; %#ok<AGROW>
                    end
                end
                stats.session_comparisons{c} = makeTest(conditionLabels{c}, ...
                    'RM ANOVA + Bonferroni paired t-tests', p, size(M, 1), pw, options.alpha);
            end
        else
            stats.session_comparisons{c} = independentTest(values(c, sIdx), sIdx, conditionLabels{c}, options.alpha);
        end
    end

    % 2. Conditions within each session (different mice -> independent tests)
    for s = 1:nS
        cIdx = find(cellfun(@(v) any(~isnan(v)), values(:, s)));
        if numel(cIdx) < 2, continue; end
        stats.condition_comparisons{s} = independentTest(values(cIdx, s), cIdx, sessionLabels{s}, options.alpha);
    end
end

function T = independentTest(vals, idxMap, name, alpha)
    idxMap = idxMap(:);
    data = [];
    grp = [];
    for i = 1:numel(vals)
        v = vals{i};
        v = v(~isnan(v));
        data = [data; v(:)];                        %#ok<AGROW>
        grp = [grp; i * ones(numel(v), 1)];         %#ok<AGROW>
    end

    if numel(vals) == 2
        [~, p] = ttest2(data(grp == 1), data(grp == 2));
        T = makeTest(name, 'two-sample t-test', p, numel(data), [idxMap(1), idxMap(2), p], alpha);
    else
        [p, ~, st] = anova1(data, grp, 'off');
        pw = zeros(0, 3);
        if p < alpha
            cm = multcompare(st, 'Alpha', alpha, 'Display', 'off');
            pw = [idxMap(cm(:, 1)), idxMap(cm(:, 2)), cm(:, 6)];
        end
        T = makeTest(name, 'one-way ANOVA + Tukey', p, numel(data), pw, alpha);
    end
end

function M = buildPairedMatrix(valsCell, idsCell)
    % Rows = mice present (non-NaN) in every session, columns = sessions
    all_ids = unique(vertcat(idsCell{:}));
    M = nan(numel(all_ids), numel(valsCell));
    for j = 1:numel(valsCell)
        [tf, loc] = ismember(idsCell{j}, all_ids);
        v = valsCell{j};
        M(loc(tf), j) = v(tf);
    end
    M = M(all(~isnan(M), 2), :);
end

function p = rmAnova1(M)
    % One-way repeated-measures ANOVA (rows = subjects, columns = sessions).
    % No sphericity correction.
    [n, k] = size(M);
    gm = mean(M(:));
    ss_cond = n * sum((mean(M, 1) - gm).^2);
    ss_subj = k * sum((mean(M, 2) - gm).^2);
    ss_err = sum((M(:) - gm).^2) - ss_cond - ss_subj;
    df1 = k - 1;
    df2 = (k - 1) * (n - 1);
    F = (ss_cond / df1) / (ss_err / df2);
    p = 1 - fcdf(F, df1, df2);
end

function T = makeTest(name, type, p, n, pw, alpha)
    T = struct('name', name, 'test_type', type, 'p_value', p, 'n', n, ...
               'significant', p < alpha, 'pairwise', pw);
end

%% ========================================================================
%  Significance brackets
%  ========================================================================
function y_top = addBehaviorSignificance(stats, pos, groupedData, errorData, values, options, offset)
    [nC, nS] = size(values);
    condTop = zeros(nC, 1);

    % Starting height per condition = highest bar/SEM/point in that group
    for c = 1:nC
        pts = vertcat(values{c, :});
        tmp = [groupedData(c, :) + errorData(c, :), pts(:)'];
        tmp = tmp(~isnan(tmp));
        if ~isempty(tmp), condTop(c) = max(tmp); end
    end

    % 1. Session comparisons within each condition (black brackets)
    for c = 1:nC
        T = stats.session_comparisons{c};
        if isempty(T) || isempty(T.pairwise), continue; end
        sig = T.pairwise(T.pairwise(:, 3) < options.alpha, :);
        [~, ord] = sort(abs(sig(:, 2) - sig(:, 1)));   % shortest brackets lowest
        sig = sig(ord, :);
        for i = 1:size(sig, 1)
            y = condTop(c) + offset;
            drawBracket(pos(c, sig(i, 1)), pos(c, sig(i, 2)), y, offset * 0.3, ...
                        getStarSymbol(sig(i, 3)), [0 0 0], options.font_size);
            condTop(c) = y + offset;
        end
    end

    % 2. Condition comparisons within each session (brackets in session colour)
    y_top = max(condTop);
    for s = 1:nS
        T = stats.condition_comparisons{s};
        if isempty(T) || isempty(T.pairwise), continue; end
        sig = T.pairwise(T.pairwise(:, 3) < options.alpha, :);
        col = options.session_colors{mod(s-1, numel(options.session_colors)) + 1};
        for i = 1:size(sig, 1)
            y = y_top + offset;
            drawBracket(pos(sig(i, 1), s), pos(sig(i, 2), s), y, offset * 0.3, ...
                        getStarSymbol(sig(i, 3)), col, options.font_size);
            y_top = y + offset;
        end
    end
end

function drawBracket(x1, x2, y, tick, label, color, font_size)
    plot([x1 x1 x2 x2], [y - tick, y, y, y - tick], '-', 'Color', color, ...
         'LineWidth', 1.3, 'HandleVisibility', 'off');
    text(mean([x1 x2]), y, label, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', font_size, 'FontWeight', 'bold', 'Color', color);
end

%% ========================================================================
%  Console summary
%  ========================================================================
function displayBehaviorStatsSummary(stats, groupedData, errorData, nData, metric, options, sessionLabels, conditionLabels)
    [~, ttl] = getMetricLabels(metric, options);
    [nC, nS] = size(groupedData);
    fprintf('\n===== %s =====\n', ttl);

    fprintf('Mean +/- SEM (n):\n');
    for c = 1:nC
        fprintf('  %s\n', conditionLabels{c});
        for s = 1:nS
            fprintf('    %-10s %8.3f +/- %.3f (n = %d)\n', sessionLabels{s}, ...
                    groupedData(c, s), errorData(c, s), nData(c, s));
        end
    end

    if ~options.show_stats
        fprintf('==================================================\n');
        return;
    end

    fprintf('\nSession comparisons within conditions:\n');
    printTests(stats.session_comparisons, sessionLabels, options.alpha);
    fprintf('\nCondition comparisons within sessions:\n');
    printTests(stats.condition_comparisons, conditionLabels, options.alpha);
    fprintf('==================================================\n');
end

function printTests(tests, levelLabels, alpha)
    any_test = false;
    for i = 1:numel(tests)
        T = tests{i};
        if isempty(T), continue; end
        any_test = true;
        sigText = '';
        if T.significant, sigText = ', significant'; end
        fprintf('  %s: p = %.4f (%s, n = %d%s)\n', T.name, T.p_value, T.test_type, T.n, sigText);
        is_posthoc = contains(T.test_type, '+');   % omnibus test with pairwise follow-ups
        if is_posthoc && ~isempty(T.pairwise)
            for r = 1:size(T.pairwise, 1)
                mark = '';
                if T.pairwise(r, 3) < alpha, mark = ' *'; end
                fprintf('      %s vs %s: p = %.4f%s\n', levelLabels{T.pairwise(r, 1)}, ...
                        levelLabels{T.pairwise(r, 2)}, T.pairwise(r, 3), mark);
            end
        end
    end
    if ~any_test, fprintf('  (no tests run)\n'); end
end

%% ========================================================================
%  Excel export
%  ========================================================================
function exportResultsToExcel(results, options, sessionLabels)
    % One sheet per metric: rows = mice grouped by condition (+ Mean, SEM, n),
    % columns = selected sessions.
    % Resolve to a full path so it is clear where the file goes
    [fdir, fbase, fext] = fileparts(options.excel_file);
    if isempty(fdir), fdir = pwd; end
    if isempty(fext), fext = '.xlsx'; end
    if ~exist(fdir, 'dir'), mkdir(fdir); end
    fname = fullfile(fdir, [fbase, fext]);

    % Start clean (no stale sheets); fails if the file is open in Excel
    if exist(fname, 'file')
        try
            delete(fname);
        catch
            error('Could not overwrite %s - is it open in Excel? Close it and run again.', fname);
        end
        if exist(fname, 'file')
            error('Could not overwrite %s - is it open in Excel? Close it and run again.', fname);
        end
    end

    order = {'speed', 'distance', 'speed_xy'};
    sheetNames = containers.Map(order, {'Speed', 'Total distance', 'Speed from xy'});
    nS = numel(options.sessions);

    for k = 1:numel(order)
        metric = order{k};
        if ~isfield(results, metric), continue; end
        R = results.(metric);
        yl = getMetricLabels(metric, options);

        C = cell(0, 2 + nS);
        C(end+1, :) = [{yl}, cell(1, 1 + nS)];
        C(end+1, :) = [{'Condition', 'Mouse'}, sessionLabels];

        for c = 1:numel(options.conditions)
            cond = formatGroupName(options.conditions{c});
            mice = unique(vertcat(R.mouse_ids{c, :}));

            for m = 1:numel(mice)
                row = [{cond, mice{m}}, cell(1, nS)];
                for s = 1:nS
                    idx = find(strcmp(R.mouse_ids{c, s}, mice{m}), 1);
                    if ~isempty(idx) && ~isnan(R.values{c, s}(idx))
                        row{2 + s} = R.values{c, s}(idx);
                    end
                end
                C(end+1, :) = row; %#ok<AGROW>
            end

            nPerSess = cellfun(@(v) sum(~isnan(v)), R.values(c, :));
            C(end+1, :) = [{cond, 'Mean'}, blankNaN(R.mean(c, :))]; %#ok<AGROW>
            C(end+1, :) = [{cond, 'SEM'},  blankNaN(R.sem(c, :))];  %#ok<AGROW>
            C(end+1, :) = [{cond, 'n'},    num2cell(nPerSess)];      %#ok<AGROW>
            C(end+1, :) = cell(1, 2 + nS);                            %#ok<AGROW>
        end

        try
            writecell(C, fname, 'Sheet', sheetNames(metric));
        catch err
            error('Writing sheet "%s" to %s failed: %s', sheetNames(metric), fname, err.message);
        end
    end
    fprintf('\nExcel file saved: %s\n', fname);
end

function c = blankNaN(v)
    c = num2cell(v);
    c(isnan(v)) = {[]};
end

%% ========================================================================
%  Helpers
%  ========================================================================
function s = setDefault(s, name, val)
    if ~isfield(s, name), s.(name) = val; end
end

function labels = getSessionLabels(sessions, session_labels_map)
    labels = cell(1, numel(sessions));
    for i = 1:numel(sessions)
        if isKey(session_labels_map, sessions{i})
            labels{i} = session_labels_map(sessions{i});
        else
            labels{i} = strrep(sessions{i}, 'sess', 'Session ');
        end
    end
end

function [yl, tl] = getMetricLabels(metric, options)
    switch metric
        case 'speed'
            yl = 'Average speed (cm/s)';
            tl = 'Average Speed - Condition Comparison';
        case 'speed_xy'
            yl = 'Average speed from x/y (cm/s)';
            tl = 'Average Speed (recalculated from x/y) - Condition Comparison';
        case 'distance'
            yl = sprintf('Total distance travelled (%s)', options.distance_units);
            tl = 'Total Distance Travelled - Condition Comparison';
    end
end

function star_symbol = getStarSymbol(p_val)
    if p_val < 0.001
        star_symbol = '***';
    elseif p_val < 0.01
        star_symbol = '**';
    elseif p_val < 0.05
        star_symbol = '*';
    else
        star_symbol = 'ns';
    end
end

function formatted_name = formatGroupName(group_name)
    if isempty(group_name)
        formatted_name = '';
        return;
    end
    formatted_name = strrep(group_name, '_', ' ');
    words = strsplit(formatted_name, ' ');
    for i = 1:numel(words)
        if ~isempty(words{i})
            words{i} = [upper(words{i}(1)), lower(words{i}(2:end))];
        end
    end
    formatted_name = strjoin(words, ' ');
end