function results = compareEventResponse(mice_all, options)
% COMPAREEVENTRESPONSE Quantifies the food visit / eating response per mouse
% and compares it between two states (default: fed vs fasted).
%
% Per-mouse response = mean ΔF/F in the W s AFTER onset
%                    - mean ΔF/F in the W s BEFORE onset
% computed on the mouse's event-averaged trace (n = mice, not events).
% Negative values = decrease.
%
% The calculation is IDENTICAL to the bar plot in plotFoodEventsAligned:
%   same event detection, same per-event z-scoring and smoothing, each
%   event interpolated onto a 1000-point grid, events averaged per mouse,
%   pre/post means taken from that grid, same min_events rule (>=).
% To get the same numbers, call plotFoodEventsAligned with the SAME
%   event_type, time_window, z_score, smooth_window, min_event_duration,
%   min_events, min_inter_event, source, bar_window = options.window,
%   and state = 'fed' / 'fasted' (one call per group).
%
% Parameters:
%   mice_all - Cell array with mouse data (same format as plotFoodEventsAligned)
%              col 1: mouse ID, col 2: state, col 3: source,
%              col 4: data matrix, col 6: discovery frame
%   options  - Struct with optional fields:
%     - groups: two states to compare (default: {'fed', 'fasted'})
%     - source: 'gel', 'food', or 'all' (default: 'all')
%     - event_type: 'visit' or 'eating' (default: 'visit')
%     - window: seconds before/after onset to average (default: 10)
%     - time_window: range extracted around each onset [before after] in s
%                    (default: [-20 60])
%     - z_score: z-score each event over time_window (default: false)
%     - smooth_window: moving-mean window in frames (default: 15, 0 = off)
%     - min_event_duration: minimum event duration in s (default: 1)
%     - min_events: mouse included if it has >= min_events valid events
%                   (default: 1)
%     - min_inter_event: exclude events that have another event active in
%                   the min_inter_event s before onset (default: 0 = off)
%     - test: 'ttest' (Welch's unpaired t-test, default) or 'ranksum'
%     - trace_xlim: x-limits for the mean-trace panel (default: time_window)
%     - group_colors: 2x3 RGB matrix (default: grey / red)
%     - export: write results to Excel (default: true)
%
% Returns:
%   results - struct with per-mouse table, group stats and test results

% Constants
COL_TIME   = 1;
COL_VISIT  = 8;
COL_EATING = 9;
COL_DFF    = 11;
N_GRID     = 1000;   % must match plotFoodEventsAligned

%% Defaults
if ~exist('options', 'var'),                options = struct();                 end
if ~isfield(options, 'groups'),             options.groups = {'fed', 'fasted'}; end
if ~isfield(options, 'source'),             options.source = 'all';             end
if ~isfield(options, 'event_type'),         options.event_type = 'visit';       end
if ~isfield(options, 'window'),             options.window = 10;                end
if ~isfield(options, 'time_window'),        options.time_window = [-20, 60];    end
if ~isfield(options, 'z_score'),            options.z_score = false;            end
if ~isfield(options, 'smooth_window'),      options.smooth_window = 15;         end
if ~isfield(options, 'min_event_duration'), options.min_event_duration = 1;     end
if ~isfield(options, 'min_events'),         options.min_events = 1;             end
if ~isfield(options, 'min_inter_event'),    options.min_inter_event = 0;        end
if ~isfield(options, 'test'),               options.test = 'ttest';             end
if ~isfield(options, 'trace_xlim'),         options.trace_xlim = options.time_window; end
if ~isfield(options, 'group_colors')
    options.group_colors = [0.45 0.45 0.45;
                            0.80 0.10 0.10];
end
if ~isfield(options, 'export'),             options.export = true;              end

groups = options.groups;
if numel(groups) ~= 2
    error('options.groups must contain exactly two states to compare.');
end

W = options.window;
if W > -options.time_window(1) || W > options.time_window(2)
    error('window (%g s) must fit inside time_window [%g %g].', ...
        W, options.time_window(1), options.time_window(2));
end

if strcmpi(options.event_type, 'visit')
    event_col = COL_VISIT;
    event_name = 'Food Visit';
elseif strcmpi(options.event_type, 'eating')
    event_col = COL_EATING;
    event_name = 'Eating';
else
    error('Invalid event_type. Must be either "visit" or "eating"');
end

if options.z_score
    data_label = 'ΔF/F (z-scored)';
else
    data_label = 'ΔF/F';
end

%% Source filter
src_mask = true(size(mice_all, 1), 1);
if ~strcmpi(options.source, 'all')
    src_mask = strcmp(mice_all(:, 3), options.source);
end

%% Common grid and pre/post masks (same as plotFoodEventsAligned)
time_grid = linspace(options.time_window(1), options.time_window(2), N_GRID);
pre_mask  = time_grid >= -W & time_grid < 0;
post_mask = time_grid >  0  & time_grid <= W;

%% Per-mouse response
mouse_ids = {};
mouse_group = {};
n_events_all = [];
pre_all = [];
post_all = [];
meets_min_all = false(0, 1);
included_all = false(0, 1);

group_traces = cell(1, 2);
resp = cell(1, 2);
group_n_events = cell(1, 2);

for g = 1:2
    idx = find(src_mask & strcmpi(mice_all(:, 2), groups{g}));
    traces = nan(numel(idx), N_GRID);
    r = nan(numel(idx), 1);
    ne = zeros(numel(idx), 1);

    for k = 1:numel(idx)
        m = idx(k);
        data = mice_all{m, 4};
        time = data(:, COL_TIME);
        dff  = data(:, COL_DFF);
        event_signal = data(:, event_col);
        discovery_time_sec = data(mice_all{m, 6}, COL_TIME);

        % ── Event detection (identical to plotFoodEventsAligned) ──
        event_changes = diff([0; event_signal; 0]);
        event_onsets  = find(event_changes == 1);
        event_offsets = find(event_changes == -1) - 1;
        event_durations = time(event_offsets) - time(event_onsets);

        valid_events = event_durations >= options.min_event_duration & ...
                       time(event_onsets) > discovery_time_sec;

        if options.min_inter_event > 0
            for j = find(valid_events)'
                prior = time >= time(event_onsets(j)) - options.min_inter_event & ...
                        time <  time(event_onsets(j));
                if any(event_signal(prior) > 0)
                    valid_events(j) = false;
                end
            end
        end

        event_onsets = event_onsets(valid_events);
        ne(k) = numel(event_onsets);

        % ── Per-event traces on the grid, then mouse average ──
        if ne(k) > 0
            events_interp = nan(ne(k), N_GRID);
            for j = 1:ne(k)
                onset_time = time(event_onsets(j));
                tm = time >= (onset_time + options.time_window(1)) & ...
                     time <= (onset_time + options.time_window(2));
                t = time(tm) - onset_time;
                y = dff(tm);

                if options.z_score
                    y = (y - nanmean(y)) / nanstd(y);
                end
                if options.smooth_window > 0
                    y = movmean(y, options.smooth_window, 'omitnan');
                end

                if numel(t) < 2
                    continue;
                end
                events_interp(j, :) = interp1(t, y, time_grid, 'linear', NaN);
            end
            mouse_trace = nanmean(events_interp, 1);
        else
            mouse_trace = nan(1, N_GRID);
        end

        % ── Pre/post from the grid (identical to the bar plot) ──
        pre  = mean(mouse_trace(pre_mask),  'omitnan');
        post = mean(mouse_trace(post_mask), 'omitnan');

        meets_min = ne(k) >= options.min_events;
        keep_mouse = meets_min && ~isnan(pre) && ~isnan(post);

        if keep_mouse
            r(k) = post - pre;
            traces(k, :) = mouse_trace;
        end

        % Mouse ID
        if ischar(mice_all{m, 1}) || isstring(mice_all{m, 1})
            id = char(mice_all{m, 1});
        else
            id = sprintf('mouse_%d', m);
        end

        mouse_ids{end+1, 1}     = id;          %#ok<AGROW>
        mouse_group{end+1, 1}   = groups{g};   %#ok<AGROW>
        n_events_all(end+1, 1)  = ne(k);       %#ok<AGROW>
        pre_all(end+1, 1)       = pre;         %#ok<AGROW>
        post_all(end+1, 1)      = post;        %#ok<AGROW>
        meets_min_all(end+1, 1) = meets_min;   %#ok<AGROW>
        included_all(end+1, 1)  = keep_mouse;  %#ok<AGROW>
    end

    n_dropped = sum(isnan(r));
    if n_dropped > 0
        fprintf('%s: %d mouse/mice excluded (fewer than %d valid %s events or no data in a window).\n', ...
            groups{g}, n_dropped, options.min_events, lower(event_name));
    end

    keep = ~isnan(r);
    resp{g} = r(keep);
    group_traces{g} = traces(keep, :);
    group_n_events{g} = ne(keep);
end

per_mouse = table(mouse_ids, mouse_group, n_events_all, pre_all, post_all, ...
    post_all - pre_all, meets_min_all, included_all, ...
    'VariableNames', {'Mouse_ID', 'State', 'N_events', 'Pre_mean', 'Post_mean', ...
                      'Response_post_minus_pre', 'Meets_min_events', 'Included'});

n1 = numel(resp{1});
n2 = numel(resp{2});
if n1 < 2 || n2 < 2
    error('Need at least 2 mice per group (got %s: %d, %s: %d).', ...
        groups{1}, n1, groups{2}, n2);
end

%% Statistics (requires Statistics and Machine Learning Toolbox)
grp_mean = [mean(resp{1}), mean(resp{2})];
grp_sem  = [std(resp{1}) / sqrt(n1), std(resp{2}) / sqrt(n2)];

% Group pre/post means (should match the bars of plotFoodEventsAligned)
grp_pre  = zeros(1, 2);
grp_post = zeros(1, 2);
for g = 1:2
    in_g = strcmp(per_mouse.State, groups{g}) & per_mouse.Included;
    grp_pre(g)  = mean(per_mouse.Pre_mean(in_g));
    grp_post(g) = mean(per_mouse.Post_mean(in_g));
end

switch lower(options.test)
    case 'ttest'
        [~, p_between, ~, st] = ttest2(resp{1}, resp{2}, 'Vartype', 'unequal');
        between_str = sprintf('Welch t(%.1f) = %.2f, p = %.4g', st.df, st.tstat, p_between);
        test_name = 'Welch two-sample t-test';
        stat_value = st.tstat;
        stat_df = st.df;
    case 'ranksum'
        [p_between, ~, st] = ranksum(resp{1}, resp{2});
        between_str = sprintf('Mann-Whitney U, p = %.4g', p_between);
        test_name = 'Mann-Whitney U (ranksum)';
        stat_value = st.ranksum;
        stat_df = NaN;
    otherwise
        error('options.test must be "ttest" or "ranksum".');
end

sd_pooled = sqrt(((n1 - 1) * var(resp{1}) + (n2 - 1) * var(resp{2})) / (n1 + n2 - 2));
cohens_d = (grp_mean(2) - grp_mean(1)) / sd_pooled;

% One-sample on post-pre = paired pre vs post (same p as plotFoodEventsAligned)
p_vs0 = nan(1, 2);
for g = 1:2
    if strcmpi(options.test, 'ttest')
        [~, p_vs0(g)] = ttest(resp{g});
    else
        p_vs0(g) = signrank(resp{g});
    end
end

%% Figure
figure('Position', [100, 100, 1100, 500]);

% ── Left: mean traces per group ──────────────────────────────────────────
subplot(1, 2, 1);
hold on;

h_pre  = patch([-W 0 0 -W], [-1 -1 1 1], [0.9 0.9 0.9], 'EdgeColor', 'none');
h_post = patch([0 W W 0],   [-1 -1 1 1], [0.95 0.88 0.88], 'EdgeColor', 'none');

line_h = gobjects(1, 2);
for g = 1:2
    tr = group_traces{g};
    mu = mean(tr, 1, 'omitnan');
    se = std(tr, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tr), 1));
    ciplot(mu - se, mu + se, time_grid, options.group_colors(g, :), 0.25);
    line_h(g) = plot(time_grid, mu, 'Color', options.group_colors(g, :), 'LineWidth', 2);
end
xline(0, '--k', 'LineWidth', 1.5);

xlim(options.trace_xlim);
yl = ylim;
set(h_pre,  'YData', [yl(1) yl(1) yl(2) yl(2)]);
set(h_post, 'YData', [yl(1) yl(1) yl(2) yl(2)]);
uistack([h_pre h_post], 'bottom');

xlabel(sprintf('Time from %s onset (s)', lower(event_name)), 'FontSize', 12);
ylabel(data_label, 'FontSize', 12);
legend(line_h, ...
    {sprintf('%s (n = %d mice, %d events)', capitalize(groups{1}), n1, sum(group_n_events{1})), ...
     sprintf('%s (n = %d mice, %d events)', capitalize(groups{2}), n2, sum(group_n_events{2}))}, ...
    'Location', 'best', 'Box', 'off');
title({sprintf('Mean response aligned to %s onset', lower(event_name)), ...
       sprintf('Shaded: %g s pre / %g s post', W, W)});
box off;

% ── Right: response (post - pre) per group ───────────────────────────────
subplot(1, 2, 2);
hold on;

b = bar(1:2, grp_mean, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 1);
b.CData = options.group_colors;
b.FaceAlpha = 0.6;

for g = 1:2
    ng = numel(resp{g});
    jitter = linspace(-0.15, 0.15, ng)';
    scatter(g + jitter, resp{g}, 35, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
end

errorbar(1:2, grp_mean, grp_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 12);
yline(0, '-k', 'LineWidth', 0.75);

all_vals = [resp{1}; resp{2}; (grp_mean + grp_sem)'; (grp_mean - grp_sem)'];
y_top = max([all_vals; 0]);
y_bot = min([all_vals; 0]);
y_range = y_top - y_bot;
if y_range == 0, y_range = 1; end
y_brk = y_top + 0.08 * y_range;
plot([1 1 2 2], [y_brk - 0.03*y_range, y_brk, y_brk, y_brk - 0.03*y_range], ...
    'k', 'LineWidth', 1.2);
text(1.5, y_brk + 0.04 * y_range, sigStars(p_between), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

xlim([0.4 2.6]);
ylim([y_bot - 0.08 * y_range, y_brk + 0.18 * y_range]);

tick_labels = cell(1, 2);
for g = 1:2
    tick_labels{g} = sprintf('%s (n=%d)\\newlinevs 0: p=%.3g %s', ...
        capitalize(groups{g}), numel(resp{g}), p_vs0(g), sigStars(p_vs0(g)));
end
set(gca, 'XTick', 1:2, 'XTickLabel', tick_labels, 'TickLabelInterpreter', 'tex');

ylabel(sprintf('Response: post - pre %s', data_label), 'FontSize', 12);
title({sprintf('%s response (%g s post - %g s pre)', event_name, W, W), ...
       sprintf('%s, d = %.2f', between_str, cohens_d)});
box off;

sgtitle(sprintf('%s vs %s: response to %s onset (source: %s, post-discovery events)', ...
    capitalize(groups{1}), capitalize(groups{2}), lower(event_name), options.source), ...
    'FontWeight', 'bold');

%% Console output
fprintf('\n%s response = mean %s (0 to %g s) - mean (-%g to 0 s), per mouse\n', ...
    event_name, data_label, W, W);
fprintf(['Settings: time_window = [%g %g], z_score = %d, smooth_window = %d, ' ...
         'min_event_duration = %g, min_events = %d, min_inter_event = %g, source = %s\n'], ...
    options.time_window(1), options.time_window(2), options.z_score, options.smooth_window, ...
    options.min_event_duration, options.min_events, options.min_inter_event, options.source);
for g = 1:2
    fprintf('  %-8s n = %2d mice (%3d events)   pre = %.3f   post = %.3f   response = %.3f ± %.3f   vs 0: p = %.4g %s\n', ...
        groups{g}, numel(resp{g}), sum(group_n_events{g}), grp_pre(g), grp_post(g), ...
        grp_mean(g), grp_sem(g), p_vs0(g), sigStars(p_vs0(g)));
end
fprintf('  %s vs %s: %s %s, Cohen''s d = %.2f\n\n', ...
    groups{1}, groups{2}, between_str, sigStars(p_between), cohens_d);

%% Results struct
results = struct();
results.per_mouse   = per_mouse;
results.groups      = groups;
results.event_type  = options.event_type;
results.window_s    = W;
results.response    = resp;
results.n_events    = group_n_events;
results.pre_mean    = grp_pre;
results.post_mean   = grp_post;
results.mean        = grp_mean;
results.sem         = grp_sem;
results.p_vs_zero   = p_vs0;
results.test        = test_name;
results.p_between   = p_between;
results.stat_value  = stat_value;
results.stat_df     = stat_df;
results.cohens_d    = cohens_d;
results.time_grid   = time_grid;
results.traces      = group_traces;

%% Excel export
if options.export
    excel_filename = sprintf('event_response_%s_%s_vs_%s_%s_%gs_min%devents_%s.xlsx', ...
        options.event_type, groups{1}, groups{2}, options.source, W, ...
        options.min_events, datestr(now, 'yyyymmdd'));

    writetable(per_mouse, excel_filename, 'Sheet', 'per_mouse');

    stats_table = table( ...
        {groups{1}; groups{2}; sprintf('%s vs %s', groups{1}, groups{2})}, ...
        [n1; n2; n1 + n2], ...
        [sum(group_n_events{1}); sum(group_n_events{2}); sum(group_n_events{1}) + sum(group_n_events{2})], ...
        [grp_pre(1); grp_pre(2); NaN], ...
        [grp_post(1); grp_post(2); NaN], ...
        [grp_mean(1); grp_mean(2); grp_mean(2) - grp_mean(1)], ...
        [grp_sem(1); grp_sem(2); NaN], ...
        {'one-sample vs 0'; 'one-sample vs 0'; test_name}, ...
        [p_vs0(1); p_vs0(2); p_between], ...
        [NaN; NaN; cohens_d], ...
        'VariableNames', {'Comparison', 'N_mice', 'N_events', 'Pre_mean', 'Post_mean', ...
                          'Mean_response', 'SEM', 'Test', 'p_value', 'Cohens_d'});
    writetable(stats_table, excel_filename, 'Sheet', 'stats');

    trace_table = table(time_grid(:), 'VariableNames', {'Time_s'});
    for g = 1:2
        tr = group_traces{g};
        trace_table.(sprintf('%s_mean', groups{g})) = mean(tr, 1, 'omitnan')';
        trace_table.(sprintf('%s_SEM', groups{g})) = ...
            (std(tr, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tr), 1)))';
    end
    writetable(trace_table, excel_filename, 'Sheet', 'mean_traces');

    fprintf('Results exported to: %s\n', excel_filename);
    results.excel_filename = excel_filename;
end
end

%% ── Helper functions ────────────────────────────────────────────────────

function s = capitalize(str)
if isempty(str)
    s = str;
else
    s = [upper(str(1)) str(2:end)];
end
end

function s = sigStars(p)
if isnan(p)
    s = 'n/a';
elseif p < 0.001
    s = '***';
elseif p < 0.01
    s = '**';
elseif p < 0.05
    s = '*';
else
    s = 'n.s.';
end
end

function ciplot(lower, upper, x, color, alpha)
x_polygon = [x, fliplr(x)];
y_polygon = [upper, fliplr(lower)];
nan_indices = isnan(x_polygon) | isnan(y_polygon);
x_polygon(nan_indices) = [];
y_polygon(nan_indices) = [];
h = fill(x_polygon, y_polygon, color);
set(h, 'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off');
end