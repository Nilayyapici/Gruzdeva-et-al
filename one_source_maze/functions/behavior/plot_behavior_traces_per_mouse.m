function plot_behavior_traces_per_mouse(mice_all, opts)
% PLOT_BEHAVIOR_TRACES_PER_MOUSE  One stacked-trace PDF per session showing
%   the behavioral predictors over time (each continuous trace z-scored
%   and offset vertically), plus food-visit / eating / grooming events.
%
%   TRACES (top to bottom), with the mice_all column each one comes from:
%     Spatial dist          col 5   shortest-path distance to food
%     Temporal dist         derived min(time since food, time to food)
%     Time since food       derived time since the mouse was last at food
%     Time to food          derived time until the mouse is next at food
%     Path from last visit  col 6   path length since the last food visit
%     Time from discovery   col 18  (add_elapsed_distance_columns)
%     Total distance        col 19  (add_elapsed_distance_columns)
%     Speed                 col 4
%     Acceleration          col 15  (add_acceleration_column)
%   EVENT ROWS (bottom): Visits (col 8), Eating (col 9), Grooming (col 10)
%
%   "At food" = spatial distance <= opts.food_dist_thresh. Time since food
%   is NaN before the first time at food; time to food is NaN after the
%   last time at food.
%
%   A trace whose column is missing from a session's matrix is left empty
%   and labeled "(n/a)", with a note printed.
%
%   OPTS (all optional):
%     opts.state             'fasted' | 'fed'   (default: 'fasted')
%     opts.source            'food' | 'gel'     (default: 'food')
%     opts.session           0 | 1              (default: 1)
%     opts.food_dist_thresh  distance counted as "at food" for the temporal
%                            distance traces.   (default: 10)
%     opts.save_dir          folder for the PDFs (default: pwd)
%     opts.z_clip            z-scored traces are clipped to +-z_clip for
%                            DISPLAY ONLY, so rare spikes (acceleration in
%                            particular) don't flatten the rest of the
%                            trace. [] = no clipping; the lanes are then
%                            sized to the largest |z| in the session.
%                            (default: 5)
%     opts.lane_margin       empty space between neighbouring trace lanes,
%                            in z units (default: 1)
%
%     opts.ask_save          after all PDFs are made, ask in the command
%                            window whether to save the trace data of
%                            selected sessions to Excel, and which traces
%                            (predictors). (default: true)
%     opts.save_sessions     skip the question and save these sessions:
%                            indices into the printed list, or session
%                            names. [] = ask (if opts.ask_save).
%     opts.save_traces       traces to save: indices 1-9 (order of the
%                            TRACES list above), trace names, or 'all'.
%                            [] = ask (if opts.ask_save).
%
%   SAVED EXCEL FILE (one per chosen session, in opts.save_dir)
%     Time_s, then for every chosen trace its RAW value and its z-scored
%     value (z-scored over the session, NOT clipped -- clipping is only for
%     display), then the Visit / Eating / Grooming events (0/1).
%
%   LAYOUT: each trace gets its own lane of height 2*z_clip, so a trace
%   can never run into its neighbours or past the edge of the plot.

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'state'), opts.state = 'fasted'; end
if ~isfield(opts,'source'), opts.source = 'food'; end
if ~isfield(opts,'session'), opts.session = 1; end   % 0 or 1
if ~isfield(opts,'food_dist_thresh'), opts.food_dist_thresh = 10; end
if ~isfield(opts,'save_dir'), opts.save_dir = pwd; end
if ~isfield(opts,'z_clip'), opts.z_clip = 5; end
if ~isfield(opts,'lane_margin'), opts.lane_margin = 1; end
if ~isfield(opts,'ask_save'), opts.ask_save = true; end
if ~isfield(opts,'save_sessions'), opts.save_sessions = []; end
if ~isfield(opts,'save_traces'), opts.save_traces = []; end

COL_TIME        = 1;
COL_SPEED       = 4;
COL_DIST        = 5;
COL_PATH_LAST_V = 6;    % path length since last food visit
COL_FOOD        = 8;
COL_EAT         = 9;
COL_GROOM       = 10;
COL_ACCEL       = 15;   % acceleration = d(speed)/dt
COL_ELAPSED     = 18;   % elapsed time from food discovery (s)
COL_TRAVEL      = 19;   % total distance traveled


sess_tag = sprintf('_sess%d', opts.session);

keep = strcmp(mice_all(:,2), opts.state) & ...
       strcmp(mice_all(:,3), opts.source) & ...
       contains(mice_all(:,1), sess_tag);

mice_sub = mice_all(keep,:);

if isempty(mice_sub)
    error('No matching sessions found.');
end

if ~exist(opts.save_dir, 'dir')
    mkdir(opts.save_dir);
end

% Trace table (top to bottom). Time since / to food use a lighter and a
% darker shade of the temporal-distance orange, since they are its two
% components.
tr_labels = {'Spatial dist', 'Temporal dist', 'Time since food', 'Time to food', ...
             'Path from last visit', 'Time from discovery', 'Total distance', ...
             'Speed', 'Acceleration'};
tr_colors = {[0.0000 0.4470 0.7410], ...   % spatial dist
             [0.8500 0.3250 0.0980], ...   % temporal dist
             [0.95   0.60   0.40  ], ...   % time since food (light orange)
             [0.60   0.20   0.05  ], ...   % time to food (dark orange)
             [0.10   0.70   0.65  ], ...   % path from last visit
             [0.20   0.45   0.45  ], ...   % time from discovery
             [0.80   0.50   0.60  ], ...   % total distance
             [0.9290 0.6940 0.1250], ...   % speed
             [0.60   0.40   0.20  ]};      % acceleration
tr_lw     = [1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.0 1.0];
n_tr      = numel(tr_labels);


n_sess   = size(mice_sub,1);
sess_out = cell(n_sess, 1);   % trace data kept for the optional export

for i = 1:n_sess

    sname = mice_sub{i,1};
    data  = mice_sub{i,4};
    data(isinf(data)) = NaN;

    time_vec = data(:,COL_TIME);
    t_sec = time_vec - time_vec(1);

    spatial_distance = data(:,COL_DIST);
    speed = data(:,COL_SPEED);

    % columns added by other pipeline steps may be missing
    path_last_visit = get_col(data, COL_PATH_LAST_V);
    accel           = get_col(data, COL_ACCEL);
    elapsed         = get_col(data, COL_ELAPSED);
    travel          = get_col(data, COL_TRAVEL);
    if size(data,2) < COL_TRAVEL
        fprintf(['  %s: matrix has %d columns -- traces from missing columns are left empty ', ...
                 '(run add_acceleration_column / add_elapsed_distance_columns first).\n'], ...
                 sname, size(data,2));
    end

    food_visit = data(:,COL_FOOD) > 0;
    eating     = data(:,COL_EAT) > 0;
    grooming   = data(:,COL_GROOM) > 0;

    % temporal distance = min(time since food, time to food)
    at_food = spatial_distance <= opts.food_dist_thresh;

    time_since_food = NaN(size(time_vec));
    last_t = NaN;
    for t = 1:length(time_vec)
        if at_food(t)
            last_t = time_vec(t);
        end
        if ~isnan(last_t)
            time_since_food(t) = time_vec(t) - last_t;
        end
    end

    time_to_food = NaN(size(time_vec));
    next_t = NaN;
    for t = length(time_vec):-1:1
        if at_food(t)
            next_t = time_vec(t);
        end
        if ~isnan(next_t)
            time_to_food(t) = next_t - time_vec(t);
        end
    end

    tsf = time_since_food;
    ttf = time_to_food;
    tsf(isnan(tsf)) = Inf;
    ttf(isnan(ttf)) = Inf;

    temporal_distance = min(tsf, ttf);
    temporal_distance(isinf(temporal_distance)) = NaN;

    % z-score continuous traces (same order as tr_labels)
    tr_data = {zscore_nan(spatial_distance), ...
               zscore_nan(temporal_distance), ...
               zscore_nan(time_since_food), ...
               zscore_nan(time_to_food), ...
               zscore_nan(path_last_visit), ...
               zscore_nan(elapsed), ...
               zscore_nan(travel), ...
               zscore_nan(speed), ...
               zscore_nan(accel)};

    % keep raw + z-scored traces for the optional export (same order)
    sess_out{i} = struct('name', sname, 't_sec', t_sec, ...
        'raw', {{spatial_distance, temporal_distance, time_since_food, time_to_food, ...
                 path_last_visit, elapsed, travel, speed, accel}}, ...
        'z', {tr_data}, 'visit', double(food_visit), 'eating', double(eating), ...
        'grooming', double(grooming));

    % ── Lane layout ──────────────────────────────────────────────────
    % Each trace lives in its own lane [off - zc, off + zc]. Lanes are
    % 2*zc + lane_margin apart, so a (clipped) trace can never reach the
    % next lane or run past the axis limits.
    if isempty(opts.z_clip)
        zc = 0;
        for k = 1:n_tr
            zc = max(zc, max(abs(tr_data{k}), [], 'omitnan'));
        end
        if isempty(zc) || isnan(zc) || zc == 0, zc = 3; end
        zc = ceil(zc);
    else
        zc = opts.z_clip;
    end
    trace_gap = 2 * zc + opts.lane_margin;

    fig = figure('Name', sname, 'Position', [50 50 1300 1100]);
    hold on;

    % plot traces, stacked from the top; baseline (z = 0) of each trace is
    % marked with a faint line
    tick_y     = zeros(1, n_tr);
    tick_label = tr_labels;
    for k = 1:n_tr
        off = (n_tr - k) * trace_gap;
        z   = tr_data{k};

        % display-only clipping; NaNs stay NaN
        n_clip = sum(abs(z) > zc);
        z(z >  zc) =  zc;
        z(z < -zc) = -zc;
        frac_clip = n_clip / max(1, sum(~isnan(z)));
        if frac_clip > 0.005
            fprintf('  %s: %s clipped at |z| = %g for %.1f%% of samples (display only)\n', ...
                sname, tr_labels{k}, zc, 100 * frac_clip);
        end

        yline(off, 'Color', [0.9 0.9 0.9], 'HandleVisibility', 'off');
        plot(t_sec, z + off, 'Color', tr_colors{k}, 'LineWidth', tr_lw(k));
        tick_y(k) = off;
        if all(isnan(z))
            tick_label{k} = [tr_labels{k}, ' (n/a)'];
        end
    end

    % ── Event rows, below the lowest lane, scaled with the lane size ──
    ev_h    = 0.2 * trace_gap;          % height of one event row
    ev_step = 1.5 * ev_h;               % row-to-row spacing
    ev_top  = -(zc + opts.lane_margin); % just below the lowest lane
    ev_rows = ev_top - (0:2) * ev_step; % top edge of each row
    add_event_patches(t_sec, food_visit, [ev_rows(1)-ev_h ev_rows(1)], [0.2 0.7 0.2], 0.35);
    add_event_patches(t_sec, eating,     [ev_rows(2)-ev_h ev_rows(2)], [0.9 0.5 0.1], 0.35);
    add_event_patches(t_sec, grooming,   [ev_rows(3)-ev_h ev_rows(3)], [0.6 0.2 0.8], 0.35);
    ev_mid = ev_rows - ev_h / 2;

    % y ticks must be ascending: events (lowest) first, then the traces
    % from the bottom one up to the top one
    yticks([ev_mid(3), ev_mid(2), ev_mid(1), fliplr(tick_y)]);
    yticklabels([{'Grooming','Eating','Visits'}, fliplr(tick_label)]);
    ylim([ev_rows(3) - ev_h - 0.5 * ev_h, (n_tr-1)*trace_gap + zc + 0.5 * opts.lane_margin]);

    xlabel('Time in session (s)');
    ylabel('Behavior / predictor traces (z-scored)');
    title(strrep(sname, '_', '\_'));

    grid on;
    box off;
    xlim([0 max(t_sec)]);

    sgtitle(sprintf('%s | %s | sess%d', opts.state, opts.source, opts.session));

    pdf_name = sprintf('behavior_traces_%s_%s_%s_sess%d.pdf', ...
        clean_name(sname), opts.state, opts.source, opts.session);

    pdf_path = fullfile(opts.save_dir, pdf_name);

    set(fig, 'PaperOrientation', 'landscape', ...
             'PaperUnits', 'normalized', ...
             'PaperPosition', [0 0 1 1]);

    print(fig, pdf_path, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved: %s\n', pdf_path);

    close(fig);
end

% ------------------------------------------------------------------ %
%  OPTIONAL: SAVE TRACE DATA FOR SELECTED SESSIONS                     %
% ------------------------------------------------------------------ %
sel_s = opts.save_sessions;
sel_t = opts.save_traces;

if isempty(sel_s) && opts.ask_save
    fprintf('\n========================================\n');
    fprintf('  SAVE TRACE DATA (Excel)\n');
    fprintf('========================================\n');
    for i = 1:n_sess
        fprintf('  %2d : %s\n', i, sess_out{i}.name);
    end
    sel_s = ask_indices('Sessions to save (e.g. [1 3]), or 0 to skip: ', n_sess);
end
if ischar(sel_s) || isstring(sel_s) || iscell(sel_s)          % names given
    nm = cellfun(@(s) s.name, sess_out, 'UniformOutput', false);
    sel_s = find(ismember(nm, cellstr(sel_s)));
end
sel_s = sel_s(sel_s >= 1 & sel_s <= n_sess);

if ~isempty(sel_s)
    if isempty(sel_t) && opts.ask_save
        fprintf('\nTraces (predictors):\n');
        for k = 1:n_tr
            fprintf('  %2d : %s\n', k, tr_labels{k});
        end
        sel_t = ask_indices('Traces to save (e.g. [1 3 5]), or 0 for all: ', n_tr);
        if isempty(sel_t), sel_t = 1:n_tr; end
    end
    if isempty(sel_t) || (ischar(sel_t) && strcmpi(sel_t, 'all'))
        sel_t = 1:n_tr;
    elseif ischar(sel_t) || isstring(sel_t) || iscell(sel_t)
        sel_t = find(ismember(lower(tr_labels), lower(cellstr(sel_t))));
    end
    sel_t = sel_t(sel_t >= 1 & sel_t <= n_tr);
    if isempty(sel_t)
        fprintf('No valid traces selected -- nothing saved.\n');
        return;
    end

    for i = sel_s(:)'
        so = sess_out{i};
        T  = table(so.t_sec(:), 'VariableNames', {'Time_s'});
        for k = sel_t(:)'
            base = matlab.lang.makeValidName(strrep(tr_labels{k}, ' ', '_'));
            T.(['Raw_' base]) = so.raw{k}(:);
            T.(['Z_'   base]) = so.z{k}(:);
        end
        T.Visit    = so.visit(:);
        T.Eating   = so.eating(:);
        T.Grooming = so.grooming(:);

        xl_name = sprintf('trace_data_%s_%s_%s_sess%d.xlsx', clean_name(so.name), ...
            opts.state, opts.source, opts.session);
        xl_path = fullfile(opts.save_dir, xl_name);
        try
            writetable(T, xl_path, 'Sheet', 'traces');
            fprintf('Saved: %s  (%d traces)\n', xl_path, numel(sel_t));
        catch ME
            fprintf('ERROR saving %s: %s\n', xl_path, ME.message);
        end
    end
else
    fprintf('No trace data saved.\n');
end

end


function idx = ask_indices(prompt, n_max)
% Ask for a numeric vector in the command window; 0 / empty -> [].
    idx = [];
    while true
        try
            v = input(prompt);
        catch
            fprintf('  Please enter numbers, e.g. [1 3].\n');
            continue;
        end
        if isempty(v) || (isscalar(v) && v == 0)
            return;
        end
        v = round(v(:)');
        bad = v < 1 | v > n_max;
        if any(bad)
            fprintf('  Ignoring out-of-range value(s): %s\n', num2str(v(bad)));
        end
        idx = unique(v(~bad), 'stable');
        return;
    end
end


function x = get_col(data, c)
% Column c of data, or NaNs if the matrix does not have that column.
    if size(data,2) >= c
        x = data(:,c);
    else
        x = NaN(size(data,1), 1);
    end
end


function z = zscore_nan(x)
    mu = mean(x, 'omitnan');
    sd = std(x, 'omitnan');
    if sd > 0
        z = (x - mu) ./ sd;
    else
        z = x * NaN;
    end
end


function add_event_patches(t_sec, event_vec, y_range, color_rgb, alpha_val)

    event_vec = logical(event_vec(:));
    t_sec = t_sec(:);

    d = diff([false; event_vec; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for k = 1:length(starts)
        x1 = t_sec(starts(k));
        x2 = t_sec(stops(k));

        if x2 <= x1
            x2 = x1 + median(diff(t_sec), 'omitnan');
        end

        patch([x1 x2 x2 x1], ...
              [y_range(1) y_range(1) y_range(2) y_range(2)], ...
              color_rgb, ...
              'FaceAlpha', alpha_val, ...
              'EdgeColor', 'none');
    end
end


function out = clean_name(str)
    out = regexprep(str, '[^a-zA-Z0-9_]', '_');
end