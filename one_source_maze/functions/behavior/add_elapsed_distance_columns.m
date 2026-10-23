function [mice_all, info] = add_elapsed_distance_columns(mice_all, opts)
% ADD_ELAPSED_DISTANCE_COLUMNS  Append two columns to every session's
%   data matrix in mice_all:
%     18 - elapsed time from food discovery (s)
%     19 - total distance traveled since the start of the recording
%
%   [mice_all, info] = ADD_ELAPSED_DISTANCE_COLUMNS(mice_all)
%   [mice_all, info] = ADD_ELAPSED_DISTANCE_COLUMNS(mice_all, opts)
%
%   RECORDINGS (how sessions are grouped)
%     Rows are grouped into recordings by name, after stripping the
%     trailing '_sess<N>': 'F10I_GLM_fasted_food_sess0' ... '_sess3'
%     all belong to the recording 'F10I_GLM_fasted_food'. Within a
%     recording, sessions are processed in session-number order and
%     treated as CONSECUTIVE chunks of one continuous recording, so
%     both columns keep running across session boundaries instead of
%     restarting in every session. A row without a '_sess<N>' suffix is
%     its own one-session recording.
%
%     This needs one continuous time axis per recording. If column 1
%     (time) is already continuous across sessions it is used as is.
%     If it restarts at a session boundary (the next session starts at
%     or before the previous one ended), that session is shifted so it
%     continues from where the previous one stopped, and a note is
%     printed. Column 1 itself is never modified.
%
%   COLUMN 18: elapsed time from food discovery
%     Discovery is taken from the DISCOVERY FRAME STORED IN mice_all{row,6}
%     (default, opts.discovery_source = 'col6') -- the same frame used by
%     plotDiscoveryAlignedData and the other analyses, so all of them use
%     one identical discovery time. mice_all{row,6} is the frame index
%     within THAT row's (session's) data matrix. Every session of a
%     recording is checked; empty, NaN, non-positive or out-of-range
%     values are ignored. If several sessions of one recording carry a
%     valid discovery frame, the EARLIEST one in time is used and a note
%     lists the others, so you can check them.
%     Column 18 = time - discovery time, in seconds. By default it is
%     NEGATIVE before discovery (time remaining until the animal finds
%     food, 0 at discovery, positive afterwards); opts.pre_discovery =
%     'nan' blanks the pre-discovery frames instead. If no session of a
%     recording has a valid discovery frame, column 18 is NaN for the
%     whole recording and a warning is printed.
%     The old behaviour-based definitions are still available through
%     opts.discovery_source ('food_visit', 'eating', 'distance').
%
%   COLUMN 19: total distance traveled
%     Cumulative path length of the CENTROID (columns 12, 13) from the
%     first frame of the recording, running across all of its sessions,
%     multiplied by opts.dist_scale. The default is 75/400 = 0.1875 cm
%     per pixel, the same calibration used for the speed column
%     (column 4), so the result is in cm. Frames with missing centroid
%     are skipped when summing steps (the displacement across such a
%     gap is counted once, and the value is interpolated across it).
%     Because every frame-to-frame step is summed, tracking glitches
%     and jitter add distance -- run remove_maze_tracking_artifacts
%     first, and use the check plots to look for sudden jumps.
%
%     CROSS-CHECK AGAINST THE SPEED COLUMN: the console table prints,
%     for every recording, the average speed implied by column 19
%     (total distance / elapsed time) next to the mean of column 4. With
%     the right scale these should be close (on a 5.6 s running segment
%     they agreed to ~1%). Over a whole recording, jitter can push the
%     distance-based figure somewhat higher, but a big factor (e.g. 35x)
%     means the wrong dist_scale was applied. The scale in use is also
%     printed above the table. A WARNING is printed under a recording's
%     row if its implied average speed exceeds 30 cm/s (far beyond what
%     a mouse sustains over a whole recording) or 3x the column-4 mean.
%
%   CHECK PLOTS
%     For opts.n_examples randomly chosen recordings, one figure each,
%     with a shared time axis: (1) dF/F with shaded patches for food
%     visits (green) and eating (red), (2) elapsed time from discovery,
%     (3) total distance traveled. A dashed black line marks the
%     discovery time and dotted gray lines mark the session boundaries
%     (so you can verify the sessions were stitched together correctly).
%
%   OPTS (all optional):
%     opts.discovery_source : what defines "discovery" --
%                             'col6'       : frame in mice_all{row,6}
%                                            (default)
%                             'food_visit' : first column 8 > 0
%                             'eating'     : first column 9 > 0
%                             'distance'   : first column 5 <=
%                                            opts.food_dist_thresh
%     opts.food_dist_thresh : only used for 'distance'. (default: 10)
%     opts.pre_discovery    : 'negative' (signed time, negative before
%                             discovery) or 'nan'. (default: 'negative')
%     opts.dist_scale       : multiplier applied to the summed pixel
%                             distance. (default: 75/400 = 0.1875, i.e.
%                             cm per pixel, matching the speed column;
%                             use 1 to keep pixels)
%     opts.dist_units       : label for the distance axis only; set to
%                             'px' if you set dist_scale = 1.
%                             (default: 'cm')
%     opts.n_examples       : number of recordings to plot as a check;
%                             0 turns plotting off. (default: 3)
%                             opts.n_example is accepted as an alias.
%     opts.example          : set to false/0 to turn plotting off
%                             (same as n_examples = 0).
%     opts.example_names    : restrict which recordings the examples are
%                             drawn from; text or cell array of text,
%                             matched as substrings of the recording
%                             name, e.g. 'fasted_food' or {'F10I','F12'}.
%                             (default: {} = any recording)
%     opts.example_seed     : RNG seed for a reproducible random choice
%                             of examples. (default: [])
%
%   OUTPUT
%     mice_all : same cell array, columns 18 and 19 written to every
%                row's data matrix (missing columns before 18, if any,
%                are padded with NaN)
%     info     : struct array, one element per recording, with fields
%                name, n_sessions, n_frames, t_discovery,
%                discovery_session, total_distance, time_restarts,
%                mean_speed_from_distance, mean_speed_col4
%
%   EXAMPLE:
%     opts = struct;
%     opts.n_examples = 4;
%     [mice_all, info] = add_elapsed_distance_columns(mice_all, opts);

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_TIME    = 1;
COL_SPEED   = 4;
COL_DIST    = 5;
COL_FOOD    = 8;
COL_EAT     = 9;
COL_DFF     = 11;
COL_X_CENT  = 12;
COL_Y_CENT  = 13;
COL_ELAPSED = 18;
COL_TRAVEL  = 19;
MICE_COL_DISCOVERY = 6;   % mice_all{row,6} = discovery frame of that session

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'discovery_source'), opts.discovery_source = 'col6';       end
if ~isfield(opts,'food_dist_thresh'), opts.food_dist_thresh = 10;           end
if ~isfield(opts,'pre_discovery'),    opts.pre_discovery    = 'negative';   end
if ~isfield(opts,'dist_scale'),       opts.dist_scale       = 75/400;       end
if ~isfield(opts,'dist_units'),       opts.dist_units       = 'cm';         end
if ~isfield(opts,'example_names'),    opts.example_names    = {};           end
if ~isfield(opts,'example_seed'),     opts.example_seed     = [];           end
if ~isfield(opts,'n_examples')
    if isfield(opts,'n_example')
        opts.n_examples = opts.n_example;
    else
        opts.n_examples = 3;
    end
end

if ~ismember(lower(opts.discovery_source), {'col6','food_visit','eating','distance'})
    error('opts.discovery_source must be ''col6'', ''food_visit'', ''eating'' or ''distance''.');
end
if ~ismember(lower(opts.pre_discovery), {'negative','nan'})
    error('opts.pre_discovery must be ''negative'' or ''nan''.');
end
opts.discovery_source = lower(opts.discovery_source);
opts.pre_discovery    = lower(opts.pre_discovery);
if ischar(opts.example_names), opts.example_names = {opts.example_names}; end

if strcmp(opts.discovery_source, 'col6') && size(mice_all, 2) < MICE_COL_DISCOVERY
    error('mice_all has only %d columns -- no discovery frame in column %d.', ...
        size(mice_all, 2), MICE_COL_DISCOVERY);
end

% Column 19 is meant to be in cm. dist_scale is cm PER PIXEL, so for this
% arena (~500 px across) it has to be well below 1 (about 0.19). A value
% >= 1 cannot be right for cm -- it is almost always an inverted value
% (e.g. 1/0.17 = 5.88 instead of 0.17), which inflated distance ~35x
% before. Fail loudly instead of writing wrong numbers into the data.
if strcmpi(opts.dist_units, 'cm') && ~(isscalar(opts.dist_scale) && ...
        isfinite(opts.dist_scale) && opts.dist_scale > 0 && opts.dist_scale < 1)
    error(['opts.dist_scale = %s is not valid for cm. It must be cm PER PIXEL ', ...
        '(about 0.19, default 75/400), i.e. between 0 and 1. If you pass 1/x ', ...
        'or a stale value from an earlier call, remove it. To keep pixels, ', ...
        'set opts.dist_scale = 1 and opts.dist_units = ''px''.'], mat2str(opts.dist_scale));
end

plot_on = opts.n_examples > 0 && ~(isfield(opts,'example') && ~opts.example);

% ------------------------------------------------------------------ %
%  GROUP ROWS INTO RECORDINGS                                         %
% ------------------------------------------------------------------ %
n_rows   = size(mice_all,1);
base     = cell(n_rows,1);
sess_num = zeros(n_rows,1);
for r = 1:n_rows
    name = char(mice_all{r,1});
    tok  = regexp(name, '_sess(\d+)$', 'tokens', 'once');
    if isempty(tok)
        base{r} = name;
    else
        base{r}     = regexprep(name, '_sess\d+$', '');
        sess_num(r) = str2double(tok{1});
    end
end
bases = unique(base, 'stable');
n_rec = numel(bases);

% ------------------------------------------------------------------ %
%  CHOOSE WHICH RECORDINGS TO PLOT AS A CHECK                          %
% ------------------------------------------------------------------ %
chosen = [];
if plot_on
    cand = 1:n_rec;
    if ~isempty(opts.example_names)
        keep = false(1,n_rec);
        for i = 1:numel(opts.example_names)
            keep = keep | contains(bases(:)', char(opts.example_names{i}));
        end
        cand = find(keep);
        if isempty(cand)
            warning('No recording name contains any of opts.example_names -- no check plots.');
        end
    end
    if ~isempty(cand)
        if ~isempty(opts.example_seed), rng(opts.example_seed); end
        n_ex   = min(opts.n_examples, numel(cand));
        chosen = cand(randperm(numel(cand), n_ex));
    end
end

% ------------------------------------------------------------------ %
%  COMPUTE COLUMNS 18 AND 19, ONE RECORDING AT A TIME                  %
% ------------------------------------------------------------------ %
info = repmat(struct('name','', 'n_sessions',0, 'n_frames',0, ...
    't_discovery',NaN, 'discovery_session',NaN, 'total_distance',NaN, ...
    'time_restarts',0, 'mean_speed_from_distance',NaN, 'mean_speed_col4',NaN), n_rec, 1);
examples = {};

fprintf('Computing elapsed time from discovery (column %d) and total distance (column %d) for %d recordings...\n', ...
    COL_ELAPSED, COL_TRAVEL, n_rec);
if strcmp(opts.discovery_source, 'col6')
    fprintf('Discovery = frame stored in mice_all{row,%d}.\n', MICE_COL_DISCOVERY);
else
    fprintf('Discovery = first frame by ''%s''.\n', opts.discovery_source);
end
fprintf('Distance = cumulative CENTROID path x dist_scale (%.4g %s per pixel).\n\n', ...
    opts.dist_scale, opts.dist_units);
fprintf('%-32s %6s %9s %14s %14s %12s %12s\n', 'Recording', 'n_sess', 'n_frames', ...
    'discovery (s)', sprintf('distance (%s)', opts.dist_units), ...
    'speed: dist/t', 'speed: col 4');
fprintf('%s\n', repmat('-', 1, 106));

for g = 1:n_rec
    rows_g = find(strcmp(base, bases{g}));
    [~, ord] = sort(sess_num(rows_g));
    rows_g = rows_g(ord);
    n_s = numel(rows_g);

    % ---- stitch the sessions of this recording together -------------
    T = []; F = []; E = []; D = []; X = []; Y = []; DFF = []; SPD = [];
    lens    = zeros(n_s,1);
    start_t = NaN(n_s,1);
    disc_global = NaN(n_s,1);     % discovery frame (index into stitched T) per session
    disc_notes  = {};
    offset = 0; prev_end = -Inf; n_restart = 0;

    for k = 1:n_s
        data = mice_all{rows_g(k),4};
        if isempty(data) || size(data,2) < COL_Y_CENT
            fprintf('  Note: %s skipped (empty or fewer than %d columns).\n', ...
                char(mice_all{rows_g(k),1}), COL_Y_CENT);
            continue;
        end
        t  = data(:, COL_TIME);
        tv = t(~isnan(t));
        if ~isempty(tv)
            if tv(1) + offset <= prev_end
                % time restarted (or overlaps) at this session boundary:
                % shift this session so it continues from the previous one
                dt = median(diff(tv));
                if isnan(dt) || dt < 0, dt = 0; end
                offset    = prev_end - tv(1) + dt;
                n_restart = n_restart + 1;
            end
        end
        t_adj = t + offset;
        tv_adj = t_adj(~isnan(t_adj));
        if ~isempty(tv_adj)
            prev_end   = max(prev_end, max(tv_adj));
            start_t(k) = tv_adj(1);
        end

        % discovery frame from mice_all{row,6} (index into THIS session)
        if strcmp(opts.discovery_source, 'col6')
            d6 = mice_all{rows_g(k), MICE_COL_DISCOVERY};
            if isnumeric(d6) && isscalar(d6) && isfinite(d6)
                d6 = round(d6);
                if d6 >= 1 && d6 <= size(data,1) && ~isnan(t_adj(d6))
                    disc_global(k) = numel(T) + d6;   % position in the stitched arrays
                else
                    disc_notes{end+1} = sprintf('%s: column-6 value %d outside 1..%d (or time NaN) -- ignored', ...
                        char(mice_all{rows_g(k),1}), d6, size(data,1)); %#ok<AGROW>
                end
            end
        end

        T   = [T;   t_adj];                  %#ok<AGROW>
        F   = [F;   data(:, COL_FOOD) > 0];  %#ok<AGROW>
        E   = [E;   data(:, COL_EAT)  > 0];  %#ok<AGROW>
        D   = [D;   data(:, COL_DIST)];      %#ok<AGROW>
        X   = [X;   data(:, COL_X_CENT)];    %#ok<AGROW>
        Y   = [Y;   data(:, COL_Y_CENT)];    %#ok<AGROW>
        DFF = [DFF; data(:, COL_DFF)];       %#ok<AGROW>
        SPD = [SPD; data(:, COL_SPEED)];     %#ok<AGROW>
        lens(k) = size(data,1);
    end

    N = sum(lens);
    info(g).name          = bases{g};
    info(g).n_sessions    = n_s;
    info(g).n_frames      = N;
    info(g).time_restarts = n_restart;
    if N == 0
        fprintf('%-32s %6d %9d %14s %14s %12s %12s\n', bases{g}, n_s, 0, '-', '-', '-', '-');
        continue;
    end

    % ---- column 18: elapsed time from discovery -------------------
    disc_sess = NaN;
    multi_note = '';
    if strcmp(opts.discovery_source, 'col6')
        has = find(~isnan(disc_global));
        if isempty(has)
            idx_disc = [];
        else
            [~, j] = min(T(disc_global(has)));   % earliest valid discovery in time
            idx_disc  = disc_global(has(j));
            disc_sess = sess_num(rows_g(has(j)));
            if numel(has) > 1
                others = setdiff(has, has(j));
                multi_note = sprintf(['    note: %d sessions carry a discovery frame in column 6; ', ...
                    'using sess%d (t = %.1f s). Others: %s'], numel(has), disc_sess, T(idx_disc), ...
                    strjoin(arrayfun(@(q) sprintf('sess%d at t = %.1f s', sess_num(rows_g(q)), ...
                    T(disc_global(q))), others(:)', 'UniformOutput', false), ', '));
            end
        end
    else
        switch opts.discovery_source
            case 'food_visit', disc_mask = F;
            case 'eating',     disc_mask = E;
            case 'distance',   disc_mask = D <= opts.food_dist_thresh;
        end
        idx_disc = find(disc_mask & ~isnan(T), 1, 'first');
        if ~isempty(idx_disc)
            cs = cumsum(lens);
            disc_sess = sess_num(rows_g(find(idx_disc <= cs, 1, 'first')));
        end
    end

    if isempty(idx_disc)
        t_disc  = NaN;
        elapsed = NaN(N,1);
    else
        t_disc  = T(idx_disc);
        elapsed = T - t_disc;
        if strcmp(opts.pre_discovery, 'nan')
            elapsed(T < t_disc) = NaN;
        end
    end

    % ---- column 19: total distance traveled -------------------------
    travel = cumulative_path_length(X, Y) * opts.dist_scale;
    total_dist = max(travel, [], 'omitnan');

    % ---- cross-check against the speed column (column 4) -----------
    t_span      = max(T) - min(T);
    speed_dist  = total_dist / t_span;
    speed_col4  = mean(SPD, 'omitnan');

    info(g).t_discovery       = t_disc;
    info(g).discovery_session = disc_sess;
    info(g).total_distance    = total_dist;
    info(g).mean_speed_from_distance = speed_dist;
    info(g).mean_speed_col4          = speed_col4;

    if isnan(t_disc)
        disc_str = 'none';
    else
        disc_str = sprintf('%.1f (s%d)', t_disc, disc_sess);
    end
    fprintf('%-32s %6d %9d %14s %14.0f %12.1f %12.1f\n', bases{g}, n_s, N, ...
        disc_str, total_dist, speed_dist, speed_col4);
    if speed_dist > 30 || (speed_col4 > 0 && speed_dist > 3*speed_col4)
        fprintf(['    WARNING: column 19 implies an average speed of %.1f cm/s over the whole ', ...
            'recording (column 4 averages %.1f cm/s). That is not plausible -- check ', ...
            'dist_scale and look for tracking glitches.\n'], speed_dist, speed_col4);
    end
    if n_restart > 0
        fprintf('    note: time column restarted at %d session boundary(ies); made continuous.\n', n_restart);
    end
    if ~isempty(multi_note)
        fprintf('%s\n', multi_note);
    end
    for q = 1:numel(disc_notes)
        fprintf('    note: %s\n', disc_notes{q});
    end
    if isnan(t_disc)
        if strcmp(opts.discovery_source, 'col6')
            fprintf('    warning: no valid discovery frame in mice_all{row,%d} -- column %d is NaN for this recording.\n', ...
                MICE_COL_DISCOVERY, COL_ELAPSED);
        else
            fprintf('    warning: no discovery event found (%s) -- column %d is NaN for this recording.\n', ...
                opts.discovery_source, COL_ELAPSED);
        end
    end

    % ---- write the columns back, session by session --------------
    pos = 0;
    for k = 1:n_s
        if lens(k) == 0, continue; end
        idx  = pos + (1:lens(k));
        pos  = pos + lens(k);
        data = mice_all{rows_g(k),4};
        if size(data,2) < COL_TRAVEL
            data(:, size(data,2)+1:COL_TRAVEL) = NaN;   % pad missing columns with NaN (not 0)
        end
        data(:, COL_ELAPSED) = elapsed(idx);
        data(:, COL_TRAVEL)  = travel(idx);
        mice_all{rows_g(k),4} = data;
    end

    % ---- keep what the check plot needs --------------------------
    if plot_on && ismember(g, chosen)
        ok = lens > 0;
        sn = sess_num(rows_g);
        b_t = start_t(ok); b_n = sn(ok);
        ex = struct();
        ex.name    = bases{g};
        ex.T       = T;
        ex.dff     = DFF;
        ex.food    = F;
        ex.eat     = E;
        ex.elapsed = elapsed;
        ex.travel  = travel;
        ex.t_disc  = t_disc;
        ex.total   = total_dist;
        ex.bound_t = b_t(2:end);
        ex.bound_lab = arrayfun(@(n) sprintf('sess%d', n), b_n(2:end), 'UniformOutput', false);
        examples{end+1} = ex; %#ok<AGROW>
    end
end

fprintf('\nDone.\n');

% ------------------------------------------------------------------ %
%  CHECK PLOTS                                                         %
% ------------------------------------------------------------------ %
if plot_on
    plot_example_recordings(examples, opts);
end

end


% ====================================================================== %
%  CUMULATIVE PATH LENGTH (skips NaN positions, interpolates across)     %
% ====================================================================== %
function cum = cumulative_path_length(x, y)

n   = numel(x);
cum = NaN(n,1);
iv  = find(~isnan(x) & ~isnan(y));
if numel(iv) < 2
    return;
end
steps = hypot(diff(x(iv)), diff(y(iv)));
cv    = [0; cumsum(steps(:))];
cum   = interp1(iv(:), cv, (1:n)', 'linear');
cum(1:iv(1)-1)       = 0;          % before the first valid position
cum(iv(end)+1:end)   = cv(end);    % after the last valid position

end


% ====================================================================== %
%  CHECK PLOT: dF/F + patches, elapsed time, total distance              %
% ====================================================================== %
function plot_example_recordings(examples, opts)

if isempty(examples)
    return;
end

col_food = [0.25 0.70 0.30];
col_eat  = [0.85 0.25 0.20];

for i = 1:numel(examples)
    e = examples{i};
    T = e.T;
    t_ok = T(~isnan(T));
    if isempty(t_ok), continue; end

    figure('Name', sprintf('Elapsed/distance check: %s', e.name), ...
        'Position', [60 60 1100 700]);
    ax = gobjects(3,1);

    % ---- (1) dF/F with food-visit and eating patches ---------------
    ax(1) = subplot(3,1,1); hold on;
    yl = [min(e.dff,[],'omitnan') max(e.dff,[],'omitnan')];
    if ~all(isfinite(yl)) || yl(1) == yl(2), yl = [-1 1]; end
    pad = 0.05 * diff(yl);
    yl  = yl + [-pad pad];
    h_food = draw_run_patches(T, e.food, yl, col_food);
    h_eat  = draw_run_patches(T, e.eat,  yl, col_eat);
    plot(T, e.dff, 'k-', 'LineWidth', 0.6);
    ylim(yl);
    ylabel('dF/F');
    hh = gobjects(0); ll = {};
    if ~isempty(h_food), hh(end+1) = h_food; ll{end+1} = 'Food visit'; end %#ok<AGROW>
    if ~isempty(h_eat),  hh(end+1) = h_eat;  ll{end+1} = 'Eating';     end %#ok<AGROW>
    if ~isempty(hh), legend(hh, ll, 'Location', 'northeast', 'FontSize', 8); end
    draw_event_lines(e, true);
    grid on; box off;

    % ---- (2) elapsed time from discovery ---------------------------
    ax(2) = subplot(3,1,2); hold on;
    plot(T, e.elapsed, '-', 'Color', [0.15 0.45 0.75], 'LineWidth', 1);
    yline(0, 'k:');
    ylabel('Elapsed from discovery (s)');
    if all(isnan(e.elapsed))
        text(0.5, 0.5, 'no discovery event detected', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.7 0.1 0.1]);
    end
    draw_event_lines(e, false);
    grid on; box off;

    % ---- (3) total distance traveled -------------------------------
    ax(3) = subplot(3,1,3); hold on;
    plot(T, e.travel, '-', 'Color', [0.55 0.30 0.10], 'LineWidth', 1);
    ylabel(sprintf('Total distance (%s)', opts.dist_units));
    xlabel('Time (s)');
    draw_event_lines(e, false);
    grid on; box off;

    linkaxes(ax, 'x');
    xlim(ax(1), [t_ok(1) t_ok(end)]);
    set(ax(1:2), 'XTickLabel', []);

    if isnan(e.t_disc)
        ttl = sprintf('%s   |   no discovery   |   total distance = %.0f %s', ...
            e.name, e.total, opts.dist_units);
    else
        ttl = sprintf('%s   |   discovery at t = %.1f s   |   total distance = %.0f %s', ...
            e.name, e.t_disc, e.total, opts.dist_units);
    end
    sgtitle(ttl, 'FontSize', 11, 'Interpreter', 'none');
end

end


% ====================================================================== %
%  SHADED PATCHES FOR EACH RUN OF TRUE IN A LOGICAL VECTOR              %
% ====================================================================== %
function h = draw_run_patches(T, mask, yl, col)

h = gobjects(0);
mask = logical(mask(:)) & ~isnan(T(:));
d = diff([0; mask; 0]);
s = find(d == 1);
e = find(d == -1) - 1;
n = numel(T);
for j = 1:numel(s)
    t0 = T(s(j));
    t1 = T(min(e(j)+1, n));        % extend to the next frame so 1-frame runs are visible
    if ~(t1 > t0), continue; end
    hp = patch('XData', [t0 t1 t1 t0], 'YData', [yl(1) yl(1) yl(2) yl(2)], ...
        'FaceColor', col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    if isempty(h)
        h = hp;
        set(h, 'HandleVisibility', 'on');   % first one carries the legend entry
    end
end

end


% ====================================================================== %
%  DISCOVERY + SESSION-BOUNDARY MARKERS (drawn in the current axes)      %
% ====================================================================== %
function draw_event_lines(e, with_labels)

if ~isnan(e.t_disc)
    if with_labels
        xline(e.t_disc, '--', 'Discovery', 'Color', [0 0 0], ...
            'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    else
        xline(e.t_disc, '--', 'Color', [0 0 0]);
    end
end
for b = 1:numel(e.bound_t)
    if isnan(e.bound_t(b)), continue; end
    if with_labels
        xline(e.bound_t(b), ':', e.bound_lab{b}, 'Color', [0.5 0.5 0.5], ...
            'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
    else
        xline(e.bound_t(b), ':', 'Color', [0.5 0.5 0.5]);
    end
end

end