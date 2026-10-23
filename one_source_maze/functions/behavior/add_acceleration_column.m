function mice_all = add_acceleration_column(mice_all, opts)
% ADD_ACCELERATION_COLUMN  Append acceleration (d(speed)/dt) as column 15
%   to every session's data matrix in mice_all.
%
%   mice_all = ADD_ACCELERATION_COLUMN(mice_all)
%   mice_all = ADD_ACCELERATION_COLUMN(mice_all, opts)
%
%   DATA COLUMNS (mice_all{row,4}) BEFORE this function runs:
%     1-time          2-x nose        3-y nose        4-speed
%     5-dist to food (shortest)       6-path from last visit
%     7-door closed(0)/open(1)        8-food visit
%     9-eating(0/1)   10-grooming(0/1) 11-dff
%     12-x centroid   13-y centroid   14-zones in maze
%
%   AFTER this function runs:
%     15-acceleration  (derivative of speed w.r.t. time)
%
%   Acceleration is computed with MATLAB's gradient(speed, time), which
%   uses the actual time coordinates rather than assuming uniform frame
%   spacing -- this keeps the derivative correct even with dropped
%   frames or variable sampling rate. Positive values = speeding up,
%   negative values = slowing down.
%
%   OPTS (all optional):
%     opts.smooth_win_sec : if > 0, speed is smoothed with a moving mean
%                           (window in seconds) before differentiating,
%                           to reduce derivative noise on jittery speed
%                           traces. (default: 0 = no smoothing, raw
%                           frame-to-frame derivative)
%     opts.example        : true/false -- after computing acceleration
%                           for all mice, plot speed and acceleration
%                           over time (same panel, two y-axes) for a
%                           handful of example sessions. (default: false)
%     opts.n_example       : number of example sessions to plot when
%                           opts.example = true. (default: 5)
%     opts.example_rows   : specific row indices into mice_all to use as
%                           the example set, overriding random pick.
%                           (default: [] = auto-select)
%     opts.example_seed   : RNG seed for reproducible random example
%                           selection. (default: [] = no seeding, i.e.
%                           a different random 5 each call)
%
%   EXAMPLE:
%     mice_all = add_acceleration_column(mice_all, struct('example', true));

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_TIME  = 1;
COL_SPEED = 4;
COL_ACCEL = 15;

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'smooth_win_sec'), opts.smooth_win_sec = 0;     end
if ~isfield(opts,'example'),        opts.example        = false; end
if ~isfield(opts,'n_example'),      opts.n_example       = 5;     end
if ~isfield(opts,'example_rows'),   opts.example_rows    = [];    end
if ~isfield(opts,'example_seed'),   opts.example_seed    = [];    end

n_rows = size(mice_all,1);
fprintf('Computing acceleration (column %d) for %d sessions', COL_ACCEL, n_rows);
if opts.smooth_win_sec > 0
    fprintf(' (speed smoothed over %.1f s before differentiating)', opts.smooth_win_sec);
end
fprintf('...\n');

n_ok = 0;
for r = 1:n_rows
    data = mice_all{r,4};
    if isempty(data) || size(data,1) < 2
        continue;
    end

    time_vec  = data(:, COL_TIME);
    speed_vec = data(:, COL_SPEED);

    % Optional smoothing of speed before differentiating -- reduces the
    % noise amplification that raw frame-to-frame differentiation causes
    if opts.smooth_win_sec > 0
        dt = median(diff(time_vec), 'omitnan');
        if isnan(dt) || dt <= 0
            speed_for_diff = speed_vec;
        else
            win_frames = max(round(opts.smooth_win_sec / dt), 1);
            speed_for_diff = movmean(speed_vec, win_frames, 'omitnan');
        end
    else
        speed_for_diff = speed_vec;
    end

    % Derivative via gradient(), using actual time coordinates so gaps
    % or non-uniform sampling don't distort the estimate
    accel_vec = NaN(size(speed_vec));
    valid = ~isnan(time_vec) & ~isnan(speed_for_diff);
    if sum(valid) >= 2
        accel_vec(valid) = gradient(speed_for_diff(valid), time_vec(valid));
    end

    if size(data,2) < COL_ACCEL
        data(:, COL_ACCEL) = NaN;   % pad any missing columns with NaN
    end
    data(:, COL_ACCEL) = accel_vec;
    mice_all{r,4} = data;
    n_ok = n_ok + 1;
end

fprintf('Done. Acceleration computed for %d/%d sessions.\n', n_ok, n_rows);

% ------------------------------------------------------------------ %
%  OPTIONAL: EXAMPLE PLOT                                             %
% ------------------------------------------------------------------ %
if opts.example
    plot_example_speed_accel(mice_all, opts, COL_TIME, COL_SPEED, COL_ACCEL);
end

end


% ====================================================================== %
%  EXAMPLE PLOT: speed & acceleration over time, a handful of sessions   %
% ====================================================================== %
function plot_example_speed_accel(mice_all, opts, COL_TIME, COL_SPEED, COL_ACCEL)

valid_rows = find(cellfun(@(d) ~isempty(d) && size(d,1) > 10 && ...
    size(d,2) >= COL_ACCEL, mice_all(:,4)));

if isempty(valid_rows)
    fprintf('No valid sessions available for example plot.\n');
    return;
end

if ~isempty(opts.example_rows)
    ex_rows = opts.example_rows(:)';
    ex_rows = ex_rows(ismember(ex_rows, valid_rows));
    if isempty(ex_rows)
        fprintf('None of opts.example_rows are valid sessions -- skipping example plot.\n');
        return;
    end
else
    if ~isempty(opts.example_seed)
        rng(opts.example_seed);
    end
    n_ex    = min(opts.n_example, length(valid_rows));
    perm    = randperm(length(valid_rows), n_ex);
    ex_rows = valid_rows(perm)';
end

n_ex = length(ex_rows);
figure('Name','Example speed & acceleration traces', ...
    'Position',[50 50 1200 200*n_ex]);

for i = 1:n_ex
    r    = ex_rows(i);
    data = mice_all{r,4};
    name = mice_all{r,1};

    t   = data(:, COL_TIME);
    spd = data(:, COL_SPEED);
    acc = data(:, COL_ACCEL);

    subplot(n_ex, 1, i); hold on;

    yyaxis left
    plot(t, spd, 'Color', [0.15 0.35 0.75], 'LineWidth', 1.0);
    ylabel('Speed');

    yyaxis right
    plot(t, acc, 'Color', [0.85 0.20 0.15], 'LineWidth', 0.7);
    yline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    ylabel('Acceleration');

    xlabel('Time (s)');
    title(strrep(name, '_', ' '), 'Interpreter', 'none', 'FontSize', 9);
    grid on; box off;
end

sgtitle('Example sessions: speed (blue, left axis) vs acceleration (red, right axis)', ...
    'FontSize', 11);

end