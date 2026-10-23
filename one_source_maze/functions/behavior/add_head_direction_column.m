function [mice_all, food_pos] = add_head_direction_column(mice_all, opts)
% ADD_HEAD_DIRECTION_COLUMN  Compute egocentric head-direction-to-food
%   angle and append it as column 16 to every session's data matrix.
%
%   [mice_all, food_pos] = ADD_HEAD_DIRECTION_COLUMN(mice_all)
%   [mice_all, food_pos] = ADD_HEAD_DIRECTION_COLUMN(mice_all, opts)
%
%   DATA COLUMNS USED (mice_all{row,4}):
%     2-x nose   3-y nose   5-dist to food (shortest)
%     12-x centroid   13-y centroid
%   COLUMN WRITTEN:
%     16-head direction to food (degrees, signed, range (-180, 180])
%
%   DEFINITION
%     For each frame:
%       facing_vec  = nose_xy  - centroid_xy    (head/body orientation)
%       to_food_vec = food_xy  - nose_xy        (direction toward food,
%                     measured from the NOSE, not centroid -- keeps it
%                     consistent with how food position itself is
%                     estimated below, and with add_path_direction_column)
%       angle       = signed angle FROM facing_vec TO to_food_vec, via
%                     atan2d(cross_z, dot), in degrees:
%         0 deg    -> animal's nose is pointed directly at the food
%         +-180 deg -> animal is facing directly away from the food
%         sign     -> which side the food is on relative to the facing
%                     direction. The 0/180 meaning is convention-free,
%                     but which sign means "left" vs "right" depends on
%                     your video's x/y axis orientation (e.g. image
%                     coordinates with y increasing downward flip the
%                     sign relative to standard math coordinates).
%                     Check this against the example plot before
%                     interpreting sign in downstream analyses.
%
%   FOOD POSITION ESTIMATION
%     mice_all carries no explicit food (x,y), so for each session food
%     position is estimated as the MEDIAN NOSE position across all
%     frames where distance-to-food (column 5) is at or below
%     opts.food_dist_thresh (frames where the animal is effectively AT
%     the food). The nose is used rather than the centroid because
%     distance-to-food is presumably computed to whichever point gets
%     closest to the food, and the nose is the point actually at the
%     food when the animal eats -- the centroid stays a body-length
%     back, which would bias the estimate. Median is used rather than
%     mean for robustness to occasional tracking glitches. If no frame
%     meets the threshold in a session (e.g. a pre-discovery session
%     where the animal never reaches food), the single frame with that
%     session's minimum distance-to-food is used instead as a
%     fallback, and the function reports how many sessions needed the
%     fallback -- those estimates are shakier and worth spot-checking
%     with the example plot.
%
%     This assumes food position is stationary within a session, which
%     should hold for this behavioral paradigm; it is estimated
%     independently per session (row of mice_all), so sess0 and sess1
%     get their own (potentially different) estimated food positions.
%
%   OPTS (all optional):
%     opts.food_dist_thresh : distance threshold (same units as column
%                             5) defining "at food" frames used to
%                             estimate food position. (default: 10,
%                             matching FOOD_DIST_THRESH used elsewhere)
%     opts.example           : true/false -- after computing the angle
%                             for all mice, plot QC figures for a
%                             handful of example sessions: trajectory +
%                             sampled facing/to-food vectors + estimated
%                             food position, and the angle over time.
%                             (default: false)
%     opts.n_example         : number of example sessions to plot.
%                             (default: 5)
%     opts.example_rows      : specific row indices into mice_all to
%                             use instead of a random pick.
%                             (default: [] = auto-select)
%     opts.example_seed      : RNG seed for reproducible random example
%                             selection. (default: [] = unseeded)
%     opts.n_vectors          : number of facing/to-food vector pairs to
%                             draw on the trajectory panel.
%                             (default: 25)
%
%   OUTPUT
%     mice_all : same cell array, with column 16 = head-direction-to-
%                food angle (degrees) appended/overwritten per session
%     food_pos : n_rows x 2 matrix of estimated [x y] food positions,
%                one row per session in mice_all (NaN row for
%                empty/too-short sessions)
%
%   EXAMPLE:
%     [mice_all, food_pos] = add_head_direction_column(mice_all, ...
%         struct('example', true));

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_X_NOSE = 2;
COL_Y_NOSE = 3;
COL_DIST   = 5;
COL_DFF    = 11;
COL_X_CENT = 12;
COL_Y_CENT = 13;
COL_HD     = 16;

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'food_dist_thresh'), opts.food_dist_thresh = 10;    end
if ~isfield(opts,'example'),          opts.example          = false; end
if ~isfield(opts,'n_example'),        opts.n_example        = 5;     end
if ~isfield(opts,'example_rows'),     opts.example_rows     = [];    end
if ~isfield(opts,'example_seed'),     opts.example_seed     = [];    end
if ~isfield(opts,'n_vectors'),        opts.n_vectors        = 25;    end

n_rows   = size(mice_all,1);
food_pos = NaN(n_rows,2);

fprintf('Computing head-direction-to-food angle (column %d) for %d sessions...\n', ...
    COL_HD, n_rows);

n_ok       = 0;
n_fallback = 0;
for r = 1:n_rows
    data = mice_all{r,4};
    if isempty(data) || size(data,1) < 2
        continue;
    end

    x_nose = data(:, COL_X_NOSE);
    y_nose = data(:, COL_Y_NOSE);
    x_cent = data(:, COL_X_CENT);
    y_cent = data(:, COL_Y_CENT);
    dist   = data(:, COL_DIST);

    % ---- Estimate food position for this session --------------------
    % Uses NOSE position (not centroid) at "at food" frames, since the
    % nose is the point actually at the food when distance-to-food ~ 0.
    at_food = dist <= opts.food_dist_thresh;
    if any(at_food)
        fx = median(x_nose(at_food), 'omitnan');
        fy = median(y_nose(at_food), 'omitnan');
    else
        [dmin, i_min] = min(dist);
        if isempty(i_min) || isnan(dmin)
            continue;   % nothing usable at all -- leave food_pos(r,:) as NaN
        end
        fx = x_nose(i_min);
        fy = y_nose(i_min);
        n_fallback = n_fallback + 1;
    end
    food_pos(r,:) = [fx, fy];

    % ---- Per-frame vectors and signed angle --------------------------
    facing_x = x_nose - x_cent;
    facing_y = y_nose - y_cent;
    tofood_x = fx - x_nose;
    tofood_y = fy - y_nose;

    cross_z   = facing_x .* tofood_y - facing_y .* tofood_x;
    dot_p     = facing_x .* tofood_x + facing_y .* tofood_y;
    angle_deg = atan2d(cross_z, dot_p);   % range (-180, 180]

    % NaN out frames with a degenerate facing vector (nose==centroid,
    % or missing tracking on either point)
    facing_len = hypot(facing_x, facing_y);
    angle_deg(facing_len < eps | isnan(facing_len)) = NaN;

    if size(data,2) < COL_HD
        data(:, COL_HD) = NaN;   % pad any missing columns with NaN
    end
    data(:, COL_HD) = angle_deg;
    mice_all{r,4}   = data;
    n_ok = n_ok + 1;
end

fprintf('Done. Head direction computed for %d/%d sessions', n_ok, n_rows);
if n_fallback > 0
    fprintf(' (%d used single-frame min-distance fallback for food position -- spot-check these)', ...
        n_fallback);
end
fprintf('.\n');

% ------------------------------------------------------------------ %
%  OPTIONAL: EXAMPLE QC PLOT                                          %
% ------------------------------------------------------------------ %
if opts.example
    plot_example_head_direction(mice_all, food_pos, opts, ...
        COL_X_NOSE, COL_Y_NOSE, COL_X_CENT, COL_Y_CENT, COL_HD, COL_DFF);
end

end


% ====================================================================== %
%  EXAMPLE QC PLOT: trajectory + vectors + angle-over-time               %
% ====================================================================== %
function plot_example_head_direction(mice_all, food_pos, opts, ...
    COL_X_NOSE, COL_Y_NOSE, COL_X_CENT, COL_Y_CENT, COL_HD, COL_DFF)

has_data   = cellfun(@(d) ~isempty(d) && size(d,1) > 10 && size(d,2) >= COL_HD, ...
    mice_all(:,4));
has_food   = ~isnan(food_pos(:,1));
valid_rows = find(has_data & has_food);

if isempty(valid_rows)
    fprintf('No valid sessions available for head-direction example plot.\n');
    return;
end

if ~isempty(opts.example_rows)
    ex_rows = opts.example_rows(:)';
    ex_rows = ex_rows(ismember(ex_rows, valid_rows));
    if isempty(ex_rows)
        fprintf('None of opts.example_rows are valid -- skipping example plot.\n');
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

for i = 1:n_ex
    r    = ex_rows(i);
    data = mice_all{r,4};
    name = mice_all{r,1};

    t         = data(:,1);          % time column
    x_nose    = data(:, COL_X_NOSE);
    y_nose    = data(:, COL_Y_NOSE);
    x_cent    = data(:, COL_X_CENT);
    y_cent    = data(:, COL_Y_CENT);
    angle_deg = data(:, COL_HD);
    dff       = data(:, COL_DFF);
    fx = food_pos(r,1); fy = food_pos(r,2);

    figure('Name', sprintf('Head direction QC: %s', strrep(name,'_',' ')), ...
        'Position', [50 50 1300 520]);

    % ---- Panel 1: trajectory + sampled facing/to-food vectors --------
    subplot(1,2,1); hold on;
    h_path = plot(x_cent, y_cent, '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.6);
    h_food = plot(fx, fy, 'p', 'MarkerSize',16, 'MarkerFaceColor',[0.9 0.2 0.1], ...
        'MarkerEdgeColor','k', 'LineWidth',1);

    n_frames = size(data,1);
    samp_idx = round(linspace(1, n_frames, min(opts.n_vectors, n_frames)));
    samp_idx = samp_idx(~isnan(angle_deg(samp_idx)));

    for k = samp_idx(:)'
        nx = x_nose(k); ny = y_nose(k);
        cx = x_cent(k); cy = y_cent(k);

        % Facing vector: true length (nose - centroid)
        fvx = nx - cx; fvy = ny - cy;
        % Vectors are drawn originating from the NOSE, matching the
        % nose-referenced to-food vector computation above.
        quiver(nx, ny, fvx, fvy, 0, 'Color',[0.15 0.35 0.75], ...
            'LineWidth',1.2, 'MaxHeadSize',2, 'HandleVisibility','off');

        % Direction-to-food: drawn as a short unit vector, rescaled to
        % roughly match the facing-vector length so both are legible
        % on the same plot (its raw length is the whole-arena distance)
        tvx = fx - nx; tvy = fy - ny;
        tlen = hypot(tvx, tvy);
        flen = hypot(fvx, fvy);
        if tlen > 0 && flen > 0
            draw_len = 3 * flen;
            quiver(nx, ny, tvx/tlen*draw_len, tvy/tlen*draw_len, 0, ...
                'Color',[0.90 0.55 0.10], 'LineWidth',1.0, 'MaxHeadSize',2, ...
                'HandleVisibility','off');
        end
    end
    plot(x_nose(samp_idx), y_nose(samp_idx), 'k.', 'MarkerSize',8, ...
        'HandleVisibility','off');

    axis equal; grid on; box off;
    xlabel('x'); ylabel('y');
    title('Facing vector (blue) vs direction-to-food (orange)', 'FontSize',9);
    legend([h_path h_food], {'Centroid path','Estimated food pos'}, ...
        'Location','best', 'FontSize',7);

    % ---- Panel 2: angle (left axis) and dF/F (right axis) over time ---
    % dF/F is in % (a few units) while the angle spans +-180, so they need
    % separate y-axes -- on a shared axis dF/F would be a flat line.
    subplot(1,2,2); hold on;
    yyaxis left
    plot(t, angle_deg, 'Color',[0.3 0.3 0.3], 'LineWidth',0.8);
    yline(0,   'Color',[0.2 0.7 0.3], 'LineStyle','--', 'LineWidth',1.2);
    yline(180, 'Color',[0.8 0.2 0.2], 'LineStyle',':',  'LineWidth',0.6);
    yline(-180,'Color',[0.8 0.2 0.2], 'LineStyle',':',  'LineWidth',0.6);
    ylim([-190 190]);
    ylabel('Head direction to food (deg)');
    yyaxis right
    plot(t, dff, 'Color',[0.85 0.20 0.15], 'LineWidth',0.9);
    ylabel('dF/F');
    xlabel('Time (s)');
    title('Angle: 0 = facing food, +-180 = facing away  |  red = dF/F', 'FontSize',9);
    grid on; box off;

    sgtitle(strrep(name,'_',' '), 'FontSize',11, 'Interpreter','none');
end

end