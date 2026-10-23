function mice_all = add_path_direction_column(mice_all, food_pos, opts)
% ADD_PATH_DIRECTION_COLUMN  Compute egocentric head-direction-to-
%   shortest-path angle and append it as column 17 to every session's
%   data matrix, using simple deterministic rules per maze ZONE
%   (column 14) instead of a nearest-observed-point lookup (which was
%   noisy -- see the previous version's discussion). Since each maze
%   arm is essentially a straight corridor, the "shortest path"
%   direction within one is just its axis, not something that needs to
%   be estimated from data at all.
%
%   mice_all = ADD_PATH_DIRECTION_COLUMN(mice_all)
%   mice_all = ADD_PATH_DIRECTION_COLUMN(mice_all, opts)
%   mice_all = ADD_PATH_DIRECTION_COLUMN(mice_all, food_pos)
%   mice_all = ADD_PATH_DIRECTION_COLUMN(mice_all, food_pos, opts)
%
%   DATA COLUMNS USED (mice_all{row,4}):
%     2-x nose   3-y nose   5-dist to food (shortest)
%     12-x centroid   13-y centroid   14-zone (1-4)
%   COLUMN WRITTEN:
%     17-head direction to shortest path (degrees, signed, same
%        convention as column 16: 0 = facing the target direction,
%        +-180 = facing directly away from it)
%
%   ZONE RULES (per your maze layout: three arms plus a horizontal
%   connector joining them. WHICH arm-zone-number actually has food can
%   differ by session/mouse -- it's detected automatically per session
%   rather than assumed, see FOOD-ARM DETECTION below):
%     the two arm-zones WITHOUT food -> direction = straight UP
%                                       (0, +1) -- the only way out of
%                                       a corridor is along its axis
%     the arm-zone WITH food         -> direction = toward this
%                                       session's estimated food
%                                       position, measured from the
%                                       NOSE (not centroid) position --
%                                       the zone label is presumably
%                                       keyed to wherever the nose is,
%                                       so this keeps position and zone
%                                       consistent right at boundaries
%     zone 4 (horizontal connector) -> direction = purely HORIZONTAL,
%                                       sign(food_x - nose_x): left
%                                       if food is to the left, right if
%                                       food is to the right, no
%                                       vertical component at all --
%                                       EXCEPT in the "turn square": the
%                                       part of zone 4 already
%                                       horizontally aligned with the
%                                       food arm's entrance (within that
%                                       arm's own width of food_x),
%                                       where continuing sideways would
%                                       be pointless -- there, direction
%                                       = toward food directly (same as
%                                       the food-arm rule), which comes
%                                       out as a diagonal-down vector
%     anything else / NaN zone      -> direction undefined, angle = NaN
%
%   FOOD-ARM DETECTION
%     Rather than assuming a fixed zone number is always the food arm
%     (zone numbering isn't necessarily consistent across sessions or
%     mice), each session determines it directly: among zones 1, 2 and
%     3, whichever one contains the frame with the SMALLEST recorded
%     distance-to-food (column 5) is treated as that session's food
%     arm. If no zone gets closer than opts.food_dist_thresh anywhere
%     (e.g. a pre-discovery session where food was never reached), it
%     falls back to opts.default_food_zone and prints a note, since
%     there's no reliable evidence to detect it from in that case.
%     Zone 4 is always assumed to be the connector (this has held
%     consistently across every zone plot seen so far).
%
%   "Up" is +y in your coordinate system (matches how zone 4 sits at
%   larger y than zone 1/2/3 in your zone plot). If your video's y-axis
%   convention is flipped, flip opts.up_sign.
%
%   FOOD POSITION (food_pos ARGUMENT, OPTIONAL)
%     If you already have it from add_head_direction_column, pass it in
%     directly. If omitted or empty, it's computed here internally
%     using the exact same method (median nose position across frames
%     where distance-to-food is at or below opts.food_dist_thresh,
%     falling back to the single minimum-distance frame).
%
%   OPTS (all optional):
%     opts.food_dist_thresh : only used if food_pos is not supplied --
%                             see add_head_direction_column. (default: 10)
%     opts.food_arm_half_width : fallback half-width (same units as x)
%                             of the food arm's horizontal extent, used
%                             to detect the zone-4 "turn square" only
%                             when a session doesn't have enough of its
%                             own food-arm data to estimate the arm's
%                             real width directly. (default: 40)
%     opts.default_food_zone   : which zone number (1, 2, or 3) to
%                             assume is the food arm ONLY when a session
%                             has no evidence to detect it from (no
%                             zone ever got within opts.food_dist_thresh
%                             of food) -- see FOOD-ARM DETECTION above.
%                             (default: 3)
%     opts.up_sign            : +1 or -1 -- the sign of "up" (toward
%                             zone 4) in your y coordinate. (default:
%                             +1, i.e. zone 4 is at LARGER y)
%     opts.example              : true/false -- plot QC figures for a
%                             handful of example sessions: trajectory
%                             colored by zone, food position, sampled
%                             facing/path vectors, and the angle over
%                             time. (default: false)
%     opts.n_example              : number of example sessions to plot.
%                             (default: 5)
%     opts.example_rows           : specific row indices into mice_all
%                             to use instead of a random pick.
%                             (default: [] = auto-select)
%     opts.example_seed           : RNG seed for reproducible random
%                             example selection. (default: [])
%     opts.n_vectors                : number of vector pairs to draw on
%                             the trajectory panel. (default: 25)
%
%   OUTPUT
%     mice_all : same cell array, with column 17 = head-direction-to-
%                shortest-path angle (degrees) appended/overwritten
%
%   EXAMPLE:
%     [mice_all, food_pos] = add_head_direction_column(mice_all);
%     mice_all = add_path_direction_column(mice_all, food_pos, ...
%         struct('example', true));

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_X_NOSE = 2;
COL_Y_NOSE = 3;
COL_DIST   = 5;
COL_X_CENT = 12;
COL_Y_CENT = 13;
COL_ZONE   = 14;
COL_PATH   = 17;

% ------------------------------------------------------------------ %
%  PARSE ARGUMENTS / OPTIONS                                          %
% ------------------------------------------------------------------ %
if nargin < 2, food_pos = []; end
if nargin < 3 || isempty(opts), opts = struct(); end
% Forgiving of add_path_direction_column(mice_all, opts): a struct in
% the food_pos slot is clearly meant to be opts, not food_pos.
if isstruct(food_pos)
    opts     = food_pos;
    food_pos = [];
end
if ~isfield(opts,'food_dist_thresh'),   opts.food_dist_thresh   = 10;    end
if ~isfield(opts,'food_arm_half_width'),opts.food_arm_half_width = 40;   end
if ~isfield(opts,'default_food_zone'),  opts.default_food_zone   = 3;    end
if ~isfield(opts,'up_sign'),            opts.up_sign            = 1;     end
if ~isfield(opts,'example'),          opts.example          = false; end
if ~isfield(opts,'n_example'),        opts.n_example        = 5;     end
if ~isfield(opts,'example_rows'),     opts.example_rows     = [];    end
if ~isfield(opts,'example_seed'),     opts.example_seed     = [];    end
if ~isfield(opts,'n_vectors'),        opts.n_vectors        = 25;    end

n_rows = size(mice_all,1);

% ------------------------------------------------------------------ %
%  FOOD POSITION (reuse if given, else compute the same way as        %
%  add_head_direction_column)                                         %
% ------------------------------------------------------------------ %
if isempty(food_pos)
    fprintf('No food_pos given -- estimating it internally (same method as add_head_direction_column).\n');
    food_pos = NaN(n_rows,2);
    for r = 1:n_rows
        data = mice_all{r,4};
        if isempty(data) || size(data,1) < 2, continue; end
        x_nose = data(:, COL_X_NOSE); y_nose = data(:, COL_Y_NOSE);
        dist   = data(:, COL_DIST);
        at_food = dist <= opts.food_dist_thresh;
        if any(at_food)
            food_pos(r,:) = [median(x_nose(at_food),'omitnan'), median(y_nose(at_food),'omitnan')];
        else
            [dmin, i_min] = min(dist);
            if ~isempty(i_min) && ~isnan(dmin)
                food_pos(r,:) = [x_nose(i_min), y_nose(i_min)];
            end
        end
    end
end

% ------------------------------------------------------------------ %
%  ASSIGN COLUMN 17 PER SESSION USING THE SIMPLE ZONE RULES           %
% ------------------------------------------------------------------ %
fprintf('Computing head-direction-to-shortest-path angle (column %d) for %d sessions...\n', ...
    COL_PATH, n_rows);

n_ok = 0;
for r = 1:n_rows
    data = mice_all{r,4};
    if isempty(data) || size(data,1) < 2 || size(data,2) < COL_ZONE
        continue;
    end

    x_nose = data(:, COL_X_NOSE); y_nose = data(:, COL_Y_NOSE);
    x_cent = data(:, COL_X_CENT); y_cent = data(:, COL_Y_CENT);
    zone   = data(:, COL_ZONE);
    dist   = data(:, COL_DIST);
    n      = size(data,1);
    fx     = food_pos(r,1); fy = food_pos(r,2);

    % ---- Detect which zone number is actually the food arm ------------
    % Zone numbering isn't necessarily consistent across sessions/mice
    % (e.g. zone 1 can be the food arm for one animal, zone 3 for
    % another) -- rather than assuming a fixed number, find it directly
    % from which of zones 1/2/3 actually gets closest to food.
    best_min_dist = Inf; food_zone = NaN;
    for zk = 1:3
        d = dist(zone == zk);
        d = d(~isnan(d));
        if ~isempty(d)
            m = min(d);
            if m < best_min_dist
                best_min_dist = m;
                food_zone = zk;
            end
        end
    end
    if isnan(food_zone) || best_min_dist > opts.food_dist_thresh
        food_zone = opts.default_food_zone;
        fprintf('  Note: %s -- no zone got within opts.food_dist_thresh of food (closest=%.1f); defaulting to zone %d as the food arm.\n', ...
            mice_all{r,1}, best_min_dist, food_zone);
    end
    up_zones = setdiff(1:3, food_zone);

    path_x = NaN(n,1); path_y = NaN(n,1);

    % The two non-food arm-zones: straight "up" toward zone 4 -- a
    % corridor only has one axis to move along, no need to estimate a
    % direction from data.
    m_up = ismember(zone, up_zones);
    path_x(m_up) = 0;
    path_y(m_up) = opts.up_sign;

    % The food-arm zone: point at the actual food position (the one
    % case where the exact target matters, since food sits partway down
    % the arm). Uses NOSE position (not centroid) as "current position"
    % -- the zone label itself is presumably keyed to wherever the nose
    % is, so measuring the delta from the nose keeps position and zone
    % label consistent, especially right at zone boundaries.
    m_food = (zone == food_zone);
    if any(m_food) && ~any(isnan(food_pos(r,:)))
        path_x(m_food) = fx - x_nose(m_food);
        path_y(m_food) = fy - y_nose(m_food);
    end

    % Zone 4: purely horizontal, toward whichever side has food --
    % EXCEPT in the "turn square": the part of zone 4 directly above
    % the food arm, where the animal is already horizontally aligned
    % with the arm's entrance. There, continuing sideways would be
    % pointless (it's already there) -- the correct direction is the
    % same as the food-arm rule: straight at the food, which naturally
    % comes out as a diagonal-down vector. The food arm's horizontal
    % extent is estimated from this session's own food-arm nose
    % positions (so it adapts to that arm's real width);
    % opts.food_arm_half_width is only used as a fallback when a
    % session doesn't have enough of its own food-arm data to estimate
    % this.
    m4 = (zone == 4);
    if any(m4) && ~isnan(fx)
        food_arm_x = x_nose(zone == food_zone);
        food_arm_x = food_arm_x(~isnan(food_arm_x));
        if length(food_arm_x) >= 20
            sx  = sort(food_arm_x);
            n3  = length(sx);
            lo  = sx(max(1, round(0.05*n3)));
            hi  = sx(max(1, round(0.95*n3)));
            arm_half_width = (hi - lo) / 2;
        else
            arm_half_width = opts.food_arm_half_width;
        end

        in_turn = m4 & (abs(x_nose - fx) <= arm_half_width);
        horiz   = m4 & ~in_turn;

        path_x(horiz) = sign(fx - x_nose(horiz));
        path_y(horiz) = 0;

        path_x(in_turn) = fx - x_nose(in_turn);
        path_y(in_turn) = fy - y_nose(in_turn);
    end

    facing_x = x_nose - x_cent;
    facing_y = y_nose - y_cent;

    cross_z   = facing_x .* path_y - facing_y .* path_x;
    dot_p     = facing_x .* path_x + facing_y .* path_y;
    angle_deg = atan2d(cross_z, dot_p);

    facing_len = hypot(facing_x, facing_y);
    path_len   = hypot(path_x, path_y);
    angle_deg(facing_len < eps | isnan(facing_len) | ...
              path_len   < eps | isnan(path_len)) = NaN;

    if size(data,2) < COL_PATH
        data(:, COL_PATH) = NaN;
    end
    data(:, COL_PATH) = angle_deg;
    mice_all{r,4}      = data;
    n_ok = n_ok + 1;
end

fprintf('Done. Head-direction-to-shortest-path computed for %d/%d sessions.\n', n_ok, n_rows);

% ------------------------------------------------------------------ %
%  OPTIONAL: EXAMPLE QC PLOT                                          %
% ------------------------------------------------------------------ %
if opts.example
    plot_example_path_direction(mice_all, food_pos, opts, ...
        COL_X_NOSE, COL_Y_NOSE, COL_X_CENT, COL_Y_CENT, COL_ZONE, COL_DIST, COL_PATH);
end

end


% ====================================================================== %
%  EXAMPLE QC PLOT: zone-colored trajectory + sampled vectors + angle    %
% ====================================================================== %
function plot_example_path_direction(mice_all, food_pos, opts, ...
    COL_X_NOSE, COL_Y_NOSE, COL_X_CENT, COL_Y_CENT, COL_ZONE, COL_DIST, COL_PATH)

has_data   = cellfun(@(d) ~isempty(d) && size(d,1) > 10 && size(d,2) >= COL_PATH, ...
    mice_all(:,4));
valid_rows = find(has_data);

if isempty(valid_rows)
    fprintf('No valid sessions available for path-direction example plot.\n');
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

zone_colors = [0.55 0.80 0.45;    % zone 1 - green
               0.35 0.70 0.90;    % zone 2 - blue
               0.35 0.25 0.65;    % zone 3 - purple
               0.95 0.85 0.30];   % zone 4 - yellow

for i = 1:length(ex_rows)
    r    = ex_rows(i);
    data = mice_all{r,4};
    name = mice_all{r,1};

    t          = data(:,1);
    x_nose     = data(:, COL_X_NOSE);
    y_nose     = data(:, COL_Y_NOSE);
    x_cent     = data(:, COL_X_CENT);
    y_cent     = data(:, COL_Y_CENT);
    zone       = data(:, COL_ZONE);
    dist       = data(:, COL_DIST);
    angle_path = data(:, COL_PATH);

    % Recompute which zone was detected as the food arm for this
    % session, purely for display in the title (same logic as the main
    % loop) -- lets you visually confirm the detection is correct.
    best_min_dist = Inf; food_zone = NaN;
    for zk = 1:3
        d = dist(zone == zk);
        d = d(~isnan(d));
        if ~isempty(d)
            m = min(d);
            if m < best_min_dist
                best_min_dist = m;
                food_zone = zk;
            end
        end
    end
    if isnan(food_zone), food_zone = opts.default_food_zone; end

    figure('Name', sprintf('Path direction QC: %s', strrep(name,'_',' ')), ...
        'Position', [50 50 1300 520]);

    % ---- Panel 1: zone-colored trajectory + sampled vectors -----------
    subplot(1,2,1); hold on;
    for zk = 1:4
        mz = zone == zk;
        if any(mz)
            plot(x_nose(mz), y_nose(mz), '.', 'Color', zone_colors(zk,:), ...
                'MarkerSize', 4, 'DisplayName', sprintf('Zone %d', zk));
        end
    end
    if ~any(isnan(food_pos(r,:)))
        plot(food_pos(r,1), food_pos(r,2), 'p', 'MarkerSize',16, ...
            'MarkerFaceColor',[0.9 0.2 0.1], 'MarkerEdgeColor','k', 'LineWidth',1, ...
            'DisplayName','Food position');
    end

    n_frames = size(data,1);
    samp_idx = round(linspace(1, n_frames, min(opts.n_vectors, n_frames)));
    samp_idx = samp_idx(~isnan(angle_path(samp_idx)));

    for k = samp_idx(:)'
        nx = x_nose(k); ny = y_nose(k);
        cx = x_cent(k); cy = y_cent(k);
        fvx = nx - cx; fvy = ny - cy;   % facing direction: nose - centroid
        flen = hypot(fvx, fvy);
        % Vectors are drawn originating from the NOSE, matching the
        % nose-referenced path-vector computation above.
        quiver(nx, ny, fvx, fvy, 0, 'Color',[0.15 0.35 0.75], ...
            'LineWidth',1.2, 'MaxHeadSize',2, 'HandleVisibility','off');

        ang = angle_path(k);
        if ~isnan(ang) && flen > 0
            th = deg2rad(ang);
            pvx = fvx*cos(th) - fvy*sin(th);
            pvy = fvx*sin(th) + fvy*cos(th);
            plen = hypot(pvx, pvy);
            if plen > 0
                draw_len = 3 * flen;
                quiver(nx, ny, pvx/plen*draw_len, pvy/plen*draw_len, 0, ...
                    'Color',[0.10 0.55 0.15], 'LineWidth',1.3, 'MaxHeadSize',2, ...
                    'HandleVisibility','off');
            end
        end
    end
    plot(x_nose(samp_idx), y_nose(samp_idx), 'k.', 'MarkerSize',8, 'HandleVisibility','off');

    axis equal; grid on; box off;
    xlabel('x'); ylabel('y');
    title('Facing (blue) vs shortest-path direction (dark green)', 'FontSize',9);
    legend('Location','bestoutside', 'FontSize',7);

    % ---- Panel 2: angle over time --------------------------------------
    subplot(1,2,2); hold on;
    plot(t, angle_path, 'Color',[0.10 0.55 0.15], 'LineWidth',0.9);
    yline(0,   'Color',[0.2 0.7 0.3], 'LineStyle','--', 'LineWidth',1.0);
    yline(180, 'Color',[0.8 0.2 0.2], 'LineStyle',':',  'LineWidth',0.6);
    yline(-180,'Color',[0.8 0.2 0.2], 'LineStyle',':',  'LineWidth',0.6);
    ylim([-190 190]);
    xlabel('Time (s)'); ylabel('Angle to shortest path (deg)');
    title('0 = facing shortest path toward food', 'FontSize',9);
    grid on; box off;

    sgtitle(sprintf('%s   (detected food zone = %d)', strrep(name,'_',' '), food_zone), ...
        'FontSize',11, 'Interpreter','none');
end

end