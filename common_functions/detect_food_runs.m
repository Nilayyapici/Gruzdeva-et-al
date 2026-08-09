function [run_data] = detect_food_runs(mice_all, food_area, cfg)
% DETECT_FOOD_RUNS  Unified run detection for single-arm and 3-arm maze experiments.
%
% Identical detection logic is applied to both paradigms so that results
% are directly comparable across experiments.
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   run_data = detect_food_runs(mice_all, food_area)
%   run_data = detect_food_runs(mice_all, food_area, cfg)
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   mice_all   - Cell array of mouse data. Two layouts are supported:
%
%     SINGLE-ARM (paradigm = 'single'):
%       Column 1 : session_info string  (e.g. 'M01_sess1')
%       Column 2 : group label          (e.g. 'fasted')
%       Column 3 : food/source type     (e.g. 'gel', 'food')
%       Column 4 : data matrix          (N x 11, columns defined below)
%       Column 5 : signal_quality       (optional string)
%       Column 6 : discovery frame      (optional)
%
%       Data matrix columns (single-arm):
%         1=time, 2=x, 3=y, 4=speed, 5=dist_to_food, 6=path_from_last_visit,
%         7=door_status, 8=interaction, 9=eating, 10=grooming, 11=dff
%
%     THREE-ARM (paradigm = 'three'):
%       Column 1 : session_info string  (e.g. 'M01_sess1')
%       Column 2 : group label          (e.g. 'control' / 'CNO')
%       Column 3 : food_arm identifier
%       Column 4 : data matrix          (N x 13+, columns defined below)
%       Column 5 : (unused / reserved)
%       Column 6 : memory field         (optional string)
%
%       Data matrix columns (three-arm):
%         1=time, 2=x, 3=y, 4=signal465, 5=signal405, 6=dff, 7=speed,
%         8=? (reserved), 9=dist_food, 10=dist_nonfood1, 11=dist_nonfood2,
%         12=door_status, 13=grooming
%
%   food_area  - Scalar. Distance threshold that defines the "food zone".
%                A towards run is kept only if min(distance) <= food_area
%                (mouse reached the food zone). An away run is kept only
%                if the first points are within food_area (mouse left from
%                the food zone). See also cfg.min_run_extent.
%
%   cfg        - (optional) Struct with any of the following fields.
%                Missing fields take the defaults shown.
%
%     .paradigm        'single' | 'three'          default: 'single'
%     .group_filter    group string | 'all'         default: 'all'
%     .source_filter   source string | 'all'        default: 'all'
%                      Single-arm only; ignored for three-arm (which has
%                      one food type). Matches mice_all column 3, which
%                      contains the source label (e.g. 'gel', 'food').
%                      Equivalent to the old positional argument:
%                        analyze_single_arm_runs_with_speed_check(
%                          mice_all, food_area, 'fasted', 'gel', ...)
%                      is now:
%                        cfg.group_filter  = 'fasted';
%                        cfg.source_filter = 'gel';
%                        detect_food_runs(mice_all, food_area, cfg)
%                      The source label is stored in run_data.source and
%                      in each run's .arm field (e.g. 'gel', 'food').
%     .time_lim        session time limit (min)     default: Inf
%                      (three-arm only; ignored for single-arm)
%     .threshold       derivative threshold         default: 0.15
%     .min_run_extent  minimum spatial extent of a valid run, in the same
%                      units as the distance column.  default: food_area
%                      Applied as two additional keep criteria:
%                        TOWARDS: distance at run START > min_run_extent
%                                 (mouse genuinely began outside the food
%                                  zone, not already sitting inside it)
%                        AWAY:    distance at run END   > min_run_extent
%                                 (mouse actually left the food zone, not
%                                  just shuffled within it)
%                      Set to 0 to disable these checks and keep all runs
%                      that satisfy the basic food-zone criteria.
%     .dist_cutoff_hz  distance filter cutoff (Hz)  default: 0.05
%                      Physical frequency in Hz. Converted per-session to
%                      a normalised Butterworth Wn = cutoff/(fs/2) so that
%                      the filter has identical physical effect regardless
%                      of the dataset's sampling rate.
%     .deriv_cutoff_hz derivative filter cutoff(Hz) default: 0.05
%                      Same physical-Hz convention as dist_cutoff_hz.
%     .deriv_order     derivative filter order      default: 1
%     .dff_cutoff_hz   dF/F filter cutoff (Hz)      default: 0.05
%                      Same physical-Hz convention as dist_cutoff_hz.
%     .max_time_gap    artifact filter: max dt (s)  default: 1.0
%                      (single-arm only)
%     .max_spatial_dist artifact filter: max jump   default: 1000 cm
%                      (single-arm only)
%     .pixels_to_cm    pixel-to-cm conversion       default: 0.17
%     .validate_plots  show validation figures       default: false
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%   run_data - Struct array, one entry per mouse/session. Fields:
%
%     .mouse_id          string
%     .session           integer (0/1/2/3)
%     .group             string
%     .source            string  source label from mice_all col 3
%                                (single-arm: e.g. 'gel', 'food')
%                                (three-arm:  '' — not applicable)
%     .paradigm          'single' | 'three'
%     .speed_correlation scalar (r between recorded and computed speed)
%     .runs              Struct array, one entry per detected run. Fields:
%
%         .id            integer (sequential within mouse)
%         .arm           'food'|'nonfood1'|'nonfood2'  (three-arm)
%                        source label e.g. 'gel'|'food' (single-arm)
%         .type          'towards' | 'away'
%         .indices       row indices into the (filtered) data matrix
%         .time          vector
%         .xcoordinate   vector
%         .ycoordinate   vector
%         .distance      vector  (to the relevant arm end)
%         .dff           vector
%         .dff_slope     scalar  (dF/F vs distance linear slope)
%         .speed         vector  (computed from x,y)
%         .signal465     vector  (three-arm only; NaN for single-arm)
%         .signal405     vector  (three-arm only; NaN for single-arm)
%         .start_time    scalar
%         .end_time      scalar
%         .duration      scalar (s)
%
% -------------------------------------------------------------------------
% DETECTION ALGORITHM (identical for both paradigms)
% -------------------------------------------------------------------------
%   1. Smooth distance trace with a 1st-order Butterworth low-pass filter
%      at cfg.dist_cutoff_hz.
%   2. Compute the discrete derivative; smooth it with a
%      cfg.deriv_order-th-order Butterworth low-pass at cfg.deriv_cutoff_hz.
%   3. TOWARDS run: starts when derivative < -threshold (anywhere along
%      the arm). Back-trace to last point where derivative > -threshold.
%      Run ends when derivative >= -threshold.
%      Run is KEPT if ALL of:
%        a) min(distance) <= food_area            [reached the food zone]
%        b) distance at run start > min_run_extent [began outside it]
%   4. AWAY run: starts when derivative > +threshold. Back-trace similarly.
%      Run ends when derivative <= +threshold.
%      Run is KEPT if ALL of:
%        a) first <= 10 points have distance <= food_area [left from food]
%        b) distance at run end > min_run_extent          [went far enough]
%   5. (single-arm only) Artifact filtering: consecutive points separated
%      by > max_time_gap seconds or > max_spatial_dist cm are split; only
%      the longest clean segment is kept.
%
% -------------------------------------------------------------------------

    % ------------------------------------------------------------------
    % 0. Parse configuration
    % ------------------------------------------------------------------
    if nargin < 3, cfg = struct(); end

    cfg = set_default(cfg, 'paradigm',         'single');
    cfg = set_default(cfg, 'group_filter',      'all');
    cfg = set_default(cfg, 'source_filter',     'all');
    cfg = set_default(cfg, 'time_lim',          Inf);
    cfg = set_default(cfg, 'threshold',         0.1);
    cfg = set_default(cfg, 'min_run_extent',    food_area); % default = food_area
    cfg = set_default(cfg, 'dist_cutoff_hz',    0.2); %if 5Hz 0.2 
    cfg = set_default(cfg, 'deriv_cutoff_hz',   0.2); %if 5Hz 0.2 if 20 Hz 0.005
    cfg = set_default(cfg, 'deriv_order',       1);
    cfg = set_default(cfg, 'dff_cutoff_hz',     0.05);% same for 5Hz and 20Hz
    cfg = set_default(cfg, 'max_time_gap',      1.0);
    cfg = set_default(cfg, 'max_spatial_dist',  1000.0);
    cfg = set_default(cfg, 'pixels_to_cm',      0.17);
    cfg = set_default(cfg, 'validate_plots',    false);

    is_single = strcmp(cfg.paradigm, 'single');

    % ------------------------------------------------------------------
    % 1. Column map — both paradigms share time/x/y; rest differs
    % ------------------------------------------------------------------
    if is_single
        COL_TIME  = 1;  COL_X    = 2;  COL_Y   = 3;
        COL_SPEED = 4;  COL_DOOR = 7;
        COL_EATING   = 9;  COL_GROOMING = 10;  COL_DFF = 11;
        % Distance to food is column 5 (single arm, single target)
        dist_cols  = 5;
        % arm_names is set per-row from the source label (e.g. 'gel', 'food')
        COL_SIG465 = NaN;  COL_SIG405 = NaN;
    else
        COL_TIME  = 1;  COL_X    = 2;  COL_Y   = 3;
        COL_SIG465 = 4; COL_SIG405 = 5; COL_DFF = 6;
        COL_SPEED = 7;  COL_DOOR = 12; COL_GROOMING = 13;
        COL_EATING = NaN; % No eating column in three-arm
        % Distance columns: 9=food, 10=nonfood1, 11=nonfood2
        dist_cols  = [9, 10, 11];
        arm_names  = {'food', 'nonfood1', 'nonfood2'};
    end

    % ------------------------------------------------------------------
    % 2. Build filters
    % ------------------------------------------------------------------
    % THREE-ARM: filters built once here with normalised cutoffs, exactly
    % as the original analyze_mouse_runs — no change to that behaviour.
    %
    % SINGLE-ARM: the pre-loop filters below are placeholders; they are
    % rebuilt per-session inside the loop once the actual sampling rate is
    % known, so that the physical cutoff frequency is preserved regardless
    % of the recording rate (see "Estimate fs" block in the loop).
    [b_dist,  a_dist]  = butter(1,              cfg.dist_cutoff_hz,  'low');
    [b_deriv, a_deriv] = butter(cfg.deriv_order, cfg.deriv_cutoff_hz, 'low');
    [b_dff,   a_dff]   = butter(1,              cfg.dff_cutoff_hz,   'low');


    % ------------------------------------------------------------------
    % 3. Initialise output
    % ------------------------------------------------------------------
    run_data = struct('mouse_id', {}, 'session', {}, 'group', {}, 'source', {}, ...
                      'paradigm', {}, 'fs', {}, 'speed_correlation', {}, 'runs', {});

    all_speed_correlations = [];
    export_cache = struct('session_info', {}, 'source', {}, 'session_number', {}, ...
                          'time', {}, 'raw_dist', {}, 'dist_cols', {}, ...
                          'towards_periods', {}, 'away_periods', {}, 'arm_names', {});

    % ------------------------------------------------------------------
    % 4. Main loop over mice_all rows
    % ------------------------------------------------------------------
    for i = 1:size(mice_all, 1)

        session_info = mice_all{i, 1};
        group        = mice_all{i, 2};
        data         = mice_all{i, 4};

        % --- group / source filter -----------------------------------
        if is_single
            source = mice_all{i, 3};   % e.g. 'gel', 'food', 'HF'
            pass_group  = strcmp(cfg.group_filter,  'all') || strcmp(group,  cfg.group_filter);
            pass_source = strcmp(cfg.source_filter, 'all') || strcmp(source, cfg.source_filter);
            if ~(pass_group && pass_source), continue; end
            % The arm label for single-arm runs IS the source type
            arm_names = {source};
        else
            source = '';               % three-arm has one food type, no label needed
            pass_group = strcmp(cfg.group_filter, 'all') || strcmp(group, cfg.group_filter);
            if ~pass_group, continue; end
        end

        % --- session number ------------------------------------------
        session_number = 0 * contains(session_info,'sess0') + ...
                         1 * contains(session_info,'sess1') + ...
                         2 * contains(session_info,'sess2') + ...
                         3 * contains(session_info,'sess3');

        % --- door / time filtering -----------------------------------
        if is_single
            % Sessions 0 & 2: closed door; sessions 1 & 3: open door
            if session_number == 0 || session_number == 2
                data = data(data(:, COL_DOOR) == 0, :);
            elseif session_number == 1 || session_number == 3
                data = data(data(:, COL_DOOR) == 1, :);
                if ~isempty(data)
                    data(:, COL_TIME) = data(:, COL_TIME) - data(1, COL_TIME);
                end
            else
                continue;
            end
        else
            if session_number == 0 || session_number == 2
                data = data(data(:, COL_DOOR) == 0, :);
            elseif session_number == 1
                data = data(data(:, COL_DOOR) == 1, :);
                if ~isempty(data)
                    data(:, COL_TIME) = data(:, COL_TIME) - data(1, COL_TIME);
                end
            else
                continue;
            end

            % Apply time limit (three-arm only)
            if isfinite(cfg.time_lim) && ~isempty(data)
                time_limit_s = cfg.time_lim * 60;
                idx = findnearest(time_limit_s, data(:, COL_TIME));
                data = data(1:idx, :);
            end
        end

        if isempty(data), continue; end

        % --- behavioural exclusions ----------------------------------
        data = data(data(:, COL_GROOMING) == 0, :);
        % if is_single && ~isnan(COL_EATING)
        %     data = data(data(:, COL_EATING) == 0, :);
        % end

        if isempty(data) || size(data,1) < 4, continue; end

        % --- estimate fs (needed for mouse_entry and single-arm filters) ---
        fs = 1 / median(diff(data(:, COL_TIME)));
        fs = max(1, min(fs, 1000));

        % --- for single-arm: rebuild filters at actual sampling rate --
        % Three-arm filters were built pre-loop with normalised Wn and
        % are left unchanged. For single-arm the recording rate may differ,
        % so fs is estimated per-session and cutoffs converted accordingly.
        if is_single
            Wn_dist  = min(cfg.dist_cutoff_hz  / (fs/2), 0.999);
            Wn_deriv = min(cfg.deriv_cutoff_hz / (fs/2), 0.999);
            Wn_dff   = min(cfg.dff_cutoff_hz   / (fs/2), 0.999);
            [b_dist,  a_dist]  = butter(1,               Wn_dist,  'low');
            [b_deriv, a_deriv] = butter(cfg.deriv_order, Wn_deriv, 'low');
            [b_dff,   a_dff]   = butter(1,               Wn_dff,   'low');
        end

        % --- speed from position (before any smoothing) --------------
        [calc_speed, speed_corr] = compute_speed(data(:,COL_TIME), ...
                                                  data(:,COL_X), ...
                                                  data(:,COL_Y), ...
                                                  data(:,COL_SPEED), ...
                                                  cfg.pixels_to_cm);
        data(:, COL_SPEED) = calc_speed;
        if ~isnan(speed_corr)
            all_speed_correlations(end+1) = speed_corr; %#ok<AGROW>
        end

        % --- save raw distance for validation plots (before smoothing) ---
        raw_dist = data(:, dist_cols);  % N x n_arms, unfiltered

        % --- smooth signals ------------------------------------------
        for dc = dist_cols
            data(:, dc) = filtfilt(b_dist, a_dist, data(:, dc));
        end
        data(:, COL_DFF)   = filtfilt(b_dff, a_dff, data(:, COL_DFF));
        data(:, COL_SPEED) = filtfilt(b_dff, a_dff, data(:, COL_SPEED));
        if ~is_single
            data(:, COL_SIG465) = filtfilt(b_dff, a_dff, data(:, COL_SIG465));
            data(:, COL_SIG405) = filtfilt(b_dff, a_dff, data(:, COL_SIG405));
        end

        % --- derivative of distance(s) -------------------------------
        n_arms = length(dist_cols);
        raw_deriv = [diff(data(:, dist_cols)); zeros(1, n_arms)];
        if size(raw_deriv, 1) >= 4
            dist_derivative = filtfilt(b_deriv, a_deriv, raw_deriv);
        else
            dist_derivative = raw_deriv;
        end

        % --- mouse ID ------------------------------------------------
        mouse_id = extractBefore(session_info, '_sess');
        if isempty(mouse_id), mouse_id = session_info; end

        % --- initialise mouse entry ----------------------------------
        mouse_entry = struct('mouse_id', mouse_id, ...
                             'session',  session_number, ...
                             'group',    group, ...
                             'source',   source, ...
                             'paradigm', cfg.paradigm, ...
                             'fs',       fs, ...
                             'speed_correlation', speed_corr, ...
                             'runs', []);
        run_counter = 1;

        % For validation plots and export cache — always collected
        towards_periods = cell(1, n_arms);
        away_periods    = cell(1, n_arms);
        for ai = 1:n_arms
            towards_periods{ai} = zeros(0,2);
            away_periods{ai}    = zeros(0,2);
        end

        % ==============================================================
        % 5. Run detection — identical logic for every arm
        % ==============================================================
        for arm_idx = 1:n_arms

            dist_col  = dist_cols(arm_idx);
            arm_name  = arm_names{arm_idx};
            deriv_col = arm_idx; % dist_derivative is N x n_arms

            in_towards      = false;
            in_away         = false;
            towards_start   = 1;   % separate start index for towards runs
            away_start      = 1;   % separate start index for away runs

            ii = 2;
            while ii <= size(data,1) - 1

                d   = dist_derivative(ii, deriv_col);
                thr = cfg.threshold;
                ext = cfg.min_run_extent;

                % ---- TOWARDS: start detection -----------------------
                if ~in_towards && d < -thr
                    in_towards = true;
                    prev = find(dist_derivative(1:ii, deriv_col) > -thr, 1, 'last');
                    towards_start = max(1, prev + 1);

                elseif in_towards && d >= -thr
                    % ---- TOWARDS: end --------------------------------
                    in_towards = false;
                    run_indices = (towards_start:ii)';

                    % Artifact filter (single-arm only)
                    if is_single
                        run_indices = filter_artifacts(data, run_indices, ...
                            cfg.max_time_gap, cfg.max_spatial_dist, cfg.pixels_to_cm);
                    end

                    reached_food = ~isempty(run_indices) && ...
                                   min(data(run_indices, dist_col)) <= food_area;
                    started_far  = ~isempty(run_indices) && ...
                                   data(run_indices(1), dist_col) > ext;

                    if length(run_indices) >= 3 && reached_food && started_far
                        run_entry = build_run_entry(run_counter, arm_name, 'towards', ...
                            run_indices, data, dist_col, COL_DFF, COL_SPEED, ...
                            COL_SIG465, COL_SIG405, is_single);
                        mouse_entry.runs = [mouse_entry.runs; run_entry];
                        run_counter = run_counter + 1;
                        towards_periods{arm_idx}(end+1,:) = [run_indices(1), run_indices(end)];
                        if cfg.validate_plots
                            % (period already stored above for export cache)
                        end
                    end
                end

                % ---- AWAY: start detection --------------------------
                % Note: uses elseif so that a point that just ended a
                % towards run cannot simultaneously start an away run.
                % This prevents the away back-trace from reaching into
                % the previous towards epoch.
                if ~in_away && ~in_towards && d > thr
                    in_away = true;
                    prev = find(dist_derivative(1:ii, deriv_col) < thr, 1, 'last');
                    away_start = max(1, prev + 1);

                elseif in_away && d <= thr
                    % ---- AWAY: end ----------------------------------
                    in_away = false;
                    run_indices = (away_start:ii)';

                    % Artifact filter (single-arm only)
                    if is_single
                        run_indices = filter_artifacts(data, run_indices, ...
                            cfg.max_time_gap, cfg.max_spatial_dist, cfg.pixels_to_cm);
                    end

                    % Check first 10 points from away_start for food zone
                    check_end = min(away_start + 10, ii);
                    left_food = ~isempty(run_indices) && ...
                                any(data(away_start:check_end, dist_col) <= food_area);
                    ended_far = ~isempty(run_indices) && ...
                                data(run_indices(end), dist_col) > ext;

                    if length(run_indices) >= 3 && left_food && ended_far
                        run_entry = build_run_entry(run_counter, arm_name, 'away', ...
                            run_indices, data, dist_col, COL_DFF, COL_SPEED, ...
                            COL_SIG465, COL_SIG405, is_single);
                        mouse_entry.runs = [mouse_entry.runs; run_entry];
                        run_counter = run_counter + 1;
                        away_periods{arm_idx}(end+1,:) = [run_indices(1), run_indices(end)];
                    end
                end

                ii = ii + 1;
            end % while

            % ---- Handle runs still open at end of session -----------
            if in_towards
                run_indices = (towards_start:size(data,1))';
                if is_single
                    run_indices = filter_artifacts(data, run_indices, ...
                        cfg.max_time_gap, cfg.max_spatial_dist, cfg.pixels_to_cm);
                end
                reached_food = ~isempty(run_indices) && ...
                               min(data(run_indices, dist_col)) <= food_area;
                started_far  = ~isempty(run_indices) && ...
                               data(run_indices(1), dist_col) > cfg.min_run_extent;
                if length(run_indices) >= 3 && reached_food && started_far
                    run_entry = build_run_entry(run_counter, arm_name, 'towards', ...
                        run_indices, data, dist_col, COL_DFF, COL_SPEED, ...
                        COL_SIG465, COL_SIG405, is_single);
                    mouse_entry.runs = [mouse_entry.runs; run_entry];
                    run_counter = run_counter + 1;
                    towards_periods{arm_idx}(end+1,:) = [run_indices(1), run_indices(end)];
                end
            end

            if in_away
                run_indices = (away_start:size(data,1))';
                if is_single
                    run_indices = filter_artifacts(data, run_indices, ...
                        cfg.max_time_gap, cfg.max_spatial_dist, cfg.pixels_to_cm);
                end
                check_end = min(away_start + 10, size(data,1));
                left_food = ~isempty(run_indices) && ...
                            any(data(away_start:check_end, dist_col) <= food_area);
                ended_far = ~isempty(run_indices) && ...
                            data(run_indices(end), dist_col) > cfg.min_run_extent;
                if length(run_indices) >= 3 && left_food && ended_far
                    run_entry = build_run_entry(run_counter, arm_name, 'away', ...
                        run_indices, data, dist_col, COL_DFF, COL_SPEED, ...
                        COL_SIG465, COL_SIG405, is_single);
                    mouse_entry.runs = [mouse_entry.runs; run_entry];
                    run_counter = run_counter + 1;
                    away_periods{arm_idx}(end+1,:) = [run_indices(1), run_indices(end)];
                end
            end

        end % for arm_idx

        % --- Store mouse entry if any runs were found ----------------
        if ~isempty(mouse_entry.runs)
            run_data = [run_data; mouse_entry]; %#ok<AGROW>

            % Always cache raw trace + run periods for potential export
            cache_entry = struct( ...
                'session_info',    session_info, ...
                'source',          source, ...
                'session_number',  session_number, ...
                'time',            data(:, COL_TIME), ...
                'raw_dist',        raw_dist, ...
                'dist_cols',       dist_cols, ...
                'towards_periods', {towards_periods}, ...
                'away_periods',    {away_periods}, ...
                'arm_names',       {arm_names});
            export_cache = [export_cache; cache_entry]; %#ok<AGROW>

            if cfg.validate_plots
                plot_validation(data, raw_dist, dist_derivative, dist_cols, ...
                    towards_periods, away_periods, cfg.threshold, ...
                    food_area, cfg.min_run_extent, session_info, arm_names);
            end
        end

    end % for i (mice_all rows)

    % ------------------------------------------------------------------
    % 6. Summary
    % ------------------------------------------------------------------
    print_summary(run_data, cfg, all_speed_correlations);

    % ------------------------------------------------------------------
    % 7. Interactive export — ask user which session to export as a table
    % ------------------------------------------------------------------
    if ~isempty(export_cache)
        export_run_table(export_cache);
    end

end % main function


% ==========================================================================
%  LOCAL HELPER FUNCTIONS
% ==========================================================================

function s = set_default(s, field, val)
    if ~isfield(s, field), s.(field) = val; end
end

% --------------------------------------------------------------------------
function run_entry = build_run_entry(id, arm_name, run_type, run_indices, ...
        data, dist_col, COL_DFF, COL_SPEED, COL_SIG465, COL_SIG405, is_single)
% Build a standardised run struct. Signal columns that don't exist for a
% paradigm are filled with NaN vectors of the correct length.

    n = length(run_indices);

    % dF/F vs distance slope
    dff_slope = compute_dff_slope(data(run_indices, dist_col), ...
                                  data(run_indices, COL_DFF), run_type);

    % Signals present in both paradigms
    t   = data(run_indices, 1);
    x   = data(run_indices, 2);
    y   = data(run_indices, 3);
    dff = data(run_indices, COL_DFF);
    spd = data(run_indices, COL_SPEED);
    dst = data(run_indices, dist_col);

    % Paradigm-specific signals
    if is_single || isnan(COL_SIG465)
        sig465 = nan(n, 1);
        sig405 = nan(n, 1);
    else
        sig465 = data(run_indices, COL_SIG465);
        sig405 = data(run_indices, COL_SIG405);
    end

    run_entry = struct( ...
        'id',           id, ...
        'arm',          arm_name, ...
        'type',         run_type, ...
        'indices',      run_indices, ...
        'time',         t, ...
        'xcoordinate',  x, ...
        'ycoordinate',  y, ...
        'distance',     dst, ...
        'dff',          dff, ...
        'dff_slope',    dff_slope, ...
        'speed',        spd, ...
        'signal465',    sig465, ...
        'signal405',    sig405, ...
        'start_time',   t(1), ...
        'end_time',     t(end), ...
        'duration',     t(end) - t(1));
end

% --------------------------------------------------------------------------
function dff_slope = compute_dff_slope(dist_vec, dff_vec, run_type)
% Linear regression of dF/F on (signed) distance.
% Towards: x = -distance  => positive slope means signal rises approaching food
% Away:    x = +distance  => negative slope means signal falls leaving food

    valid = ~isnan(dist_vec) & ~isnan(dff_vec) & ...
             isfinite(dist_vec) & isfinite(dff_vec);
    if sum(valid) < 3
        dff_slope = NaN;
        return;
    end

    d = dist_vec(valid);
    f = dff_vec(valid);

    if strcmp(run_type, 'towards')
        x = -d;
    else
        x = d;
    end

    try
        p = polyfit(x, f, 1);
        dff_slope = p(1);
    catch
        dff_slope = NaN;
    end
end

% --------------------------------------------------------------------------
function [calc_speed, r] = compute_speed(t, x, y, orig_speed, px2cm)
% Calculate instantaneous speed (cm/s) from position and time.

    if length(t) < 2
        calc_speed = zeros(size(t));
        r = NaN;
        return;
    end

    dt  = diff(t);
    spd = sqrt(diff(x).^2 + diff(y).^2) * px2cm ./ dt;
    spd(dt < 1e-6) = 0;
    spd(~isfinite(spd)) = 0;

    win = min(5, length(spd));
    if win > 1, spd = movmean(spd, win); end

    calc_speed = [spd(1); spd];  % prepend first value to match length

    try
        r = corr(calc_speed, orig_speed);
    catch
        r = NaN;
    end
end

% --------------------------------------------------------------------------
function out_idx = filter_artifacts(data, idx, max_dt, max_dist_cm, px2cm)
% Remove segments separated by time gaps > max_dt or spatial jumps > max_dist_cm.
% Returns the longest clean contiguous segment (or [] if none is long enough).

    if length(idx) < 2
        out_idx = idx;
        return;
    end

    segments = {};
    seg = idx(1);

    for k = 2:length(idx)
        ci = idx(k);  pi = idx(k-1);
        dt = data(ci,1) - data(pi,1);
        ds = sqrt((data(ci,2)-data(pi,2))^2 + (data(ci,3)-data(pi,3))^2) * px2cm;

        if dt <= max_dt && ds <= max_dist_cm
            seg(end+1) = ci; %#ok<AGROW>
        else
            if length(seg) >= 3, segments{end+1} = seg(:); end %#ok<AGROW>
            seg = ci;
        end
    end
    if length(seg) >= 3, segments{end+1} = seg(:); end

    if isempty(segments)
        out_idx = [];
    else
        [~, best] = max(cellfun(@numel, segments));
        out_idx = segments{best};
    end
end

% --------------------------------------------------------------------------
function print_summary(run_data, cfg, all_speed_correlations)

    if isempty(run_data)
        fprintf('No valid runs found (paradigm="%s", group="%s").\n', ...
                cfg.paradigm, cfg.group_filter);
        return;
    end

    total = 0;  n_tow = 0;  n_away = 0;
    for i = 1:numel(run_data)
        for j = 1:numel(run_data(i).runs)
            total = total + 1;
            if strcmp(run_data(i).runs(j).type, 'towards'), n_tow  = n_tow  + 1; end
            if strcmp(run_data(i).runs(j).type, 'away'),    n_away = n_away + 1; end
        end
    end

    fprintf('\n===== detect_food_runs SUMMARY =====\n');
    fprintf('Paradigm : %s\n', cfg.paradigm);
    fprintf('Mice     : %d\n', numel(run_data));
    fprintf('Runs     : %d  (%d towards, %d away)\n', total, n_tow, n_away);
    fprintf('Filters  : threshold=%.3f\n', cfg.threshold);
    fprintf('           dist_lp=%.4f Hz | deriv_lp=%.4f Hz (order %d) | dff_lp=%.4f Hz\n', ...
            cfg.dist_cutoff_hz, cfg.deriv_cutoff_hz, cfg.deriv_order, cfg.dff_cutoff_hz);
    fprintf('           (cutoffs are physical Hz, normalised per-session to actual fs)\n');

    % Report sampling rates found
    all_fs = [run_data.fs];
    fprintf('Sampling : min=%.1f Hz, max=%.1f Hz, median=%.1f Hz across %d sessions\n', ...
            min(all_fs), max(all_fs), median(all_fs), numel(all_fs));
    fprintf('           min_run_extent=%.1f | towards: start>ext | away: end>ext\n', ...
            cfg.min_run_extent);
    if strcmp(cfg.paradigm, 'single')
        fprintf('           artifact: dt<=%.1fs, jump<=%.1fcm\n', ...
                cfg.max_time_gap, cfg.max_spatial_dist);
    end
    if ~isempty(all_speed_correlations)
        fprintf('Speed corr: %.3f ± %.3f (n=%d sessions)\n', ...
                mean(all_speed_correlations), std(all_speed_correlations), ...
                numel(all_speed_correlations));
    end
    fprintf('=====================================\n\n');
end

% --------------------------------------------------------------------------
function plot_validation(data, raw_dist, dist_derivative, dist_cols, ...
        towards_periods, away_periods, threshold, food_area, min_run_extent, ...
        session_info, arm_names)
% One subplot per arm.
% Background trace : raw (unfiltered) distance in grey.
% Overlay          : filtered distance in black (used for detection).
% Towards runs highlighted in blue, away runs in red — both on the raw trace.

    n_arms = length(dist_cols);
    t = data(:, 1);

    figure('Name', ['Run Validation: ' session_info], ...
           'Position', [100 100 1200 300*n_arms]);

    for ai = 1:n_arms
        subplot(n_arms, 1, ai);
        hold on;

        % Raw distance — light grey background
        plot(t, raw_dist(:, ai), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);

        % Filtered distance — black, used for detection
        plot(t, data(:, dist_cols(ai)), 'k-', 'LineWidth', 1.2);

        % Reference lines
        yline(food_area,      'k--', 'Food zone',      'LineWidth', 1.5);
        yline(min_run_extent, 'k:',  'Min run extent', 'LineWidth', 1.5);

        % Towards runs (blue) plotted on raw distance
        for p = 1:size(towards_periods{ai}, 1)
            s = towards_periods{ai}(p,1);  e = towards_periods{ai}(p,2);
            if s <= size(data,1) && e <= size(data,1)
                scatter(t(s:e), raw_dist(s:e, ai), 12, 'b', 'filled', ...
                        'MarkerFaceAlpha', 0.7);
            end
        end

        % Away runs (red) plotted on raw distance
        for p = 1:size(away_periods{ai}, 1)
            s = away_periods{ai}(p,1);  e = away_periods{ai}(p,2);
            if s <= size(data,1) && e <= size(data,1)
                scatter(t(s:e), raw_dist(s:e, ai), 12, 'r', 'filled', ...
                        'MarkerFaceAlpha', 0.7);
            end
        end

        title(sprintf('%s — %s', session_info, arm_names{ai}));
        xlabel('Time (s)');  ylabel('Distance');
        legend('Raw', 'Filtered', 'Location', 'best', 'FontSize', 8);
        box off;  hold off;
    end
end

% --------------------------------------------------------------------------
function export_run_table(export_cache)
% Interactive export: lists all cached sessions, asks the user to pick one,
% then builds a table with columns:
%   time | distance | towards (0/1) | away (0/1)
% (three-arm: one distance/towards/away triplet per arm, prefixed by arm name)
%
% The table is assigned to the base workspace as  run_table_<session_info>
% and written to  run_table_<session_info>.csv  in the current folder.
% A quick confirmation plot is also shown.

    % --- List available sessions ------------------------------------
    fprintf('\n--- Sessions available for export ---\n');
    for k = 1:numel(export_cache)
        e = export_cache(k);
        src = e.source;  if isempty(src), src = '(three-arm)'; end
        fprintf('  [%d]  %-40s  source: %-10s  session: %d\n', ...
                k, e.session_info, src, e.session_number);
    end
    fprintf('  [0]  Skip\n');

    % --- Get user input ---------------------------------------------
    choice = input('Enter number to export (0 to skip): ');
    if isempty(choice) || ~isnumeric(choice) || choice == 0 || ...
            floor(choice) ~= choice || choice > numel(export_cache)
        fprintf('Export skipped.\n');  return;
    end

    e      = export_cache(choice);
    t      = e.time;
    N      = length(t);
    n_arms = length(e.dist_cols);

    % --- Build output label -----------------------------------------
    label = strrep([e.session_info, '_', e.source], ' ', '_');
    label = strrep(label, '__', '_');   % clean up if source is empty

    % --- Assemble table columns -------------------------------------
    all_vars   = {'time'};
    all_arrays = {t};

    for ai = 1:n_arms
        arm      = e.arm_names{ai};
        dist_vec = e.raw_dist(:, ai);

        % Towards binary vector
        tow = zeros(N, 1);
        for p = 1:size(e.towards_periods{ai}, 1)
            s  = max(1, min(e.towards_periods{ai}(p,1), N));
            en = max(1, min(e.towards_periods{ai}(p,2), N));
            tow(s:en) = 1;
        end

        % Away binary vector
        aw = zeros(N, 1);
        for p = 1:size(e.away_periods{ai}, 1)
            s  = max(1, min(e.away_periods{ai}(p,1), N));
            en = max(1, min(e.away_periods{ai}(p,2), N));
            aw(s:en) = 1;
        end

        if n_arms == 1
            % Single-arm: plain names
            all_vars   = [all_vars,   {'distance', 'towards', 'away'}];
            all_arrays = [all_arrays, {dist_vec,    tow,       aw}];
        else
            % Three-arm: prefix with arm name
            all_vars   = [all_vars,   {['dist_' arm], ['towards_' arm], ['away_' arm]}];
            all_arrays = [all_arrays, {dist_vec,        tow,              aw}];
        end
    end

    T = table(all_arrays{:}, 'VariableNames', all_vars);

    % --- Save -------------------------------------------------------
    ws_name  = ['run_table_' label];
    csv_name = [ws_name '.csv'];
    assignin('base', ws_name, T);
    writetable(T, csv_name);
    fprintf('\nExported %d rows → workspace variable "%s"  and  "%s"\n\n', ...
            height(T), ws_name, csv_name);

    % --- Confirmation plot ------------------------------------------
    figure('Name', ['Export: ' label], 'Position', [150 150 1100 250*n_arms]);
    for ai = 1:n_arms
        subplot(n_arms, 1, ai);  hold on;
        arm = e.arm_names{ai};
        if n_arms == 1
            dcol = 'distance'; tcol = 'towards'; acol = 'away';
        else
            dcol = ['dist_' arm];  tcol = ['towards_' arm];  acol = ['away_' arm];
        end
        plot(T.time, T.(dcol), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        scatter(T.time(logical(T.(tcol))), T.(dcol)(logical(T.(tcol))), 10, 'b', 'filled');
        scatter(T.time(logical(T.(acol))), T.(dcol)(logical(T.(acol))), 10, 'r', 'filled');
        title(sprintf('%s — %s', label, arm));
        xlabel('Time (s)');  ylabel('Raw distance');
        legend('Distance', 'Towards', 'Away', 'Location', 'best', 'FontSize', 8);
        box off;  hold off;
    end
end