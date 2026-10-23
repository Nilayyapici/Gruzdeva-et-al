function glm_results = glm_dff_model_comparison(mice_all, opts)
%   Per-animal GLM analysis with cross-validation and
%   shuffle control.  Fits models separately for each animal and plots group-level summaries.
%
%   CALL:
%     glm_results = glm_dff_model_comparison(mice_all, opts)
%
%   OPTS FIELDS (all optional):
%     opts.state            : 'fasted' | 'fed'       (default: 'fasted')
%     opts.source           : 'food'   | 'gel'       (default: 'food')
%     opts.excl_grooming    : true | false            (default: false)
%     opts.excl_food_events : true | false            (default: false)
%     opts.excl_abs_time    : true | false            (default: false)
%     opts.n_cv_iter        : number of CV iterations (default: 100)
%     opts.cv_train_frac    : train fraction for CV   (default: 0.80)
%     opts.lag_sec          : include lagged dF/F as predictor, e.g. 5
%                             means dF/F(t - 5s) is added. (default: []  = off)
%     opts.excl_lag         : true | false  exclude lag predictor even if
%                             opts.lag_sec is set. (default: false)
%     opts.excl_elapsed     : true | false  exclude elapsed time from food
%                             discovery (column 18) from all models.
%                             NOTE: within one session this is time minus a
%                             constant, i.e. IDENTICAL to the time predictor
%                             after z-scoring, so whenever time is in the
%                             model it is dropped automatically (zero
%                             variance -> contribution NaN, and the
%                             NO_ELAPSED model equals FULL) with a summary
%                             note printed. It only acts as a separate
%                             predictor if opts.excl_abs_time = true, where
%                             it simply stands in for time. (default: false)
%     opts.excl_travel      : true | false  exclude total distance traveled
%                             (column 19) from all models. Nearly collinear
%                             with time, so its unique contribution is
%                             small and hard to interpret; the correlation
%                             with time is printed. (default: false)
%
%   Columns 18 and 19 are written by add_elapsed_distance_columns -- run
%   that first (or set opts.excl_elapsed = true and opts.excl_travel = true).
%
%     opts.headdirection    : true | false  add cos(head direction to food)
%                             (column 16, written by add_head_direction_column)
%                             to the FULL model, plus two models:
%                             NO_HEADDIR (FULL minus it) and HEADDIR_ONLY
%                             (cos head direction as the only predictor).
%                             cos = +1 head pointing at the food, -1 pointing
%                             away; it carries no left/right information.
%                             (default: false)
%     opts.min_food_dist    : [] | number.  Use ONLY frames whose distance to
%                             food (column 5, its own units) is >= this value;
%                             closer frames are removed from training AND
%                             prediction (CV, shuffles, R2, everything).
%                             E.g. 100 to study how the models work far from
%                             the food.  The time-based variables (temporal
%                             distance, time since/to food, the lagged dF/F)
%                             and the baseline correction are computed on the
%                             full session BEFORE the filter, so they stay
%                             correct.  dF/F and the continuous predictors are
%                             re-standardised on the retained frames.
%                             Predictors with no variance among them (e.g.
%                             food_visit and eating, which are 0 far from
%                             food) are dropped automatically.  Compare the
%                             models with each other and with their shuffles
%                             WITHIN a run -- not R2 values against an
%                             unfiltered run, which has a different range of
%                             distances.  (default: [] = use all frames)
%
%   CORRELATED PREDICTORS (e.g. spatial distance and path from last visit)
%     The "unique" contribution of a predictor (R2 lost when it alone is
%     removed) is small for each of two redundant predictors, because the
%     other covers for it.  Three tools quantify this:
%       1. printed list of predictor pairs with median |r| >= 0.5
%       2. commonality table for spatial distance vs path: unique, shared
%          and joint R2, in-sample and cross-validated (uses the
%          NO_SPAT_PATH model)
%       3. Shapley (dominance) shares of the FULL in-sample R2: each
%          predictor's average gain over all orders of adding predictors;
%          sums to R2.  Figure 3b + Shapley_*.xlsx.
%
%   OUTPUT STRUCT fields (glm_results):
%     .shap_r2, .shap_total, .commonality, .pair_r_keys/.pair_r_vals
%     are the results of the three tools above.
%     .per_animal   - per-animal R2, betas, MSE (real and shuffled)
%     .group        - mean/SEM across animals, t-test vs shuffled
%
%   CROSS-VALIDATION:
%     For each animal x model x iteration:
%       1. Randomly assign 80% of frames to train, 20% to test
%       2. Fit GLM on train set
%       3. Predict on test set, compute MSE
%     Repeated 100 times; mean MSE reported per animal.
%     Same procedure applied to shuffled dF/F (circular shift within session).
%
%   SHUFFLE CONTROL:
%     dF/F is circularly shifted by a random offset (10-90% of session
%     length), preserving autocorrelation while destroying alignment with
%     predictors.  Used consistently for both the full-dataset shuffle R2
%     (averaged over N_SHUF_R2=20 shifts) and the CV shuffle MSE.

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_TIME  = 1;
COL_SPEED = 4;
COL_DIST  = 5;
COL_PATH_LAST_V = 6;   % path length traveled since last food visit
COL_FOOD  = 8;
COL_EAT   = 9;
COL_GROOM = 10;
COL_DFF   = 11;
COL_ACCEL = 15;         % acceleration = d(speed)/dt
COL_HD    = 16;         % head direction to food (deg; 0 = facing food)
COL_ELAPSED = 18;       % elapsed time from food discovery (s)
COL_TRAVEL  = 19;       % total distance traveled since recording start
FOOD_DIST_THRESH = 5;

% Number of independent circular shifts used to estimate shuffle R2
% on the full dataset.  Averaging reduces single-sample noise and is
% consistent with the CV shuffle which already averages over n_cv_iter.
N_SHUF_R2 = 20;

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'state'),            opts.state            = 'fasted'; end
if ~isfield(opts,'source'),           opts.source           = 'food';   end
if ~isfield(opts,'excl_grooming'),    opts.excl_grooming    = false;    end
if ~isfield(opts,'excl_food_events'), opts.excl_food_events = false;    end
if ~isfield(opts,'excl_abs_time'),    opts.excl_abs_time    = false;    end
if ~isfield(opts,'n_cv_iter'),        opts.n_cv_iter        = 100;      end
if ~isfield(opts,'cv_train_frac'),    opts.cv_train_frac    = 0.80;     end
if ~isfield(opts,'lag_sec'),          opts.lag_sec          = [];        end
if ~isfield(opts,'excl_lag'),         opts.excl_lag         = false;     end
if ~isfield(opts,'excl_elapsed'),     opts.excl_elapsed     = false;     end
if ~isfield(opts,'excl_travel'),      opts.excl_travel      = false;     end
if ~isfield(opts,'headdirection'),    opts.headdirection    = false;     end
if ~isfield(opts,'min_food_dist'),    opts.min_food_dist    = [];        end
if ~isempty(opts.min_food_dist) && ...
        (~isnumeric(opts.min_food_dist) || ~isscalar(opts.min_food_dist) || ...
         isnan(opts.min_food_dist) || opts.min_food_dist < 0)
    error('opts.min_food_dist must be [] (off) or a single number >= 0.');
end
if ~isfield(opts,'tsf_type'), opts.tsf_type = 'temporal_distance'; end
% opts.tsf_type controls which temporal food variable is used:
%   'temporal_distance' : min(time_since_food, time_to_food)  (default)
%                         symmetric -- seconds to nearest food visit
%   'time_since_food'   : seconds since last time dist <= threshold
%                         backward-looking only
%   'time_to_food'      : seconds until next time dist <= threshold
%                         forward-looking only
if ~isfield(opts,'baseline_correct'), opts.baseline_correct = false;    end
% opts.baseline_correct can be:
%   false             : no correction (default)
%   'sliding_mean'    : subtract sliding window mean
%   'sliding_pct'     : subtract sliding window percentile (default 8th)
%                       most common in photometry field
%   'poly'            : subtract polynomial trend (degree 2)
% opts.baseline_win_sec : window length in seconds for sliding methods
%                         (default: 60 s -- preserves slow hunger dynamics)
% opts.baseline_pct     : percentile for 'sliding_pct' method (default: 8)
if ~isfield(opts,'baseline_win_sec'), opts.baseline_win_sec = 60;      end
if ~isfield(opts,'baseline_pct'),     opts.baseline_pct     = 8;       end
if ~isfield(opts,'cv_type'),          opts.cv_type          = 'random'; end
% opts.cv_type:
%   'random'   : frames randomly split 80/20 each iteration (default)
%   'blocked'  : session divided into n_cv_iter equal time-blocks;
%                within each block the first 80% is train, last 20%
%                is test -- temporal order preserved.
if ~isfield(opts,'shuf_type'),        opts.shuf_type        = 'circular'; end
% opts.shuf_type:
%   'circular' : roll dF/F by random offset, preserves autocorrelation
%                and temporal structure -- only destroys alignment with
%                predictors. More conservative null. (default)
%   'random'   : randomly permute dF/F -- destroys all temporal structure

fprintf('\n=== Per-animal GLM: state=%s  source=%s ===\n\n', ...
    opts.state, opts.source);
if ~isempty(opts.lag_sec) && opts.lag_sec > 0
    if opts.excl_lag
        fprintf('NOTE: lag_sec=%.1f set but excl_lag=true -- lag predictor excluded.\n', opts.lag_sec);
    else
        fprintf('NOTE: lagged dF/F predictor at %.1f s added to all models.\n', opts.lag_sec);
    end
end

inc_new = {};
if ~opts.excl_elapsed, inc_new{end+1} = 'elapsed_discovery (col 18)'; end
if ~opts.excl_travel,  inc_new{end+1} = 'distance_traveled (col 19)'; end
if opts.headdirection, inc_new{end+1} = 'cos_head_dir (col 16)'; end
if ~isempty(inc_new)
    fprintf('NOTE: added to the FULL model: %s\n', strjoin(inc_new, ', '));
end
if opts.headdirection
    fprintf('NOTE: head direction on -- extra models NO_HEADDIR and HEADDIR_ONLY.\n');
end
if ~isempty(opts.min_food_dist)
    fprintf(['NOTE: distance filter on -- only frames with distance to food >= %g ', ...
        '(units of column 5) are used for training AND prediction.\n'], opts.min_food_dist);
end

if ~isequal(opts.baseline_correct, false)
    fprintf('NOTE: baseline correction = %s', ...
        char(string(opts.baseline_correct)));
    if ismember(opts.baseline_correct, {'sliding_mean','sliding_pct'})
        fprintf('  (window=%.0f s', opts.baseline_win_sec);
        if strcmp(opts.baseline_correct,'sliding_pct')
            fprintf(', percentile=%.0f', opts.baseline_pct);
        end
        fprintf(')');
    end
    fprintf('\n');
end

fprintf('Temporal variable: %s\n', opts.tsf_type);
if strcmp(opts.cv_type, 'blocked')
    fprintf('CV type:    blocked (%d time-blocks, 80%% train / 20%% test within each)\n', opts.n_cv_iter);
else
    fprintf('CV type:    random (80/20 random split, %d iterations)\n', opts.n_cv_iter);
end
fprintf('Shuffle:    %s  (R2 shuffle averaged over %d circular shifts)\n', ...
    opts.shuf_type, N_SHUF_R2);

% Suppress rank-deficiency warning: expected in sess0 where food_visit
% and eating are all-zero (mouse hasn't found food yet).
warning('off', 'stats:LinearModel:RankDefDesignMat');
warning('off', 'MATLAB:rankDeficientMatrix');
cleanupObj = onCleanup(@() warning('on', 'stats:LinearModel:RankDefDesignMat'));

% ------------------------------------------------------------------ %
%  BUILD BASE FORMULA                                                  %
% ------------------------------------------------------------------ %
if opts.excl_abs_time
    base_preds = {'speed'};
else
    base_preds = {'time', 'speed'};
end
if ~opts.excl_food_events
    base_preds{end+1} = 'food_visit';
    base_preds{end+1} = 'eating';
end
if ~opts.excl_grooming
    base_preds{end+1} = 'grooming';
end
base_preds_orig = base_preds;
base = ['dff ~ ', strjoin(base_preds, ' + ')];

% 24 models (26 with opts.headdirection): full model and systematic
% ablations/isolations
% Model naming convention:
%   FULL          = base + spatial_distance + temporal_distance + path_last_visit + acceleration
%                   + elapsed_discovery + distance_traveled
%   NO_SPAT       = FULL - spatial_distance
%   NO_TEMP       = FULL - temporal_distance
%   NO_DIST       = FULL - (spatial_distance + temporal_distance)
%   NO_FOOD_VIS   = FULL - food_visit
%   NO_EATING     = FULL - eating
%   NO_FOOD_EAT   = FULL - (food_visit + eating)
%   NO_GROOM      = FULL - grooming
%   NO_SPEED      = FULL - speed
%   DIST_ONLY     = dff ~ spatial_distance + temporal_distance
%   FOOD_ONLY     = dff ~ food_visit + eating
%   GROOM_ONLY    = dff ~ grooming
%   SPEED_ONLY    = dff ~ speed
%   NO_PATH_LAST_V = FULL - path_last_visit
%   NO_ACCEL      = FULL - acceleration
%   NO_ELAPSED    = FULL - elapsed_discovery   (equals FULL whenever time is
%                   in the model: elapsed_discovery is then dropped, see above)
%   NO_TRAVEL     = FULL - distance_traveled
%   PATH_ONLY     = dff ~ path_last_visit
%   ACCEL_ONLY    = dff ~ acceleration
%   ELAPSED_ONLY  = dff ~ elapsed time from discovery (column 18).  Within a
%                   session this is time minus a constant, so it is the same
%                   model as "time only": identical R2 and CV R2.
%   TRAVEL_ONLY   = dff ~ distance_traveled (column 19)
%   NO_SPAT_PATH  = FULL - (spatial_distance + path_last_visit)
%                   spatial distance and path since the last visit are
%                   correlated, so each one's UNIQUE contribution
%                   (FULL - NO_SPAT, FULL - NO_PATH_LAST_V) under-states
%                   it; this model measures what the two lose TOGETHER.
%                   See the commonality table printed after the group
%                   statistics (unique / shared / joint).
model_keys = {'FULL','NO_SPAT','NO_TEMP','NO_DIST', ...
    'NO_FOOD_VIS','NO_EATING','NO_FOOD_EAT', ...
    'NO_GROOM','NO_SPEED', ...
    'SPAT_ONLY','TEMP_ONLY','DIST_ONLY','FOOD_ONLY','GROOM_ONLY','SPEED_ONLY', ...
    'PATH_ONLY','ACCEL_ONLY','ELAPSED_ONLY','TRAVEL_ONLY', ...
    'NO_PATH_LAST_V','NO_SPAT_PATH','NO_ACCEL','NO_ELAPSED','NO_TRAVEL'};
model_labs = {'Full model', ...
    'No spat dist','No temp dist','No spat+temp dist', ...
    'No food visit','No eating','No food+eating', ...
    'No grooming','No speed', ...
    'Spat dist only','Temp dist only','Spat+Temp dist only', ...
    'Food+Eating only','Grooming only','Speed only', ...
    'Path from last visit only','Acceleration only', ...
    'Time from discovery only','Total distance only', ...
    'No path from last visit','No spat dist + path','No acceleration', ...
    'No elapsed from discovery','No distance traveled'};

% Head-direction models (opts.headdirection) are appended at the END so the
% positions of all other models are unchanged.
%   NO_HEADDIR    = FULL - cos_head_dir
%   HEADDIR_ONLY  = dff ~ cos_head_dir
if opts.headdirection
    model_keys{end+1} = 'NO_HEADDIR';
    model_labs{end+1} = 'No head direction (cos)';
    model_keys{end+1} = 'HEADDIR_ONLY';
    model_labs{end+1} = 'Head direction (cos) only';
end

% Determine which temporal variable column to use
switch opts.tsf_type
    case 'temporal_distance'
        tp_var = 'temporal_distance';
    case 'time_to_food'
        tp_var = 'time_to_food';
    otherwise  % 'time_since_food'
        tp_var = 'time_since_food';
end

% NOTE: formulas are built per-session inside the loop because the
% full predictor list depends on which predictors are available (e.g.
% temporal_distance may be all-NaN in sess0) and on opts flags.
% The variable 'formulas' is a placeholder here; it is overwritten
% inside the loop before any fitting occurs.
formulas = struct();   % placeholder, rebuilt per session below

n_models = length(model_keys);

% ------------------------------------------------------------------ %
%  FILTER DATA                                                         %
% ------------------------------------------------------------------ %
is_state  = strcmp(mice_all(:,2), opts.state);
is_source = strcmp(mice_all(:,3), opts.source);
is_sess0  = contains(mice_all(:,1), '_sess0');
is_sess1  = contains(mice_all(:,1), '_sess1');
keep      = is_state & is_source & (is_sess0 | is_sess1);
mice_sub  = mice_all(keep, :);
if sum(keep)==0
    error('No sessions found for state="%s" source="%s".', opts.state, opts.source);
end

% Columns 18/19 must exist if either new predictor is requested
if ~opts.excl_elapsed || ~opts.excl_travel
    cols_ok = cellfun(@(d) isempty(d) || size(d,2) >= COL_TRAVEL, mice_sub(:,4));
    if ~all(cols_ok)
        error(['%d of %d selected sessions have fewer than %d columns. Run ', ...
               'add_elapsed_distance_columns first, or set opts.excl_elapsed = true ', ...
               'and opts.excl_travel = true.'], sum(~cols_ok), numel(cols_ok), COL_TRAVEL);
    end
end

% Column 16 must exist if head direction is requested
if opts.headdirection
    hd_ok = cellfun(@(d) isempty(d) || size(d,2) >= COL_HD, mice_sub(:,4));
    if ~all(hd_ok)
        error(['%d of %d selected sessions have fewer than %d columns. Run ', ...
               'add_head_direction_column first, or set opts.headdirection = false.'], ...
               sum(~hd_ok), numel(hd_ok), COL_HD);
    end
end

all_names   = mice_sub(:,1);
mouse_ids   = regexprep(all_names, '_sess\d$', '');
unique_mice = unique(mouse_ids, 'stable');
n_mice      = length(unique_mice);

fprintf('%d sessions, %d animals\n\n', sum(keep), n_mice);

% ------------------------------------------------------------------ %
%  PRE-ALLOCATE RESULTS                                               %
% ------------------------------------------------------------------ %
r2_real   = NaN(n_mice, n_models, 2);   % in-sample adjusted R2
r2_shuf   = NaN(n_mice, n_models, 2);   % shuffle in-sample R2
mse_real  = NaN(n_mice, n_models, 2);   % CV MSE (real)
mse_shuf  = NaN(n_mice, n_models, 2);   % CV MSE (shuffled)
cv_r2_real = NaN(n_mice, n_models, 2);  % CV R2 on held-out test set (real)
cv_r2_shuf = NaN(n_mice, n_models, 2);  % CV R2 on held-out test set (shuffled)

n_elapsed_dropped = 0;    % sessions where elapsed_discovery duplicated time
r_travel_time     = [];   % corr(distance_traveled, time) per session fit

% Predictor contribution list
pred_names = {};
if ~opts.excl_abs_time,    pred_names{end+1} = 'time'; end
pred_names{end+1} = 'speed';
pred_names{end+1} = 'spatial_distance';
pred_names{end+1} = tp_var;
pred_names{end+1} = 'path_last_visit';
pred_names{end+1} = 'acceleration';
if ~opts.excl_elapsed, pred_names{end+1} = 'elapsed_discovery'; end
if ~opts.excl_travel,  pred_names{end+1} = 'distance_traveled'; end
if opts.headdirection, pred_names{end+1} = 'cos_head_dir'; end
if ~opts.excl_food_events
    pred_names{end+1} = 'food_visit';
    pred_names{end+1} = 'eating';
end
if ~opts.excl_grooming,    pred_names{end+1} = 'grooming'; end
if ~isempty(opts.lag_sec) && opts.lag_sec > 0 && ~opts.excl_lag
    pred_names{end+1} = 'dff_lag';
end
n_preds    = length(pred_names);
r2_contrib = NaN(n_mice, n_preds, 2);

% Shapley (dominance) decomposition of the in-sample R2 of the FULL model:
% each predictor's share, averaged over every order in which predictors
% could be added.  Unlike the unique contribution above, the shares SUM to
% the full model's R2, so variance shared by correlated predictors is
% divided between them instead of being credited to nobody.
shap_r2    = NaN(n_mice, n_preds, 2);
shap_total = NaN(n_mice, 2);       % plain in-sample R2 of FULL (= sum of shares)
MAX_SHAP_PREDS = 13;               % 2^13 subsets per session is still fast
has_path   = false(n_mice, 2);     % path_last_visit usable in this session (for the commonality table)

% pairwise correlations among the continuous predictors, one map per
% session type; key = 'a x b' (names sorted), value = r in each session
pair_r = {containers.Map('KeyType','char','ValueType','any'), ...
          containers.Map('KeyType','char','ValueType','any')};

% Trace candidates for sess1 — 6 models stored per animal:
%   FULL, NO_DIST, NO_TEMP, NO_SPAT (=NO_DIST alias in plot), SPAT_ONLY, DIST_ONLY
% Each entry: struct with fields time, dff, preds{6}, r2, animal
trace_model_plot_keys = {'FULL','NO_SPAT','NO_TEMP','NO_DIST','SPAT_ONLY','TEMP_ONLY'};

trace_plot_labs = {'Full model', ...
    'No spatial distance', ...
    'No temporal distance', ...
    'No spatial + temporal distance', ...
    'Only spatial distance', ...
    'Only temporal distance'};
n_trace_models  = length(trace_model_plot_keys);
trace_cands_s1  = {};   % filled during main loop, one entry per animal
for m = 1:n_mice
    mid = unique_mice{m};
    fprintf('Animal %d/%d: %s\n', m, n_mice, mid);

    for sess_idx = 0:1
        sess_str = sprintf('_sess%d', sess_idx);
        row = find(strcmp(mice_sub(:,1), [mid, sess_str]));
        if isempty(row), continue; end

        data = mice_sub{row, 4};
        data(isinf(data)) = NaN;

        time_vec = data(:, COL_TIME);
        dist_vec = data(:, COL_DIST);

        % time_since_food
        at_food       = dist_vec <= FOOD_DIST_THRESH;
        time_since_fd = NaN(size(time_vec));
        last_t        = NaN;
        for t = 1:length(time_vec)
            if at_food(t), last_t = time_vec(t); end
            if ~isnan(last_t), time_since_fd(t) = time_vec(t) - last_t; end
        end

        % time_to_food
        time_to_fd = NaN(size(time_vec));
        next_t     = NaN;
        for t = length(time_vec):-1:1
            if at_food(t), next_t = time_vec(t); end
            if ~isnan(next_t), time_to_fd(t) = next_t - time_vec(t); end
        end

        % temporal_distance = min(tsf, ttf), NaN-safe
        % Equals 0 at food visits and grows with time away from food.
        tsf_for_min = time_since_fd;
        ttf_for_min = time_to_fd;
        tsf_for_min(isnan(tsf_for_min)) = Inf;
        ttf_for_min(isnan(ttf_for_min)) = Inf;
        temp_dist = min(tsf_for_min, ttf_for_min);
        temp_dist(isinf(temp_dist)) = NaN;  % NaN only where BOTH were NaN

        T = table();
        T.time             = time_vec;
        T.speed            = data(:, COL_SPEED);
        T.spatial_distance = dist_vec;          % renamed from 'distance'
        T.path_last_visit  = data(:, COL_PATH_LAST_V);   % path length since last food visit
        T.acceleration     = data(:, COL_ACCEL);          % d(speed)/dt
        T.food_visit       = double(data(:, COL_FOOD) > 0);
        T.eating           = double(data(:, COL_EAT)  > 0);
        T.grooming         = double(data(:, COL_GROOM) > 0);
        T.time_since_food  = time_since_fd;
        T.time_to_food     = time_to_fd;
        T.temporal_distance = temp_dist;        % renamed from 'temporal_proximity'
        T.dff              = data(:, COL_DFF);

        if ~opts.excl_elapsed, T.elapsed_discovery = data(:, COL_ELAPSED); end
        if ~opts.excl_travel,  T.distance_traveled = data(:, COL_TRAVEL);  end

        tsf_available    = sum(~isnan(T.time_since_food)) > 0;
        tp_var_local     = tp_var;
        tp_col_available = sum(~isnan(T.(tp_var_local))) > 0;
        path_available   = sum(~isnan(T.path_last_visit)) > 0;
        accel_available  = sum(~isnan(T.acceleration)) > 0;
        elapsed_available = ~opts.excl_elapsed && sum(~isnan(T.elapsed_discovery)) > 0;
        travel_available  = ~opts.excl_travel  && sum(~isnan(T.distance_traveled)) > 0;

        % cos of the angle between head direction and the vector to food:
        % +1 = pointing at the food, -1 = pointing away
        if opts.headdirection, T.cos_head_dir = cosd(data(:, COL_HD)); end
        hd_available = opts.headdirection && sum(~isnan(T.cos_head_dir)) > 0;

        % Drop NaN rows on core predictors.
        % If temporal variable is unavailable, simply omit it from the
        % NaN-check -- do NOT substitute another variable.
        core = {'speed','spatial_distance','dff'};
        if tp_col_available
            core{end+1} = tp_var_local;
        end
        if path_available
            core{end+1} = 'path_last_visit';
        end
        if accel_available
            core{end+1} = 'acceleration';
        end
        if elapsed_available, core{end+1} = 'elapsed_discovery'; end
        if travel_available,  core{end+1} = 'distance_traveled'; end
        if hd_available,      core{end+1} = 'cos_head_dir'; end
        if ~opts.excl_abs_time, core{end+1} = 'time'; end
        bad = false(height(T),1);
        for v = 1:length(core), bad = bad | isnan(T.(core{v})); end
        T = T(~bad,:);

        % Raw distance to food, kept un-standardised for the optional
        % distance filter below (spatial_distance is z-scored later).
        % Not a predictor, so it never enters a formula.
        T.dist_raw = T.spatial_distance;

        % Build session-specific formulas from the actual available predictors.
        % All 24 models are derived from the same full predictor set so
        % that ablations are meaningful and consistent.
        if ~tp_col_available
            fprintf('  Note: %s all-NaN in %s%s -- temporal predictor absent from models\n', ...
                tp_var_local, mid, sess_str);
        end
        if ~path_available
            fprintf('  Note: path_last_visit all-NaN in %s%s -- predictor absent from models\n', ...
                mid, sess_str);
        end
        if ~accel_available
            fprintf('  Note: acceleration all-NaN in %s%s -- predictor absent from models\n', ...
                mid, sess_str);
        end
        if ~opts.excl_elapsed && ~elapsed_available
            fprintf('  Note: elapsed_discovery all-NaN in %s%s (no food discovery?) -- predictor absent from models\n', ...
                mid, sess_str);
        end
        if ~opts.excl_travel && ~travel_available
            fprintf('  Note: distance_traveled all-NaN in %s%s -- predictor absent from models\n', ...
                mid, sess_str);
        end
        if opts.headdirection && ~hd_available
            fprintf('  Note: head direction (col %d) all-NaN in %s%s -- predictor absent from models\n', ...
                COL_HD, mid, sess_str);
        end

        % Helper: build formula string from predictor cell array
        mk = @(preds) ['dff ~ ', strjoin(preds, ' + ')];
        % Helper: remove one or more predictors
        rm = @(preds, excl) preds(~ismember(preds, excl));

        % Full predictor list for this session
        full_preds = base_preds;   % time/speed/food_visit/eating/grooming from opts
        full_preds{end+1} = 'spatial_distance';
        if tp_col_available
            full_preds{end+1} = tp_var_local;
        end
        if path_available
            full_preds{end+1} = 'path_last_visit';
        end
        if accel_available
            full_preds{end+1} = 'acceleration';
        end
        if elapsed_available
            full_preds{end+1} = 'elapsed_discovery';
        end
        if travel_available
            full_preds{end+1} = 'distance_traveled';
        end
        if hd_available
            full_preds{end+1} = 'cos_head_dir';
        end
        % (lag added below after lagged predictor is computed)

        sess_formulas = struct();
        sess_formulas.FULL        = mk(full_preds);
        sess_formulas.NO_SPAT     = mk(rm(full_preds, {'spatial_distance'}));
        sess_formulas.NO_TEMP     = mk(rm(full_preds, {tp_var_local}));
        sess_formulas.NO_DIST     = mk(rm(full_preds, {'spatial_distance', tp_var_local}));
        sess_formulas.NO_FOOD_VIS = mk(rm(full_preds, {'food_visit'}));
        sess_formulas.NO_EATING   = mk(rm(full_preds, {'eating'}));
        sess_formulas.NO_FOOD_EAT = mk(rm(full_preds, {'food_visit','eating'}));
        sess_formulas.NO_GROOM    = mk(rm(full_preds, {'grooming'}));
        sess_formulas.NO_SPEED    = mk(rm(full_preds, {'speed'}));
        sess_formulas.NO_PATH_LAST_V = mk(rm(full_preds, {'path_last_visit'}));
        sess_formulas.NO_ACCEL    = mk(rm(full_preds, {'acceleration'}));
        sess_formulas.NO_ELAPSED  = mk(rm(full_preds, {'elapsed_discovery'}));
        sess_formulas.NO_TRAVEL   = mk(rm(full_preds, {'distance_traveled'}));
        sess_formulas.NO_SPAT_PATH = mk(rm(full_preds, {'spatial_distance','path_last_visit'}));
        sess_formulas.PATH_ONLY   = 'dff ~ path_last_visit';
        sess_formulas.ACCEL_ONLY  = 'dff ~ acceleration';
        % 'elapsed_only' is an untouched copy of elapsed_discovery made
        % before the collinearity safeguard below zeroes the latter (when
        % time is in the model); it is missing if elapsed time is excluded
        % or unavailable, and the fit then fails and gives NaN.
        sess_formulas.ELAPSED_ONLY = 'dff ~ elapsed_only';
        sess_formulas.TRAVEL_ONLY = 'dff ~ distance_traveled';
        if opts.headdirection
            sess_formulas.NO_HEADDIR   = mk(rm(full_preds, {'cos_head_dir'}));
            sess_formulas.HEADDIR_ONLY = 'dff ~ cos_head_dir';
        end
        % Isolation models (ignore base, use only named predictors)
        dist_preds = {'spatial_distance'};
        if tp_col_available, dist_preds{end+1} = tp_var_local; end
        sess_formulas.SPAT_ONLY   = 'dff ~ spatial_distance';
        sess_formulas.TEMP_ONLY   = ['dff ~ ', tp_var_local];
        sess_formulas.DIST_ONLY   = mk(dist_preds);
        sess_formulas.FOOD_ONLY   = 'dff ~ food_visit + eating';
        sess_formulas.GROOM_ONLY  = 'dff ~ grooming';
        sess_formulas.SPEED_ONLY  = 'dff ~ speed';

        if height(T) < 50
            fprintf('  Skipping %s%s: only %d rows\n', mid, sess_str, height(T));
            continue;
        end


        % Baseline correction
        if ~isequal(opts.baseline_correct, false)
            T.dff = apply_baseline_correction(T.dff, T.time, opts);
        end

        % Z-score dF/F and continuous predictors
        dff_std = std(T.dff,'omitnan');
        if dff_std > 0
            T.dff = (T.dff - mean(T.dff,'omitnan')) / dff_std;
        end
        cont_vars = {'spatial_distance', 'speed'};
        if tp_col_available
            cont_vars{end+1} = tp_var_local;
        end
        if path_available
            cont_vars{end+1} = 'path_last_visit';
        end
        if accel_available
            cont_vars{end+1} = 'acceleration';
        end
        if elapsed_available, cont_vars{end+1} = 'elapsed_discovery'; end
        if travel_available,  cont_vars{end+1} = 'distance_traveled'; end
        if hd_available,      cont_vars{end+1} = 'cos_head_dir'; end
        if ~opts.excl_abs_time, cont_vars{end+1} = 'time'; end
        for v = 1:length(cont_vars)
            c = T.(cont_vars{v}); s = std(c,'omitnan');
            if s > 0, T.(cont_vars{v}) = (c - mean(c,'omitnan')) / s; end
        end

        % Collinearity with time (checked AFTER z-scoring so scale is
        % irrelevant).  Within one session elapsed_discovery = time -
        % constant, so its z-score is IDENTICAL to time's: keeping both
        % makes the design rank-deficient and every contribution ~0.
        % Replace it by a constant -> zero variance -> safe_fitlm drops
        % it, its contribution is reported as NaN (not estimable) instead
        % of a misleading 0, and NO_ELAPSED then equals FULL.  Skipped
        % when time is excluded (opts.excl_abs_time), where it just
        % stands in for time.
        if elapsed_available
            T.elapsed_only = T.elapsed_discovery;   % copy for ELAPSED_ONLY (z-scored, never zeroed)
        end
        if elapsed_available && ~opts.excl_abs_time
            r_te = corr(T.elapsed_discovery, T.time, 'rows', 'complete');
            if abs(r_te) > 0.999
                T.elapsed_discovery = zeros(height(T), 1);
                n_elapsed_dropped = n_elapsed_dropped + 1;
            end
        end
        % distance_traveled is monotone in time, so usually very highly
        % (but not perfectly) correlated with it -- kept, r reported.
        if travel_available && ~opts.excl_abs_time
            r_tt = corr(T.distance_traveled, T.time, 'rows', 'complete');
            r_travel_time(end+1) = r_tt; %#ok<AGROW>
        end

        % Lagged dF/F predictor (optional)
        sess_formulas_lag = sess_formulas;
        if ~isempty(opts.lag_sec) && opts.lag_sec > 0 && ~opts.excl_lag
            lag_frames = round(opts.lag_sec / median(diff(T.time)));
            lag_frames = max(lag_frames, 1);
            dff_lag    = NaN(height(T), 1);
            dff_lag(lag_frames+1:end) = T.dff(1:end-lag_frames);
            T.dff_lag = dff_lag;
            T = T(~isnan(T.dff_lag), :);
            lg = T.dff_lag; s_lg = std(lg,'omitnan');
            if s_lg > 0, T.dff_lag = (lg - mean(lg,'omitnan')) / s_lg; end
            % Append dff_lag to every model formula
            fn_lag = fieldnames(sess_formulas_lag);
            for klag = 1:length(fn_lag)
                sess_formulas_lag.(fn_lag{klag}) = ...
                    [sess_formulas_lag.(fn_lag{klag}), ' + dff_lag'];
            end
        end

        % ----- DISTANCE FILTER (opts.min_food_dist) --------------------
        % Keep only frames at least opts.min_food_dist from the food; all
        % others are removed from training AND prediction, because every
        % fit, CV fold and shuffle below uses T.  This happens HERE, after
        % the baseline correction, the temporal-distance / time-since-food
        % variables and the lag were computed on the full session, so they
        % still refer to the real timeline.  dF/F, the continuous
        % predictors and the lag are then re-standardised on the retained
        % frames.
        if ~isempty(opts.min_food_dist)
            n_before = height(T);
            T = T(T.dist_raw >= opts.min_food_dist, :);
            fprintf('  %s%s: distance filter (>= %g) kept %d of %d frames (%.0f%%)\n', ...
                mid, sess_str, opts.min_food_dist, height(T), n_before, ...
                100 * height(T) / max(n_before, 1));
            if height(T) < 50
                fprintf('  Skipping %s%s: only %d frames left after the distance filter\n', ...
                    mid, sess_str, height(T));
                continue;
            end
            restd_vars = [{'dff'}, cont_vars];
            if ismember('dff_lag', T.Properties.VariableNames)
                restd_vars{end+1} = 'dff_lag';
            end
            for v = 1:numel(restd_vars)
                c = T.(restd_vars{v}); s = std(c, 'omitnan');
                if s > 0, T.(restd_vars{v}) = (c - mean(c, 'omitnan')) / s; end
            end
        end

        % Pairwise correlations among the continuous predictors that are
        % actually in the models (zero-variance ones, e.g. a dropped
        % elapsed_discovery, are skipped).  Computed on the frames the
        % models actually use (after the distance filter, if any).
        % Summarised after the main loop: this is the direct check of how
        % redundant the predictors are.
        has_path(m, sess_idx+1) = path_available;
        cv_names = cont_vars(cellfun(@(vn) std(T.(vn), 'omitnan') > 0, cont_vars));
        if numel(cv_names) >= 2
            Rc = corr(T{:, cv_names}, 'rows', 'complete');
            for ia = 1:numel(cv_names)-1
                for ib = ia+1:numel(cv_names)
                    pkey = strjoin(sort({cv_names{ia}, cv_names{ib}}), ' x ');
                    pm   = pair_r{sess_idx+1};
                    if isKey(pm, pkey)
                        pm(pkey) = [pm(pkey), Rc(ia, ib)]; %#ok<NASGU>
                    else
                        pm(pkey) = Rc(ia, ib);
                    end
                end
            end
        end

        if height(T) < 20, continue; end

        % ----- FIT ALL MODELS (real data, full dataset R2) ----------
        for ki = 1:n_models
            key = model_keys{ki};
            f   = sess_formulas_lag.(key);
            try
                mdl = safe_fitlm(T, f);
                r2_real(m, ki, sess_idx+1) = mdl.Rsquared.Adjusted;
            catch
            end
        end

        % ----- STORE TRACE CANDIDATES (sess1 only) ------------------
        % Fit the 6 display models and store time-sorted predictions
        % so Figure 5 can plot actual vs predicted for every animal.
        if sess_idx == 1
            try
                [~, tidx] = sort(T.time);
                preds_store = cell(1, n_trace_models);
                for pk = 1:n_trace_models
                    key_k = trace_model_plot_keys{pk};
                    mdl_k = safe_fitlm(T, sess_formulas_lag.(key_k));
                    preds_store{pk} = predict(mdl_k, T);
                    preds_store{pk} = preds_store{pk}(tidx);
                end
                cand        = struct();
                cand.time   = T.time(tidx);
                cand.dff    = T.dff(tidx);
                cand.preds  = preds_store;
                cand.r2     = r2_real(m, 1, sess_idx+1);
                cand.animal = mid;
                trace_cands_s1{end+1} = cand;
            catch ME_tr
                fprintf('  Warning: trace storage failed for %s sess1: %s\n', ...
                    mid, ME_tr.message);
            end
        end

        % ----- PREDICTOR CONTRIBUTIONS (R2 drop when removed) -------
        r2_full = r2_real(m, 1, sess_idx+1);   % index 1 = FULL model
        if ~isnan(r2_full)
            try
                mdl_full_check = safe_fitlm(T, sess_formulas_lag.FULL);
                fitted_preds = mdl_full_check.CoefficientNames;
                fitted_preds = fitted_preds(~strcmp(fitted_preds, '(Intercept)'));
            catch
                fitted_preds = {};
            end
            for pi = 1:n_preds
                pn = pred_names{pi};
                if ~any(strcmp(fitted_preds, pn))
                    r2_contrib(m, pi, sess_idx+1) = NaN; continue;
                end
                reduced_preds = fitted_preds(~strcmp(fitted_preds, pn));
                if isempty(reduced_preds)
                    r2_contrib(m, pi, sess_idx+1) = NaN; continue;
                end
                reduced_f = ['dff ~ ', strjoin(reduced_preds, ' + ')];
                try
                    mdl_red = safe_fitlm(T, reduced_f);
                    r2_contrib(m, pi, sess_idx+1) = ...
                        r2_full - mdl_red.Rsquared.Adjusted;
                catch
                    r2_contrib(m, pi, sess_idx+1) = NaN;
                end
            end

            % ----- SHAPLEY (DOMINANCE) DECOMPOSITION of the FULL R2 ---
            % Same predictors as the FULL fit above.  In-sample, plain R2.
            if numel(fitted_preds) >= 2 && numel(fitted_preds) <= MAX_SHAP_PREDS
                try
                    [sv, r2_plain] = shapley_r2(T, fitted_preds);
                    for pi = 1:n_preds
                        jj = find(strcmp(fitted_preds, pred_names{pi}), 1);
                        if ~isempty(jj), shap_r2(m, pi, sess_idx+1) = sv(jj); end
                    end
                    shap_total(m, sess_idx+1) = r2_plain;
                catch ME_shap
                    fprintf('  Note: Shapley decomposition failed for %s%s: %s\n', ...
                        mid, sess_str, ME_shap.message);
                end
            elseif numel(fitted_preds) > MAX_SHAP_PREDS
                fprintf('  Note: %d predictors in %s%s -- Shapley decomposition skipped (max %d)\n', ...
                    numel(fitted_preds), mid, sess_str, MAX_SHAP_PREDS);
            end
        end

        % ----- SHUFFLE R2: average over N_SHUF_R2 circular shifts ---
        n = height(T);
        min_shift = max(1, round(0.10 * n));
        max_shift = round(0.90 * n);
        for ki = 1:n_models
            f = sess_formulas_lag.(model_keys{ki});
            r2_shuf_acc = NaN(N_SHUF_R2, 1);
            for si_shuf = 1:N_SHUF_R2
                shift_amt = min_shift + randi(max_shift - min_shift + 1) - 1;
                T_shuf     = T;
                T_shuf.dff = circshift(T.dff, shift_amt);
                try
                    mdl_s = safe_fitlm(T_shuf, f);
                    r2_shuf_acc(si_shuf) = mdl_s.Rsquared.Adjusted;
                catch
                end
            end
            r2_shuf(m, ki, sess_idx+1) = nanmean(r2_shuf_acc);
        end

        % ----- CROSS-VALIDATION (real + shuffled) ------------------
        n_train        = round(n * opts.cv_train_frac);
        mse_real_iter  = NaN(opts.n_cv_iter, n_models);
        mse_shuf_iter  = NaN(opts.n_cv_iter, n_models);
        cv_r2_real_iter = NaN(opts.n_cv_iter, n_models);
        cv_r2_shuf_iter = NaN(opts.n_cv_iter, n_models);

        [~, t_order] = sort(T.time);

        if strcmp(opts.cv_type, 'blocked')
            block_size  = floor(n / opts.n_cv_iter);
            n_test_blk  = block_size - max(1, round(block_size * opts.cv_train_frac));
            if n_test_blk < 1
                warning('Block size too small for 80/20 split -- using random CV');
                cv_type_eff = 'random';
            else
                cv_type_eff = 'blocked';
            end
        else
            cv_type_eff = 'random';
        end

        for it = 1:opts.n_cv_iter
            if strcmp(cv_type_eff, 'blocked')
                blk_start = (it-1)*block_size + 1;
                blk_end   = min(it*block_size, n);
                blk_idx   = t_order(blk_start:blk_end);
                n_tr_blk  = round(length(blk_idx) * opts.cv_train_frac);
                idx_tr    = blk_idx(1:n_tr_blk);
                idx_te    = blk_idx(n_tr_blk+1:end);
            else
                perm   = randperm(n);
                idx_tr = perm(1:n_train);
                idx_te = perm(n_train+1:end);
            end
            T_tr = T(idx_tr, :);
            T_te = T(idx_te, :);

            % Circular shift for shuffled CV (same method as R2 shuffle)
            shift_amt        = min_shift + randi(max_shift - min_shift + 1) - 1;
            dff_full_shifted = circshift(T.dff, shift_amt);
            T_tr_sh          = T_tr;
            T_te_sh          = T_te;
            T_tr_sh.dff      = dff_full_shifted(idx_tr);
            T_te_sh.dff      = dff_full_shifted(idx_te);

            % Precompute test-set variance for CV R2 denominator
            var_te    = var(T_te.dff,    'omitnan');
            var_te_sh = var(T_te_sh.dff, 'omitnan');

            for ki = 1:n_models
                f = sess_formulas_lag.(model_keys{ki});
                try
                    mdl_tr  = safe_fitlm(T_tr, f);
                    pred_te = predict(mdl_tr, T_te);
                    resid   = T_te.dff - pred_te;
                    mse_val = mean(resid.^2, 'omitnan');
                    mse_real_iter(it,ki)   = mse_val;
                    % CV R2 = 1 - MSE / var(y_test)
                    % Clamped at -1 to avoid extreme negative values from
                    % very bad predictions on small test folds
                    if var_te > 0
                        cv_r2_real_iter(it,ki) = max(1 - mse_val/var_te, -1);
                    end
                catch
                end
                try
                    mdl_sh  = safe_fitlm(T_tr_sh, f);
                    pred_sh = predict(mdl_sh, T_te_sh);
                    resid_sh  = T_te_sh.dff - pred_sh;
                    mse_sh_val = mean(resid_sh.^2, 'omitnan');
                    mse_shuf_iter(it,ki)   = mse_sh_val;
                    if var_te_sh > 0
                        cv_r2_shuf_iter(it,ki) = max(1 - mse_sh_val/var_te_sh, -1);
                    end
                catch
                end
            end
        end

        mse_real(m,:,sess_idx+1)   = nanmean(mse_real_iter,  1);
        mse_shuf(m,:,sess_idx+1)   = nanmean(mse_shuf_iter,  1);
        cv_r2_real(m,:,sess_idx+1) = nanmean(cv_r2_real_iter, 1);
        cv_r2_shuf(m,:,sess_idx+1) = nanmean(cv_r2_shuf_iter, 1);

        fprintf('  %s%s: %d rows, R2(FULL)=%.3f, CV-R2=%.3f, MSE=%.4f\n', ...
            mid, sess_str, height(T), ...
            r2_real(m,1,sess_idx+1), cv_r2_real(m,1,sess_idx+1), mse_real(m,1,sess_idx+1));
    end
end

% ------------------------------------------------------------------ %
%  NOTES ON THE COLUMN 18 / 19 PREDICTORS                             %
% ------------------------------------------------------------------ %
if n_elapsed_dropped > 0
    fprintf(['\nNOTE: elapsed_discovery was dropped from %d session fit(s) because within ', ...
             'a session it equals time minus a constant (r > 0.999 after z-scoring), so it ', ...
             'cannot add anything to the time predictor; its contribution is NaN there and ', ...
             'NO_ELAPSED equals FULL. Set opts.excl_abs_time = true if you want it to stand ', ...
             'in for time, or opts.excl_elapsed = true to silence this.\n'], n_elapsed_dropped);
end
if ~isempty(r_travel_time)
    fprintf(['NOTE: distance_traveled vs time correlation across session fits: median r = %.3f ', ...
             '(range %.3f to %.3f). The closer to 1, the less unique variance it can explain ', ...
             'beyond time, so read its contribution cautiously.\n'], ...
             median(r_travel_time, 'omitnan'), min(r_travel_time), max(r_travel_time));
end

% ------------------------------------------------------------------ %
%  REDUNDANCY AMONG PREDICTORS: correlated pairs                      %
% ------------------------------------------------------------------ %
% Two correlated predictors can each show a small UNIQUE contribution
% even when each is informative on its own, because removing one leaves
% the other to cover for it.  This lists which pairs are that redundant.
fprintf('\n=== Correlated predictors (median r across animals; pairs with |r| >= 0.5) ===\n');
for sess_idx = 0:1
    pm = pair_r{sess_idx+1};
    ks = keys(pm);
    med = cellfun(@(kk) median(pm(kk), 'omitnan'), ks);
    n_s = cellfun(@(kk) numel(pm(kk)), ks);
    [~, o_med] = sort(abs(med), 'descend');
    fprintf('--- sess%d ---\n', sess_idx);
    n_shown = 0;
    for oi = o_med(:)'
        if abs(med(oi)) < 0.5, break; end
        fprintf('  r = %+.2f   %-45s (n=%d)\n', med(oi), ks{oi}, n_s(oi));
        n_shown = n_shown + 1;
    end
    if n_shown == 0, fprintf('  none\n'); end
    k_sp = 'path_last_visit x spatial_distance';
    if isKey(pm, k_sp)
        v_sp = pm(k_sp);
        fprintf('  [spatial_distance vs path_last_visit: median r = %+.2f, range %+.2f to %+.2f]\n', ...
            median(v_sp, 'omitnan'), min(v_sp), max(v_sp));
    end
end

% ------------------------------------------------------------------ %
%  GROUP STATISTICS                                                   %
% ------------------------------------------------------------------ %
fprintf('\n=== Group Statistics ===\n');
group = struct();

for sess_idx = 0:1
    sl = sprintf('sess%d', sess_idx);
    fprintf('\n--- %s ---\n', sl);
    fprintf('%-16s  R2_real(mean+SEM)  R2_shuf(mean+SEM)  p(real vs shuf)  MSE_real  MSE_shuf\n', 'Model');
    fprintf('%s\n', repmat('-',1,95));

    for ki = 1:n_models
        r2r = r2_real(:,ki,sess_idx+1);
        r2s = r2_shuf(:,ki,sess_idx+1);
        mr  = mse_real(:,ki,sess_idx+1);
        ms  = mse_shuf(:,ki,sess_idx+1);
        cvr = cv_r2_real(:,ki,sess_idx+1);
        cvs = cv_r2_shuf(:,ki,sess_idx+1);

        valid = ~isnan(r2r) & ~isnan(r2s);
        if sum(valid) < 3, continue; end

        r2r_v = r2r(valid); r2s_v = r2s(valid);
        mr_v  = mr(valid);  ms_v  = ms(valid);
        cvr_v = cvr(valid); cvs_v = cvs(valid);

        [~,p_r2]   = ttest(r2r_v, r2s_v);
        [~,p_mse]  = ttest(mr_v,  ms_v);
        [~,p_cvr2] = ttest(cvr_v, cvs_v);

        group.(sl).(model_keys{ki}).r2_real_mean    = mean(r2r_v);
        group.(sl).(model_keys{ki}).r2_real_sem     = std(r2r_v)/sqrt(length(r2r_v));
        group.(sl).(model_keys{ki}).r2_shuf_mean    = mean(r2s_v);
        group.(sl).(model_keys{ki}).r2_shuf_sem     = std(r2s_v)/sqrt(length(r2s_v));
        group.(sl).(model_keys{ki}).mse_real_mean   = mean(mr_v);
        group.(sl).(model_keys{ki}).mse_shuf_mean   = mean(ms_v);
        group.(sl).(model_keys{ki}).cv_r2_real_mean = mean(cvr_v);
        group.(sl).(model_keys{ki}).cv_r2_real_sem  = std(cvr_v)/sqrt(length(cvr_v));
        group.(sl).(model_keys{ki}).cv_r2_shuf_mean = mean(cvs_v);
        group.(sl).(model_keys{ki}).cv_r2_shuf_sem  = std(cvs_v)/sqrt(length(cvs_v));
        group.(sl).(model_keys{ki}).p_r2    = p_r2;
        group.(sl).(model_keys{ki}).p_mse   = p_mse;
        group.(sl).(model_keys{ki}).p_cvr2  = p_cvr2;
        group.(sl).(model_keys{ki}).n       = sum(valid);

        fprintf('%-22s  R2=%.3f+/-%.3f  CV-R2=%.3f+/-%.3f  MSE=%.4f\n', ...
            model_labs{ki}, mean(r2r_v), std(r2r_v)/sqrt(length(r2r_v)), ...
            mean(cvr_v), std(cvr_v)/sqrt(length(cvr_v)), mean(mr_v));
    end
end

% ------------------------------------------------------------------ %
%  COMMONALITY: spatial distance vs path from last visit              %
% ------------------------------------------------------------------ %
% If two predictors carry overlapping information, the drop in R2 when
% ONE is removed ("unique") under-states both, because the other covers
% for it.  Using the model without BOTH (NO_SPAT_PATH):
%   unique_spat = FULL - NO_SPAT
%   unique_path = FULL - NO_PATH_LAST_V
%   joint       = FULL - NO_SPAT_PATH      what the pair adds on top of
%                                          everything else in the model
%   shared      = joint - unique_spat - unique_path
%                 explained by either one, credited to neither by the
%                 unique contributions.  shared >> unique means the two
%                 are largely interchangeable.  (shared < 0 would mean
%                 suppression: each is more useful with the other.)
% Computed for in-sample adjusted R2 and for CV R2 (the latter is the
% one to trust: it asks what the pair adds on held-out data).
ik_model = @(kk) find(strcmp(model_keys, kk));
commonality  = struct();
comm_metric  = {'insample_adjR2', 'cv_R2'};
comm_data    = {r2_real, cv_r2_real};
for sess_idx = 0:1
    sl = sprintf('sess%d', sess_idx);
    fprintf('\n--- Commonality: spatial distance vs path from last visit, %s ---\n', sl);
    fprintf('%-15s %-14s %-14s %-14s %-14s %s\n', 'metric', 'unique spat', ...
        'unique path', 'shared', 'joint (pair)', 'n');
    for mt = 1:2
        M  = comm_data{mt};
        F  = M(:, ik_model('FULL'),           sess_idx+1);
        NS = M(:, ik_model('NO_SPAT'),        sess_idx+1);
        NP = M(:, ik_model('NO_PATH_LAST_V'), sess_idx+1);
        NB = M(:, ik_model('NO_SPAT_PATH'),   sess_idx+1);
        uS = F - NS;  uP = F - NP;  jt = F - NB;  sh = jt - uS - uP;
        ok = ~isnan(uS) & ~isnan(uP) & ~isnan(jt) & has_path(:, sess_idx+1);
        if sum(ok) < 1
            fprintf('%-15s not available (no usable path_last_visit)\n', comm_metric{mt});
            continue;
        end
        ms = @(v) sprintf('%+.3f+/-%.3f', mean(v(ok)), std(v(ok))/sqrt(sum(ok)));
        fprintf('%-15s %-14s %-14s %-14s %-14s %d\n', comm_metric{mt}, ...
            ms(uS), ms(uP), ms(sh), ms(jt), sum(ok));
        commonality.(sl).(comm_metric{mt}) = struct('unique_spat', uS, ...
            'unique_path', uP, 'shared', sh, 'joint', jt, 'valid', ok);
    end
end

% ------------------------------------------------------------------ %
%  STORE RESULTS                                                       %
% ------------------------------------------------------------------ %
glm_results.r2_real        = r2_real;
glm_results.r2_contrib     = r2_contrib;
glm_results.pred_names     = pred_names;
glm_results.r2_shuf        = r2_shuf;
glm_results.mse_real       = mse_real;
glm_results.mse_shuf       = mse_shuf;
glm_results.cv_r2_real     = cv_r2_real;
glm_results.cv_r2_shuf     = cv_r2_shuf;
glm_results.unique_mice    = unique_mice;
glm_results.model_keys     = model_keys;
glm_results.model_labs     = model_labs;
glm_results.group          = group;
glm_results.trace_cands_s1 = trace_cands_s1;
glm_results.trace_plot_labs = trace_plot_labs;
glm_results.shap_r2        = shap_r2;       % animals x predictors x session
glm_results.shap_total     = shap_total;    % in-sample R2 of FULL = sum of shares
glm_results.commonality    = commonality;   % spat vs path: unique/shared/joint
glm_results.pair_r_keys    = {keys(pair_r{1}), keys(pair_r{2})};
glm_results.pair_r_vals    = {values(pair_r{1}), values(pair_r{2})};

% ------------------------------------------------------------------ %
%  CONDITION LABEL FOR PLOTS                                          %
% ------------------------------------------------------------------ %
cond_parts = {opts.state, opts.source};
if opts.excl_abs_time,    cond_parts{end+1} = 'no-time'; end
if opts.excl_grooming,    cond_parts{end+1} = 'no-groom'; end
if opts.excl_food_events, cond_parts{end+1} = 'no-food-ev'; end
% Tag runs that include the column 18/19 predictors so their output
% files are not confused with (or overwrite) earlier runs without them.
if ~opts.excl_elapsed,    cond_parts{end+1} = 'elapsed'; end
if ~opts.excl_travel,     cond_parts{end+1} = 'travel'; end
if opts.headdirection,    cond_parts{end+1} = 'headdir'; end
if ~isempty(opts.min_food_dist)
    cond_parts{end+1} = sprintf('dist-ge%g', opts.min_food_dist);
end
if ~isempty(opts.lag_sec) && opts.lag_sec > 0 && ~opts.excl_lag
    cond_parts{end+1} = sprintf('lag%.0fs', opts.lag_sec);
end
cond_parts{end+1} = strrep(opts.tsf_type, '_', '-');
if ~isequal(opts.baseline_correct, false)
    cond_parts{end+1} = sprintf('bl-%s', char(string(opts.baseline_correct)));
end
if strcmp(opts.cv_type, 'blocked'),    cond_parts{end+1} = 'CV-time-blocks'; end
if strcmp(opts.shuf_type, 'circular'), cond_parts{end+1} = 'shuf-circ';      end
condition_label = strjoin(cond_parts, ' | ');

plot_model_comparison(glm_results, condition_label, model_labs, model_keys, ...
    pred_names, r2_contrib, unique_mice, ...
    trace_cands_s1, trace_plot_labs, condition_label);

fprintf('\nDone.\n');
end


% ====================================================================== %
%  INTERACTIVE PDF SAVING                                                %
% ====================================================================== %
function dff_out = apply_baseline_correction(dff, time_vec, opts)
dff_out = dff;
n       = length(dff);
method  = opts.baseline_correct;

if strcmp(method, 'poly')
    t_norm  = (time_vec - mean(time_vec)) / std(time_vec);
    p       = polyfit(t_norm, dff, 2);
    dff_out = dff - polyval(p, t_norm);
    return;
end

dt       = median(diff(time_vec), 'omitnan');
half_win = max(round((opts.baseline_win_sec / 2) / dt), 1);
baseline = NaN(n, 1);

for t = 1:n
    i1  = max(1, t - half_win);
    i2  = min(n, t + half_win);
    seg = dff(i1:i2);
    seg = seg(~isnan(seg));
    if isempty(seg), continue; end
    if strcmp(method, 'sliding_mean')
        baseline(t) = mean(seg);
    elseif strcmp(method, 'sliding_pct')
        baseline(t) = prctile(seg, opts.baseline_pct);
    end
end

valid = ~isnan(baseline);
dff_out(valid) = dff(valid) - baseline(valid);
end


% ====================================================================== %
%  SAFE FITLM                                                            %
% ====================================================================== %
function mdl = safe_fitlm(T, formula)
parts   = strtrim(strsplit(formula, '~'));
rhs     = strtrim(strsplit(parts{2}, '+'));
outcome = strtrim(parts{1});
keep_preds = {};
for k = 1:length(rhs)
    pn = strtrim(rhs{k});
    if ismember(pn, T.Properties.VariableNames)
        col = T.(pn);
        if isnumeric(col) && std(col, 'omitnan') > 0
            keep_preds{end+1} = pn;
        end
    end
end
if isempty(keep_preds)
    error('No valid predictors remain after removing zero-variance columns.');
end
mdl = fitlm(T, [outcome, ' ~ ', strjoin(keep_preds, ' + ')]);
end


% ====================================================================== %
%  SHAPLEY (DOMINANCE) DECOMPOSITION OF R2                               %
% ====================================================================== %
function [sv, r2_full] = shapley_r2(T, pred_list)
% Shapley value of each predictor's contribution to the in-sample R2 of
% the model that contains all of pred_list: the average, over all subsets
% S of the OTHER predictors (weighted by the number of orderings that put
% exactly S before it), of the R2 gained by adding the predictor to S.
% The shares sum exactly to the R2 of the full model, so variance shared
% by correlated predictors is split between them rather than credited to
% neither.  R2 of every subset comes from the correlation matrix, so no
% model is refit.
%
% NOTE: splitting shared variance equally is a convention, not a finding
% -- it says how much of the R2 each predictor can claim on average, not
% which one is the true driver.
p = numel(pred_list);
X = zeros(height(T), p);
for j = 1:p, X(:, j) = T.(pred_list{j}); end
y = T.dff;

ok = all(isfinite([X, y]), 2);               % complete rows only
X = X(ok, :);  y = y(ok);
n = size(X, 1);
X = (X - mean(X, 1)) ./ std(X, 0, 1);
y = (y - mean(y)) / std(y);
C = (X' * X) / (n - 1);                      % predictor correlations
c = (X' * y) / (n - 1);                      % predictor-outcome correlations

n_sub = 2^p;
R2 = zeros(n_sub, 1);                        % R2(mask+1); mask 0 = empty set
pc = zeros(n_sub, 1);                        % number of predictors in each subset
for mask = 1:n_sub-1
    idx = find(bitget(mask, 1:p));
    pc(mask+1) = numel(idx);
    R2(mask+1) = c(idx)' * (pinv(C(idx, idx)) * c(idx));
end

sv = zeros(1, p);
fp = factorial(p);
for j = 1:p
    bj = 2^(j-1);
    for mask = 0:n_sub-1
        if bitand(mask, bj) == 0             % subsets that do not contain j
            k = pc(mask+1);
            w = factorial(k) * factorial(p-k-1) / fp;
            sv(j) = sv(j) + w * (R2(mask+bj+1) - R2(mask+1));
        end
    end
end
r2_full = R2(n_sub);
end


% ====================================================================== %
%  PLOTS                                                                  %
% ====================================================================== %
function cmap = rdbu_colormap(n)
% Returns an n-row RGB colormap going blue -> white -> red.
% Mimics the matplotlib/ColorBrewer RdBu diverging colormap.
% Blue end:  [0.02 0.44 0.69]
% White mid: [1.00 1.00 1.00]
% Red end:   [0.84 0.19 0.15]
if nargin < 1, n = 256; end
half = ceil(n/2);
% Blue to white
blue_end = [0.02 0.44 0.69];
white    = [1.00 1.00 1.00];
red_end  = [0.84 0.19 0.15];
t1 = linspace(0, 1, half)';
t2 = linspace(0, 1, n - half)';
top    = blue_end + t1 .* (white - blue_end);   % blue -> white
bottom = white    + t2 .* (red_end - white);     % white -> red
cmap   = [top; bottom];
end


% ====================================================================== %
%  PLOTS                                                                  %
% ====================================================================== %
function plot_model_comparison(res, condition_label, model_labs, model_keys, ...
    pred_names, r2_contrib, unique_mice, ...
    trace_cands_s1, trace_plot_labs, cond_label)

n_models    = length(model_keys);
n_mice      = length(unique_mice);
n_preds     = length(pred_names);
shuf_alpha  = 0.35;
bw          = 0.28;
sess_labels = {'Sess 0 (pre-discovery)','Sess 1 (post-discovery)'};
short_labels = cellfun(@(s) s(max(1,end-6):end), unique_mice, 'UniformOutput', false);

% Color scheme: group models visually
% 1=FULL(blue) 2-4=reds(spat/temp) 5-7=greens(food) 8=purple(groom) 9=orange(speed)
% 10-12=isolation dist, 13=food, 14=groom, 15=speed,
% 16=path-only, 17=accel-only, 18=elapsed-only, 19=travel-only,
% 20=path-last-visit knockout, 21=spat+path knockout,
% 22=acceleration knockout, 23=elapsed-from-discovery knockout,
% 24=distance-traveled knockout, 25=head-direction knockout and
% 26=head-direction-only (rows 25-26 are used only if opts.headdirection)
bar_colors = [0.15 0.35 0.75;   % FULL             - dark blue
    0.85 0.20 0.15;   % NO_SPAT          - red
    0.95 0.45 0.40;   % NO_TEMP          - light red
    0.65 0.10 0.10;   % NO_DIST          - dark red
    0.15 0.65 0.30;   % NO_FOOD_VIS      - green
    0.45 0.80 0.50;   % NO_EATING        - light green
    0.05 0.45 0.15;   % NO_FOOD_EAT      - dark green
    0.55 0.15 0.65;   % NO_GROOM         - purple
    0.95 0.55 0.10;   % NO_SPEED         - orange
    0.90 0.25 0.20;   % SPAT_ONLY        - red (matches spat dist color)
    0.20 0.65 0.85;   % TEMP_ONLY        - teal (matches temp dist color)
    0.10 0.45 0.60;   % DIST_ONLY        - dark teal
    0.30 0.75 0.40;   % FOOD_ONLY        - lime
    0.70 0.35 0.80;   % GROOM_ONLY       - violet
    0.85 0.70 0.05;   % SPEED_ONLY       - gold
    0.15 0.80 0.72;   % PATH_ONLY        - light teal-green (matches path color)
    0.78 0.58 0.38;   % ACCEL_ONLY       - light brown (matches accel color)
    0.35 0.60 0.60;   % ELAPSED_ONLY     - light teal-gray
    0.88 0.65 0.72;   % TRAVEL_ONLY      - light mauve
    0.10 0.70 0.65;   % NO_PATH_LAST_V   - teal-green
    0.50 0.35 0.45;   % NO_SPAT_PATH     - slate
    0.60 0.40 0.20;   % NO_ACCEL         - brown
    0.20 0.45 0.45;   % NO_ELAPSED       - dark teal-gray
    0.80 0.50 0.60;   % NO_TRAVEL        - mauve
    0.85 0.25 0.60;   % NO_HEADDIR       - magenta   (only with opts.headdirection)
    0.93 0.55 0.78];  % HEADDIR_ONLY     - pink      (only with opts.headdirection)

% ---- Figure 1: R2 real vs shuffled — all models -------------------- %
% Use two subplots (sess0, sess1) with horizontal bars for readability
% since the models don't fit well as vertical bars.
r2_all = [];
for sp_pre = 1:2
    r2r_pre = res.r2_real(:,:,sp_pre);
    r2s_pre = res.r2_shuf(:,:,sp_pre);
    r2r_m = nanmean(r2r_pre,1) + nanstd(r2r_pre,0,1)./sqrt(sum(~isnan(r2r_pre),1));
    r2s_m = nanmean(r2s_pre,1) + nanstd(r2s_pre,0,1)./sqrt(sum(~isnan(r2s_pre),1));
    r2_all = [r2_all, r2r_m, r2s_m];
end
shared_xlim = [0, max(r2_all,[],'omitnan') * 1.30];

figure('Name','Model comparison R2: real vs shuffled', ...
    'Position',[50 50 1300 720]);
for sp = 1:2
    subplot(1,2,sp); hold on;
    r2r = res.r2_real(:,:,sp);
    r2s = res.r2_shuf(:,:,sp);
    r2r_mean = nanmean(r2r,1); r2r_sem = nanstd(r2r,0,1)./sqrt(sum(~isnan(r2r),1));
    r2s_mean = nanmean(r2s,1); r2s_sem = nanstd(r2s,0,1)./sqrt(sum(~isnan(r2s),1));

    % Plot from bottom (model 1=FULL at top, last model at bottom)
    % y positions: n_models down to 1
    for ki = 1:n_models
        ypos  = n_models - ki + 1;
        col_r = bar_colors(ki,:);
        col_s = col_r + (1-col_r)*shuf_alpha;

        % Real bar
        barh(ypos+bw/2, r2r_mean(ki), bw, 'FaceColor', col_r, 'EdgeColor','none');
        % Shuffled bar
        barh(ypos-bw/2, r2s_mean(ki), bw, 'FaceColor', col_s, 'EdgeColor', col_r, ...
            'LineWidth',0.8,'LineStyle','--');
        % Error bars
        errorbar(r2r_mean(ki), ypos+bw/2, r2r_sem(ki), 'horizontal', 'k.', ...
            'LineWidth',1.2,'CapSize',3);
        errorbar(r2s_mean(ki), ypos-bw/2, r2s_sem(ki), 'horizontal', ...
            'Color',col_r*0.6,'LineStyle','none','LineWidth',1,'CapSize',3);
        % Individual animal dots
        valid = ~isnan(r2r(:,ki));
        jitter = (rand(sum(valid),1)-0.5)*0.12;
        scatter(r2r(valid,ki), ypos+bw/2+jitter, 14, col_r*0.7, 'filled', ...
            'MarkerFaceAlpha',0.5);

        % Significance star
        sl_key = sprintf('sess%d',sp-1);
        if isfield(res.group,sl_key) && isfield(res.group.(sl_key),model_keys{ki})
            p_val = res.group.(sl_key).(model_keys{ki}).p_r2;
            if     p_val < 0.001, sig = '***';
            elseif p_val < 0.01,  sig = '**';
            elseif p_val < 0.05,  sig = '*';
            else,                 sig = ''; end
            if ~isempty(sig)
                xstar = max(r2r_mean(ki)+r2r_sem(ki), r2s_mean(ki)+r2s_sem(ki)) + ...
                    0.015*shared_xlim(2);
                text(xstar, ypos+bw/2, sig, 'VerticalAlignment','middle', ...
                    'FontSize',9,'FontWeight','bold','Color',col_r*0.8);
            end
        end
    end

    % Group separators
    yline(n_models-0.5,'k:','LineWidth',0.5);   % after FULL
    yline(n_models-3.5,'k:','LineWidth',0.5);   % after NO_DIST group
    yline(n_models-6.5,'k:','LineWidth',0.5);   % after food group
    yline(n_models-8.5,'k:','LineWidth',0.5);   % after groom/speed

    set(gca,'YTick',1:n_models,'YTickLabel',flip(model_labs),'FontSize',8.5);
    xlabel('Adjusted R^2  (mean \pm SEM)');
    title(sess_labels{sp},'FontSize',11);
    xlim(shared_xlim); grid on; box off;
end
subplot(1,2,1);
h1 = barh(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
h2 = barh(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
legend([h1 h2],{'Real dF/F','Shuffled (circular)'},'Location','southeast','FontSize',8);
sgtitle(sprintf('Model comparison: Adj R^2 real vs shuffled  [%s]\n(* p<0.05, ** p<0.01, *** p<0.001)', ...
    condition_label),'FontSize',11);

% ---- Figure 2: CV MSE: real vs shuffled --------------------------- %
mse_all = [res.mse_real(:); res.mse_shuf(:)];
xlim_mse = [0, max(mse_all(isfinite(mse_all)),[],'omitnan') * 1.15];

figure('Name','Model comparison MSE: real vs shuffled', ...
    'Position',[50 750 1300 720]);
for sp = 1:2
    subplot(1,2,sp); hold on;
    mr = res.mse_real(:,:,sp); ms = res.mse_shuf(:,:,sp);
    mr_mean = nanmean(mr,1); mr_sem = nanstd(mr,0,1)./sqrt(sum(~isnan(mr),1));
    ms_mean = nanmean(ms,1); ms_sem = nanstd(ms,0,1)./sqrt(sum(~isnan(ms),1));

    for ki = 1:n_models
        ypos  = n_models - ki + 1;
        col_r = bar_colors(ki,:);
        col_s = col_r + (1-col_r)*shuf_alpha;
        barh(ypos+bw/2, mr_mean(ki), bw, 'FaceColor',col_r,'EdgeColor','none');
        barh(ypos-bw/2, ms_mean(ki), bw, 'FaceColor',col_s,'EdgeColor',col_r, ...
            'LineWidth',0.8,'LineStyle','--');
        errorbar(mr_mean(ki),ypos+bw/2,mr_sem(ki),'horizontal','k.','LineWidth',1.2,'CapSize',3);
        errorbar(ms_mean(ki),ypos-bw/2,ms_sem(ki),'horizontal', ...
            'Color',col_r*0.6,'LineStyle','none','LineWidth',1,'CapSize',3);
        valid = ~isnan(mr(:,ki));
        jitter = (rand(sum(valid),1)-0.5)*0.12;
        scatter(mr(valid,ki), ypos+bw/2+jitter, 14, col_r*0.7,'filled','MarkerFaceAlpha',0.5);
    end
    yline(n_models-0.5,'k:','LineWidth',0.5);
    yline(n_models-3.5,'k:','LineWidth',0.5);
    yline(n_models-6.5,'k:','LineWidth',0.5);
    yline(n_models-8.5,'k:','LineWidth',0.5);
    set(gca,'YTick',1:n_models,'YTickLabel',flip(model_labs),'FontSize',8.5);
    xlabel('CV MSE  (mean \pm SEM, lower = better)');
    title(sess_labels{sp},'FontSize',11);
    xlim(xlim_mse); grid on; box off;
end
subplot(1,2,1);
h1 = barh(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
h2 = barh(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
legend([h1 h2],{'Real dF/F','Shuffled (circular)'},'Location','southeast','FontSize',8);
sgtitle(sprintf('Model comparison: CV MSE  [%s]',condition_label),'FontSize',11);

% ---- Figure 2b: CV R2 (held-out test set) -------------------------- %
% CV R2 = 1 - MSE/var(y_test), averaged across folds.
% This is the out-of-sample R2: how much variance is explained in
% held-out data.  Directly comparable to in-sample R2 (Figure 1).
cvr2_all = [res.cv_r2_real(:); res.cv_r2_shuf(:)];
xlim_cvr2 = [min(cvr2_all(isfinite(cvr2_all)),[],'omitnan') * 1.1, ...
    max(cvr2_all(isfinite(cvr2_all)),[],'omitnan') * 1.30];
xlim_cvr2(1) = min(xlim_cvr2(1), 0);   % always show 0

figure('Name','Model comparison CV-R2: real vs shuffled', ...
    'Position',[50 50 1300 720]);
for sp = 1:2
    subplot(1,2,sp); hold on;
    cvr = res.cv_r2_real(:,:,sp); cvs = res.cv_r2_shuf(:,:,sp);
    cvr_mean = nanmean(cvr,1); cvr_sem = nanstd(cvr,0,1)./sqrt(sum(~isnan(cvr),1));
    cvs_mean = nanmean(cvs,1); cvs_sem = nanstd(cvs,0,1)./sqrt(sum(~isnan(cvs),1));
    for ki = 1:n_models
        ypos = n_models-ki+1; col_r = bar_colors(ki,:); col_s = col_r+(1-col_r)*shuf_alpha;
        barh(ypos+bw/2, cvr_mean(ki), bw, 'FaceColor',col_r,'EdgeColor','none');
        barh(ypos-bw/2, cvs_mean(ki), bw, 'FaceColor',col_s,'EdgeColor',col_r,'LineWidth',0.8,'LineStyle','--');
        errorbar(cvr_mean(ki),ypos+bw/2,cvr_sem(ki),'horizontal','k.','LineWidth',1.2,'CapSize',3);
        errorbar(cvs_mean(ki),ypos-bw/2,cvs_sem(ki),'horizontal','Color',col_r*0.6,'LineStyle','none','LineWidth',1,'CapSize',3);
        valid = ~isnan(cvr(:,ki));
        jitter = (rand(sum(valid),1)-0.5)*0.12;
        scatter(cvr(valid,ki),ypos+bw/2+jitter,14,col_r*0.7,'filled','MarkerFaceAlpha',0.5);
        sl_key = sprintf('sess%d',sp-1);
        if isfield(res.group,sl_key) && isfield(res.group.(sl_key),model_keys{ki})
            if isfield(res.group.(sl_key).(model_keys{ki}),'p_cvr2')
                p_val = res.group.(sl_key).(model_keys{ki}).p_cvr2;
                if p_val<0.001,sig='***'; elseif p_val<0.01,sig='**'; elseif p_val<0.05,sig='*'; else,sig=''; end
                if ~isempty(sig)
                    text(max(cvr_mean(ki)+cvr_sem(ki),cvs_mean(ki)+cvs_sem(ki))+0.015*range(xlim_cvr2), ...
                        ypos+bw/2, sig,'VerticalAlignment','middle','FontSize',9,'FontWeight','bold','Color',col_r*0.8);
                end
            end
        end
    end
    xline(0,'k--','LineWidth',0.8);
    yline(n_models-0.5,'k:','LineWidth',0.5); yline(n_models-3.5,'k:','LineWidth',0.5);
    yline(n_models-6.5,'k:','LineWidth',0.5); yline(n_models-8.5,'k:','LineWidth',0.5);
    set(gca,'YTick',1:n_models,'YTickLabel',flip(model_labs),'FontSize',8.5);
    xlabel('CV R^2  (held-out test set, mean \pm SEM)');
    title(sess_labels{sp},'FontSize',11);
    xlim(xlim_cvr2); grid on; box off;
end
subplot(1,2,1);
legend([barh(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none'), ...
    barh(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--')], ...
    {'Real dF/F','Shuffled (circular)'},'Location','southeast','FontSize',8);
sgtitle(sprintf('Model comparison: CV R^2  [%s]\n(* p<0.05, ** p<0.01, *** p<0.001)',condition_label),'FontSize',11);

% ---- Save combined Excel: in-sample R2, CV R2, MSE ---------------- %
try
    % Three sheets: one per metric
    xl_fname_all = sprintf('R2_CV_MSE_model_comparison_%s.xlsx', ...
        strrep(strrep(condition_label,' | ','_'),' ','_'));

    metrics = {'InSample_R2', 'CV_R2', 'CV_MSE'};
    data_real = {res.r2_real, res.cv_r2_real, res.mse_real};
    data_shuf = {res.r2_shuf, res.cv_r2_shuf, res.mse_shuf};

    for met_idx = 1:3
        xl_headers = {'Animal'};
        for sp_xl=1:2, for ki_xl=1:n_models
                xl_headers{end+1} = sprintf('S%d_%s_real', sp_xl-1, model_keys{ki_xl});
                xl_headers{end+1} = sprintf('S%d_%s_shuf', sp_xl-1, model_keys{ki_xl});
        end; end
    xl_data = {};
    d_r = data_real{met_idx}; d_s = data_shuf{met_idx};
    for mi=1:n_mice
        row = {unique_mice{mi}};
        for sp_xl=1:2, for ki_xl=1:n_models
                row{end+1}=d_r(mi,ki_xl,sp_xl); row{end+1}=d_s(mi,ki_xl,sp_xl);
        end; end
    xl_data(end+1,:) = row;
    end
    mean_row={'Mean'}; sem_row={'SEM'};
    for sp_xl=1:2, for ki_xl=1:n_models
            vr=d_r(:,ki_xl,sp_xl); vs=d_s(:,ki_xl,sp_xl);
            mean_row{end+1}=nanmean(vr); mean_row{end+1}=nanmean(vs);
            sem_row{end+1}=nanstd(vr)/sqrt(sum(~isnan(vr)));
            sem_row{end+1}=nanstd(vs)/sqrt(sum(~isnan(vs)));
    end; end
xl_data(end+1,:)=mean_row; xl_data(end+1,:)=sem_row;
writetable(cell2table(xl_data,'VariableNames',xl_headers), ...
    xl_fname_all, 'Sheet', metrics{met_idx});
    end
    fprintf('\nCombined Excel saved: %s  (sheets: InSample_R2, CV_R2, CV_MSE)\n', xl_fname_all);
catch ME_xl
    fprintf('\nCould not save combined Excel: %s\n', ME_xl.message);
end

% ---- Figure 3: Predictor contributions ---------------------------- %
if n_preds > 0 && ~all(isnan(r2_contrib(:)))
    pred_disp = pred_names;
    pred_disp = strrep(pred_disp,'spatial_distance','spatial dist');
    pred_disp = strrep(pred_disp,'temporal_distance','temporal dist');
    pred_disp = strrep(pred_disp,'time_since_food','tsf');
    pred_disp = strrep(pred_disp,'food_visit','food visit');
    pred_disp = strrep(pred_disp,'path_last_visit','path last visit');
    pred_disp = strrep(pred_disp,'acceleration','accel');
    pred_disp = strrep(pred_disp,'_',' ');

    pred_colors = repmat([0.55 0.65 0.85], n_preds, 1);
    for pi = 1:n_preds
        if strcmp(pred_names{pi},'spatial_distance'),  pred_colors(pi,:) = [0.9  0.25 0.2 ]; end
        if strcmp(pred_names{pi},'temporal_distance'), pred_colors(pi,:) = [0.2  0.65 0.85]; end
        if strcmp(pred_names{pi},'time_since_food'),   pred_colors(pi,:) = [0.2  0.7  0.4 ]; end
        if strcmp(pred_names{pi},'time'),              pred_colors(pi,:) = [0.4  0.4  0.4 ]; end
        if strcmp(pred_names{pi},'grooming'),          pred_colors(pi,:) = [0.55 0.15 0.65]; end
        if strcmp(pred_names{pi},'speed'),             pred_colors(pi,:) = [0.95 0.55 0.10]; end
        if strcmp(pred_names{pi},'food_visit'),        pred_colors(pi,:) = [0.15 0.65 0.30]; end
        if strcmp(pred_names{pi},'eating'),            pred_colors(pi,:) = [0.45 0.80 0.50]; end
        if strcmp(pred_names{pi},'path_last_visit'),   pred_colors(pi,:) = [0.10 0.70 0.65]; end
        if strcmp(pred_names{pi},'acceleration'),      pred_colors(pi,:) = [0.60 0.40 0.20]; end
        if strcmp(pred_names{pi},'elapsed_discovery'), pred_colors(pi,:) = [0.20 0.45 0.45]; end
        if strcmp(pred_names{pi},'distance_traveled'), pred_colors(pi,:) = [0.80 0.50 0.60]; end
        if strcmp(pred_names{pi},'cos_head_dir'),      pred_colors(pi,:) = [0.85 0.25 0.60]; end
    end

    figure('Name','Predictor contributions to R2','Position',[50 50 900 480]);
    for sp = 1:2
        subplot(1,2,sp); hold on;
        contrib = r2_contrib(:,:,sp);
        c_mean  = nanmean(contrib,1);
        c_sem   = nanstd(contrib,0,1)./sqrt(sum(~isnan(contrib),1));
        % NaN-contribution predictors (e.g. elapsed_discovery when it
        % duplicates time) go last instead of first
        [~,sort_idx] = sort(c_mean,'descend','MissingPlacement','last');

        for pi = 1:n_preds
            oi  = sort_idx(pi);
            col = pred_colors(oi,:);
            bar(pi, c_mean(oi), 0.6, 'FaceColor',col,'EdgeColor','none');
            errorbar(pi, c_mean(oi), c_sem(oi), 'k.','LineWidth',1.2,'CapSize',4);
            valid = ~isnan(contrib(:,oi));
            jitter = (rand(sum(valid),1)-0.5)*0.18;
            scatter(pi+jitter, contrib(valid,oi), 18, col*0.7,'filled','MarkerFaceAlpha',0.5);
            v = contrib(valid,oi);
            if length(v) >= 3
                [~,p_v] = ttest(v);
                if     p_v < 0.001, sig = '***';
                elseif p_v < 0.01,  sig = '**';
                elseif p_v < 0.05,  sig = '*';
                else,               sig = ''; end
                if ~isempty(sig)
                    text(pi, c_mean(oi)+c_sem(oi)+0.003, sig, ...
                        'HorizontalAlignment','center','FontSize',11, ...
                        'FontWeight','bold','Color',col*0.75);
                end
            end
        end
        yline(0,'k--','LineWidth',0.8);
        set(gca,'XTick',1:n_preds,'XTickLabel',pred_disp(sort_idx), ...
            'XTickLabelRotation',35,'FontSize',9);
        ylabel('Delta R^2  (unique contribution to Full model)');
        title(sess_labels{sp},'FontSize',11);
        grid on; box off;
        ymax = max([c_mean+c_sem, 0.01],[],'omitnan');
        ymin = min([c_mean-c_sem, 0  ],[],'omitnan');
        ylim([ymin*1.2, ymax*1.35]);
    end
    sgtitle(sprintf('Unique predictor contributions  [%s]  (* p<0.05 vs 0)', ...
        condition_label),'FontSize',11);

    % Save Excel for contributions
    try
        xl_headers_c = {'Animal'};
        for sp_xl = 1:2
            for pi_xl = 1:n_preds
                xl_headers_c{end+1} = sprintf('S%d_%s_dR2',sp_xl-1,pred_disp{pi_xl});
            end
        end
        xl_data_c = {};
        for mi = 1:n_mice
            row = {unique_mice{mi}};
            for sp_xl = 1:2
                for pi_xl = 1:n_preds
                    row{end+1} = r2_contrib(mi,pi_xl,sp_xl);
                end
            end
            xl_data_c(end+1,:) = row;
        end
        mean_row_c = {'Mean'}; sem_row_c = {'SEM'}; p_row_c = {'p (vs 0)'};
        for sp_xl = 1:2
            for pi_xl = 1:n_preds
                v = r2_contrib(:,pi_xl,sp_xl); v = v(~isnan(v));
                mean_row_c{end+1} = nanmean(v);
                sem_row_c{end+1}  = std(v)/sqrt(length(v));
                if length(v) >= 3, [~,pv] = ttest(v); p_row_c{end+1} = pv;
                else, p_row_c{end+1} = NaN; end
            end
        end
        xl_data_c(end+1,:) = mean_row_c;
        xl_data_c(end+1,:) = sem_row_c;
        xl_data_c(end+1,:) = p_row_c;
        xl_fname_c = sprintf('Contributions_model_comparison_%s.xlsx', ...
            strrep(strrep(condition_label,' | ','_'),' ','_'));
        writetable(cell2table(xl_data_c,'VariableNames',xl_headers_c), xl_fname_c);
        fprintf('Contributions Excel saved: %s\n', xl_fname_c);
    catch ME_xl
        fprintf('Could not save contributions Excel: %s\n', ME_xl.message);
    end

    % ---- Figure 3b: Shapley share vs unique contribution ------------- %
    % The unique dR2 above ignores variance shared with correlated
    % predictors, so two redundant predictors can both look unimportant.
    % The Shapley share divides shared variance between them; shares sum
    % to the FULL model's in-sample R2.  Compare the two bars per predictor:
    % a big gap means most of that predictor's information is shared.
    if isfield(res, 'shap_r2') && ~all(isnan(res.shap_r2(:)))
        figure('Name','Shapley share vs unique contribution','Position',[50 50 1000 480]);
        for sp = 1:2
            subplot(1,2,sp); hold on;
            sh = res.shap_r2(:,:,sp);
            un = r2_contrib(:,:,sp);
            sh_mean = nanmean(sh,1);
            sh_sem  = nanstd(sh,0,1)./sqrt(sum(~isnan(sh),1));
            un_mean = nanmean(un,1);
            un_sem  = nanstd(un,0,1)./sqrt(sum(~isnan(un),1));
            [~, so] = sort(sh_mean, 'descend', 'MissingPlacement', 'last');
            for k = 1:n_preds
                oi  = so(k);
                col = pred_colors(oi,:);
                bar(k-0.18, un_mean(oi), 0.32, 'FaceColor', col + (1-col)*0.55, 'EdgeColor', col);
                bar(k+0.18, sh_mean(oi), 0.32, 'FaceColor', col, 'EdgeColor', 'none');
                errorbar(k-0.18, un_mean(oi), un_sem(oi), 'k.', 'LineWidth',1, 'CapSize',3);
                errorbar(k+0.18, sh_mean(oi), sh_sem(oi), 'k.', 'LineWidth',1, 'CapSize',3);
                v = sh(~isnan(sh(:,oi)), oi);
                scatter(k+0.18+(rand(numel(v),1)-0.5)*0.12, v, 12, col*0.7, 'filled', ...
                    'MarkerFaceAlpha',0.5);
            end
            yline(0,'k--','LineWidth',0.8);
            set(gca,'XTick',1:n_preds,'XTickLabel',pred_disp(so), ...
                'XTickLabelRotation',35,'FontSize',9);
            ylabel('R^2');
            title(sprintf('%s   (Shapley shares sum to R^2 = %.3f)', sess_labels{sp}, ...
                nanmean(res.shap_total(:,sp))), 'FontSize',10);
            grid on; box off;
        end
        subplot(1,2,1);
        h_u = bar(nan, nan, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor',[0.4 0.4 0.4]);
        h_s = bar(nan, nan, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor','none');
        legend([h_u h_s], {'Unique (adj. R^2 lost if removed)', ...
            'Shapley share (shared variance split)'}, 'Location','northeast','FontSize',8);
        sgtitle(sprintf('Unique contribution vs Shapley share  [%s]', condition_label), ...
            'FontSize', 11);

        % Excel: Shapley shares per animal
        try
            xl_h_s = {'Animal'};
            for sp_xl = 1:2
                for pi_xl = 1:n_preds
                    xl_h_s{end+1} = sprintf('S%d_%s_shapley_R2', sp_xl-1, pred_disp{pi_xl}); %#ok<AGROW>
                end
            end
            xl_d_s = {};
            for mi = 1:n_mice
                row = {unique_mice{mi}};
                for sp_xl = 1:2
                    for pi_xl = 1:n_preds
                        row{end+1} = res.shap_r2(mi,pi_xl,sp_xl); %#ok<AGROW>
                    end
                end
                xl_d_s(end+1,:) = row; %#ok<AGROW>
            end
            mean_s = {'Mean'}; sem_s = {'SEM'};
            for sp_xl = 1:2
                for pi_xl = 1:n_preds
                    v = res.shap_r2(:,pi_xl,sp_xl); v = v(~isnan(v));
                    mean_s{end+1} = mean(v);                  %#ok<AGROW>
                    sem_s{end+1}  = std(v)/sqrt(numel(v));    %#ok<AGROW>
                end
            end
            xl_d_s(end+1,:) = mean_s;
            xl_d_s(end+1,:) = sem_s;
            xl_fname_s = sprintf('Shapley_model_comparison_%s.xlsx', ...
                strrep(strrep(condition_label,' | ','_'),' ','_'));
            writetable(cell2table(xl_d_s, 'VariableNames', matlab.lang.makeValidName(xl_h_s)), ...
                xl_fname_s);
            fprintf('Shapley Excel saved: %s\n', xl_fname_s);
        catch ME_xl_s
            fprintf('Could not save Shapley Excel: %s\n', ME_xl_s.message);
        end
    end
end

% ---- Figure 4: Per-animal R2 heatmap (Full model) ----------------- %
figure('Name','Per-animal R2 heatmap (Full model)','Position',[700 50 700 500]);
r2_full_heat = squeeze(res.r2_real(:,1,:));
imagesc(r2_full_heat); colormap(flipud(hot)); cb = colorbar;
cb.Label.String = 'Adjusted R^2';
set(gca,'XTick',1:2,'XTickLabel',{'Sess 0','Sess 1'}, ...
    'YTick',1:n_mice,'YTickLabel',short_labels,'FontSize',9);
xlabel('Session'); ylabel('Animal');
title(sprintf('Full model  Adj R^2 per animal  [%s]',condition_label),'FontSize',11);
for mi = 1:n_mice
    for si = 1:2
        v = r2_full_heat(mi,si);
        if ~isnan(v)
            text(si,mi,sprintf('%.3f',v),'HorizontalAlignment','center', ...
                'FontSize',8,'Color','w','FontWeight','bold');
        end
    end
end

% ---- Figure 5: Predicted vs actual traces — sess1, all animals ---- %
% 6 columns (one per model), rows = actual + predicted overlay.
% Sorted by sess1 R2 descending.  Each animal auto-saved as a PDF.
if isempty(trace_cands_s1)
    fprintf('No sess1 trace candidates available for plotting.\n');
    return;
end

n_trace_models = length(trace_plot_labs);

% Colors matching the model identity
trace_colors = {[0.15 0.35 0.75], ...  % FULL
                [0.85 0.20 0.15], ...  % NO_SPAT
                [0.95 0.45 0.40], ...  % NO_TEMP
                [0.65 0.10 0.10], ...  % NO_DIST
                [0.90 0.25 0.20], ...  % SPAT_ONLY
                [0.20 0.65 0.85]};     % TEMP_ONLY
col_actual = [0.15 0.15 0.15];

% Sort by R2 descending
r2_vals    = cellfun(@(c) c.r2, trace_cands_s1);
[~, sord]  = sort(r2_vals, 'descend', 'MissingPlacement', 'last');

cond_safe = regexprep(strrep(cond_label,' | ','_'), '[^a-zA-Z0-9_]', '_');
cond_safe = regexprep(cond_safe, '_+', '_');

fprintf('\nGenerating and saving sess1 trace PDFs (%d animals)...\n', ...
    length(trace_cands_s1));

for ci = 1:length(trace_cands_s1)
    cand  = trace_cands_s1{sord(ci)};
    mid   = cand.animal;
    t_sec = cand.time - cand.time(1);
    dff_r = cand.dff;

    % Shared y-limits across all panels
    all_y = dff_r;
    for pk = 1:n_trace_models
        if pk <= length(cand.preds) && ~isempty(cand.preds{pk})
            all_y = [all_y; cand.preds{pk}];
        end
    end
    y_pad = 0.15 * max(range(all_y), 0.1);
    ylims = [min(all_y)-y_pad, max(all_y)+y_pad];

    fig = figure('Name', sprintf('Traces sess1: %s', strrep(mid,'_',' ')), ...
        'Position', [30 30 1300 750]);

    for pk = 1:n_trace_models
        subplot(2, 3, pk); hold on;
        col = trace_colors{pk};

        plot(t_sec, dff_r, 'Color', col_actual, 'LineWidth', 0.6, ...
            'DisplayName', 'Actual');

        if pk <= length(cand.preds) && ~isempty(cand.preds{pk})
            pred_k  = cand.preds{pk};
            resid_k = dff_r - pred_k;
            r2_k    = 1 - var(resid_k) / max(var(dff_r), 1e-10);

            % Shaded residual band
            fill([t_sec; flipud(t_sec)], ...
                [pred_k + abs(resid_k); flipud(pred_k - abs(resid_k))], ...
                col, 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
            plot(t_sec, pred_k, 'Color', col, 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%s  r^2=%.3f', trace_plot_labs{pk}, r2_k));
        end

        yline(0, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
        ylim(ylims); xlim([0 max(t_sec)]);
        grid on; box off;
        legend('Location', 'best', 'FontSize', 6);
        xlabel('Time (s)', 'FontSize', 8);
        if pk == 1, ylabel('z-scored dF/F', 'FontSize', 8); end
        title(trace_plot_labs{pk}, 'FontSize', 8, 'FontWeight', 'bold');
    end

    sgtitle(sprintf('%s  —  Sess 1  |  R^2(Full)=%.3f  [%s]', ...
        strrep(mid,'_',' '), cand.r2, cond_label), 'FontSize', 9);

    % Auto-save as PDF (painters renderer = true vector output)
    mid_safe  = regexprep(strrep(mid,'-','_'), '[^a-zA-Z0-9_]', '');
    pdf_fname = sprintf('traces_sess1__%s__%s.pdf', mid_safe, cond_safe);
    if length(pdf_fname) > 200, pdf_fname = [pdf_fname(1:190), '.pdf']; end
    try
        set(fig, 'PaperOrientation', 'landscape', ...
            'PaperUnits',       'normalized', ...
            'PaperPosition',    [0 0 1 1]);
        print(fig, pdf_fname, '-dpdf', '-painters', '-bestfit');
        fprintf('  Saved: %s\n', pdf_fname);
    catch ME_pdf
        fprintf('  ERROR saving PDF for %s: %s\n', mid, ME_pdf.message);
    end
end
fprintf('Trace PDFs complete.\n');

% ---- Interactive: save trace data to Excel for selected animals ---- %
% Lists all animals with their sess1 R2, lets user pick which ones to
% save.  Each animal gets one Excel file with columns:
%   time | dff_actual | pred_<model1> | pred_<model2> | ...
if isempty(trace_cands_s1), return; end

fprintf('\n========================================\n');
fprintf('  TRACE DATA EXPORT (Excel)\n');
fprintf('========================================\n');
fprintf('Available animals (sorted by sess1 R2):\n');
r2_vals_disp = cellfun(@(c) c.r2, trace_cands_s1);
[~, sord_disp] = sort(r2_vals_disp, 'descend', 'MissingPlacement','last');
for ci = 1:length(trace_cands_s1)
    cand_d = trace_cands_s1{sord_disp(ci)};
    fprintf('  %2d : %-35s  R2(Full)=%.3f\n', ci, ...
            strrep(cand_d.animal,'_',' '), cand_d.r2);
end
fprintf('\nEnter animal indices to save (e.g. [1 3 5]), or 0 to skip: ');
sel = input('');
if isempty(sel) || (isscalar(sel) && sel == 0)
    fprintf('No trace data saved.\n');
    return;
end
sel = sel(sel >= 1 & sel <= length(trace_cands_s1));
if isempty(sel)
    fprintf('No valid indices selected.\n');
    return;
end

% Build column headers from model labels
col_headers = [{'time_s', 'dff_actual'}, ...
               cellfun(@(l) strrep(strrep(l,' ','_'),'+','plus'), ...
                       trace_plot_labs, 'UniformOutput', false)];

for si_sel = sel(:)'
    cand_s = trace_cands_s1{sord_disp(si_sel)};
    mid_s  = cand_s.animal;
    t_sec  = cand_s.time - cand_s.time(1);

    % Build numeric matrix: time | dff | pred1 | pred2 | ...
    n_rows   = length(t_sec);
    mat      = NaN(n_rows, 2 + length(trace_plot_labs));
    mat(:,1) = t_sec;
    mat(:,2) = cand_s.dff;
    for pk = 1:length(trace_plot_labs)
        if pk <= length(cand_s.preds) && ~isempty(cand_s.preds{pk})
            mat(:, 2+pk) = cand_s.preds{pk};
        end
    end

    xl_tbl = array2table(mat, 'VariableNames', col_headers);

    mid_safe_s = regexprep(strrep(mid_s,'-','_'), '[^a-zA-Z0-9_]','');
    xl_trace_fname = sprintf('trace_data_sess1__%s__%s.xlsx', mid_safe_s, cond_safe);
    if length(xl_trace_fname) > 200, xl_trace_fname = [xl_trace_fname(1:190),'.xlsx']; end

    try
        writetable(xl_tbl, xl_trace_fname);
        fprintf('  Saved: %s\n', xl_trace_fname);
    catch ME_xl_tr
        fprintf('  ERROR saving trace data for %s: %s\n', mid_s, ME_xl_tr.message);
    end
end
fprintf('Trace data export complete.\n');
end