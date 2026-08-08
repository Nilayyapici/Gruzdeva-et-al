function glm_results = glm_dff_analysis(mice_all, opts)
% GLM_DFF_ANALYSIS  Fits GLMs predicting z-scored dF/F from behavioral
%   predictors for a chosen state and source.
%
%   CALL (defaults: fasted, food):
%     glm_results = glm_dff_analysis(mice_all);
%     glm_results = glm_dff_analysis(mice_all, opts);
%
%   OPTS FIELDS (all optional):
%     opts.state            : 'fasted' | 'fed'         (default: 'fasted')
%     opts.source           : 'food'   | 'gel'         (default: 'food')
%     opts.excl_grooming    : true | false              (default: false)
%     opts.excl_food_events : true | false              (default: false)
%     opts.excl_abs_time    : true | false              (default: false)
%     opts.lag_sec          : include dF/F(t - N s) as predictor (default: [])
%     opts.excl_lag         : true | false              (default: false)
%     opts.tsf_type         : 'temporal_distance' | 'time_since_food' |
%                             'time_to_food'            (default: 'temporal_distance')
%     opts.n_cv_iter        : number of CV iterations   (default: 100)
%     opts.cv_train_frac    : train fraction for CV     (default: 0.80)
%     opts.baseline_correct : false | 'sliding_mean' | 'sliding_pct' | 'poly'
%                             (default: false)
%     opts.baseline_win_sec : window for sliding methods in seconds (default: 60)
%     opts.baseline_pct     : percentile for 'sliding_pct'          (default: 8)
%
%   RENAMED PREDICTORS (consistent with glm_dff_per_animal):
%     'spatial_distance'  -- distance to food (was 'distance')
%     'temporal_distance' -- min(time_since_food, time_to_food) (was 'temporal_proximity')
%
%   CROSS-VALIDATION:
%     For each model x iteration:
%       1. Randomly assign 80% of rows to train, 20% to test
%       2. Fit on train, predict on test, compute MSE
%     Repeated n_cv_iter times; mean MSE reported.
%
%   SHUFFLE CONTROL:
%     dF/F is circularly shifted by a random offset (10-90% of session length).
%     Shuffle R2 averaged over N_SHUF_R2=20 independent shifts.
%     CV shuffle uses one fresh circular shift per iteration.
%
%   DATA COLUMNS:
%     1-time  2-x  3-y  4-speed  5-dist_to_food  6-path_from_last_visit
%     7-door  8-food_interaction  9-eating  10-grooming  11-dff

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_TIME  = 1;
COL_SPEED = 4;
COL_DIST  = 5;
COL_FOOD  = 8;
COL_EAT   = 9;
COL_GROOM = 10;
COL_DFF   = 11;

FOOD_DIST_THRESH = 10;   % cm -- "at food" threshold

% Number of independent circular shifts for full-dataset shuffle R2.
% Averaging reduces single-sample noise.
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
if ~isfield(opts,'lag_sec'),          opts.lag_sec          = [];       end
if ~isfield(opts,'excl_lag'),         opts.excl_lag         = false;    end
if ~isfield(opts,'tsf_type'),         opts.tsf_type         = 'temporal_distance'; end
% opts.tsf_type:
%   'temporal_distance' : min(time_since_food, time_to_food)  (default)
%   'time_since_food'   : seconds since last approach
%   'time_to_food'      : seconds until next approach
if ~isfield(opts,'n_cv_iter'),        opts.n_cv_iter        = 100;      end
if ~isfield(opts,'cv_train_frac'),    opts.cv_train_frac    = 0.80;     end
if ~isfield(opts,'baseline_correct'), opts.baseline_correct = false;    end
if ~isfield(opts,'baseline_win_sec'), opts.baseline_win_sec = 60;       end
if ~isfield(opts,'baseline_pct'),     opts.baseline_pct     = 8;        end

valid_states  = {'fasted','fed'};
valid_sources = {'food','gel'};
if ~ismember(opts.state,  valid_states)
    error('opts.state must be "fasted" or "fed", got "%s"', opts.state);
end
if ~ismember(opts.source, valid_sources)
    error('opts.source must be "food" or "gel", got "%s"', opts.source);
end

use_lag = ~isempty(opts.lag_sec) && opts.lag_sec > 0 && ~opts.excl_lag;

% Determine temporal variable column name
switch opts.tsf_type
    case 'temporal_distance'
        tp_var = 'temporal_distance';
    case 'time_to_food'
        tp_var = 'time_to_food';
    otherwise  % 'time_since_food'
        tp_var = 'time_since_food';
end

fprintf('\n=== GLM analysis: state=%s  source=%s ===\n\n', opts.state, opts.source);
fprintf('Temporal variable:  %s\n', opts.tsf_type);
if use_lag
    fprintf('NOTE: dff_lag (%.0f s) included as predictor in all models.\n', opts.lag_sec);
end
if ~isequal(opts.baseline_correct, false)
    fprintf('NOTE: baseline correction = %s', char(string(opts.baseline_correct)));
    if ismember(opts.baseline_correct, {'sliding_mean','sliding_pct'})
        fprintf('  (window=%.0f s', opts.baseline_win_sec);
        if strcmp(opts.baseline_correct,'sliding_pct')
            fprintf(', percentile=%.0f', opts.baseline_pct);
        end
        fprintf(')');
    end
    fprintf('\n');
end
fprintf('CV:    random 80/20, %d iterations\n', opts.n_cv_iter);
fprintf('Shuffle: circular, R2 averaged over %d shifts\n', N_SHUF_R2);

warning('off', 'stats:LinearModel:RankDefDesignMat');
warning('off', 'MATLAB:rankDeficientMatrix');
cleanupObj = onCleanup(@() warning('on', 'stats:LinearModel:RankDefDesignMat'));

% ------------------------------------------------------------------ %
%  FILTER                                                              %
% ------------------------------------------------------------------ %
is_state  = strcmp(mice_all(:,2), opts.state);
is_source = strcmp(mice_all(:,3), opts.source);
is_sess0  = contains(mice_all(:,1), '_sess0');
is_sess1  = contains(mice_all(:,1), '_sess1');
keep      = is_state & is_source & (is_sess0 | is_sess1);
mice_sub  = mice_all(keep, :);
if sum(keep) == 0
    error('No sessions found for state="%s" source="%s".', opts.state, opts.source);
end
fprintf('%d sessions selected (%s + %s, sess0+sess1)\n\n', ...
    sum(keep), opts.state, opts.source);

% ------------------------------------------------------------------ %
%  BUILD PREDICTOR TABLES                                             %
% ------------------------------------------------------------------ %
tables_s0 = {};  tables_s1 = {};  mouse_ids = {};

for i = 1:size(mice_sub, 1)
    sname = mice_sub{i,1};
    data  = mice_sub{i,4};
    is_s0 = contains(sname, '_sess0');

    data(isinf(data)) = NaN;

    time_vec = data(:, COL_TIME);
    dist_vec = data(:, COL_DIST);

    % time_since_food
    at_food       = dist_vec <= FOOD_DIST_THRESH;
    time_since_fd = NaN(size(time_vec));
    last_food_t   = NaN;
    for t = 1:length(time_vec)
        if at_food(t),          last_food_t = time_vec(t); end
        if ~isnan(last_food_t), time_since_fd(t) = time_vec(t) - last_food_t; end
    end

    % time_to_food
    time_to_fd  = NaN(size(time_vec));
    next_food_t = NaN;
    for t = length(time_vec):-1:1
        if at_food(t),          next_food_t = time_vec(t); end
        if ~isnan(next_food_t), time_to_fd(t) = next_food_t - time_vec(t); end
    end

    % temporal_distance = min(tsf, ttf), NaN-safe
    % 0 at food visits, grows with time away -- symmetric.
    tsf_for_min = time_since_fd; tsf_for_min(isnan(tsf_for_min)) = Inf;
    ttf_for_min = time_to_fd;   ttf_for_min(isnan(ttf_for_min)) = Inf;
    temp_dist   = min(tsf_for_min, ttf_for_min);
    temp_dist(isinf(temp_dist)) = NaN;   % NaN only where BOTH were NaN

    T = table();
    T.time             = time_vec;
    T.speed            = data(:, COL_SPEED);
    T.spatial_distance = dist_vec;
    T.food_visit       = double(data(:, COL_FOOD) > 0);
    T.eating           = double(data(:, COL_EAT)  > 0);
    T.grooming         = double(data(:, COL_GROOM) > 0);
    T.time_since_food  = time_since_fd;
    T.time_to_food     = time_to_fd;
    T.temporal_distance = temp_dist;
    T.dff              = data(:, COL_DFF);
    T.session          = double(~is_s0) * ones(height(T), 1);
    mouse_id           = regexprep(sname, '_sess\d$', '');
    T.mouse_id         = repmat({mouse_id}, height(T), 1);

    % Availability checks
    tsf_available    = sum(~isnan(T.time_since_food)) > 0;
    tp_col_available = sum(~isnan(T.(tp_var))) > 0;
    if ~tsf_available
        fprintf('  Note: time_since_food all-NaN for %s\n', sname);
    end

    % Drop NaN rows on core predictors
    core = {'spatial_distance','speed','dff'};
    if ~opts.excl_abs_time,    core{end+1} = 'time'; end
    if ~opts.excl_food_events
        core{end+1} = 'food_visit';
        core{end+1} = 'eating';
    end
    if ~opts.excl_grooming,    core{end+1} = 'grooming'; end
    if tp_col_available
        core{end+1} = tp_var;
    elseif tsf_available && ~strcmp(tp_var,'time_since_food')
        core{end+1} = 'time_since_food';
    end
    bad = false(height(T),1);
    for v = 1:length(core), bad = bad | isnan(T.(core{v})); end
    T = T(~bad,:);

    if height(T) < 20
        fprintf('  Skipping %s: %d rows\n', sname, height(T));
        continue;
    end

    % ---- BASELINE CORRECTION (before z-scoring) --------------------
    if ~isequal(opts.baseline_correct, false)
        T.dff = apply_baseline_correction(T.dff, T.time, opts);
    end

    % z-score dF/F
    dff_std = std(T.dff,'omitnan');
    if dff_std > 0
        T.dff = (T.dff - mean(T.dff,'omitnan')) / dff_std;
    end

    % Compute dt before z-scoring time
    dt_raw = median(diff(T.time), 'omitnan');

    % z-score continuous predictors
    cvars = {'spatial_distance','speed'};
    if ~opts.excl_abs_time, cvars{end+1} = 'time'; end
    if tp_col_available
        cvars{end+1} = tp_var;
    elseif tsf_available && ~strcmp(tp_var,'time_since_food')
        cvars{end+1} = 'time_since_food';
    end
    for v = 1:length(cvars)
        c = T.(cvars{v}); s = std(c,'omitnan');
        if s > 0, T.(cvars{v}) = (c - mean(c,'omitnan')) / s; end
    end

    % Lagged dF/F predictor
    if use_lag
        n_lag   = max(round(opts.lag_sec / dt_raw), 1);
        dff_lag = NaN(height(T), 1);
        dff_lag(n_lag+1:end) = T.dff(1:end-n_lag);
        T.dff_lag = dff_lag;
        T = T(~isnan(T.dff_lag), :);
        lg = T.dff_lag; s_lg = std(lg,'omitnan');
        if s_lg > 0, T.dff_lag = (lg - mean(lg,'omitnan')) / s_lg; end
    end

    fprintf('  %s: %d valid rows\n', sname, height(T));
    if is_s0, tables_s0{end+1} = T; else, tables_s1{end+1} = T; end
    mouse_ids{end+1} = mouse_id;
end

T_s0 = []; T_s1 = [];
if ~isempty(tables_s0), T_s0 = vertcat(tables_s0{:}); end
if ~isempty(tables_s1), T_s1 = vertcat(tables_s1{:}); end
fprintf('\nPooled: sess0=%d rows, sess1=%d rows\n', height(T_s0), height(T_s1));

% ================================================================== %
%  MODEL FORMULAS                                                      %
% ================================================================== %
base_preds = {};
if ~opts.excl_abs_time,    base_preds{end+1} = 'time'; end
base_preds{end+1} = 'speed';
if ~opts.excl_food_events
    base_preds{end+1} = 'food_visit';
    base_preds{end+1} = 'eating';
end
if ~opts.excl_grooming,    base_preds{end+1} = 'grooming'; end
if use_lag,                base_preds{end+1} = 'dff_lag';  end
base = ['dff ~ ', strjoin(base_preds, ' + ')];

% Six models -- same structure as glm_dff_per_animal:
%   ALL_PRED      = base + spatial_distance + tp_var
%   NO_DIST       = base + tp_var
%   NO_TP         = base + spatial_distance
%   NO_DIST_NO_TP = base
%   DIST_ALONE    = dff ~ spatial_distance
%   TP_ALONE      = dff ~ tp_var
model_keys = {'ALL_PRED','NO_DIST','NO_TP','NO_DIST_NO_TP','DIST_ALONE','TP_ALONE'};
model_labs = {'ALL PRED','NO SPAT DIST','NO TEMP DIST','NO DIST+TD', ...
    'SPAT DIST ALONE','TEMP DIST ALONE'};
n_models   = length(model_keys);

% Build per-session formulas (tsf fallback if tp_var all-NaN in pooled table)
f_s0 = build_formulas(base, tp_var, T_s0, model_keys);
f_s1 = build_formulas(base, tp_var, T_s1, model_keys);

% ================================================================== %
%  FIT ALL MODELS (full dataset)                                       %
% ================================================================== %
fprintf('\n--- Fitting session 0 models ---\n');
mdl_s0 = fit_all_models(T_s0, f_s0, model_keys);
glm_results.sess0 = mdl_s0;

fprintf('\n--- Fitting session 1 models ---\n');
mdl_s1 = fit_all_models(T_s1, f_s1, model_keys);
glm_results.sess1 = mdl_s1;

% ================================================================== %
%  PRINT SUMMARIES                                                     %
% ================================================================== %
print_session_summary(mdl_s0, 'SESSION 0 (pre-discovery)', model_keys, model_labs);
print_session_summary(mdl_s1, 'SESSION 1 (post-discovery)', model_keys, model_labs);

% ================================================================== %
%  VIF CHECK                                                           %
% ================================================================== %
fprintf('\n[sess0] VIF:\n'); compute_vif(T_s0, f_s0.ALL_PRED);
fprintf('\n[sess1] VIF:\n'); compute_vif(T_s1, f_s1.ALL_PRED);

% ================================================================== %
%  SHUFFLE R2 (full dataset, circular, averaged over N_SHUF_R2)       %
% ================================================================== %
fprintf('\n--- Computing shuffle R2 (circular, %d shifts per model) ---\n', N_SHUF_R2);
r2_real  = NaN(2, n_models);
r2_shuf  = NaN(2, n_models);

for si = 1:2
    if si == 1, T_sp = T_s0; f_sp = f_s0; else, T_sp = T_s1; f_sp = f_s1; end
    if isempty(T_sp), continue; end
    n_sp      = height(T_sp);
    min_shift = max(1, round(0.10 * n_sp));
    max_shift = round(0.90 * n_sp);

    for ki = 1:n_models
        key = model_keys{ki};
        % Real R2
        if isfield(mdl_s0,'ALL_PRED') || si==2  % already fitted above
            sess_mdl = get_model(mdl_s0, mdl_s1, si);
            if isfield(sess_mdl, key) && isobject(sess_mdl.(key))
                r2_real(si, ki) = sess_mdl.(key).Rsquared.Adjusted;
            end
        end
        % Shuffle R2: average over N_SHUF_R2 circular shifts
        r2_acc = NaN(N_SHUF_R2, 1);
        for sh = 1:N_SHUF_R2
            shift_amt    = min_shift + randi(max_shift - min_shift + 1) - 1;
            T_sh         = T_sp;
            T_sh.dff     = circshift(T_sp.dff, shift_amt);
            try
                mdl_sh       = safe_fitlm(T_sh, f_sp.(key));
                r2_acc(sh)   = mdl_sh.Rsquared.Adjusted;
            catch
            end
        end
        r2_shuf(si, ki) = nanmean(r2_acc);
    end
end

fprintf('\n%-18s  %s\n', '', 'Sess0                          Sess1');
fprintf('%-18s  %-10s  %-10s  %-10s  %-10s\n', 'Model', 'R2_real','R2_shuf','R2_real','R2_shuf');
fprintf('%s\n', repmat('-',1,65));
for ki = 1:n_models
    fprintf('%-18s  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        model_labs{ki}, r2_real(1,ki), r2_shuf(1,ki), r2_real(2,ki), r2_shuf(2,ki));
end

% ================================================================== %
%  CROSS-VALIDATION (real + circular shuffle)                         %
% ================================================================== %
fprintf('\n--- Cross-validation (%d iterations, 80/20 random split) ---\n', opts.n_cv_iter);

mse_real = NaN(2, n_models);
mse_shuf = NaN(2, n_models);

for si = 1:2
    if si == 1, T_sp = T_s0; f_sp = f_s0; else, T_sp = T_s1; f_sp = f_s1; end
    if isempty(T_sp), continue; end

    n_sp      = height(T_sp);
    n_train   = round(n_sp * opts.cv_train_frac);
    min_shift = max(1, round(0.10 * n_sp));
    max_shift = round(0.90 * n_sp);

    mse_real_iter = NaN(opts.n_cv_iter, n_models);
    mse_shuf_iter = NaN(opts.n_cv_iter, n_models);

    for it = 1:opts.n_cv_iter
        % Random 80/20 split
        perm   = randperm(n_sp);
        idx_tr = perm(1:n_train);
        idx_te = perm(n_train+1:end);
        T_tr   = T_sp(idx_tr, :);
        T_te   = T_sp(idx_te, :);

        % Circular shift for shuffle
        shift_amt        = min_shift + randi(max_shift - min_shift + 1) - 1;
        dff_shifted      = circshift(T_sp.dff, shift_amt);
        T_tr_sh          = T_tr;
        T_te_sh          = T_te;
        T_tr_sh.dff      = dff_shifted(idx_tr);
        T_te_sh.dff      = dff_shifted(idx_te);

        for ki = 1:n_models
            f = f_sp.(model_keys{ki});
            % Real CV
            try
                mdl_tr  = safe_fitlm(T_tr, f);
                pred_te = predict(mdl_tr, T_te);
                mse_real_iter(it,ki) = mean((T_te.dff - pred_te).^2);
            catch
            end
            % Shuffle CV
            try
                mdl_sh  = safe_fitlm(T_tr_sh, f);
                pred_sh = predict(mdl_sh, T_te_sh);
                mse_shuf_iter(it,ki) = mean((T_te_sh.dff - pred_sh).^2);
            catch
            end
        end
    end

    mse_real(si,:) = nanmean(mse_real_iter, 1);
    mse_shuf(si,:) = nanmean(mse_shuf_iter, 1);
    fprintf('  Sess%d done.\n', si-1);
end

fprintf('\n%-18s  %s\n', '', 'Sess0                          Sess1');
fprintf('%-18s  %-10s  %-10s  %-10s  %-10s\n', 'Model','MSE_real','MSE_shuf','MSE_real','MSE_shuf');
fprintf('%s\n', repmat('-',1,65));
for ki = 1:n_models
    fprintf('%-18s  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n', ...
        model_labs{ki}, mse_real(1,ki), mse_shuf(1,ki), mse_real(2,ki), mse_shuf(2,ki));
end

% ================================================================== %
%  PER-MOUSE GLMs                                                      %
% ================================================================== %
fprintf('\n============================================================\n');
fprintf('  PER-MOUSE GLMs\n');
fprintf('============================================================\n');

unique_mice = unique(mouse_ids,'stable');
n_mice      = length(unique_mice);

fields = {'r2_full','r2_nodist','delta_r2','beta_dist','pval_dist'};
pm_s0 = struct(); pm_s1 = struct();
for fld = 1:length(fields)
    pm_s0.(fields{fld}) = NaN(n_mice,1);
    pm_s1.(fields{fld}) = NaN(n_mice,1);
end

for m = 1:n_mice
    mid = unique_mice{m};
    for si = 1:2
        if si==1, T_sp=T_s0; f_sp=f_s0; pm=pm_s0;
        else,     T_sp=T_s1; f_sp=f_s1; pm=pm_s1; end
        if isempty(T_sp), continue; end
        Tm = T_sp(strcmp(T_sp.mouse_id, mid), :);
        if height(Tm) < 20, continue; end
        try
            mf = safe_fitlm(Tm, f_sp.ALL_PRED);
            mr = safe_fitlm(Tm, f_sp.NO_DIST_NO_TP);
            pm.r2_full(m)   = mf.Rsquared.Adjusted;
            pm.r2_nodist(m) = mr.Rsquared.Adjusted;
            pm.delta_r2(m)  = pm.r2_full(m) - pm.r2_nodist(m);
            ct = mf.Coefficients;
            dr = strcmp(ct.Properties.RowNames, 'spatial_distance');
            if any(dr)
                pm.beta_dist(m) = ct.Estimate(dr);
                pm.pval_dist(m) = ct.pValue(dr);
            end
        catch
        end
        if si==1, pm_s0=pm; else, pm_s1=pm; end
    end
end

short_labels = cellfun(@(s) s(max(1,end-6):end), unique_mice, 'UniformOutput', false);
fprintf('\nMouse                    | S0 beta | S0 dR2  | S0 p   | S1 beta | S1 dR2  | S1 p\n');
fprintf('-------------------------+---------+---------+--------+---------+---------+------\n');
for m = 1:n_mice
    fprintf('%-25s| %7.4f | %7.4f | %.3f  | %7.4f | %7.4f | %.3f\n', ...
        unique_mice{m}, pm_s0.beta_dist(m), pm_s0.delta_r2(m), pm_s0.pval_dist(m), ...
        pm_s1.beta_dist(m), pm_s1.delta_r2(m), pm_s1.pval_dist(m));
end

for sp = 1:2
    if sp==1, sl='sess0'; pm=pm_s0; else, sl='sess1'; pm=pm_s1; end
    v  = pm.delta_r2(~isnan(pm.delta_r2));
    bv = pm.beta_dist(~isnan(pm.beta_dist));
    if length(v) >= 3
        [~,p_dr2,~,st_dr2] = ttest(v);
        [~,p_b,~,st_b]     = ttest(bv);
        fprintf('\n[%s] dR2: mean=%.4f SEM=%.4f t(%d)=%.3f p=%.4f\n', ...
            sl, mean(v), std(v)/sqrt(length(v)), st_dr2.df, st_dr2.tstat, p_dr2);
        fprintf('[%s] beta(spatial_dist): mean=%.4f SEM=%.4f t(%d)=%.3f p=%.4f\n', ...
            sl, mean(bv), std(bv)/sqrt(length(bv)), st_b.df, st_b.tstat, p_b);
    end
end

% ================================================================== %
%  STORE RESULTS                                                       %
% ================================================================== %
glm_results.r2_real        = r2_real;
glm_results.r2_shuf        = r2_shuf;
glm_results.mse_real       = mse_real;
glm_results.mse_shuf       = mse_shuf;
glm_results.model_keys     = model_keys;
glm_results.model_labs     = model_labs;
glm_results.per_mouse.mouse_ids = unique_mice;
glm_results.per_mouse.sess0     = pm_s0;
glm_results.per_mouse.sess1     = pm_s1;

% ================================================================== %
%  PLOTS                                                               %
% ================================================================== %
cond_parts = {opts.state, opts.source};
if opts.excl_abs_time,    cond_parts{end+1} = 'no-time'; end
if opts.excl_grooming,    cond_parts{end+1} = 'no-groom'; end
if opts.excl_food_events, cond_parts{end+1} = 'no-food-ev'; end
if use_lag,               cond_parts{end+1} = sprintf('lag%.0fs', opts.lag_sec); end
cond_parts{end+1} = strrep(opts.tsf_type,'_','-');
if ~isequal(opts.baseline_correct, false)
    cond_parts{end+1} = sprintf('bl-%s', char(string(opts.baseline_correct)));
end
condition_label = strjoin(cond_parts, ' | ');

plot_all(glm_results, unique_mice, short_labels, T_s0, T_s1, f_s0, f_s1, ...
    mdl_s0, mdl_s1, model_keys, model_labs, r2_real, r2_shuf, ...
    mse_real, mse_shuf, condition_label, n_mice);

fprintf('\nDone.\n');
end


% ====================================================================== %
%  BUILD FORMULAS (per session, handles tp_var fallback)                 %
% ====================================================================== %
function f = build_formulas(base, tp_var, T_sp, model_keys)
% Determine if chosen tp_var is available in this pooled table
if isempty(T_sp)
    tp_available = false;
    tsf_available = false;
else
    tp_available  = ismember(tp_var, T_sp.Properties.VariableNames) && ...
        sum(~isnan(T_sp.(tp_var))) > 0;
    tsf_available = ismember('time_since_food', T_sp.Properties.VariableNames) && ...
        sum(~isnan(T_sp.time_since_food)) > 0;
end

% Choose effective temporal variable for formulas
if tp_available
    tv = tp_var;
elseif tsf_available && ~strcmp(tp_var,'time_since_food')
    tv = 'time_since_food';
    fprintf('  Note: %s unavailable -- falling back to time_since_food\n', tp_var);
else
    tv = '';
end

if ~isempty(tv)
    f.ALL_PRED      = [base, ' + spatial_distance + ', tv];
    f.NO_DIST       = [base, ' + ', tv];
    f.NO_TP         = [base, ' + spatial_distance'];
    f.NO_DIST_NO_TP =  base;
    f.DIST_ALONE    = 'dff ~ spatial_distance';
    f.TP_ALONE      = ['dff ~ ', tv];
else
    f.ALL_PRED      = [base, ' + spatial_distance'];
    f.NO_DIST       =  base;
    f.NO_TP         = [base, ' + spatial_distance'];
    f.NO_DIST_NO_TP =  base;
    f.DIST_ALONE    = 'dff ~ spatial_distance';
    f.TP_ALONE      =  base;
end
% Aliases used in per-mouse and VIF sections
f.FULL        = f.ALL_PRED;
f.BASE        = base;
f.NO_DIST_OLD = f.NO_DIST;
end


% ====================================================================== %
%  FIT ALL MODELS                                                         %
% ====================================================================== %
function mdl = fit_all_models(T, f, model_keys)
mdl = struct();
if isempty(T), return; end
for k = 1:length(model_keys)
    key = model_keys{k};
    if ~isfield(f, key), continue; end
    try
        mdl.(key) = safe_fitlm(T, f.(key));
        fprintf('  %-18s  Adj R2 = %.4f\n', key, mdl.(key).Rsquared.Adjusted);
    catch ME
        fprintf('  Could not fit %s: %s\n', key, ME.message);
    end
end
% Aliases
if isfield(mdl,'ALL_PRED'), mdl.FULL = mdl.ALL_PRED; end
end


% ====================================================================== %
%  SAFE FITLM                                                             %
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
        if isnumeric(col) && std(col,'omitnan') > 0
            keep_preds{end+1} = pn;
        end
    end
end
if isempty(keep_preds)
    error('No valid predictors remain.');
end
mdl = fitlm(T, [outcome, ' ~ ', strjoin(keep_preds, ' + ')]);
end


% ====================================================================== %
%  HELPER: get model struct for session index                             %
% ====================================================================== %
function sess_mdl = get_model(mdl_s0, mdl_s1, si)
if si == 1, sess_mdl = mdl_s0; else, sess_mdl = mdl_s1; end
end


% ====================================================================== %
%  PRINT SESSION SUMMARY                                                  %
% ====================================================================== %
function print_session_summary(mdl, label, model_keys, model_labs)
fprintf('\n============================================================\n');
fprintf('  %s\n', label);
fprintf('============================================================\n');
fprintf('  %-18s  Adj R2\n','Model');
fprintf('  %-18s  ------\n','-----');
for k = 1:length(model_keys)
    key = model_keys{k};
    if isfield(mdl,key) && isobject(mdl.(key))
        fprintf('  %-18s  %.4f\n', model_labs{k}, mdl.(key).Rsquared.Adjusted);
    end
end
% Unique contribution of spatial_distance
if isfield(mdl,'ALL_PRED') && isfield(mdl,'NO_DIST') && ...
        isobject(mdl.ALL_PRED) && isobject(mdl.NO_DIST)
    dr2   = mdl.ALL_PRED.Rsquared.Adjusted - mdl.NO_DIST.Rsquared.Adjusted;
    rss_f = sum(mdl.ALL_PRED.Residuals.Raw.^2);
    rss_r = sum(mdl.NO_DIST.Residuals.Raw.^2);
    df1   = mdl.ALL_PRED.NumCoefficients - mdl.NO_DIST.NumCoefficients;
    df2   = mdl.ALL_PRED.NumObservations  - mdl.ALL_PRED.NumCoefficients;
    if df1 > 0
        F = ((rss_r-rss_f)/df1) / (rss_f/df2);
        p = 1 - fcdf(F, df1, df2);
        fprintf('\n  dR2(spatial_distance, unique) = %.4f (%.2f%%)  F(%d,%d)=%.2f  p=%.4g\n', ...
            dr2, dr2*100, df1, df2, F, p);
    end
end
% Unique contribution of temporal_distance
if isfield(mdl,'ALL_PRED') && isfield(mdl,'NO_TP') && ...
        isobject(mdl.ALL_PRED) && isobject(mdl.NO_TP)
    dr2_tp = mdl.ALL_PRED.Rsquared.Adjusted - mdl.NO_TP.Rsquared.Adjusted;
    fprintf('  dR2(temporal_distance, unique) = %.4f (%.2f%%)\n', dr2_tp, dr2_tp*100);
end
if isfield(mdl,'ALL_PRED') && isobject(mdl.ALL_PRED)
    fprintf('\nALL_PRED model coefficients:\n');
    disp(mdl.ALL_PRED.Coefficients);
end
end


% ====================================================================== %
%  VIF CHECK                                                              %
% ====================================================================== %
function compute_vif(T, formula)
if isempty(T), return; end
rhs   = strtrim(strsplit(formula,'~'));
preds = strtrim(strsplit(rhs{2},' + '));
preds(strcmp(preds,'dff')) = [];
numeric_preds = {};
for k = 1:length(preds)
    pn = strtrim(preds{k});
    if ismember(pn, T.Properties.VariableNames)
        col = T.(pn);
        if isnumeric(col) && std(col,'omitnan') > 0
            numeric_preds{end+1} = pn;
        end
    end
end
fprintf('  %-22s  VIF\n','Predictor');
for k = 1:length(numeric_preds)
    target = numeric_preds{k};
    others = numeric_preds(~strcmp(numeric_preds, target));
    if isempty(others), continue; end
    try
        mv  = fitlm(T, [target,' ~ ',strjoin(others,' + ')]);
        vif = 1 / max(1 - mv.Rsquared.Ordinary, 1e-10);
        flag = '';
        if vif > 10, flag = '***HIGH'; elseif vif > 5, flag = '**MOD'; end
        fprintf('  %-22s  %.2f  %s\n', target, vif, flag);
    catch
        fprintf('  %-22s  (could not compute)\n', target);
    end
end
end


% ====================================================================== %
%  BASELINE CORRECTION                                                    %
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
    if strcmp(method,'sliding_mean')
        baseline(t) = mean(seg);
    elseif strcmp(method,'sliding_pct')
        baseline(t) = prctile(seg, opts.baseline_pct);
    end
end
valid = ~isnan(baseline);
dff_out(valid) = dff(valid) - baseline(valid);
end


% ====================================================================== %
%  PLOTS                                                                  %
% ====================================================================== %
function plot_all(glm_results, unique_mice, short_labels, T_s0, T_s1, f_s0, f_s1, ...
    mdl_s0, mdl_s1, model_keys, model_labs, r2_real, r2_shuf, ...
    mse_real, mse_shuf, condition_label, n_mice)

has_s0 = ~isempty(T_s0);
has_s1 = ~isempty(T_s1);
n_models = length(model_keys);

bar_colors = [0.15 0.35 0.75;   % ALL PRED          - dark blue
    0.95 0.55 0.1;    % NO SPAT DIST      - orange
    0.55 0.15 0.65;   % NO TEMP DIST      - purple
    0.65 0.65 0.65;   % NO DIST+TD        - grey
    0.9  0.25 0.2;    % SPAT DIST ALONE   - red
    0.2  0.65 0.85];  % TEMP DIST ALONE   - teal
shuf_alpha = 0.35;
bw         = 0.3;

% ---- Figure 1: Coefficient plot ----------------------------------- %
figure('Name','GLM Coefficients (ALL PRED model)','Position',[50 50 1000 480]);
sess_list_coef = {};
if has_s0 && isfield(mdl_s0,'ALL_PRED') && isobject(mdl_s0.ALL_PRED)
    sess_list_coef{end+1} = {'Sess 0 (pre-discovery)', mdl_s0.ALL_PRED};
end
if has_s1 && isfield(mdl_s1,'ALL_PRED') && isobject(mdl_s1.ALL_PRED)
    sess_list_coef{end+1} = {'Sess 1 (post-discovery)', mdl_s1.ALL_PRED};
end

for sp = 1:length(sess_list_coef)
    ttl  = sess_list_coef{sp}{1};
    mdlF = sess_list_coef{sp}{2};
    subplot(1, length(sess_list_coef), sp);
    coef_tbl  = mdlF.Coefficients(2:end,:);
    valid_rows = coef_tbl.SE > 0;
    coef_tbl  = coef_tbl(valid_rows,:);
    n_coef    = height(coef_tbl);
    names     = coef_tbl.Properties.RowNames;

    colors = repmat([0.4 0.6 0.9], n_coef, 1);
    for k = 1:n_coef
        if strcmp(names{k},'spatial_distance'),  colors(k,:) = [0.9  0.3  0.2 ]; end
        if strcmp(names{k},'temporal_distance'), colors(k,:) = [0.2  0.65 0.85]; end
        if strcmp(names{k},'time_since_food'),   colors(k,:) = [0.3  0.75 0.45]; end
        if strcmp(names{k},'dff_lag'),           colors(k,:) = [0.6  0.2  0.7 ]; end
    end

    hb = barh(1:n_coef, coef_tbl.Estimate); hb.FaceColor = 'flat';
    for k = 1:n_coef, hb.CData(k,:) = colors(k,:); end
    hold on;
    ci = coef_tbl.SE * 1.96;
    errorbar(coef_tbl.Estimate, 1:n_coef, ci, 'horizontal', 'k.', 'LineWidth', 1.2);
    for k = 1:n_coef
        pv = coef_tbl.pValue(k);
        if     pv < 0.001, sig = '***'; elseif pv < 0.01, sig = '**';
        elseif pv < 0.05,  sig = '*';   else,              sig = ''; end
        if ~isempty(sig)
            text(coef_tbl.Estimate(k)+ci(k)+0.01, k, sig, 'FontSize',9,'VerticalAlignment','middle');
        end
    end
    xline(0,'k--');
    set(gca,'YTick',1:n_coef,'YTickLabel',names,'FontSize',9);
    xlabel('Std \beta (per SD of predictor)');
    dr2 = NaN;
    if isfield(mdl_s0,'NO_DIST') && sp==1 && isobject(mdl_s0.NO_DIST)
        dr2 = mdl_s0.ALL_PRED.Rsquared.Adjusted - mdl_s0.NO_DIST.Rsquared.Adjusted;
    elseif isfield(mdl_s1,'NO_DIST') && sp==2 && isobject(mdl_s1.NO_DIST)
        dr2 = mdl_s1.ALL_PRED.Rsquared.Adjusted - mdl_s1.NO_DIST.Rsquared.Adjusted;
    end
    if ~isnan(dr2)
        title(sprintf('%s\nAdj R^2=%.3f  dR^2(spat dist)=%.4f', ...
            ttl, mdlF.Rsquared.Adjusted, dr2), 'FontSize',10);
    else
        title(sprintf('%s\nAdj R^2=%.3f', ttl, mdlF.Rsquared.Adjusted), 'FontSize',10);
    end
    grid on; box off;
end
sgtitle(sprintf('ALL PRED model coefficients  [%s]\n(red=spat dist  teal=temp dist  purple=lag  blue=others)', ...
    condition_label), 'FontSize',11);

% ---- Figure 2: R2 real vs shuffled -------------------------------- %
% Shared y-axis
all_r2 = [r2_real(:); r2_shuf(:)];
ymax_r2 = max(all_r2(isfinite(all_r2)), [], 'omitnan') * 1.35;
if isempty(ymax_r2) || ymax_r2 <= 0, ymax_r2 = 0.1; end

figure('Name','R2: real vs shuffled (circular)','Position',[50 580 1200 500]);
sess_titles = {'Sess 0 (pre-discovery)','Sess 1 (post-discovery)'};
for sp = 1:2
    subplot(1,2,sp); hold on;
    r2r = r2_real(sp,:);
    r2s = r2_shuf(sp,:);

    for ki = 1:n_models
        col_r = bar_colors(ki,:);
        col_s = col_r + (1 - col_r) * shuf_alpha;
        if ~isnan(r2r(ki))
            bar(ki-bw/2, r2r(ki), bw, 'FaceColor', col_r, 'EdgeColor','none');
        end
        if ~isnan(r2s(ki))
            bar(ki+bw/2, r2s(ki), bw, 'FaceColor', col_s, 'EdgeColor', col_r, ...
                'LineWidth',0.8,'LineStyle','--');
        end
        % Bracket: dR2 for ALL_PRED vs NO_DIST (ki=1 vs ki=2)
        if ki == 1 && ~isnan(r2r(1)) && ~isnan(r2r(2))
            yb = max(r2r(1), r2r(2)) + 0.05*ymax_r2;
            line([1 2],[yb yb],'Color','k','LineWidth',1.2);
            line([1 1],[r2r(1) yb],'Color','k','LineWidth',1);
            line([2 2],[r2r(2) yb],'Color','k','LineWidth',1);
            text(1.5, yb+0.02*ymax_r2, sprintf('dR2(spat dist)=%.4f',r2r(1)-r2r(2)), ...
                'HorizontalAlignment','center','FontSize',8,'FontWeight','bold', ...
                'Color',[0.9 0.25 0.2]);
        end
        % dR2 for temporal distance (ALL_PRED vs NO_TP, ki=1 vs ki=3)
        if ki == 1 && ~isnan(r2r(1)) && ~isnan(r2r(3))
            yb2 = max(r2r(1), r2r(3)) + 0.14*ymax_r2;
            line([1 3],[yb2 yb2],'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',1.2);
            text(2, yb2+0.02*ymax_r2, sprintf('dR2(temp dist)=%.4f',r2r(1)-r2r(3)), ...
                'HorizontalAlignment','center','FontSize',8,'Color',[0.2 0.65 0.85]);
        end
    end
    set(gca,'XTick',1:n_models,'XTickLabel',model_labs,'XTickLabelRotation',30,'FontSize',9);
    ylabel('Adjusted R^2');
    title(sess_titles{sp},'FontSize',11);
    ylim([0 ymax_r2]); grid on; box off;
end
subplot(1,2,1);
h1 = bar(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
h2 = bar(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
legend([h1 h2],{'Real dF/F',sprintf('Shuffled (circ, N=%d)',20)},'Location','northwest','FontSize',8);
sgtitle(sprintf('Model R^2 — real vs shuffled  [%s]',condition_label),'FontSize',12);

% ---- Save Excel: R2 real/shuf and MSE real/shuf across models ----- %
try
    % Headers: Model | S0_R2_real | S0_R2_shuf | S1_R2_real | S1_R2_shuf
    %                | S0_MSE_real | S0_MSE_shuf | S1_MSE_real | S1_MSE_shuf
    xl_r2_headers = {'Model', ...
        'S0_R2_real','S0_R2_shuf','S1_R2_real','S1_R2_shuf', ...
        'S0_MSE_real','S0_MSE_shuf','S1_MSE_real','S1_MSE_shuf'};
    xl_r2_data = {};
    for ki = 1:n_models
        xl_r2_data(end+1,:) = { ...
            model_labs{ki}, ...
            r2_real(1,ki), r2_shuf(1,ki), r2_real(2,ki), r2_shuf(2,ki), ...
            mse_real(1,ki), mse_shuf(1,ki), mse_real(2,ki), mse_shuf(2,ki) };
    end
    xl_fname_r2 = sprintf('R2_models_%s.xlsx', ...
        strrep(strrep(condition_label,' | ','_'),' ','_'));
    writetable(cell2table(xl_r2_data, 'VariableNames', xl_r2_headers), xl_fname_r2);
    fprintf('\nR2/MSE Excel saved: %s\n', xl_fname_r2);
catch ME_xl
    fprintf('\nCould not save R2 Excel: %s\n', ME_xl.message);
end

% ---- Figure 3: CV MSE: real vs shuffled --------------------------- %
all_mse = [mse_real(:); mse_shuf(:)];
ymax_mse = max(all_mse(isfinite(all_mse)),[],'omitnan') * 1.25;
if isempty(ymax_mse) || ymax_mse <= 0, ymax_mse = 1; end

figure('Name','CV MSE: real vs shuffled','Position',[50 50 1200 480]);
for sp = 1:2
    subplot(1,2,sp); hold on;
    mr = mse_real(sp,:);
    ms = mse_shuf(sp,:);

    for ki = 1:n_models
        col_r = bar_colors(ki,:);
        col_s = col_r + (1 - col_r) * shuf_alpha;
        if ~isnan(mr(ki))
            bar(ki-bw/2, mr(ki), bw, 'FaceColor', col_r, 'EdgeColor','none');
            text(ki-bw/2, mr(ki)+0.01*ymax_mse, sprintf('%.4f',mr(ki)), ...
                'HorizontalAlignment','center','FontSize',7);
        end
        if ~isnan(ms(ki))
            bar(ki+bw/2, ms(ki), bw, 'FaceColor', col_s, 'EdgeColor', col_r, ...
                'LineWidth',0.8,'LineStyle','--');
        end
    end
    set(gca,'XTick',1:n_models,'XTickLabel',model_labs,'XTickLabelRotation',30,'FontSize',9);
    ylabel(sprintf('CV MSE  (mean over %d iterations)',size(mse_real,1)*0+100));
    title(sess_titles{sp},'FontSize',11);
    ylim([0 ymax_mse]); grid on; box off;
end
subplot(1,2,1);
h1 = bar(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
h2 = bar(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
legend([h1 h2],{'Real dF/F','Shuffled (circular)'},'Location','northeast','FontSize',8);
sgtitle(sprintf('Cross-validation MSE — real vs shuffled  [%s]\n(lower = better prediction)', ...
    condition_label),'FontSize',11);

% ---- Figure 4: Predicted vs observed dF/F for all 6 models -------- %
figure('Name','Predicted vs Observed dF/F','Position',[50 50 1300 700]);
sess_list = {};
if has_s0, sess_list{end+1} = {T_s0, mdl_s0, 'Sess 0'}; end
if has_s1, sess_list{end+1} = {T_s1, mdl_s1, 'Sess 1'}; end

for row = 1:length(sess_list)
    T_sp   = sess_list{row}{1};
    mdl_sp = sess_list{row}{2};
    s_lbl  = sess_list{row}{3};

    n_sub = min(3000, height(T_sp));
    rng(42);
    idx   = sort(randperm(height(T_sp), n_sub));
    T_sub = T_sp(idx,:);
    [~, si] = sort(T_sub.dff);
    T_sub   = T_sub(si,:);

    for col = 1:n_models
        key = model_keys{col};
        subplot(length(sess_list), n_models, (row-1)*n_models + col);
        if ~isfield(mdl_sp, key) || ~isobject(mdl_sp.(key)), continue; end
        try
            pred = predict(mdl_sp.(key), T_sub);
        catch, continue; end
        scatter(T_sub.dff, pred, 4, bar_colors(col,:), 'filled','MarkerFaceAlpha',0.3);
        hold on;
        lims = [min(T_sub.dff) max(T_sub.dff)];
        plot(lims, lims, 'k--', 'LineWidth',1);
        r_val  = corr(T_sub.dff, pred, 'rows','complete');
        r2_val = mdl_sp.(key).Rsquared.Adjusted;
        xlabel('Observed dF/F (z)','FontSize',8);
        if col==1, ylabel(sprintf('%s\nPredicted dF/F (z)',s_lbl),'FontSize',8); end
        title(sprintf('%s\nAdj R^2=%.3f  r=%.3f',model_labs{col},r2_val,r_val),'FontSize',8);
        grid on; box off; axis equal tight;
    end
end
sgtitle(sprintf('Predicted vs Observed dF/F  [%s]',condition_label),'FontSize',11);

% ---- Figure 5: Per-mouse dR2 and beta ----------------------------- %
figure('Name','Per-mouse dR2 and beta','Position',[50 50 1100 550]);

all_bv  = [glm_results.per_mouse.sess0.beta_dist; glm_results.per_mouse.sess1.beta_dist];
all_dr2 = [glm_results.per_mouse.sess0.delta_r2;  glm_results.per_mouse.sess1.delta_r2];
all_bv  = all_bv(~isnan(all_bv));
all_dr2 = all_dr2(~isnan(all_dr2));
if isempty(all_bv)  || range(all_bv)==0,  all_bv  = [-1 1]; end
if isempty(all_dr2) || range(all_dr2)==0, all_dr2 = [ 0 1]; end
pad     = 0.10;
bv_lo   = min(min(all_bv)  - pad*range(all_bv),  0);
bv_hi   = max(all_bv)  + pad*range(all_bv);
dr2_lo  = min(min(all_dr2) - pad*range(all_dr2), 0);
dr2_hi  = max(all_dr2) + pad*range(all_dr2);
if bv_lo  >= bv_hi,  bv_hi  = bv_lo  + 0.1; end
if dr2_lo >= dr2_hi, dr2_hi = dr2_lo + 0.1; end
bv_lim  = [bv_lo,  bv_hi];
dr2_lim = [dr2_lo, dr2_hi];

for sp = 1:2
    if sp==1, sl='sess0'; pm=glm_results.per_mouse.sess0;
    else,     sl='sess1'; pm=glm_results.per_mouse.sess1; end
    dr2 = pm.delta_r2; pv = pm.pval_dist; bv = pm.beta_dist;
    cols = zeros(n_mice,3);
    for m = 1:n_mice
        if ~isnan(pv(m)) && pv(m)<0.05, cols(m,:)=[0.85 0.25 0.2];
        else, cols(m,:)=[0.7 0.7 0.7]; end
    end

    subplot(2,2,(sp-1)*2+1); hold on;
    for m = 1:n_mice
        if ~isnan(bv(m)), bar(m,bv(m),'FaceColor',cols(m,:)); end
    end
    yline(0,'k--'); ylim(bv_lim);
    set(gca,'XTick',1:n_mice,'XTickLabel',short_labels,'XTickLabelRotation',45,'FontSize',8);
    ylabel('beta(spatial\_distance)');
    valid_bv = bv(~isnan(bv));
    if length(valid_bv)>=3
        [~,p_b,~,st_b] = ttest(valid_bv);
        title(sprintf('%s  beta(spat dist)\nmean=%.3f  p=%.3f',sl,mean(valid_bv),p_b),'FontSize',9);
    end
    grid on; box off;

    subplot(2,2,(sp-1)*2+2); hold on;
    for m = 1:n_mice
        if ~isnan(dr2(m)), bar(m,dr2(m),'FaceColor',cols(m,:)); end
    end
    yline(0,'k--'); ylim(dr2_lim);
    set(gca,'XTick',1:n_mice,'XTickLabel',short_labels,'XTickLabelRotation',45,'FontSize',8);
    ylabel('\DeltaR^2');
    valid_dr2 = dr2(~isnan(dr2));
    if length(valid_dr2)>=3
        [~,p_d,~,st_d] = ttest(valid_dr2);
        title(sprintf('%s  dR^2(spat dist)\nmean=%.4f  p=%.3f',sl,mean(valid_dr2),p_d),'FontSize',9);
    end
    grid on; box off;
end
sgtitle(sprintf('Per-mouse spatial distance effect  [%s]  (red=p<0.05)',condition_label),'FontSize',11);

% ---- Save Excel: per-mouse beta and dR2 --------------------------- %
try
    % Sheet 1: coefficients from ALL_PRED model for each session
    % Rows: predictor names | Cols: mouse x session
    coef_headers = {'Predictor'};
    for m = 1:n_mice
        coef_headers{end+1} = sprintf('%s_S0', short_labels{m});
        coef_headers{end+1} = sprintf('%s_S1', short_labels{m});
    end

    % Collect all predictor names across both sessions
    all_coef_names = {};
    for si = 1:2
        if si==1, mdl_sp=mdl_s0; else, mdl_sp=mdl_s1; end
        if isfield(mdl_sp,'ALL_PRED') && isobject(mdl_sp.ALL_PRED)
            ct = mdl_sp.ALL_PRED.Coefficients;
            rn = ct.Properties.RowNames;
            rn = rn(~strcmp(rn,'(Intercept)'));
            all_coef_names = union(all_coef_names, rn, 'stable');
        end
    end

    % Build coefficient table: rows=predictors, cols=estimate/SE/p per session
    coef_xl_headers = {'Predictor', ...
        'S0_Estimate','S0_SE','S0_p', ...
        'S1_Estimate','S1_SE','S1_p'};
    coef_xl_data = {};
    for pi = 1:length(all_coef_names)
        pn = all_coef_names{pi};
        row = {pn, NaN, NaN, NaN, NaN, NaN, NaN};
        for si = 1:2
            if si==1, mdl_sp=mdl_s0; else, mdl_sp=mdl_s1; end
            if isfield(mdl_sp,'ALL_PRED') && isobject(mdl_sp.ALL_PRED)
                ct = mdl_sp.ALL_PRED.Coefficients;
                idx_pn = strcmp(ct.Properties.RowNames, pn);
                if any(idx_pn)
                    row{(si-1)*3 + 2} = ct.Estimate(idx_pn);
                    row{(si-1)*3 + 3} = ct.SE(idx_pn);
                    row{(si-1)*3 + 4} = ct.pValue(idx_pn);
                end
            end
        end
        coef_xl_data(end+1,:) = row;
    end

    % Per-mouse beta(spatial_distance) and dR2 table
    pm_headers = {'Mouse', ...
        'S0_beta_spat_dist','S0_dR2','S0_pval', ...
        'S1_beta_spat_dist','S1_dR2','S1_pval'};
    pm_data = {};
    for m = 1:n_mice
        pm_data(end+1,:) = { unique_mice{m}, ...
            glm_results.per_mouse.sess0.beta_dist(m), ...
            glm_results.per_mouse.sess0.delta_r2(m), ...
            glm_results.per_mouse.sess0.pval_dist(m), ...
            glm_results.per_mouse.sess1.beta_dist(m), ...
            glm_results.per_mouse.sess1.delta_r2(m), ...
            glm_results.per_mouse.sess1.pval_dist(m) };
    end
    % Add mean and SEM rows
    for sl_idx = 1:2
        if sl_idx==1, pm=glm_results.per_mouse.sess0; tag='Mean';
        else,          pm=glm_results.per_mouse.sess0; tag='SEM'; end
    end
    % Mean row
    bv0=glm_results.per_mouse.sess0.beta_dist; dr0=glm_results.per_mouse.sess0.delta_r2;
    bv1=glm_results.per_mouse.sess1.beta_dist; dr1=glm_results.per_mouse.sess1.delta_r2;
    pm_data(end+1,:) = {'Mean', ...
        nanmean(bv0), nanmean(dr0), NaN, nanmean(bv1), nanmean(dr1), NaN};
    pm_data(end+1,:) = {'SEM', ...
        nanstd(bv0)/sqrt(sum(~isnan(bv0))), nanstd(dr0)/sqrt(sum(~isnan(dr0))), NaN, ...
        nanstd(bv1)/sqrt(sum(~isnan(bv1))), nanstd(dr1)/sqrt(sum(~isnan(dr1))), NaN};

    xl_fname_beta = sprintf('Coefficients_%s.xlsx', ...
        strrep(strrep(condition_label,' | ','_'),' ','_'));
    writetable(cell2table(coef_xl_data, 'VariableNames', coef_xl_headers), ...
        xl_fname_beta, 'Sheet', 'ALL_PRED_coefficients');
    writetable(cell2table(pm_data, 'VariableNames', pm_headers), ...
        xl_fname_beta, 'Sheet', 'per_mouse_beta_dR2');
    fprintf('Coefficients Excel saved: %s\n', xl_fname_beta);
catch ME_xl
    fprintf('\nCould not save coefficients Excel: %s\n', ME_xl.message);
end

% ---- Figure 6: Partial regression --------------------------------- %
figure('Name','Partial Regression: spatial_distance','Position',[50 580 900 400]);
preg_data = {};
if has_s0, preg_data{end+1} = {T_s0, f_s0.ALL_PRED, 'Sess 0'}; end
if has_s1, preg_data{end+1} = {T_s1, f_s1.ALL_PRED, 'Sess 1'}; end

for sp = 1:length(preg_data)
    T_sp   = preg_data{sp}{1};
    ff     = preg_data{sp}{2};
    ptitle = preg_data{sp}{3};

    ff_nodist       = strrep(ff,' + spatial_distance','');
    ff_dist_on_rest = strrep(ff,'dff ~','spatial_distance ~');
    ff_dist_on_rest = strrep(ff_dist_on_rest,' + spatial_distance','');

    try
        m1 = fitlm(T_sp, ff_nodist);
        m2 = fitlm(T_sp, ff_dist_on_rest);
        res_dff  = m1.Residuals.Raw;
        res_dist = m2.Residuals.Raw;
        subplot(1, length(preg_data), sp);
        n_sub = min(4000, length(res_dff));
        idx   = randperm(length(res_dff), n_sub);
        scatter(res_dist(idx), res_dff(idx), 4, [0.5 0.5 0.7], 'filled','MarkerFaceAlpha',0.25);
        hold on;
        pf = polyfit(res_dist, res_dff, 1);
        xl = linspace(min(res_dist), max(res_dist), 100);
        plot(xl, polyval(pf,xl), 'r-', 'LineWidth', 2);
        [r_p, p_p] = corr(res_dist, res_dff);
        xlabel('spatial\_distance residuals'); ylabel('dF/F residuals');
        title(sprintf('%s\npartial r=%.3f  p=%.4f', ptitle, r_p, p_p),'FontSize',10);
        grid on; box off;
    catch ME
        fprintf('  Partial regression failed for %s: %s\n', ptitle, ME.message);
    end
end
sgtitle(sprintf('Partial regression: spatial\_distance vs dF/F  [%s]',condition_label),'FontSize',11);

% ---- Figure 7: Spatial/temporal distance traces for ALL mice -------- %
col_dist = [0.9 0.25 0.2];
col_td   = [0.2 0.65 0.85];
col_dff  = [0.15 0.15 0.15];

cond_safe = regexprep(strrep(condition_label,' | ','_'), '[^a-zA-Z0-9_]', '_');
cond_safe = regexprep(cond_safe, '_+', '_');

for m_ex = 1:n_mice
    mid_ex = unique_mice{m_ex};

    has_mouse_s0 = has_s0 && any(strcmp(T_s0.mouse_id, mid_ex));
    has_mouse_s1 = has_s1 && any(strcmp(T_s1.mouse_id, mid_ex));

    if ~has_mouse_s0 && ~has_mouse_s1
        continue;
    end

    ex_T = struct();
    if has_mouse_s0
        ex_T.sess0 = T_s0(strcmp(T_s0.mouse_id, mid_ex), :);
    else
        ex_T.sess0 = [];
    end

    if has_mouse_s1
        ex_T.sess1 = T_s1(strcmp(T_s1.mouse_id, mid_ex), :);
    else
        ex_T.sess1 = [];
    end

    fig = figure('Name', sprintf('Distance traces: %s', strrep(mid_ex,'_',' ')), ...
        'Position', [50 50 1300 680]);

    sess_names  = {'sess0','sess1'};
    sess_titles = {'Sess 0 (pre-discovery)','Sess 1 (post-discovery)'};

    for si = 1:2
        Tex = ex_T.(sess_names{si});
        if isempty(Tex), continue; end

        [~, tidx] = sort(Tex.time);
        Tex = Tex(tidx,:);
        t_sec = Tex.time - Tex.time(1);

        subplot(3,2,(si-1)+1); hold on;
        plot(t_sec, Tex.spatial_distance, 'Color', col_dist, 'LineWidth',1);
        yline(0,'k--','LineWidth',0.5);
        ylabel('Spatial distance (z)');
        title(sess_titles{si},'FontSize',10,'FontWeight','bold');
        grid on; box off; xlim([0 max(t_sec)]);
        if si==1
            legend('Spatial distance','FontSize',7,'Location','best');
        end

        subplot(3,2,(si-1)+3); hold on;
        td_col = 'temporal_distance';
        if ismember(td_col, Tex.Properties.VariableNames) && sum(~isnan(Tex.(td_col))) > 0
            plot(t_sec, Tex.(td_col), 'Color', col_td, 'LineWidth',1);
            yline(0,'k--','LineWidth',0.5);
            ylabel('Temporal distance (z)');
            if si==1
                legend('Temporal distance','FontSize',7,'Location','best');
            end
        else
            text(0.5,0.5,'temporal\_distance not available', ...
                'HorizontalAlignment','center','Units','normalized','FontSize',9);
            ylabel('Temporal distance');
        end
        grid on; box off; xlim([0 max(t_sec)]);

        subplot(3,2,(si-1)+5); hold on;
        plot(t_sec, Tex.dff, 'Color', col_dff, 'LineWidth',0.8);
        yline(0,'Color',[0.6 0.6 0.6],'LineWidth',0.5);
        ylabel('dF/F (z)');
        xlabel('Time in session (s)');
        grid on; box off; xlim([0 max(t_sec)]);
        if si==1
            legend('dF/F','FontSize',7,'Location','best');
        end
    end

    sgtitle(sprintf('Distance and dF/F traces: %s  [%s]\nred=spatial dist  teal=temporal dist  grey=dF/F', ...
        strrep(mid_ex,'_',' '), condition_label), 'FontSize',11);

    mid_safe = regexprep(strrep(mid_ex,'-','_'), '[^a-zA-Z0-9_]', '');
    pdf_fname = sprintf('distance_traces__%s__%s.pdf', mid_safe, cond_safe);

    try
        set(fig, 'PaperOrientation', 'landscape', ...
            'PaperUnits', 'normalized', ...
            'PaperPosition', [0 0 1 1]);
        print(fig, pdf_fname, '-dpdf', '-painters', '-bestfit');
        fprintf('  Saved distance trace PDF: %s\n', pdf_fname);
    catch ME_pdf
        fprintf('  ERROR saving distance trace PDF for %s: %s\n', mid_ex, ME_pdf.message);
    end
end

end