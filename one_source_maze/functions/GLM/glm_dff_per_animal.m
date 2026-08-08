function glm_results = glm_dff_per_animal(mice_all, opts)
% GLM_DFF_PER_ANIMAL  Per-animal GLM analysis with cross-validation and
%   shuffle control.  Mirrors glm_dff_analysis() but fits models
%   separately for each animal and plots group-level summaries.
%
%   CALL:
%     glm_results = glm_dff_per_animal(mice_all);
%     glm_results = glm_dff_per_animal(mice_all, opts);
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
%
%   OUTPUT STRUCT fields (glm_results):
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
    COL_FOOD  = 8;
    COL_EAT   = 9;
    COL_GROOM = 10;
    COL_DFF   = 11;
    FOOD_DIST_THRESH = 10;

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

    % Six models:
    %   ALL_PRED      = base + spatial_distance + temporal_distance
    %   NO_DIST       = base + temporal_distance
    %   NO_TP         = base + spatial_distance
    %   NO_DIST_NO_TP = base
    %   DIST_ALONE    = dff ~ spatial_distance
    %   TP_ALONE      = dff ~ temporal_distance
    model_keys = {'ALL_PRED','NO_DIST','NO_TP','NO_DIST_NO_TP','DIST_ALONE','TP_ALONE'};
    model_labs = {'ALL PRED','NO SPAT DIST','NO TEMP DIST','NO DIST+TD','SPAT DIST ALONE','TEMP DIST ALONE'};

    % Determine which temporal variable column to use
    switch opts.tsf_type
        case 'temporal_distance'
            tp_var = 'temporal_distance';
        case 'time_to_food'
            tp_var = 'time_to_food';
        otherwise  % 'time_since_food'
            tp_var = 'time_since_food';
    end

    formulas.ALL_PRED      = [base, ' + spatial_distance + ', tp_var];
    formulas.NO_DIST       = [base, ' + ', tp_var];
    formulas.NO_TP         = [base, ' + spatial_distance'];
    formulas.NO_DIST_NO_TP =  base;
    formulas.DIST_ALONE    = 'dff ~ spatial_distance';
    formulas.TP_ALONE      = ['dff ~ ', tp_var];

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

    all_names   = mice_sub(:,1);
    mouse_ids   = regexprep(all_names, '_sess\d$', '');
    unique_mice = unique(mouse_ids, 'stable');
    n_mice      = length(unique_mice);

    fprintf('%d sessions, %d animals\n\n', sum(keep), n_mice);

    % ------------------------------------------------------------------ %
    %  PRE-ALLOCATE RESULTS                                               %
    % ------------------------------------------------------------------ %
    r2_real           = NaN(n_mice, n_models, 2);
    r2_shuf           = NaN(n_mice, n_models, 2);
    mse_real          = NaN(n_mice, n_models, 2);
    mse_shuf          = NaN(n_mice, n_models, 2);
    beta_dist         = NaN(n_mice, 2);
    n_rows_per_animal = NaN(n_mice, 2);

    example_traces = struct();
    for si = 1:2
        example_traces(si).time     = [];
        example_traces(si).dff_real = [];
        example_traces(si).dff_pred = [];
        example_traces(si).r2       = NaN;
        example_traces(si).animal   = '';
        example_traces(si).r2_all   = [];
        example_traces(si).cands    = {};
    end

    % Predictor contribution list
    pred_names = {};
    if ~opts.excl_abs_time,    pred_names{end+1} = 'time'; end
    pred_names{end+1} = 'speed';
    pred_names{end+1} = 'spatial_distance';
    pred_names{end+1} = tp_var;
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

    % ------------------------------------------------------------------ %
    %  MAIN LOOP: per animal x session                                    %
    % ------------------------------------------------------------------ %
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
            T.food_visit       = double(data(:, COL_FOOD) > 0);
            T.eating           = double(data(:, COL_EAT)  > 0);
            T.grooming         = double(data(:, COL_GROOM) > 0);
            T.time_since_food  = time_since_fd;
            T.time_to_food     = time_to_fd;
            T.temporal_distance = temp_dist;        % renamed from 'temporal_proximity'
            T.dff              = data(:, COL_DFF);

            tsf_available    = sum(~isnan(T.time_since_food)) > 0;
            tp_var_local     = tp_var;
            tp_col_available = sum(~isnan(T.(tp_var_local))) > 0;

            % Drop NaN rows on core predictors.
            % If temporal variable is unavailable, simply omit it from the
            % NaN-check -- do NOT substitute another variable.
            core = {'speed','spatial_distance','dff'};
            if tp_col_available
                core{end+1} = tp_var_local;
            end
            if ~opts.excl_abs_time, core{end+1} = 'time'; end
            bad = false(height(T),1);
            for v = 1:length(core), bad = bad | isnan(T.(core{v})); end
            T = T(~bad,:);

            % Session-specific formulas.
            % If the temporal variable is all-NaN, use the default formulas
            % unchanged -- safe_fitlm will silently drop zero-variance or
            % absent predictors.  We do NOT substitute time_since_food or
            % any other variable.
            sess_formulas = formulas;
            if ~tp_col_available
                fprintf('  Note: %s all-NaN in %s%s -- temporal predictor absent from models\n', ...
                        tp_var_local, mid, sess_str);
            end

            if height(T) < 50
                fprintf('  Skipping %s%s: only %d rows\n', mid, sess_str, height(T));
                continue;
            end

            n_rows_per_animal(m, sess_idx+1) = height(T);

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
            if ~opts.excl_abs_time, cont_vars{end+1} = 'time'; end
            for v = 1:length(cont_vars)
                c = T.(cont_vars{v}); s = std(c,'omitnan');
                if s > 0, T.(cont_vars{v}) = (c - mean(c,'omitnan')) / s; end
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
                fn_lag = fieldnames(sess_formulas_lag);
                for klag = 1:length(fn_lag)
                    sess_formulas_lag.(fn_lag{klag}) = ...
                        [sess_formulas_lag.(fn_lag{klag}), ' + dff_lag'];
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
                    if strcmp(key,'ALL_PRED')
                        ct = mdl.Coefficients;
                        dr = strcmp(ct.Properties.RowNames,'spatial_distance');
                        if any(dr)
                            beta_dist(m, sess_idx+1) = ct.Estimate(dr);
                        end
                        try
                            pred_full   = predict(mdl, T);
                            pred_nodist = [];
                            pred_notp   = [];
                            try
                                mdl_nd = safe_fitlm(T, sess_formulas_lag.NO_DIST);
                                pred_nodist = predict(mdl_nd, T);
                            catch, end
                            try
                                mdl_nt = safe_fitlm(T, sess_formulas_lag.NO_TP);
                                pred_notp = predict(mdl_nt, T);
                            catch, end
                            cand = {T.time, T.dff, pred_full, ...
                                    mdl.Rsquared.Adjusted, mid, ...
                                    pred_nodist, pred_notp};
                            example_traces(sess_idx+1).cands{end+1} = cand;
                        catch
                        end
                    end
                catch
                end
            end

            % ----- PREDICTOR CONTRIBUTIONS (R2 drop when removed) -------
            % Contribution of predictor X = R2(ALL_PRED) - R2(ALL_PRED \ X).
            %
            % IMPORTANT: we use the actual predictor names from the FITTED
            % model coefficients, NOT from the formula string.  safe_fitlm
            % may silently drop zero-variance predictors (e.g. temporal_
            % distance has near-zero variance in some sessions), so the
            % formula string and the fitted model can differ.  If we parsed
            % the formula string, r2_full would already exclude the dropped
            % predictor, making its contribution spuriously ~0 instead of NaN.
            r2_full = r2_real(m, 1, sess_idx+1);
            if ~isnan(r2_full)
                % Re-fit ALL_PRED to get the actual coefficient names
                try
                    mdl_full_check = safe_fitlm(T, sess_formulas_lag.ALL_PRED);
                    % Get predictor names actually in the model
                    % (excludes Intercept, excludes anything dropped by safe_fitlm)
                    fitted_preds = mdl_full_check.CoefficientNames;
                    fitted_preds = fitted_preds(~strcmp(fitted_preds, '(Intercept)'));
                catch
                    fitted_preds = {};
                end

                for pi = 1:n_preds
                    pn = pred_names{pi};
                    % If predictor was not in the FITTED model, store NaN
                    if ~any(strcmp(fitted_preds, pn))
                        r2_contrib(m, pi, sess_idx+1) = NaN;
                        continue;
                    end
                    % Remove pn from the fitted predictor list
                    reduced_preds = fitted_preds(~strcmp(fitted_preds, pn));
                    if isempty(reduced_preds)
                        r2_contrib(m, pi, sess_idx+1) = NaN;
                        continue;
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
            end

            % ----- SHUFFLE R2: average over N_SHUF_R2 circular shifts ---
            % Consistent with the CV shuffle (also circular).
            % Averaged to reduce noise from a single shift sample.
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
            n_train       = round(n * opts.cv_train_frac);
            mse_real_iter = NaN(opts.n_cv_iter, n_models);
            mse_shuf_iter = NaN(opts.n_cv_iter, n_models);

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

                for ki = 1:n_models
                    f = sess_formulas_lag.(model_keys{ki});
                    try
                        mdl_tr  = safe_fitlm(T_tr, f);
                        pred_te = predict(mdl_tr, T_te);
                        mse_real_iter(it,ki) = mean((T_te.dff - pred_te).^2);
                    catch
                    end
                    try
                        mdl_sh  = safe_fitlm(T_tr_sh, f);
                        pred_sh = predict(mdl_sh, T_te_sh);
                        mse_shuf_iter(it,ki) = mean((T_te_sh.dff - pred_sh).^2);
                    catch
                    end
                end
            end

            mse_real(m,:,sess_idx+1) = nanmean(mse_real_iter, 1);
            mse_shuf(m,:,sess_idx+1) = nanmean(mse_shuf_iter, 1);

            fprintf('  %s%s: %d rows, R2(ALL_PRED)=%.3f, MSE(ALL_PRED)=%.4f\n', ...
                mid, sess_str, height(T), ...
                r2_real(m,1,sess_idx+1), mse_real(m,1,sess_idx+1));
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

            valid = ~isnan(r2r) & ~isnan(r2s);
            if sum(valid) < 3, continue; end

            r2r_v = r2r(valid); r2s_v = r2s(valid);
            mr_v  = mr(valid);  ms_v  = ms(valid);

            [~,p_r2]  = ttest(r2r_v, r2s_v);
            [~,p_mse] = ttest(mr_v,  ms_v);

            group.(sl).(model_keys{ki}).r2_real_mean  = mean(r2r_v);
            group.(sl).(model_keys{ki}).r2_real_sem   = std(r2r_v)/sqrt(length(r2r_v));
            group.(sl).(model_keys{ki}).r2_shuf_mean  = mean(r2s_v);
            group.(sl).(model_keys{ki}).r2_shuf_sem   = std(r2s_v)/sqrt(length(r2s_v));
            group.(sl).(model_keys{ki}).mse_real_mean = mean(mr_v);
            group.(sl).(model_keys{ki}).mse_shuf_mean = mean(ms_v);
            group.(sl).(model_keys{ki}).p_r2          = p_r2;
            group.(sl).(model_keys{ki}).p_mse         = p_mse;
            group.(sl).(model_keys{ki}).n             = sum(valid);

            fprintf('%-16s  %.3f +/- %.3f      %.3f +/- %.3f      p=%.4f          %.4f    %.4f\n', ...
                model_labs{ki}, mean(r2r_v), std(r2r_v)/sqrt(length(r2r_v)), ...
                mean(r2s_v), std(r2s_v)/sqrt(length(r2s_v)), p_r2, ...
                mean(mr_v), mean(ms_v));
        end
    end

    % ------------------------------------------------------------------ %
    %  PICK MEDIAN-R2 EXAMPLE ANIMAL FOR TRACE PLOT                       %
    % ------------------------------------------------------------------ %
    for si = 1:2
        cands = example_traces(si).cands;
        if isempty(cands), continue; end
        r2_cands = cellfun(@(c) c{4}, cands);
        med_r2   = median(r2_cands, 'omitnan');
        [~, best_idx] = min(abs(r2_cands - med_r2));
        best = cands{best_idx};
        example_traces(si).time     = best{1};
        example_traces(si).dff_real = best{2};
        example_traces(si).dff_pred = best{3};
        example_traces(si).r2       = best{4};
        example_traces(si).animal   = best{5};
        if length(best) >= 6, example_traces(si).dff_nodist    = best{6}; end
        if length(best) >= 7, example_traces(si).dff_nodisttsf = best{7}; end
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
    glm_results.beta_dist      = beta_dist;
    glm_results.unique_mice    = unique_mice;
    glm_results.model_keys     = model_keys;
    glm_results.model_labs     = model_labs;
    glm_results.group          = group;
    glm_results.n_rows         = n_rows_per_animal;
    glm_results.example_traces = example_traces;

    % ------------------------------------------------------------------ %
    %  CONDITION LABEL FOR PLOTS                                          %
    % ------------------------------------------------------------------ %
    cond_parts = {opts.state, opts.source};
    if opts.excl_abs_time,    cond_parts{end+1} = 'no-time'; end
    if opts.excl_grooming,    cond_parts{end+1} = 'no-groom'; end
    if opts.excl_food_events, cond_parts{end+1} = 'no-food-ev'; end
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

    plot_per_animal(glm_results, condition_label, model_labs, model_keys, pred_names, r2_contrib, example_traces, mice_sub);

    fprintf('\nDone.\n');

    % ------------------------------------------------------------------ %
    %  INTERACTIVE PDF SAVING                                             %
    % ------------------------------------------------------------------ %
    save_pdfs_interactive(unique_mice, condition_label);
end


% ====================================================================== %
%  INTERACTIVE PDF SAVING                                                %
% ====================================================================== %
function save_pdfs_interactive(unique_mice, condition_label)
% Prompts the user to choose which figures to save as PDF.
% Two groups of figures are offered:
%   (A) GLM prediction traces  -- figures named "Traces: ..."
%   (B) Distance traces        -- figures named "Distance traces: ..."

    cond_safe = strrep(strrep(condition_label, ' | ', '_'), ' ', '_');

    % Collect all open figures and classify them by name
    all_figs = findall(0, 'Type', 'figure');
    if isempty(all_figs)
        fprintf('No open figures to save.\n');
        return;
    end

    % Build cell arrays of {handle, name} for each group
    trace_list = {};   % GLM prediction traces
    dist_list  = {};   % Distance time traces

    for fi = 1:length(all_figs)
        try
            fname = get(all_figs(fi), 'Name');
        catch
            continue;  % skip invalid handles
        end
        if strncmp(fname, 'Traces:', 7)
            trace_list{end+1} = all_figs(fi);
        elseif strncmp(fname, 'Distance traces:', 16)
            dist_list{end+1}  = all_figs(fi);
        end
    end

    n_trace = length(trace_list);
    n_dist  = length(dist_list);

    fprintf('\n========================================\n');
    fprintf('  PDF EXPORT\n');
    fprintf('========================================\n');
    fprintf('Found %d GLM prediction trace figure(s).\n', n_trace);
    fprintf('Found %d distance trace figure(s).\n', n_dist);
    fprintf('\nOptions:\n');
    fprintf('  1 - Save ALL GLM prediction trace figures\n');
    fprintf('  2 - Save SELECTED GLM prediction trace figures\n');
    fprintf('  3 - Save ALL distance trace figures\n');
    fprintf('  4 - Save SELECTED distance trace figures\n');
    fprintf('  5 - Save ALL of the above\n');
    fprintf('  0 - Skip / done\n');

    while true
        answer = input('\nEnter option (0-5, or multiple e.g. [1 3]): ');
        if isempty(answer), answer = 0; end
        answer = answer(:)';

        if any(answer == 0)
            fprintf('No figures saved.\n');
            return;
        end

        save_trace_all    = any(answer == 1) || any(answer == 5);
        save_trace_select = any(answer == 2);
        save_dist_all     = any(answer == 3) || any(answer == 5);
        save_dist_select  = any(answer == 4);

        % --- GLM prediction traces ---
        if save_trace_all && n_trace > 0
            fprintf('Saving %d GLM prediction trace PDFs...\n', n_trace);
            for fi = 1:n_trace
                save_fig_pdf(trace_list{fi}, 'traces', cond_safe);
            end
        end

        if save_trace_select && n_trace > 0
            fprintf('\nGLM prediction trace figures:\n');
            for fi = 1:n_trace
                try
                    fprintf('  %2d : %s\n', fi, get(trace_list{fi}, 'Name'));
                catch
                    fprintf('  %2d : (figure unavailable)\n', fi);
                end
            end
            sel = input(sprintf('Enter indices to save (e.g. [1 3 5], max %d): ', n_trace));
            if ~isempty(sel)
                sel = sel(sel >= 1 & sel <= n_trace);
                for fi = sel(:)'
                    save_fig_pdf(trace_list{fi}, 'traces', cond_safe);
                end
            end
        end

        % --- Distance traces ---
        if save_dist_all && n_dist > 0
            fprintf('Saving %d distance trace PDFs...\n', n_dist);
            for fi = 1:n_dist
                save_fig_pdf(dist_list{fi}, 'dist_traces', cond_safe);
            end
        end

        if save_dist_select && n_dist > 0
            fprintf('\nDistance trace figures:\n');
            for fi = 1:n_dist
                try
                    fprintf('  %2d : %s\n', fi, get(dist_list{fi}, 'Name'));
                catch
                    fprintf('  %2d : (figure unavailable)\n', fi);
                end
            end
            sel = input(sprintf('Enter indices to save (e.g. [1 2], max %d): ', n_dist));
            if ~isempty(sel)
                sel = sel(sel >= 1 & sel <= n_dist);
                for fi = sel(:)'
                    save_fig_pdf(dist_list{fi}, 'dist_traces', cond_safe);
                end
            end
        end

        more = input('\nSave more figures? (1=yes, 0=no): ');
        if isempty(more) || more == 0, break; end
    end

    fprintf('PDF export complete.\n');
end


% ====================================================================== %
%  SAVE ONE FIGURE AS PDF                                                %
% ====================================================================== %
function save_fig_pdf(fig_handle, prefix, cond_safe)
% Save fig_handle to a PDF in the current directory.
% Filename is built from prefix + sanitised figure name + condition.

    raw_name = get(fig_handle, 'Name');

    % Sanitise: keep alphanumeric, spaces→underscore, remove special chars
    safe_name = regexprep(raw_name, '[^a-zA-Z0-9_ ]', '');
    safe_name = strtrim(strrep(safe_name, ' ', '_'));
    safe_name = regexprep(safe_name, '_+', '_');  % collapse multiple underscores

    fname = sprintf('%s__%s__%s.pdf', prefix, safe_name, cond_safe);
    % Truncate if very long (filesystem limit)
    if length(fname) > 200
        fname = [fname(1:190), '.pdf'];
    end

    try
        % Use painters renderer for vector output
        set(fig_handle, 'PaperOrientation', 'landscape');
        set(fig_handle, 'PaperUnits',       'normalized');
        set(fig_handle, 'PaperPosition',    [0 0 1 1]);
        print(fig_handle, fname, '-dpdf', '-painters', '-bestfit');
        fprintf('  Saved: %s\n', fname);
    catch ME
        fprintf('  ERROR saving %s: %s\n', fname, ME.message);
    end
end


% ====================================================================== %
%  BASELINE CORRECTION HELPER                                            %
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
%  PLOTS                                                                  %
% ====================================================================== %
function plot_per_animal(res, condition_label, model_labs, model_keys, pred_names, r2_contrib, example_traces, mice_sub)

    % Column indices and threshold -- must redeclare here since subfunctions
    % do not share the main function's workspace in MATLAB.
    COL_TIME         = 1;
    COL_DIST         = 5;
    COL_DFF          = 11;
    FOOD_DIST_THRESH = 10;

    n_models = length(model_keys);
    n_mice   = length(res.unique_mice);
    n_preds  = length(pred_names);
    bw       = 0.3;   % bar half-width (used in Figs 1 and 2)
    sess_labels = {'Sess 0 (pre-discovery)', 'Sess 1 (post-discovery)'};

    bar_colors = [0.15 0.35 0.75;   % ALL PRED          - dark blue
                  0.95 0.55 0.1;    % NO SPAT DIST      - orange
                  0.55 0.15 0.65;   % NO TEMP DIST      - purple
                  0.65 0.65 0.65;   % NO DIST+TD        - grey
                  0.9  0.25 0.2;    % SPAT DIST ALONE   - red
                  0.2  0.65 0.85];  % TEMP DIST ALONE   - teal
    shuf_alpha = 0.35;

    short_labels = cellfun(@(s) s(max(1,end-6):end), res.unique_mice, 'UniformOutput', false);

    % ---- Figure 1: Adjusted R2 per model, real vs shuffled ------------ %
    r2_all_sessions = [];
    for sp_pre = 1:2
        r2r_pre = res.r2_real(:,:,sp_pre);
        r2s_pre = res.r2_shuf(:,:,sp_pre);
        r2r_m   = nanmean(r2r_pre,1) + nanstd(r2r_pre,0,1)./sqrt(sum(~isnan(r2r_pre),1));
        r2s_m   = nanmean(r2s_pre,1) + nanstd(r2s_pre,0,1)./sqrt(sum(~isnan(r2s_pre),1));
        r2_all_sessions = [r2_all_sessions, r2r_m, r2s_m];
    end
    shared_r2_ylim = [0, max(r2_all_sessions,[],'omitnan') * 1.30];

    figure('Name','Per-animal R2: real vs shuffled','Position',[50 50 1200 600]);
    for sp = 1:2
        subplot(1,2,sp); hold on;
        r2r = res.r2_real(:,:,sp);
        r2s = res.r2_shuf(:,:,sp);
        r2r_mean = nanmean(r2r,1);  r2r_sem = nanstd(r2r,0,1)./sqrt(sum(~isnan(r2r),1));
        r2s_mean = nanmean(r2s,1);  r2s_sem = nanstd(r2s,0,1)./sqrt(sum(~isnan(r2s),1));

        for ki = 1:n_models
            col_r = bar_colors(ki,:);
            col_s = col_r + (1-col_r)*shuf_alpha;
            bar(ki-bw/2, r2r_mean(ki), bw, 'FaceColor', col_r, 'EdgeColor', 'none');
            bar(ki+bw/2, r2s_mean(ki), bw, 'FaceColor', col_s, 'EdgeColor', col_r, ...
                'LineWidth', 0.8, 'LineStyle', '--');
            errorbar(ki-bw/2, r2r_mean(ki), r2r_sem(ki), 'k.', 'LineWidth', 1.2, 'CapSize', 4);
            errorbar(ki+bw/2, r2s_mean(ki), r2s_sem(ki), 'Color', col_r*0.6, ...
                     'LineStyle','none', 'LineWidth', 1, 'CapSize', 4);
            valid = ~isnan(r2r(:,ki));
            jitter = (rand(sum(valid),1)-0.5)*0.12;
            scatter(ki-bw/2 + jitter, r2r(valid,ki), 18, col_r*0.7, 'filled', ...
                    'MarkerFaceAlpha', 0.5);
            % Significance stars
            sl_key = sprintf('sess%d', sp-1);
            if isfield(res.group, sl_key) && isfield(res.group.(sl_key), model_keys{ki})
                p_val = res.group.(sl_key).(model_keys{ki}).p_r2;
                if     p_val < 0.001, sig = '***';
                elseif p_val < 0.01,  sig = '**';
                elseif p_val < 0.05,  sig = '*';
                else,                 sig = ''; end
                if ~isempty(sig)
                    ystar = max(r2r_mean(ki)+r2r_sem(ki), r2s_mean(ki)+r2s_sem(ki)) + 0.01;
                    text(ki, ystar, sig, 'HorizontalAlignment','center', ...
                         'FontSize',11,'FontWeight','bold','Color',col_r*0.8);
                end
            end
        end
        set(gca,'XTick',1:n_models,'XTickLabel',model_labs,'XTickLabelRotation',30,'FontSize',9);
        ylabel('Adjusted R^2  (mean \pm SEM)');
        title(sess_labels{sp},'FontSize',11);
        ylim(shared_r2_ylim); grid on; box off;
    end
    subplot(1,2,1);
    h1 = bar(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
    h2 = bar(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
    legend([h1 h2],{'Real dF/F','Shuffled dF/F (circular)'},'Location','northwest','FontSize',8);
    sgtitle(sprintf('Per-animal Adj R2 - real vs shuffled  [%s]\n(* p<0.05, ** p<0.01, *** p<0.001 paired t-test)', ...
            condition_label),'FontSize',11);

    % Save Excel for Fig 1
    try
        xl_headers = {'Animal'};
        for sp_xl = 1:2
            for ki_xl = 1:n_models
                xl_headers{end+1} = sprintf('S%d_%s_real', sp_xl-1, model_labs{ki_xl});
                xl_headers{end+1} = sprintf('S%d_%s_shuf', sp_xl-1, model_labs{ki_xl});
            end
        end
        xl_data = {};
        for mi_xl = 1:n_mice
            row = {res.unique_mice{mi_xl}};
            for sp_xl = 1:2
                for ki_xl = 1:n_models
                    row{end+1} = res.r2_real(mi_xl,ki_xl,sp_xl);
                    row{end+1} = res.r2_shuf(mi_xl,ki_xl,sp_xl);
                end
            end
            xl_data(end+1,:) = row;
        end
        mean_row = {'Mean'}; sem_row = {'SEM'};
        for sp_xl = 1:2
            for ki_xl = 1:n_models
                v_r = res.r2_real(:,ki_xl,sp_xl); v_s = res.r2_shuf(:,ki_xl,sp_xl);
                mean_row{end+1} = nanmean(v_r); mean_row{end+1} = nanmean(v_s);
                sem_row{end+1}  = nanstd(v_r)/sqrt(sum(~isnan(v_r)));
                sem_row{end+1}  = nanstd(v_s)/sqrt(sum(~isnan(v_s)));
            end
        end
        xl_data(end+1,:) = mean_row; xl_data(end+1,:) = sem_row;
        xl_fname = sprintf('R2_per_animal_%s.xlsx', strrep(strrep(condition_label,' | ','_'),' ','_'));
        writetable(cell2table(xl_data,'VariableNames',xl_headers), xl_fname);
        fprintf('\nExcel saved: %s\n', xl_fname);
    catch ME_xl
        fprintf('\nCould not save Excel: %s\n', ME_xl.message);
    end

    % ---- Figure 2: CV MSE: real vs shuffled --------------------------- %
    figure('Name','Cross-validation MSE: real vs shuffled','Position',[50 680 1200 500]);
    for sp = 1:2
        subplot(1,2,sp); hold on;
        mr = res.mse_real(:,:,sp);
        ms = res.mse_shuf(:,:,sp);
        mr_mean = nanmean(mr,1); mr_sem = nanstd(mr,0,1)./sqrt(sum(~isnan(mr),1));
        ms_mean = nanmean(ms,1); ms_sem = nanstd(ms,0,1)./sqrt(sum(~isnan(ms),1));

        for ki = 1:n_models
            col_r = bar_colors(ki,:);
            col_s = col_r + (1-col_r)*shuf_alpha;
            bar(ki-bw/2, mr_mean(ki), bw, 'FaceColor', col_r, 'EdgeColor', 'none');
            bar(ki+bw/2, ms_mean(ki), bw, 'FaceColor', col_s, 'EdgeColor', col_r, ...
                'LineWidth',0.8,'LineStyle','--');
            errorbar(ki-bw/2, mr_mean(ki), mr_sem(ki), 'k.', 'LineWidth',1.2,'CapSize',4);
            errorbar(ki+bw/2, ms_mean(ki), ms_sem(ki), 'Color',col_r*0.6, ...
                     'LineStyle','none','LineWidth',1,'CapSize',4);
            valid = ~isnan(mr(:,ki));
            jitter = (rand(sum(valid),1)-0.5)*0.12;
            scatter(ki-bw/2 + jitter, mr(valid,ki), 18, col_r*0.7, 'filled','MarkerFaceAlpha',0.5);
            sl_key = sprintf('sess%d',sp-1);
            if isfield(res.group,sl_key) && isfield(res.group.(sl_key),model_keys{ki})
                p_val = res.group.(sl_key).(model_keys{ki}).p_mse;
                if     p_val < 0.001, sig = '***';
                elseif p_val < 0.01,  sig = '**';
                elseif p_val < 0.05,  sig = '*';
                else,                 sig = ''; end
                if ~isempty(sig)
                    ystar = max(mr_mean(ki)+mr_sem(ki), ms_mean(ki)+ms_sem(ki)) + 0.02;
                    text(ki,ystar,sig,'HorizontalAlignment','center', ...
                         'FontSize',11,'FontWeight','bold','Color',col_r*0.8);
                end
            end
        end
        set(gca,'XTick',1:n_models,'XTickLabel',model_labs,'XTickLabelRotation',30,'FontSize',9);
        ylabel('CV MSE  (mean \pm SEM)');
        title(sess_labels{sp},'FontSize',11);
        grid on; box off;
        ylim([0, max([mr_mean+mr_sem, ms_mean+ms_sem])*1.25]);
    end
    subplot(1,2,1);
    h1 = bar(nan,nan,'FaceColor',[0.3 0.3 0.3],'EdgeColor','none');
    h2 = bar(nan,nan,'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0.3 0.3 0.3],'LineStyle','--');
    legend([h1 h2],{'Real dF/F','Shuffled dF/F (circular)'},'Location','northeast','FontSize',8);
    sgtitle(sprintf('Cross-validation MSE  [%s]  (lower = better prediction)',condition_label),'FontSize',11);

    % ---- Figure 3: Per-animal R2 heatmap ------------------------------ %
    figure('Name','Per-animal R2 heatmap','Position',[700 50 700 500]);
    r2_allpred = squeeze(res.r2_real(:,1,:));
    imagesc(r2_allpred); colormap(flipud(hot)); cb = colorbar;
    cb.Label.String = 'Adjusted R^2';
    set(gca,'XTick',1:2,'XTickLabel',{'Sess 0','Sess 1'}, ...
            'YTick',1:n_mice,'YTickLabel',short_labels,'FontSize',9);
    xlabel('Session'); ylabel('Animal');
    title(sprintf('ALL PRED model  Adj R^2 per animal  [%s]',condition_label),'FontSize',11);
    for mi = 1:n_mice
        for si = 1:2
            v = r2_allpred(mi,si);
            if ~isnan(v)
                text(si,mi,sprintf('%.3f',v),'HorizontalAlignment','center', ...
                     'FontSize',8,'Color','w','FontWeight','bold');
            end
        end
    end

    % ---- Figure 4: R2 sess0 vs sess1 scatter -------------------------- %
    figure('Name','R2 Sess0 vs Sess1 scatter','Position',[700 620 500 450]);
    hold on;
    r2r0 = res.r2_real(:,1,1); r2r1 = res.r2_real(:,1,2);
    r2s0 = res.r2_shuf(:,1,1); r2s1 = res.r2_shuf(:,1,2);
    scatter(r2r0,r2r1,60,[0.15 0.35 0.75],'filled','DisplayName','Real');
    scatter(r2s0,r2s1,40,[0.7 0.7 0.7],'^','DisplayName','Shuffled');
    for mi = 1:n_mice
        if ~isnan(r2r0(mi)) && ~isnan(r2r1(mi))
            text(r2r0(mi)+0.005,r2r1(mi),short_labels{mi},'FontSize',7,'Color',[0.2 0.2 0.4]);
        end
    end
    ax_lim = [0, max([r2r0; r2r1; 0.1],[],'omitnan')*1.1];
    plot(ax_lim,ax_lim,'k--','LineWidth',1,'DisplayName','Unity');
    valid = ~isnan(r2r0) & ~isnan(r2r1);
    if sum(valid) >= 3
        [~,p_sess,~,st] = ttest(r2r0(valid),r2r1(valid));
        text(ax_lim(1)+0.01,ax_lim(2)*0.95, ...
             sprintf('t(%d)=%.2f, p=%.3f',st.df,st.tstat,p_sess), ...
             'FontSize',9,'Color',[0.3 0.3 0.3]);
    end
    xlim(ax_lim); ylim(ax_lim); axis square;
    xlabel('Adj R^2  Sess 0'); ylabel('Adj R^2  Sess 1');
    legend('Location','southeast','FontSize',8);
    title(sprintf('ALL PRED model  Sess0 vs Sess1  [%s]',condition_label),'FontSize',11);
    grid on; box off;

    % ---- Figure 5: Predictor contributions ---------------------------- %
    if n_preds > 0 && ~all(isnan(r2_contrib(:)))
        pred_disp = pred_names;
        pred_disp = strrep(pred_disp, 'spatial_distance',  'spatial dist');
        pred_disp = strrep(pred_disp, 'temporal_distance', 'temporal dist');
        pred_disp = strrep(pred_disp, 'time_since_food',   'tsf');
        pred_disp = strrep(pred_disp, 'food_visit',        'food visit');
        pred_disp = strrep(pred_disp, '_', ' ');

        pred_colors = repmat([0.55 0.65 0.85], n_preds, 1);
        for pi = 1:n_preds
            if strcmp(pred_names{pi},'spatial_distance'),  pred_colors(pi,:) = [0.9  0.25 0.2 ]; end
            if strcmp(pred_names{pi},'temporal_distance'), pred_colors(pi,:) = [0.2  0.65 0.85]; end
            if strcmp(pred_names{pi},'time_since_food'),   pred_colors(pi,:) = [0.2  0.7  0.4 ]; end
            if strcmp(pred_names{pi},'time'),              pred_colors(pi,:) = [0.4  0.4  0.4 ]; end
        end

        figure('Name','Predictor contributions to R2','Position',[50 50 900 480]);
        for sp = 1:2
            subplot(1,2,sp); hold on;
            contrib = r2_contrib(:,:,sp);
            c_mean  = nanmean(contrib,1);
            c_sem   = nanstd(contrib,0,1)./sqrt(sum(~isnan(contrib),1));
            [~, sort_idx] = sort(c_mean,'descend');

            for pi = 1:n_preds
                oi  = sort_idx(pi);
                col = pred_colors(oi,:);
                bar(pi, c_mean(oi), 0.6, 'FaceColor', col, 'EdgeColor', 'none');
                errorbar(pi, c_mean(oi), c_sem(oi), 'k.', 'LineWidth',1.2,'CapSize',4);
                valid = ~isnan(contrib(:,oi));
                jitter = (rand(sum(valid),1)-0.5)*0.18;
                scatter(pi+jitter, contrib(valid,oi), 18, col*0.7, 'filled','MarkerFaceAlpha',0.5);
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
            ylabel('Delta R^2  (unique contribution to ALL PRED)');
            title(sess_labels{sp},'FontSize',11);
            grid on; box off;
            ymax = max([c_mean+c_sem, 0.01],[],'omitnan');
            ymin = min([c_mean-c_sem, 0  ],[],'omitnan');
            ylim([ymin*1.2, ymax*1.35]);
        end
        sgtitle(sprintf('Unique predictor contributions to R^2  [%s]  (* p<0.05 vs 0, sorted by mean)', ...
                condition_label),'FontSize',11);

        % Save Excel for contributions
        try
            xl_headers_c = {'Animal'};
            for sp_xl = 1:2
                for pi_xl = 1:n_preds
                    xl_headers_c{end+1} = sprintf('S%d_%s_dR2', sp_xl-1, pred_disp{pi_xl});
                end
            end
            xl_data_c = {};
            for mi_xl = 1:n_mice
                row = {res.unique_mice{mi_xl}};
                for sp_xl = 1:2
                    for pi_xl = 1:n_preds
                        row{end+1} = r2_contrib(mi_xl,pi_xl,sp_xl);
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
            xl_fname_c = sprintf('Contributions_per_animal_%s.xlsx', ...
                strrep(strrep(condition_label,' | ','_'),' ','_'));
            writetable(cell2table(xl_data_c,'VariableNames',xl_headers_c), xl_fname_c);
            fprintf('\nContributions Excel saved: %s\n', xl_fname_c);
        catch ME_xl
            fprintf('\nCould not save contributions Excel: %s\n', ME_xl.message);
        end
    end

    % ---- Figure 6: Predicted vs actual traces — ALL animals ------------ %
    col_actual = [0.15 0.15 0.15];
    col_full   = [0.15 0.35 0.75];
    col_nodist = [0.95 0.55 0.10];
    col_notp   = [0.55 0.15 0.65];
    model_trace_labels = {'ALL PRED','NO SPAT DIST','NO TEMP DIST'};
    model_trace_colors = {col_full, col_nodist, col_notp};

    r2_s1 = res.r2_real(:,1,2);
    [~, sort_mice] = sort(r2_s1,'descend','MissingPlacement','last');

    trace_lookup = struct();
    for si = 1:2
        cands = example_traces(si).cands;
        for ci = 1:length(cands)
            cand = cands{ci};
            aname_safe = strrep(cand{5},'-','_');
            if ~isfield(trace_lookup, aname_safe)
                trace_lookup.(aname_safe) = cell(1,2);
            end
            trace_lookup.(aname_safe){si} = cand;
        end
    end

    for mi_plot = 1:n_mice
        m_idx    = sort_mice(mi_plot);
        mid      = res.unique_mice{m_idx};
        mid_safe = strrep(mid,'-','_');
        if ~isfield(trace_lookup, mid_safe), continue; end
        animal_cands = trace_lookup.(mid_safe);
        if isempty(animal_cands{1}) && isempty(animal_cands{2}), continue; end

        r2_s0_v = res.r2_real(m_idx,1,1);
        r2_s1_v = res.r2_real(m_idx,1,2);
        figure('Name', sprintf('Traces: %s  (S0 R2=%.3f, S1 R2=%.3f)', ...
               strrep(mid,'_',' '), r2_s0_v, r2_s1_v), ...
               'Position',[50 50 1400 560]);

        for si = 1:2
            cand = animal_cands{si};
            if isempty(cand), continue; end
            [t_sort, tidx] = sort(cand{1});
            t_sec = t_sort - t_sort(1);
            dff_r = cand{2}(tidx);
            preds = {};
            preds{1} = cand{3}(tidx);
            if length(cand) >= 6 && ~isempty(cand{6}), preds{2} = cand{6}(tidx); else, preds{2} = []; end
            if length(cand) >= 7 && ~isempty(cand{7}), preds{3} = cand{7}(tidx); else, preds{3} = []; end

            all_y = dff_r;
            for pk = 1:3
                if ~isempty(preds{pk}), all_y = [all_y; preds{pk}]; end
            end
            y_pad = 0.15 * max(range(all_y), 0.1);
            ylims = [min(all_y)-y_pad, max(all_y)+y_pad];

            for pk = 1:3
                subplot(2,3,(si-1)*3+pk); hold on;
                col = model_trace_colors{pk};
                plot(t_sec, dff_r, 'Color', col_actual, 'LineWidth', 0.6, 'DisplayName', 'Actual');
                if ~isempty(preds{pk})
                    pred_k  = preds{pk};
                    resid_k = dff_r - pred_k;
                    r2_k    = 1 - var(resid_k)/max(var(dff_r),1e-10);
                    fill([t_sec; flipud(t_sec)], ...
                         [pred_k+abs(resid_k); flipud(pred_k-abs(resid_k))], ...
                         col,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
                    plot(t_sec, pred_k, 'Color', col, 'LineWidth', 1.4, ...
                         'DisplayName', sprintf('%s  r^2=%.3f', model_trace_labels{pk}, r2_k));
                end
                yline(0,'Color',[0.75 0.75 0.75],'LineWidth',0.5,'HandleVisibility','off');
                ylim(ylims); xlim([0 max(t_sec)]);
                grid on; box off;
                legend('Location','best','FontSize',7);
                xlabel('Time (s)');
                if pk == 1
                    ylabel(sprintf('z-scored dF/F\n%s',sess_labels{si}),'FontSize',9);
                end
                if si == 1
                    title(model_trace_labels{pk},'FontSize',10,'FontWeight','bold');
                end
            end
        end
        sgtitle(sprintf('%s  [%s]\nDark=actual | Blue=ALL PRED | Orange=NO SPAT DIST | Purple=NO TEMP DIST', ...
                strrep(mid,'_',' '), condition_label),'FontSize',10);
    end

    % ---- Figure 7: Spatial distance vs Temporal distance scatter ------- %
    % For up to 6 animals (highest sess1 R2), both sessions side by side.
    % Purpose: verify the two predictors are correlated as expected and
    % diagnose pathological cases (e.g. temporal_distance all-NaN in sess0).
    % Points are colored by raw dF/F so we can see if both predictors track
    % the neural signal in the same direction.
    n_scatter = min(6, n_mice);
    scatter_mice_idx = sort_mice(1:n_scatter);

    col_sd = [0.9 0.25 0.2];
    col_td = [0.2 0.65 0.85];

    figure('Name','Spatial vs Temporal distance scatter', ...
           'Position',[50 50 1100 260*n_scatter]);

    for mi_sc = 1:n_scatter
        m_idx = scatter_mice_idx(mi_sc);
        mid   = res.unique_mice{m_idx};

        for si = 1:2
            sess_str_sc = sprintf('_sess%d', si-1);
            row_sc = find(strcmp(mice_sub(:,1), [mid, sess_str_sc]));
            if isempty(row_sc), continue; end

            data_sc = mice_sub{row_sc, 4};
            data_sc(isinf(data_sc)) = NaN;
            tv_sc = data_sc(:, COL_TIME);
            dv_sc = data_sc(:, COL_DIST);

            at_f = dv_sc <= FOOD_DIST_THRESH;
            tsf_sc = NaN(size(tv_sc)); lt = NaN;
            for tt = 1:length(tv_sc)
                if at_f(tt), lt = tv_sc(tt); end
                if ~isnan(lt), tsf_sc(tt) = tv_sc(tt) - lt; end
            end
            ttf_sc = NaN(size(tv_sc)); nt_v = NaN;
            for tt = length(tv_sc):-1:1
                if at_f(tt), nt_v = tv_sc(tt); end
                if ~isnan(nt_v), ttf_sc(tt) = nt_v - tv_sc(tt); end
            end
            tsf_m = tsf_sc; tsf_m(isnan(tsf_m)) = Inf;
            ttf_m = ttf_sc; ttf_m(isnan(ttf_m)) = Inf;
            td_sc = min(tsf_m, ttf_m);
            td_sc(isinf(td_sc)) = NaN;

            dff_sc  = data_sc(:, COL_DFF);
            valid_sc = ~isnan(dv_sc) & ~isnan(td_sc) & ~isnan(dff_sc);

            sp_idx = (mi_sc-1)*2 + si;
            ax = subplot(n_scatter, 2, sp_idx);

            if sum(valid_sc) < 10
                text(0.5,0.5, sprintf('%s Sess%d\ntemporal\_distance all-NaN', ...
                     strrep(mid,'_',' '), si-1), ...
                     'HorizontalAlignment','center','Units','normalized', ...
                     'FontSize',8,'Color',[0.5 0.5 0.5]);
                axis off; continue;
            end

            sd_v  = dv_sc(valid_sc);
            td_v  = td_sc(valid_sc);
            dff_v = dff_sc(valid_sc);

            n_pts  = min(2000, sum(valid_sc));
            rng(42 + mi_sc*10 + si);
            idx_sc = randperm(sum(valid_sc), n_pts);
            sd_p   = sd_v(idx_sc);
            td_p   = td_v(idx_sc);
            dff_p  = dff_v(idx_sc);

            r_sdtd = corr(sd_v, td_v, 'rows','complete');

            dff_z    = (dff_p - mean(dff_p,'omitnan')) / max(std(dff_p,'omitnan'),1e-6);
            clim_val = max(abs(dff_z));

            scatter(sd_p, td_p, 8, dff_z, 'filled', 'MarkerFaceAlpha', 0.5);
            colormap(ax, rdbu_colormap(256));
            if clim_val > 0, clim(ax, [-clim_val clim_val]); end
            hold on;

            % OLS fit line
            cf = polyfit(sd_v, td_v, 1);
            xl_sc = linspace(min(sd_v), max(sd_v), 100);
            plot(xl_sc, polyval(cf, xl_sc), 'k-', 'LineWidth', 1.5);

            cb = colorbar; cb.Label.String = 'dF/F'; cb.FontSize = 7;

            xlabel('Spatial distance (cm)', 'Color', col_sd, 'FontSize', 8);
            ylabel('Temporal distance (s)', 'Color', col_td, 'FontSize', 8);
            title(sprintf('%s  Sess%d  |  r(sd,td)=%.3f  n=%d', ...
                  strrep(mid,'_',' '), si-1, r_sdtd, sum(valid_sc)), ...
                  'FontSize', 8, 'FontWeight','bold');
            grid on; box off;
        end
    end
    sgtitle(sprintf('Spatial vs Temporal distance  [%s]\nColor = raw dF/F (red=high, blue=low) | r = Pearson correlation', ...
            condition_label), 'FontSize', 11);

    % ---- Figure 8: Time traces of spatial & temporal distance per mouse - %
    % One figure per mouse.  Each figure has 2 rows (sess0, sess1) x 3 cols:
    %   Col 1 : spatial_distance  (cm, raw)
    %   Col 2 : temporal_distance (s, raw)
    %   Col 3 : dF/F overlay -- both predictors z-scored on same axis so
    %            their relationship to the neural signal is visible.
    % Food-zone threshold marked as a horizontal dashed line on col 1.
    % Food visits (dist <= threshold) shown as grey shaded regions.

    col_sd  = [0.9  0.25 0.2 ];   % red   -- spatial distance
    col_td  = [0.2  0.65 0.85];   % teal  -- temporal distance
    col_dff = [0.15 0.15 0.15];   % dark  -- dF/F

    for mi_tr = 1:n_mice
        mid_tr = res.unique_mice{mi_tr};

        figure('Name', sprintf('Distance traces: %s', strrep(mid_tr,'_',' ')), ...
               'Position', [50 50 1400 520]);

        for si = 1:2
            sess_str_tr = sprintf('_sess%d', si-1);
            row_tr = find(strcmp(mice_sub(:,1), [mid_tr, sess_str_tr]));
            if isempty(row_tr), continue; end

            data_tr = mice_sub{row_tr, 4};
            data_tr(isinf(data_tr)) = NaN;
            tv_tr = data_tr(:, COL_TIME);
            dv_tr = data_tr(:, COL_DIST);   % spatial distance, raw cm
            dff_tr = data_tr(:, COL_DFF);

            % Recompute temporal_distance (NaN-safe min)
            at_f_tr = dv_tr <= FOOD_DIST_THRESH;
            tsf_tr  = NaN(size(tv_tr)); lt_tr = NaN;
            for tt = 1:length(tv_tr)
                if at_f_tr(tt), lt_tr = tv_tr(tt); end
                if ~isnan(lt_tr), tsf_tr(tt) = tv_tr(tt) - lt_tr; end
            end
            ttf_tr = NaN(size(tv_tr)); nt_tr = NaN;
            for tt = length(tv_tr):-1:1
                if at_f_tr(tt), nt_tr = tv_tr(tt); end
                if ~isnan(nt_tr), ttf_tr(tt) = nt_tr - tv_tr(tt); end
            end
            tsf_m_tr = tsf_tr; tsf_m_tr(isnan(tsf_m_tr)) = Inf;
            ttf_m_tr = ttf_tr; ttf_m_tr(isnan(ttf_m_tr)) = Inf;
            td_tr = min(tsf_m_tr, ttf_m_tr);
            td_tr(isinf(td_tr)) = NaN;

            % Time axis in seconds from session start
            valid_t = ~isnan(tv_tr);
            t_sec_tr = tv_tr - tv_tr(find(valid_t,1));

            sess_lbl = sprintf('Sess %d', si-1);

            % ---- Col 1: Spatial distance --------------------------------
            subplot(2, 3, (si-1)*3 + 1); hold on;

            % Shade food-zone visits
            in_zone = dv_tr <= FOOD_DIST_THRESH;
            % Find contiguous runs
            starts = find(diff([0; in_zone]) == 1);
            ends   = find(diff([in_zone; 0]) == -1);
            for k = 1:length(starts)
                t1 = t_sec_tr(starts(k));
                t2 = t_sec_tr(ends(k));
                patch([t1 t2 t2 t1], [0 0 FOOD_DIST_THRESH*1.5 FOOD_DIST_THRESH*1.5], ...
                      [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.5, ...
                      'HandleVisibility','off');
            end

            plot(t_sec_tr, dv_tr, 'Color', col_sd, 'LineWidth', 0.8);
            yline(FOOD_DIST_THRESH, 'k--', 'LineWidth', 0.8, ...
                  'DisplayName', sprintf('threshold (%d cm)', FOOD_DIST_THRESH));
            xlabel('Time (s)'); ylabel('Spatial distance (cm)');
            title(sprintf('%s  |  %s', sess_lbl, 'Spatial distance'), 'FontSize', 9);
            legend('FontSize', 7, 'Location', 'best');
            grid on; box off; xlim([0 max(t_sec_tr,[],'omitnan')]);

            % ---- Col 2: Temporal distance --------------------------------
            subplot(2, 3, (si-1)*3 + 2); hold on;

            % Same food-zone shading
            for k = 1:length(starts)
                t1 = t_sec_tr(starts(k));
                t2 = t_sec_tr(ends(k));
                td_max = max(td_tr,[],'omitnan');
                if isnan(td_max), td_max = 1; end
                patch([t1 t2 t2 t1], [0 0 td_max*1.1 td_max*1.1], ...
                      [0.85 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.5, ...
                      'HandleVisibility','off');
            end

            if sum(~isnan(td_tr)) > 0
                plot(t_sec_tr, td_tr, 'Color', col_td, 'LineWidth', 0.8);
                ylabel('Temporal distance (s)');
            else
                text(0.5, 0.5, 'temporal\_distance all-NaN', ...
                     'HorizontalAlignment','center','Units','normalized', ...
                     'FontSize', 9, 'Color', [0.5 0.5 0.5]);
                ylabel('Temporal distance (s)');
            end
            xlabel('Time (s)');
            title(sprintf('%s  |  %s', sess_lbl, 'Temporal distance'), 'FontSize', 9);
            grid on; box off; xlim([0 max(t_sec_tr,[],'omitnan')]);

            % ---- Col 3: Overlay -- z-scored sd, td, dF/F ---------------
            subplot(2, 3, (si-1)*3 + 3); hold on;

            % Z-score each for overlay
            zs = @(x) (x - mean(x,'omitnan')) / max(std(x,'omitnan'), 1e-6);
            sd_z  = zs(dv_tr);
            td_z  = zs(td_tr);
            dff_z = zs(dff_tr);

            plot(t_sec_tr, dff_z,  'Color', col_dff, 'LineWidth', 0.6, ...
                 'DisplayName', 'dF/F (z)');
            plot(t_sec_tr, sd_z,   'Color', col_sd,  'LineWidth', 0.9, ...
                 'DisplayName', 'Spatial dist (z)');
            if sum(~isnan(td_z)) > 0
                plot(t_sec_tr, td_z, 'Color', col_td, 'LineWidth', 0.9, ...
                     'DisplayName', 'Temporal dist (z)');
            end
            yline(0, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, ...
                  'HandleVisibility', 'off');
            xlabel('Time (s)'); ylabel('z-scored value');
            title(sprintf('%s  |  Overlay (z-scored)', sess_lbl), 'FontSize', 9);
            legend('FontSize', 7, 'Location', 'best');
            grid on; box off; xlim([0 max(t_sec_tr,[],'omitnan')]);
        end

        sgtitle(sprintf('%s  [%s]\nGrey shading = at food zone (dist <= %d cm)', ...
                strrep(mid_tr,'_',' '), condition_label, FOOD_DIST_THRESH), ...
                'FontSize', 10);
    end
end


% ====================================================================== %
%  RED-WHITE-BLUE COLORMAP (MATLAB-native replacement for RdBu)          %
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