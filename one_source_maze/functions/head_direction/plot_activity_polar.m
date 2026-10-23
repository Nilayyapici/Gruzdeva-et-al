function stats = plot_activity_polar(mice_all, opts)
% PLOT_ACTIVITY_POLAR  Plot mean dF/F and occupancy as a function of an
%   angle column (e.g. column 16 = head-direction-to-food, or column 17 =
%   head-direction-to-shortest-path) as polar "tuning curves", for three
%   spatial conditions, and TEST whether activity depends on direction:
%       1) Full maze           all frames
%       2) Near food           distance to food <  opts.dist_threshold
%       3) Far from food       distance to food >= opts.dist_threshold
%
%   stats = PLOT_ACTIVITY_POLAR(mice_all)
%   stats = PLOT_ACTIVITY_POLAR(mice_all, opts)
%
%   FIGURES (per angle column)
%     Comparison figure, one column per condition:
%       row 1 : group-average dF/F tuning curve (mean +- SEM across mice,
%               blue) with the group harmonic fit overlaid (red), SAME
%               radial scale in all three panels; harmonic p-values in
%               the title
%       row 2 : occupancy (% of that condition's time per angle bin,
%               averaged across mice), SAME radial scale
%       row 3 : sector test -- mean dF/F in front / +90 side / back /
%               -90 side (or 3 sectors with opts.merge_sides = true), one line
%               per mouse, Friedman test, stars vs chance, brackets for
%               sector-vs-sector differences
%     Individual mice (opts.individual): one figure per condition.
%
%   STATISTICS (the mouse is the unit -- frames are never pooled across
%   mice for a test, since frames within a mouse are autocorrelated)
%
%   1) HARMONIC REGRESSION (main test). For each mouse and condition:
%        dF/F ~ c0 + a1*cos(th) + b1*sin(th) + a2*cos(2th) + b2*sin(2th)
%      fitted on every frame (no binning).
%        1st harmonic (a1,b1): ONE preferred direction.
%            a1 > 0  -> higher when facing the target (0 deg)
%            b1 != 0 -> left/right asymmetry
%            preferred direction = atan2(b1, a1)
%        2nd harmonic (a2,b2): TWO opposite lobes, e.g. higher when
%            facing sideways on either side vs facing toward/away.
%            preferred axis = atan2(b2, a2)/2 (and that + 180 deg).
%      Group test per harmonic: one-sample Hotelling's T^2 on the per-mouse
%      (a,b) vectors -- "do the mice share a consistent directional
%      tuning?" Plus a t-test on each coefficient separately.
%      Mice need >= opts.min_frames_for_fit frames AND at least
%      opts.min_coverage of the angular bins sampled, otherwise their fit
%      would be extrapolated from a narrow range of directions.
%
%   2) SECTOR TEST (easy to read). Mean dF/F per mouse in fixed sectors:
%        opts.merge_sides = false (default): 4 sectors
%          Front |th| <= 45, +90 side 45 < th <= 135, Back |th| > 135,
%          -90 side -135 <= th < -45
%        opts.merge_sides = true: 3 sectors
%          Front, Sides 45 < |th| <= 135 (both sides pooled), Back
%      Friedman test across mice (mice with all sectors sampled), then
%      pairwise comparisons between sectors (opts.sector_posthoc, default
%      Wilcoxon signed-rank, paired across mice), Holm-corrected; the
%      significant pairs are drawn as brackets on the sector panel.
%      Sectors are labeled by SIGN (+90 / -90), not left/right, so the
%      labels can't be wrong if the sign convention is mixed up.
%      EACH SECTOR vs CHANCE (stars on the bars): chance is estimated per
%      mouse by circularly shifting dF/F in time relative to the angle
%      within each session (random offset 10-90% of the session length,
%      opts.n_shuffle_sector times). This keeps occupancy and the
%      autocorrelation of the signal and destroys only their alignment.
%      Per mouse and sector: deviation = real mean - mean shuffled mean.
%      Deviations are tested vs 0 across mice with a one-sample t-test
%      (opts.sector_chance_test; 'signrank' also available), WITHOUT
%      multiple-comparison correction by default
%      (opts.sector_chance_correction = 'none'; 'holm' to correct over the
%      sectors). Stars above a bar = higher than chance, below a bar =
%      lower than chance; dashed line = chance.
%      NOTE: the shift also breaks the relation of dF/F with distance and
%      approach/departure, so "above chance" means activity is related to
%      facing direction in that sector, not that direction is the cause.
%
%   3) NEAR vs FAR (paired, mice with valid fits in both): Hotelling's T^2
%      on the per-mouse difference of (a,b) per harmonic -- "is the tuning
%      different near vs far from food?" -- plus a paired t-test on the
%      modulation amplitude sqrt(a^2+b^2).
%
%   DATA ASSUMED (mice_all{row,4}):
%     5-distance to food   9-eating   10-grooming   11-dff
%     16 or 17 - an angle column (degrees, signed, (-180,180])
%   ROW NAMING ASSUMED: '<mouseID>_..._sess<N>', e.g.
%     'F10I_GLM_fasted_food_sess0' -> mouse 'F10I', session number 0
%
%   OPTS (all optional):
%     opts.sessions   : session numbers to pool per mouse. (default: 0:3)
%     opts.state      : 'fasted' | 'fed'. (default: 'fasted')
%     opts.source     : 'food' | 'gel'. (default: 'food')
%     opts.angle_col  : 16, 17, or [16 17]. (default: 16)
%     opts.dff_col    : activity column. (default: 11)
%     opts.n_bins     : angular bins for plotting. (default: 24 = 15 deg)
%     opts.zscore_dff : z-score dff per session on all kept frames,
%                       BEFORE the near/far split. (default: true)
%     opts.min_frames_per_bin : minimum frames for a bin to be plotted /
%                       counted as sampled. (default: 5)
%     opts.exclude_eating   : drop eating frames. (default: false)
%     opts.exclude_grooming : drop grooming frames. (default: false)
%     opts.speed_threshold  : drop frames slower than this. (default: [])
%     opts.speed_col        : speed column. (default: 4)
%     opts.dist_threshold   : near/far split distance. (default: 10)
%     opts.dist_col         : distance-to-food column. (default: 5)
%     opts.dist_remove      : remove frames where distance to food is
%                       BELOW this value, from ALL conditions, before
%                       z-scoring and before any statistics (e.g. to drop
%                       frames where the mouse is so close to the food
%                       that the angle to it is unreliable). [] = off.
%                       With it set, "Near food" becomes
%                       dist_remove <= dist < dist_threshold. (default: [])
%     opts.theta_dir        : 'counterclockwise' | 'clockwise'.
%                       (default: 'counterclockwise')
%     opts.weight_by_occupancy : occupancy-weight the plotted group curve.
%                       (default: true)
%     opts.individual       : per-mouse figures. (default: true)
%     opts.n_harmonics      : 1 or 2 harmonics in the regression.
%                       (default: 2)
%     opts.min_frames_for_fit : (default: 50)
%     opts.min_coverage     : fraction of the n_bins angular bins that
%                       must be sampled for a mouse's fit. (default: 0.5)
%     opts.min_frames_per_sector : (default: 10)
%     opts.merge_sides      : pool +90 and -90 sides into one "Sides"
%                       sector. (default: false = keep +90 and -90 separate)
%     opts.sector_posthoc   : 'signrank' | 'ttest' for the pairwise
%                       sector comparisons. (default: 'signrank')
%     opts.sector_ylim      : [ymin ymax] y-limits for the sector bar
%                       plots (same for all three conditions); [] = auto,
%                       with room above the bars for stars and brackets.
%                       Leave space above the bars yourself when setting
%                       it, otherwise brackets can be cut off. (default: [])
%     opts.n_shuffle_sector : number of circular shifts for the sector
%                       chance level; 0 = off. (default: 200)
%     opts.sector_chance_test : 'ttest' | 'signrank' for each sector vs
%                       chance. (default: 'ttest')
%     opts.sector_chance_correction : 'none' | 'holm' for the sector-vs-
%                       chance p-values. (default: 'none')
%     opts.shuffle_seed     : RNG seed for reproducible shuffles; [] =
%                       unseeded. (default: 0)
%     opts.export           : save the plotted data and statistics to an
%                       Excel file, one per angle column. (default: true)
%     opts.export_dir       : folder for the Excel files. (default: pwd)
%
%   EXCEL SHEETS (opts.export)
%     polar_group        : per angle bin and condition: group mean, SEM,
%                          n mice, occupancy % (rows 1-2 of the figure)
%     harmonic_fit_curve : the red group fit curve per condition
%     polar_per_mouse    : per mouse, condition and angle bin: mean dF/F
%                          and frame count (individual-mouse figures)
%     harmonic_per_mouse : per mouse and condition: fit coefficients
%     harmonic_stats     : Hotelling / coefficient tests, and near vs far
%     sectors_per_mouse  : per mouse, condition and sector: real, chance,
%                          real - chance, included in sector tests
%     sector_stats       : per condition and sector: plotted bar (mean,
%                          SEM), chance line, tests vs chance
%     sector_pairs       : sector-vs-sector tests (brackets)
%     friedman           : Friedman test per condition
%
%   OUTPUT
%     stats.col16 (etc.) with .full, .near, .far, each containing
%       .harmonic(k) : .a, .b (per mouse), .c0, .valid, .n, .mean_a,
%                      .mean_b, .amp, .pref_deg, .p_hotelling, .T2,
%                      .p_a, .p_b
%       .sector      : .names, .means (mice x 4), .complete, .n,
%                      .p_friedman, .pairs, .p_pairs, .p_pairs_holm,
%                      .chance (mice x 4), .dev (real - chance),
%                      .p_vs_chance, .p_vs_chance_holm,
%                      .p_vs_chance_used (the one the stars use),
%                      .chance_test, .chance_correction
%     and .near_vs_far.harmonic(k) : .n, .p_hotelling, .mean_amp_near,
%                      .mean_amp_far, .p_amp
%     and .excel_filename (if opts.export)
%
%   EXAMPLE:
%     stats = plot_activity_polar(mice_all, struct('angle_col', [16 17], ...
%         'dist_threshold', 10, 'exclude_eating', true));

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'sessions'),            opts.sessions            = 0:3;      end
if ~isfield(opts,'state'),               opts.state               = 'fasted'; end
if ~isfield(opts,'source'),              opts.source              = 'food';   end
if ~isfield(opts,'angle_col'),           opts.angle_col           = 16;       end
if ~isfield(opts,'dff_col'),             opts.dff_col             = 11;       end
if ~isfield(opts,'n_bins'),              opts.n_bins              = 24;       end
if ~isfield(opts,'zscore_dff'),          opts.zscore_dff          = true;     end
if ~isfield(opts,'min_frames_per_bin'),  opts.min_frames_per_bin  = 5;        end
if ~isfield(opts,'exclude_eating'),      opts.exclude_eating      = false;    end
if ~isfield(opts,'exclude_grooming'),    opts.exclude_grooming    = false;    end
if ~isfield(opts,'speed_threshold'),     opts.speed_threshold     = [];       end
if ~isfield(opts,'speed_col'),           opts.speed_col           = 4;        end
if ~isfield(opts,'dist_threshold') || isempty(opts.dist_threshold)
    opts.dist_threshold = 10;
end
if ~isfield(opts,'dist_col'),            opts.dist_col            = 5;        end
if ~isfield(opts,'dist_remove'),         opts.dist_remove         = [];       end
if ~isempty(opts.dist_remove) && opts.dist_remove >= opts.dist_threshold
    warning(['dist_remove (%g) >= dist_threshold (%g): the "Near food" ', ...
             'condition will be empty.'], opts.dist_remove, opts.dist_threshold);
end
if ~isfield(opts,'theta_dir'),           opts.theta_dir           = 'counterclockwise'; end
if ~isfield(opts,'weight_by_occupancy'), opts.weight_by_occupancy = true;     end
if ~isfield(opts,'individual'),          opts.individual          = true;     end
if ~isfield(opts,'n_harmonics'),         opts.n_harmonics         = 2;        end
if ~isfield(opts,'min_frames_for_fit'),  opts.min_frames_for_fit  = 50;       end
if ~isfield(opts,'min_coverage'),        opts.min_coverage        = 0.5;      end
if ~isfield(opts,'min_frames_per_sector'), opts.min_frames_per_sector = 10;   end
if ~isfield(opts,'n_shuffle_sector'),    opts.n_shuffle_sector    = 200;      end
if ~isfield(opts,'merge_sides'),         opts.merge_sides         = false;    end
if ~isfield(opts,'sector_posthoc'),      opts.sector_posthoc      = 'signrank'; end
if ~isfield(opts,'sector_ylim'),         opts.sector_ylim         = [];       end
if ~isempty(opts.sector_ylim) && ...
        (numel(opts.sector_ylim) ~= 2 || opts.sector_ylim(2) <= opts.sector_ylim(1))
    error('opts.sector_ylim must be [ymin ymax] with ymax > ymin.');
end
if ~isfield(opts,'sector_chance_test'),  opts.sector_chance_test  = 'ttest';  end
if ~isfield(opts,'sector_chance_correction'), opts.sector_chance_correction = 'none'; end
if ~isfield(opts,'shuffle_seed'),        opts.shuffle_seed        = 0;        end
if ~isfield(opts,'export'),              opts.export              = true;     end
if ~isfield(opts,'export_dir'),          opts.export_dir          = pwd;      end
if ~isempty(opts.shuffle_seed), rng(opts.shuffle_seed); end
opts.n_harmonics = max(1, min(2, round(opts.n_harmonics)));
if isfield(opts,'show_occupancy') || isfield(opts,'average')
    fprintf(['NOTE: opts.show_occupancy / opts.average are no longer used -- ', ...
             'the comparison figure always shows dF/F and occupancy.\n']);
end

n_rows = size(mice_all,1);

% ------------------------------------------------------------------ %
%  PARSE SESSION NUMBER + MOUSE ID FROM EACH ROW'S NAME                %
% ------------------------------------------------------------------ %
sess_num = NaN(n_rows,1);
mouse_id = cell(n_rows,1);
for r = 1:n_rows
    name = mice_all{r,1};
    tok  = regexp(name, '_sess(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        sess_num(r) = str2double(tok{1});
    end
    mid = regexp(name, '^[^_]+', 'match', 'once');
    if isempty(mid), mid = name; end
    mouse_id{r} = mid;
end

% ------------------------------------------------------------------ %
%  FILTER ROWS BY state / source / sessions                            %
% ------------------------------------------------------------------ %
keep = strcmp(mice_all(:,2), opts.state) & strcmp(mice_all(:,3), opts.source) & ...
       ismember(sess_num, opts.sessions);
rows = find(keep);
if isempty(rows)
    error('No sessions match state="%s", source="%s", sessions=[%s].', ...
        opts.state, opts.source, num2str(opts.sessions));
end

fprintf('Found %d matching sessions across %d unique mice [%s | %s | sessions %s].\n', ...
    length(rows), length(unique(mouse_id(rows))), opts.state, opts.source, num2str(opts.sessions));

stats = struct();
for ac = 1:length(opts.angle_col)
    stats.(sprintf('col%d', opts.angle_col(ac))) = ...
        process_one_angle(mice_all, rows, mouse_id, opts.angle_col(ac), opts);
end

end


% ====================================================================== %
%  BUILD + TEST + PLOT FOR ONE ANGLE COLUMN                             %
% ====================================================================== %
function out = process_one_angle(mice_all, rows, mouse_id, angle_col, opts)

edges       = linspace(-180, 180, opts.n_bins+1);
centers_deg = edges(1:end-1) + diff(edges)/2;
centers_rad = deg2rad(centers_deg);

COL_EATING   = 9;
COL_GROOMING = 10;
thr = opts.dist_threshold;
K   = opts.n_harmonics;

cond_keys  = {'full', 'near', 'far'};
if isempty(opts.dist_remove)
    cond_names = {'Full maze', ...
                  sprintf('Near food (dist < %g)', thr), ...
                  sprintf('Far from food (dist >= %g)', thr)};
else
    dr = opts.dist_remove;
    cond_names = {sprintf('Full maze (dist >= %g)', dr), ...
                  sprintf('Near food (%g <= dist < %g)', dr, thr), ...
                  sprintf('Far from food (dist >= %g)', thr)};
end
n_cond = numel(cond_keys);

if opts.merge_sides
    sector_names = {'Front', 'Sides', 'Back'};
else
    sector_names = {'Front', '+90 side', 'Back', '-90 side'};
end
n_sec = numel(sector_names);

mice_this = unique(mouse_id(rows), 'stable');
n_mice    = length(mice_this);

curves = cell(1, n_cond); occ = cell(1, n_cond); n_used = cell(1, n_cond);
coef   = cell(1, n_cond);   % n_mice x (1 + 2K): [c0 a1 b1 a2 b2]
fit_ok = cell(1, n_cond);   % n_mice x 1 logical
sec_m  = cell(1, n_cond);   % n_mice x 4 sector means
sec_ch = cell(1, n_cond);   % n_mice x 4 chance (shuffled) sector means
for c = 1:n_cond
    curves{c} = NaN(n_mice, opts.n_bins);
    occ{c}    = zeros(n_mice, opts.n_bins);
    n_used{c} = zeros(n_mice, 1);
    coef{c}   = NaN(n_mice, 1 + 2*K);
    fit_ok{c} = false(n_mice, 1);
    sec_m{c}  = NaN(n_mice, n_sec);
    sec_ch{c} = NaN(n_mice, n_sec);   % shuffle chance level per sector
end

for m = 1:n_mice
    mid        = mice_this{m};
    mouse_rows = rows(strcmp(mouse_id(rows), mid));

    angle_all = []; dff_all = []; dist_all = []; sess_all = []; k_sess = 0;
    for ri = mouse_rows(:)'
        data = mice_all{ri,4};
        if isempty(data) || size(data,2) < max([angle_col, opts.dff_col, opts.dist_col])
            continue;
        end
        ang  = data(:, angle_col);
        dff  = data(:, opts.dff_col);
        dist = data(:, opts.dist_col);
        dist(isinf(dist)) = NaN;

        % behavioural exclusions (apply to all three conditions)
        exclude = false(size(ang));
        if opts.exclude_eating && size(data,2) >= COL_EATING
            exclude = exclude | (data(:,COL_EATING) > 0);
        end
        if opts.exclude_grooming && size(data,2) >= COL_GROOMING
            exclude = exclude | (data(:,COL_GROOMING) > 0);
        end
        if ~isempty(opts.speed_threshold) && size(data,2) >= opts.speed_col
            exclude = exclude | (data(:,opts.speed_col) < opts.speed_threshold);
        end
        if ~isempty(opts.dist_remove)
            exclude = exclude | (dist < opts.dist_remove);   % too close to food
        end
        ang(exclude) = NaN;
        dff(exclude) = NaN;

        % z-score on all kept frames, BEFORE the near/far split
        if opts.zscore_dff
            s = std(dff, 'omitnan');
            if s > 0
                dff = (dff - mean(dff,'omitnan')) / s;
            end
        end
        k_sess    = k_sess + 1;
        angle_all = [angle_all; ang];  %#ok<AGROW>
        dff_all   = [dff_all;   dff];  %#ok<AGROW>
        dist_all  = [dist_all;  dist]; %#ok<AGROW>
        sess_all  = [sess_all;  k_sess * ones(size(ang))]; %#ok<AGROW>
    end

    if isempty(angle_all), continue; end
    valid  = ~isnan(angle_all) & ~isnan(dff_all);
    region = {valid, ...
              valid & dist_all <  thr, ...   % NaN distance -> false
              valid & dist_all >= thr};

    bin_idx = discretize(angle_all, edges);
    sec_idx = sector_index(angle_all, opts.merge_sides);

    for c = 1:n_cond
        sel = region{c};
        n_used{c}(m) = sum(sel);

        % binned curve (for plotting)
        for b = 1:opts.n_bins
            in_bin = sel & bin_idx == b;
            occ{c}(m,b) = sum(in_bin);
            if occ{c}(m,b) >= opts.min_frames_per_bin
                curves{c}(m,b) = mean(dff_all(in_bin));
            end
        end

        % harmonic regression on all frames of this condition
        coverage = mean(occ{c}(m,:) >= opts.min_frames_per_bin);
        if n_used{c}(m) >= opts.min_frames_for_fit && coverage >= opts.min_coverage
            th = deg2rad(angle_all(sel));
            coef{c}(m,:) = fit_harmonics(th, dff_all(sel), K);
            fit_ok{c}(m) = true;
        end

        % sector means
        for s = 1:n_sec
            in_s = sel & sec_idx == s;
            if sum(in_s) >= opts.min_frames_per_sector
                sec_m{c}(m,s) = mean(dff_all(in_s));
            end
        end
    end

    % ---- sector chance level: circular shift of dF/F vs angle ----------
    if opts.n_shuffle_sector > 0
        acc = zeros(n_cond, n_sec);
        cnt = zeros(n_cond, n_sec);
        sess_ids = unique(sess_all)';
        for i_sh = 1:opts.n_shuffle_sector
            dff_sh = dff_all;
            for sid = sess_ids
                ii = find(sess_all == sid);
                n  = numel(ii);
                if n < 10, continue; end
                lo_sh = max(1, round(0.1 * n));
                hi_sh = max(lo_sh, round(0.9 * n));
                dff_sh(ii) = circshift(dff_all(ii), randi([lo_sh, hi_sh]));
            end
            ok_sh  = ~isnan(angle_all) & ~isnan(dff_sh);
            reg_sh = {ok_sh, ok_sh & dist_all < thr, ok_sh & dist_all >= thr};
            for c = 1:n_cond
                for s = 1:n_sec
                    if isnan(sec_m{c}(m,s)), continue; end
                    in_s = reg_sh{c} & sec_idx == s;
                    if any(in_s)
                        acc(c,s) = acc(c,s) + mean(dff_sh(in_s));
                        cnt(c,s) = cnt(c,s) + 1;
                    end
                end
            end
        end
        for c = 1:n_cond
            ch = acc(c,:) ./ cnt(c,:);
            ch(cnt(c,:) == 0 | isnan(sec_m{c}(m,:))) = NaN;
            sec_ch{c}(m,:) = ch;
        end
    end
end

col_label = sprintf('column %d', angle_col);
if angle_col == 16, col_label = 'head direction to food'; end
if angle_col == 17, col_label = 'head direction to shortest path'; end
cond_label = sprintf('%s | %s | sessions %s', opts.state, opts.source, num2str(opts.sessions));
if ~isempty(opts.dist_remove)
    cond_label = sprintf('%s | frames with dist < %g removed', cond_label, opts.dist_remove);
end

fprintf('\n=== %s  [%s] ===\n', col_label, cond_label);
fprintf('Frames used per condition (total across mice):\n');
for c = 1:n_cond
    fprintf('  %-28s %8d frames  (%d mice with data, %d with a valid harmonic fit)\n', ...
        cond_names{c}, sum(n_used{c}), sum(n_used{c} > 0), sum(fit_ok{c}));
end
n_nan_dist = sum(n_used{1}) - sum(n_used{2}) - sum(n_used{3});
if n_nan_dist > 0
    fprintf('  (%d frames with missing distance are in "Full maze" only)\n', n_nan_dist);
end

% ------------------------------------------------------------------ %
%  STATISTICS                                                          %
% ------------------------------------------------------------------ %
out = struct();
for c = 1:n_cond
    res = struct();

    % ---- harmonic tests ----
    for k = 1:K
        A  = coef{c}(:, 2*k);
        B  = coef{c}(:, 2*k+1);
        ok = fit_ok{c};
        h  = struct('a', A, 'b', B, 'c0', coef{c}(:,1), 'valid', ok, 'n', sum(ok), ...
                    'mean_a', NaN, 'mean_b', NaN, 'amp', NaN, 'pref_deg', NaN, ...
                    'p_hotelling', NaN, 'T2', NaN, 'p_a', NaN, 'p_b', NaN);
        if h.n >= 1
            h.mean_a = mean(A(ok));
            h.mean_b = mean(B(ok));
            h.amp    = hypot(h.mean_a, h.mean_b);
            h.pref_deg = atan2d(h.mean_b, h.mean_a) / k;   % k = 2 -> axis
        end
        if h.n >= 3
            [h.p_hotelling, h.T2] = hotelling_one_sample([A(ok), B(ok)]);
        end
        if h.n >= 2
            [~, h.p_a] = ttest(A(ok));
            [~, h.p_b] = ttest(B(ok));
        end
        res.harmonic(k) = h;
    end

    % ---- sector test ----
    S        = sec_m{c};
    complete = all(~isnan(S), 2);
    sec = struct('names', {sector_names}, 'means', S, 'complete', complete, ...
                 'n', sum(complete), 'p_friedman', NaN, 'pairs', [], ...
                 'p_pairs', [], 'p_pairs_holm', [], 'posthoc', opts.sector_posthoc);
    if sec.n >= 3
        try
            sec.p_friedman = friedman(S(complete,:), 1, 'off');
        catch
        end
        pairs = nchoosek(1:n_sec, 2);
        p_t   = NaN(size(pairs,1), 1);
        p_sr  = NaN(size(pairs,1), 1);
        for pp = 1:size(pairs,1)
            xa = S(complete, pairs(pp,1));
            xb = S(complete, pairs(pp,2));
            [~, p_t(pp)] = ttest(xa, xb);
            try, p_sr(pp) = signrank(xa, xb); catch, end
        end
        if strcmpi(opts.sector_posthoc, 'ttest'), p_pairs = p_t; else, p_pairs = p_sr; end
        sec.pairs          = pairs;
        sec.pair_names     = arrayfun(@(i) sprintf('%s vs %s', sector_names{pairs(i,1)}, ...
                                 sector_names{pairs(i,2)}), (1:size(pairs,1))', 'UniformOutput', false);
        sec.pair_mean_diff = mean(S(complete, pairs(:,1)) - S(complete, pairs(:,2)), 1)';
        sec.posthoc        = opts.sector_posthoc;
        sec.p_pairs        = p_pairs;
        sec.p_pairs_ttest  = p_t;
        sec.p_pairs_signrank = p_sr;
        sec.p_pairs_holm   = holm(p_pairs);
    end

    % ---- each sector vs chance (shuffle) ----
    sec.chance           = sec_ch{c};
    sec.dev              = S - sec_ch{c};
    sec.chance_test      = opts.sector_chance_test;
    sec.p_vs_chance      = NaN(1, n_sec);
    sec.p_vs_chance_t    = NaN(1, n_sec);
    sec.p_vs_chance_sr   = NaN(1, n_sec);
    sec.p_vs_chance_holm = NaN(1, n_sec);
    sec.p_vs_chance_used = NaN(1, n_sec);
    sec.chance_correction = opts.sector_chance_correction;
    if opts.n_shuffle_sector > 0 && sec.n >= 3
        for s = 1:n_sec
            d = sec.dev(complete, s);
            d = d(~isnan(d));
            if numel(d) < 3, continue; end
            [~, sec.p_vs_chance_t(s)] = ttest(d);
            try, sec.p_vs_chance_sr(s) = signrank(d); catch, end
        end
        if strcmpi(opts.sector_chance_test, 'ttest')
            sec.p_vs_chance = sec.p_vs_chance_t;
        else
            sec.p_vs_chance = sec.p_vs_chance_sr;
        end
        sec.p_vs_chance_holm = holm(sec.p_vs_chance(:))';
        if strcmpi(opts.sector_chance_correction, 'holm')
            sec.p_vs_chance_used = sec.p_vs_chance_holm;
        else
            sec.p_vs_chance_used = sec.p_vs_chance;
        end
    end
    res.sector = sec;

    out.(cond_keys{c}) = res;
end

% ---- near vs far, paired ----
both = fit_ok{2} & fit_ok{3};
for k = 1:K
    nv = struct('n', sum(both), 'p_hotelling', NaN, 'T2', NaN, ...
                'mean_amp_near', NaN, 'mean_amp_far', NaN, 'p_amp', NaN);
    if nv.n >= 1
        dA = coef{2}(both, 2*k)   - coef{3}(both, 2*k);
        dB = coef{2}(both, 2*k+1) - coef{3}(both, 2*k+1);
        amp_n = hypot(coef{2}(both,2*k), coef{2}(both,2*k+1));
        amp_f = hypot(coef{3}(both,2*k), coef{3}(both,2*k+1));
        nv.mean_amp_near = mean(amp_n);
        nv.mean_amp_far  = mean(amp_f);
        if nv.n >= 3
            [nv.p_hotelling, nv.T2] = hotelling_one_sample([dA, dB]);
        end
        if nv.n >= 2
            [~, nv.p_amp] = ttest(amp_n, amp_f);
        end
    end
    out.near_vs_far.harmonic(k) = nv;
end

% ---- console report ----
harm_name = {'1st harmonic (one preferred direction)', ...
             '2nd harmonic (two opposite lobes)'};
for c = 1:n_cond
    res = out.(cond_keys{c});
    fprintf('\n--- %s ---\n', cond_names{c});
    for k = 1:K
        h = res.harmonic(k);
        if k == 1
            dir_txt = sprintf('preferred direction = %+.0f deg', h.pref_deg);
        else
            dir_txt = sprintf('preferred axis = %+.0f / %+.0f deg', h.pref_deg, wrap180(h.pref_deg + 180));
        end
        fprintf('  %-40s n=%2d  Hotelling p = %.4g %s   (a: mean %+.3f, p = %.3g; b: mean %+.3f, p = %.3g)\n', ...
            harm_name{k}, h.n, h.p_hotelling, sig_stars(h.p_hotelling), ...
            h.mean_a, h.p_a, h.mean_b, h.p_b);
        fprintf('  %-40s amplitude = %.3f, %s\n', '', h.amp, dir_txt);
    end
    s = res.sector;
    fprintf('  Sector test (n=%d mice with all 4 sectors): Friedman p = %.4g %s\n', ...
        s.n, s.p_friedman, sig_stars(s.p_friedman));
    if opts.n_shuffle_sector > 0 && ~all(isnan(s.p_vs_chance))
        fprintf('  Each sector vs chance (%d circular shifts, %s, correction: %s; stars use %s):\n', ...
            opts.n_shuffle_sector, s.chance_test, s.chance_correction, ...
            ternary(strcmpi(s.chance_correction,'holm'), 'Holm p', 'uncorrected p'));
        for ss = 1:n_sec
            dv = mean(s.dev(s.complete, ss), 'omitnan');
            if dv > 0, dtxt = 'higher'; else, dtxt = 'lower'; end
            fprintf('      %-9s real %+.3f  chance %+.3f  diff %+.3f (%s)  p = %.4g  Holm p = %.4g %s\n', ...
                sector_names{ss}, mean(s.means(s.complete, ss), 'omitnan'), ...
                mean(s.chance(s.complete, ss), 'omitnan'), dv, dtxt, ...
                s.p_vs_chance(ss), s.p_vs_chance_holm(ss), sig_stars(s.p_vs_chance_used(ss)));
        end
    end
    if ~isempty(s.pairs)
        fprintf('  Pairwise sector comparisons (%s, paired across mice, Holm over %d pairs):\n', ...
            s.posthoc, size(s.pairs,1));
        for pp = 1:size(s.pairs,1)
            fprintf('      %-9s vs %-9s  diff %+.3f  p = %.4g  Holm p = %.4g %s\n', ...
                sector_names{s.pairs(pp,1)}, sector_names{s.pairs(pp,2)}, ...
                s.pair_mean_diff(pp), s.p_pairs(pp), s.p_pairs_holm(pp), sig_stars(s.p_pairs_holm(pp)));
        end
    end
end
fprintf('\n--- Near vs far (paired across mice with valid fits in both) ---\n');
for k = 1:K
    nv = out.near_vs_far.harmonic(k);
    fprintf('  %-40s n=%2d  tuning differs: Hotelling p = %.4g %s   amplitude near %.3f vs far %.3f, p = %.4g\n', ...
        harm_name{k}, nv.n, nv.p_hotelling, sig_stars(nv.p_hotelling), ...
        nv.mean_amp_near, nv.mean_amp_far, nv.p_amp);
end

% ------------------------------------------------------------------ %
%  GROUP CURVES, FITS, OCCUPANCY                                       %
% ------------------------------------------------------------------ %
theta_fit = linspace(-pi, pi, 181);
mean_c = cell(1, n_cond); sem_c = cell(1, n_cond);
occ_pct = cell(1, n_cond); n_mice_c = zeros(1, n_cond);
fit_c = cell(1, n_cond);
for c = 1:n_cond
    [mean_c{c}, sem_c{c}] = group_curve(curves{c}, occ{c}, opts.weight_by_occupancy);
    tot = sum(occ{c}, 2);
    has = tot > 0;
    occ_pct{c} = NaN(1, opts.n_bins);
    if any(has)
        occ_pct{c} = mean(occ{c}(has,:) ./ tot(has), 1) * 100;
    end
    n_mice_c(c) = sum(has);

    % group fit = curve from the mean coefficients of mice with a valid fit
    fit_c{c} = [];
    if any(fit_ok{c})
        mc = mean(coef{c}(fit_ok{c}, :), 1);
        f  = mc(1) * ones(size(theta_fit));
        for k = 1:K
            f = f + mc(2*k) * cos(k*theta_fit) + mc(2*k+1) * sin(k*theta_fit);
        end
        fit_c{c} = f;
    end
end

% Shared radial scales (include the fitted curves)
lo = Inf; hi = -Inf;
for c = 1:n_cond
    s = sem_c{c}; s(isnan(s)) = 0;
    lo = min([lo, min(mean_c{c} - s, [], 'omitnan')]);
    hi = max([hi, max(mean_c{c} + s, [], 'omitnan')]);
    if ~isempty(fit_c{c})
        lo = min(lo, min(fit_c{c})); hi = max(hi, max(fit_c{c}));
    end
end
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = 0; hi = 1;
end
pad = 0.05 * (hi - lo);
dff_rlim = [lo - pad, hi + pad];
occ_max  = max(cellfun(@(v) max([v, 0], [], 'omitnan'), occ_pct));
if occ_max <= 0, occ_max = 1; end

% shared y for sector panels
all_sec = [];
for c = 1:n_cond, all_sec = [all_sec; sec_m{c}(:); sec_ch{c}(:)]; end %#ok<AGROW>
all_sec = all_sec(~isnan(all_sec));
if isempty(all_sec), all_sec = [-1 1]; end
sec_ylim = [min([all_sec; 0]), max([all_sec; 0])];
sec_rng  = max(diff(sec_ylim), 0.1);
n_pairs_max = nchoosek(n_sec, 2);
sec_ylim = [sec_ylim(1) - 0.18*sec_rng, ...
            sec_ylim(2) + (0.30 + 0.12*n_pairs_max) * sec_rng];   % room for stars + brackets
if ~isempty(opts.sector_ylim)
    sec_ylim = opts.sector_ylim(:)';   % manual limits
end

% ------------------------------------------------------------------ %
%  COMPARISON FIGURE                                                   %
% ------------------------------------------------------------------ %
weight_txt = '';
if opts.weight_by_occupancy, weight_txt = ', occ.-weighted'; end
if opts.zscore_dff, dff_unit = 'z-scored dF/F'; else, dff_unit = 'dF/F'; end

figure('Name', sprintf('Polar tuning (%s) -- full / near / far', col_label), ...
    'Position', [50 30 1400 1150]);
for c = 1:n_cond
    res = out.(cond_keys{c});

    % row 1: tuning curve + harmonic fit
    subplot(3, n_cond, c);
    h_txt = sprintf('H1 p = %.3g %s', res.harmonic(1).p_hotelling, sig_stars(res.harmonic(1).p_hotelling));
    if K >= 2
        h_txt = sprintf('%s | H2 p = %.3g %s', h_txt, res.harmonic(2).p_hotelling, ...
            sig_stars(res.harmonic(2).p_hotelling));
    end
    plot_one_polar_curve(centers_rad, mean_c{c}, sem_c{c}, ...
        sprintf('%s\n%s, n = %d mice%s\n%s (n fit = %d)', cond_names{c}, dff_unit, ...
        n_mice_c(c), weight_txt, h_txt, res.harmonic(1).n), ...
        opts.theta_dir, dff_rlim, theta_fit, fit_c{c});

    % row 2: occupancy
    subplot(3, n_cond, n_cond + c);
    plot_occupancy(centers_rad, occ_pct{c}, ...
        sprintf('Occupancy (%% of time, mean across mice)'), opts.theta_dir, occ_max);

    % row 3: sector test
    subplot(3, n_cond, 2*n_cond + c); hold on;
    S  = res.sector.means;
    cm = res.sector.complete;
    if any(cm)
        plot(1:n_sec, S(cm,:)', '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.7);
    end
    mu = mean(S(cm,:), 1);
    se = std(S(cm,:), 0, 1) / sqrt(max(sum(cm), 1));
    if any(cm)
        bar(1:n_sec, mu, 0.55, 'FaceColor', [0.15 0.45 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        errorbar(1:n_sec, mu, se, 'k.', 'LineWidth', 1.2, 'CapSize', 6);

        % chance level (dashed) and stars vs chance
        if opts.n_shuffle_sector > 0
            ch_mu = mean(res.sector.chance(cm,:), 1, 'omitnan');
            yr    = diff(sec_ylim);
            for ss = 1:n_sec
                plot(ss + [-0.32 0.32], [ch_mu(ss) ch_mu(ss)], 'k--', 'LineWidth', 1.1);
                st = sig_stars(res.sector.p_vs_chance_used(ss));
                if isempty(st) || strcmp(st, '(n.s.)'), continue; end
                dv = mean(res.sector.dev(cm, ss), 'omitnan');
                se_ss = se(ss); if isnan(se_ss), se_ss = 0; end
                if dv > 0
                    text(ss, max(mu(ss), ch_mu(ss)) + se_ss + 0.04*yr, st, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                        'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.8 0.1 0.1]);
                else
                    text(ss, min(mu(ss), ch_mu(ss)) - se_ss - 0.04*yr, st, ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                        'FontWeight', 'bold', 'FontSize', 11, 'Color', [0.1 0.3 0.8]);
                end
            end
        end

        % brackets for significant pairwise differences between sectors
        if ~isempty(res.sector.pairs)
            yr  = diff(sec_ylim);
            top = mu + se;  top(isnan(top)) = mu(isnan(top));
            if opts.n_shuffle_sector > 0
                top = max(top, mean(res.sector.chance(cm,:), 1, 'omitnan'));
            end
            base = max(top) + 0.14 * yr;   % above the chance stars
            lvl  = 0;
            for pp = 1:size(res.sector.pairs, 1)
                st = sig_stars(res.sector.p_pairs_holm(pp));
                if isempty(st) || strcmp(st, '(n.s.)'), continue; end
                x1 = res.sector.pairs(pp,1); x2 = res.sector.pairs(pp,2);
                yb = base + lvl * 0.10 * yr;
                plot([x1 x1 x2 x2], [yb - 0.02*yr, yb, yb, yb - 0.02*yr], 'k-', 'LineWidth', 1);
                text(mean([x1 x2]), yb, st, 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 10);
                lvl = lvl + 1;
            end
        end
    end
    yline(0, 'k:', 'LineWidth', 0.6);
    set(gca, 'XTick', 1:n_sec, 'XTickLabel', sector_names, 'FontSize', 8);
    xlim([0.4 n_sec+0.6]); ylim(sec_ylim); grid on; box off;
    if c == 1, ylabel(['Mean ', dff_unit]); end
    if opts.n_shuffle_sector > 0
        title(sprintf(['Sectors (n = %d): Friedman p = %.3g %s\n', ...
            '-- = chance; red/blue * = vs chance (%s, %s)\nbrackets = sector vs sector (%s, Holm)'], ...
            res.sector.n, res.sector.p_friedman, sig_stars(res.sector.p_friedman), ...
            res.sector.chance_test, ternary(strcmpi(res.sector.chance_correction,'holm'), ...
            'Holm', 'uncorrected'), res.sector.posthoc), 'FontSize', 8);
    else
        title(sprintf('Sectors (n = %d): Friedman p = %.3g %s', res.sector.n, ...
            res.sector.p_friedman, sig_stars(res.sector.p_friedman)), 'FontSize', 9);
    end
end
nv1 = out.near_vs_far.harmonic(1);
sg = sprintf('%s  [%s]\nNear vs far, H1 tuning: p = %.3g', col_label, cond_label, nv1.p_hotelling);
if K >= 2
    sg = sprintf('%s,  H2 tuning: p = %.3g', sg, out.near_vs_far.harmonic(2).p_hotelling);
end
sgtitle(sg, 'FontSize', 11, 'Interpreter', 'none');

% ------------------------------------------------------------------ %
%  INDIVIDUAL MICE (one figure per condition)                          %
% ------------------------------------------------------------------ %
if opts.individual
    n_cols = ceil(sqrt(n_mice));
    n_rows_fig = ceil(n_mice / n_cols);
    for c = 1:n_cond
        figure('Name', sprintf('Polar tuning (%s) -- %s -- individual mice', col_label, cond_names{c}), ...
            'Position', [50 50 min(220*n_cols,1400) min(220*n_rows_fig,1000)]);
        for m = 1:n_mice
            subplot(n_rows_fig, n_cols, m);
            fm = [];
            if fit_ok{c}(m)
                cf = coef{c}(m,:);
                fm = cf(1) * ones(size(theta_fit));
                for k = 1:K
                    fm = fm + cf(2*k) * cos(k*theta_fit) + cf(2*k+1) * sin(k*theta_fit);
                end
            end
            plot_one_polar_curve(centers_rad, curves{c}(m,:), [], ...
                sprintf('%s (n=%d)', mice_this{m}, n_used{c}(m)), opts.theta_dir, [], ...
                theta_fit, fm);
        end
        sgtitle(sprintf('%s -- %s  [%s]  (red = harmonic fit)', col_label, cond_names{c}, cond_label), ...
            'FontSize', 11, 'Interpreter', 'none');
    end
end

% ------------------------------------------------------------------ %
%  EXPORT PLOTTED DATA + STATISTICS                                    %
% ------------------------------------------------------------------ %
if opts.export
    try
        if ~exist(opts.export_dir, 'dir'), mkdir(opts.export_dir); end
        tag = sprintf('col%d_%s_%s_sess%s_thr%g', angle_col, opts.state, opts.source, ...
            regexprep(strtrim(num2str(opts.sessions)), '\s+', '-'), thr);
        if ~isempty(opts.dist_remove), tag = sprintf('%s_rm%g', tag, opts.dist_remove); end
        if opts.merge_sides, tag = [tag '_3sectors']; end
        fname = fullfile(opts.export_dir, sprintf('polar_%s_%s.xlsx', tag, datestr(now, 'yyyymmdd')));

        % ---- polar_group (rows 1-2 of the figure) ----
        T = table(centers_deg(:), 'VariableNames', {'Angle_deg'});
        for c = 1:n_cond
            k_ = cond_keys{c};
            T.([k_ '_mean'])          = mean_c{c}(:);
            T.([k_ '_SEM'])           = sem_c{c}(:);
            T.([k_ '_n_mice'])        = sum(~isnan(curves{c}), 1)';
            T.([k_ '_occupancy_pct']) = occ_pct{c}(:);
        end
        writetable(T, fname, 'Sheet', 'polar_group');

        % ---- harmonic_fit_curve (red line) ----
        T = table(rad2deg(theta_fit(:)), 'VariableNames', {'Angle_deg'});
        for c = 1:n_cond
            if isempty(fit_c{c}), f = NaN(numel(theta_fit), 1); else, f = fit_c{c}(:); end
            T.([cond_keys{c} '_fit']) = f;
        end
        writetable(T, fname, 'Sheet', 'harmonic_fit_curve');

        % ---- polar_per_mouse (long format) ----
        Ts = cell(n_cond, 1);
        for c = 1:n_cond
            Ts{c} = table(repmat(mice_this(:), opts.n_bins, 1), ...
                repmat(cond_keys(c), n_mice * opts.n_bins, 1), ...
                kron(centers_deg(:), ones(n_mice, 1)), curves{c}(:), occ{c}(:), ...
                'VariableNames', {'Mouse_ID', 'Condition', 'Angle_deg', 'Mean_dFF', 'N_frames'});
        end
        writetable(vertcat(Ts{:}), fname, 'Sheet', 'polar_per_mouse');

        % ---- harmonic_per_mouse ----
        coef_names = {'c0'};
        for k = 1:K, coef_names = [coef_names, {sprintf('a%d', k), sprintf('b%d', k)}]; end %#ok<AGROW>
        Ts = cell(n_cond, 1);
        for c = 1:n_cond
            Tc = table(mice_this(:), repmat(cond_keys(c), n_mice, 1), fit_ok{c}, n_used{c}, ...
                'VariableNames', {'Mouse_ID', 'Condition', 'Valid_fit', 'N_frames'});
            for j = 1:numel(coef_names)
                Tc.(coef_names{j}) = coef{c}(:, j);
            end
            for k = 1:K
                Tc.(sprintf('amp%d', k)) = hypot(coef{c}(:, 2*k), coef{c}(:, 2*k+1));
            end
            Ts{c} = Tc;
        end
        writetable(vertcat(Ts{:}), fname, 'Sheet', 'harmonic_per_mouse');

        % ---- harmonic_stats (+ near vs far) ----
        rc = {}; rk = []; rn = []; ra = []; rb = []; ramp = []; rpref = [];
        rph = []; rpa = []; rpb = [];
        for c = 1:n_cond
            for k = 1:K
                h = out.(cond_keys{c}).harmonic(k);
                rc{end+1,1} = cond_keys{c}; rk(end+1,1) = k; rn(end+1,1) = h.n; %#ok<AGROW>
                ra(end+1,1) = h.mean_a; rb(end+1,1) = h.mean_b; ramp(end+1,1) = h.amp; %#ok<AGROW>
                rpref(end+1,1) = h.pref_deg; rph(end+1,1) = h.p_hotelling; %#ok<AGROW>
                rpa(end+1,1) = h.p_a; rpb(end+1,1) = h.p_b; %#ok<AGROW>
            end
        end
        for k = 1:K
            nv = out.near_vs_far.harmonic(k);
            rc{end+1,1} = 'near_vs_far'; rk(end+1,1) = k; rn(end+1,1) = nv.n; %#ok<AGROW>
            ra(end+1,1) = nv.mean_amp_near; rb(end+1,1) = nv.mean_amp_far; %#ok<AGROW>
            ramp(end+1,1) = NaN; rpref(end+1,1) = NaN; rph(end+1,1) = nv.p_hotelling; %#ok<AGROW>
            rpa(end+1,1) = nv.p_amp; rpb(end+1,1) = NaN; %#ok<AGROW>
        end
        T = table(rc, rk, rn, ra, rb, ramp, rpref, rph, rpa, rpb, 'VariableNames', ...
            {'Condition', 'Harmonic', 'N_mice', 'Mean_a_or_amp_near', 'Mean_b_or_amp_far', ...
             'Group_amplitude', 'Preferred_deg', 'p_Hotelling', 'p_a_or_p_amp_near_vs_far', 'p_b'});
        writetable(T, fname, 'Sheet', 'harmonic_stats');

        % ---- sectors_per_mouse ----
        Ts = cell(n_cond, 1);
        for c = 1:n_cond
            sc = out.(cond_keys{c}).sector;
            Ts{c} = table(repmat(mice_this(:), n_sec, 1), repmat(cond_keys(c), n_mice * n_sec, 1), ...
                repelem(sector_names(:), n_mice), sc.means(:), sc.chance(:), sc.dev(:), ...
                repmat(sc.complete, n_sec, 1), ...
                'VariableNames', {'Mouse_ID', 'Condition', 'Sector', 'Mean_dFF', ...
                                  'Chance', 'Real_minus_chance', 'In_sector_tests'});
        end
        writetable(vertcat(Ts{:}), fname, 'Sheet', 'sectors_per_mouse');

        % ---- sector_stats (bars, chance lines, stars) ----
        Ts = cell(n_cond, 1);
        for c = 1:n_cond
            sc = out.(cond_keys{c}).sector;
            cm = sc.complete;
            nc = sum(cm);
            mu_ = mean(sc.means(cm,:), 1);
            se_ = std(sc.means(cm,:), 0, 1) / sqrt(max(nc, 1));
            Ts{c} = table(repmat(cond_keys(c), n_sec, 1), sector_names(:), repmat(nc, n_sec, 1), ...
                mu_(:), se_(:), mean(sc.chance(cm,:), 1, 'omitnan')', ...
                mean(sc.dev(cm,:), 1, 'omitnan')', sc.p_vs_chance_t(:), sc.p_vs_chance_sr(:), ...
                sc.p_vs_chance_holm(:), sc.p_vs_chance_used(:), ...
                'VariableNames', {'Condition', 'Sector', 'N_mice', 'Bar_mean', 'Bar_SEM', ...
                'Chance_mean', 'Real_minus_chance', 'p_ttest_vs_chance', ...
                'p_signrank_vs_chance', 'p_Holm_vs_chance', 'p_used_for_stars'});
        end
        writetable(vertcat(Ts{:}), fname, 'Sheet', 'sector_stats');

        % ---- sector_pairs (brackets) ----
        Ts = {};
        for c = 1:n_cond
            sc = out.(cond_keys{c}).sector;
            if isempty(sc.pairs), continue; end
            np = size(sc.pairs, 1);
            Ts{end+1, 1} = table(repmat(cond_keys(c), np, 1), sc.pair_names, ...
                sc.pair_mean_diff, sc.p_pairs_ttest, sc.p_pairs_signrank, sc.p_pairs, ...
                sc.p_pairs_holm, repmat({sc.posthoc}, np, 1), ...
                'VariableNames', {'Condition', 'Pair', 'Mean_difference', 'p_ttest', ...
                'p_signrank', 'p_used', 'p_Holm', 'Test_used'}); %#ok<AGROW>
        end
        if ~isempty(Ts)
            writetable(vertcat(Ts{:}), fname, 'Sheet', 'sector_pairs');
        end

        % ---- friedman ----
        pf = NaN(n_cond, 1); nf = NaN(n_cond, 1);
        for c = 1:n_cond
            pf(c) = out.(cond_keys{c}).sector.p_friedman;
            nf(c) = out.(cond_keys{c}).sector.n;
        end
        writetable(table(cond_keys(:), cond_names(:), nf, pf, 'VariableNames', ...
            {'Condition', 'Description', 'N_mice', 'p_Friedman'}), fname, 'Sheet', 'friedman');

        fprintf('\nPlotted data and statistics saved to: %s\n', fname);
        out.excel_filename = fname;
    catch ME
        fprintf('\nCould not save Excel file: %s\n', ME.message);
    end
end

end


% ====================================================================== %
%  HARMONIC REGRESSION                                                   %
% ====================================================================== %
function b = fit_harmonics(th, y, K)
X = ones(numel(th), 1 + 2*K);
for k = 1:K
    X(:, 2*k)   = cos(k*th);
    X(:, 2*k+1) = sin(k*th);
end
b = (X \ y(:))';
end


% ====================================================================== %
%  ONE-SAMPLE HOTELLING'S T^2 (H0: mean vector = 0)                      %
% ====================================================================== %
function [p, T2] = hotelling_one_sample(X)
[n, d] = size(X);
mu = mean(X, 1)';
S  = cov(X);
T2 = n * (mu' * pinv(S) * mu);
F  = (n - d) / (d * (n - 1)) * T2;
p  = 1 - fcdf(F, d, n - d);
end


% ====================================================================== %
%  HOLM-BONFERRONI CORRECTION                                            %
% ====================================================================== %
function p_adj = holm(p)
p_adj = NaN(size(p));
ok = ~isnan(p);
pv = p(ok);
m  = numel(pv);
[ps, o] = sort(pv);
adj = zeros(m, 1);
running = 0;
for i = 1:m
    running = max(running, min(1, (m - i + 1) * ps(i)));
    adj(i) = running;
end
tmp = zeros(m, 1);
tmp(o) = adj;
p_adj(ok) = tmp;
end


% ====================================================================== %
%  SECTOR INDEX                                                          %
%    merged : 1 front, 2 sides (both), 3 back                            %
%    split  : 1 front, 2 +90 side, 3 back, 4 -90 side                    %
% ====================================================================== %
function s = sector_index(ang, merge_sides)
s = NaN(size(ang));
s(abs(ang) <= 45) = 1;
s(abs(ang) > 135) = 3;
if merge_sides
    s(abs(ang) > 45 & abs(ang) <= 135) = 2;
else
    s(ang > 45 & ang <= 135)   = 2;
    s(ang >= -135 & ang < -45) = 4;
end
end


% ====================================================================== %
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


% ====================================================================== %
function a = wrap180(a)
a = mod(a + 180, 360) - 180;
end


% ====================================================================== %
function s = sig_stars(p)
if isnan(p),      s = '';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
else,             s = '(n.s.)';
end
end


% ====================================================================== %
%  GROUP CURVE (mean +- SEM across mice)                                 %
% ====================================================================== %
function [mean_curve, sem_curve] = group_curve(curves, occ, weight_by_occupancy)
if weight_by_occupancy
    w = occ;
    w(isnan(curves)) = 0;
    cc = curves;
    cc(isnan(curves)) = 0;
    sum_w = sum(w, 1);
    mean_curve = sum(cc .* w, 1) ./ sum_w;
    mean_curve(sum_w == 0) = NaN;
else
    mean_curve = mean(curves, 1, 'omitnan');
end
n_per_bin = sum(~isnan(curves), 1);
sem_curve = std(curves, 0, 1, 'omitnan') ./ sqrt(max(n_per_bin, 1));
sem_curve(n_per_bin < 2) = NaN;
end


% ====================================================================== %
%  DRAW ONE POLAR TUNING CURVE (+ optional fitted curve in red)          %
%  Radius is offset so it is never negative on screen; radial ticks are  %
%  relabeled with the true values. rlim_true = [lo hi] fixes the scale   %
%  (shared across panels); [] = auto per panel.                          %
% ====================================================================== %
function plot_one_polar_curve(theta, rho, rho_sem, name, theta_dir, rlim_true, fit_theta, fit_rho)

if nargin < 5 || isempty(theta_dir), theta_dir = 'counterclockwise'; end
if nargin < 6, rlim_true = []; end
if nargin < 8, fit_theta = []; fit_rho = []; end

valid = ~isnan(rho);
if ~any(valid)
    ax0 = gca; pa = polaraxes('Position', ax0.Position); delete(ax0);
    title(pa, strrep(name,'_',' '), 'FontSize', 9, 'Interpreter', 'none');
    return;
end

has_sem = ~isempty(rho_sem) && any(~isnan(rho_sem));

if isempty(rlim_true)
    min_val = min(rho(valid));
    if has_sem
        lb = rho(valid) - rho_sem(valid);
        lb = lb(~isnan(lb));
        if ~isempty(lb), min_val = min(min_val, min(lb)); end
    end
    if ~isempty(fit_rho), min_val = min(min_val, min(fit_rho)); end
    offset = min_val;
else
    offset = rlim_true(1);
end
rho_disp = rho - offset;

theta_c = [theta, theta(1)];
rho_c   = [rho_disp, rho_disp(1)];

polarplot(theta_c, rho_c, '-o', 'LineWidth', 1.3, 'MarkerSize', 3, ...
    'Color', [0.15 0.45 0.75]);
hold on;
if has_sem
    upper = rho_disp + rho_sem;
    lower = rho_disp - rho_sem;
    polarplot(theta_c, [upper, upper(1)], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    polarplot(theta_c, [lower, lower(1)], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
end
if ~isempty(fit_rho)
    polarplot(fit_theta, fit_rho - offset, '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.6);
end

ax = gca;
ax.ThetaZeroLocation = 'top';   % 0 deg (facing target) at 12 o'clock
ax.ThetaDir = theta_dir;
ax.ThetaTick = 0:30:330;
if ~isempty(rlim_true)
    rlim(ax, [0, rlim_true(2) - rlim_true(1)]);
end

rticks_disp = ax.RTick;
ax.RTickLabel = arrayfun(@(v) sprintf('%.2f', v), rticks_disp + offset, ...
    'UniformOutput', false);

title(strrep(name,'_',' '), 'FontSize', 9, 'Interpreter', 'none');
end


% ====================================================================== %
%  DRAW ONE OCCUPANCY PLOT                                               %
% ====================================================================== %
function plot_occupancy(theta, occ_pct, name, theta_dir, rmax)
if all(isnan(occ_pct))
    ax0 = gca; pa = polaraxes('Position', ax0.Position); delete(ax0);
    title(pa, name, 'FontSize', 9, 'Interpreter', 'none');
    return;
end
theta_c = [theta, theta(1)];
occ_c   = [occ_pct, occ_pct(1)];
polarplot(theta_c, occ_c, '-', 'LineWidth', 1.3, 'Color', [0.5 0.5 0.5]);
hold on;
polarplot(theta_c, occ_c, 'o', 'MarkerSize', 3, 'Color', [0.5 0.5 0.5], ...
    'MarkerFaceColor', [0.5 0.5 0.5]);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = theta_dir;
ax.ThetaTick = 0:30:330;
rlim(ax, [0, rmax * 1.05]);
title(name, 'FontSize', 9, 'Interpreter', 'none');
end