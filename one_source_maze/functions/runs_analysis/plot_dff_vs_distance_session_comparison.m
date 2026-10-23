function results = plot_dff_vs_distance_session_comparison(run_data, options)
    % Plot dF/F vs distance for food runs, comparing session 0 and session 1
    % - Both towards and away runs on the same subplot
    % - Session 0 block on top, Session 1 block on bottom
    % - Under every run trace: bar plot of per-mouse slopes (towards vs away),
    %   mean +/- SEM, dots = mice, stars = one-sample t-test vs 0
    % - Extra figures (one per session x direction): slopes across runs 1..n
    %   with ANOVA across runs and a linear trend test over run number
    %
    % Slope method (same as analyze_single_arm_slopes):
    %   For each mouse and run, (z-scored) dF/F is averaged in distance bins of
    %   width .bin_width within .slope_range of the food. Towards bins get
    %   negative x (-distance), away bins positive x (+distance). A line is fit
    %   to the bin means (>= .slope_min_bins bins). So for towards runs a
    %   positive slope means activity rises while approaching food, for away
    %   runs a positive slope means activity rises while leaving food.
    %
    % Inputs:
    %   run_data - Structure with run data from analyze_single_arm_runs_with_speed_check
    %   options - Structure with the following fields (all optional):
    %       .food_source - Food source to analyze (default: 'food')
    %       .smoothing - Window size for smoothing the traces (default: 5)
    %       .max_distance - Maximum distance range to plot (default: 200)
    %       .ylim_dff - Y-axis limits for dF/F (default: automatic)
    %       .distance_points - Number of distance points for interpolation (default: 200)
    %       .min_mice - Minimum number of mice required for a run (default: 1)
    %       .max_runs - Maximum number of runs to plot (default: 10)
    %       .use_zscore - Whether to use z-scored dF/F (default: true)
    %       .gap - Gap between towards and away runs (default: 10)
    %       .scalebar_distance - Length of distance scale bar in units (default: 50)
    %       .scalebar_dff - Height of dF/F scale bar (default: 1)
    %       .run_category - 'food', 'not_food', 'all' (default: 'food')
    %       .slope_range - Distance from food used for slopes (default: 50)
    %       .bin_width - Distance bin width for slopes (default: 1)
    %       .slope_min_bins - Min. filled bins for a mouse slope (default: 3)
    %       .slope_per_distance - Slopes reported per this many cm (default: 1)
    %       .slope_alpha - Significance level (default: 0.05)
    %       .ylim_slope - Y-axis limits for slope bar plots (default: automatic)
    %       .plot_run_anova - Make the slope-across-runs figures (default: true)
    %       .save_data - Save all plotted data (traces + bar plots) (default: false)
    %       .save_path - Folder to save into (created if missing). If it ends in
    %                    .mat/.xlsx/.csv it is used as the file name instead.
    %                    (default: current folder)
    %       .save_name - Base file name, no extension (default:
    %                    dff_vs_distance_<food_source>_<run_category>_<timestamp>)
    %       .save_format - 'xlsx', 'mat', 'csv', 'both' (= mat + xlsx), 'all',
    %                      or a cell array e.g. {'mat','xlsx'} (default: 'both')
    %
    %   Excel file tabs:
    %       Sess0 traces / Sess1 traces - mean and SEM across mice of the plotted
    %           traces; rows = distance (towards then away), columns = runs
    %       Sess0 slopes / Sess1 slopes - mean, SEM, n, p (vs 0) of the slopes,
    %           columns = runs; below: every mouse's slope per run
    %       Runs ANOVA - ANOVA / trend across runs for each session x direction
    %
    % Output (optional):
    %   results.slopes - per-mouse slopes for every session/direction/run
    %   results.anova  - ANOVA and trend results for every session/direction
    %   results.plotted_data - everything that was plotted (see save_plotted_data)
    
    if nargin < 2
        options = struct();
    end
    results = struct();
    
    defaults = struct( ...
        'smoothing', 5, ...
        'food_source', 'food', ...
        'max_distance', 200, ...
        'distance_points', 200, ...
        'min_mice', 1, ...
        'max_runs', 10, ...
        'use_zscore', true, ...
        'gap', 10, ...
        'scalebar_distance', 50, ...
        'scalebar_dff', 1, ...
        'run_category', 'food', ...
        'slope_range', 50, ...
        'bin_width', 1, ...
        'slope_min_bins', 3, ...
        'slope_per_distance', 1, ...
        'slope_alpha', 0.05, ...
        'plot_run_anova', true, ...
        'save_data', false, ...
        'save_path', '', ...
        'save_name', '', ...
        'save_format', 'both');
    default_names = fieldnames(defaults);
    for i = 1:numel(default_names)
        nm = default_names{i};
        if ~isfield(options, nm) || isempty(options.(nm))
            options.(nm) = defaults.(nm);
        end
    end
    if ~isfield(options, 'ylim_dff')
        options.ylim_dff = [];
    end
    if ~isfield(options, 'ylim_slope')
        options.ylim_slope = [];
    end
    
    switch options.run_category
        case 'food'
            valid_towards_types = {'towards'};
            valid_away_types = {'away'};
            category_str = ' (Food Runs)';
        case 'not_food'
            valid_towards_types = {'not_food_towards'};
            valid_away_types = {'not_food_away'};
            category_str = ' (Non-Food Runs)';
        case 'all'
            valid_towards_types = {'towards', 'not_food_towards'};
            valid_away_types = {'away', 'not_food_away'};
            category_str = ' (All Runs)';
        otherwise
            error('Invalid run_category option. Must be ''food'', ''not_food'', or ''all''');
    end
    
    session_nums = [0, 1];
    
    fprintf('Analyzing runs for sessions 0 and 1, food source: %s%s\n', options.food_source, category_str);
    fprintf('Maximum runs to plot: %d, Minimum mice per run: %d\n', options.max_runs, options.min_mice);
    
    % Organize mice by session
    session_data = struct('session', {}, 'mice', {});
    for m = 1:length(run_data)
        sess = run_data(m).session;
        if ~ismember(sess, session_nums)
            continue;
        end
        sess_idx = find([session_data.session] == sess, 1);
        if isempty(sess_idx)
            sess_idx = length(session_data) + 1;
            session_data(sess_idx).session = sess;
            session_data(sess_idx).mice = [];
        end
        session_data(sess_idx).mice = [session_data(sess_idx).mice, m];
    end
    
    session_indices = [session_data.session];
    if ~ismember(0, session_indices)
        warning('No data found for session 0');
        return;
    end
    if ~ismember(1, session_indices)
        warning('No data found for session 1');
        return;
    end
    
    % Gather run information per session / direction
    session_run_data = struct();
    for s = 1:length(session_data)
        sess = session_data(s).session;
        mice_indices = session_data(s).mice;
        session_run_data.(sprintf('sess%d', sess)).session = sess;
        
        for direction_idx = 1:2
            if direction_idx == 1
                valid_types = valid_towards_types;
                type_label = 'towards';
            else
                valid_types = valid_away_types;
                type_label = 'away';
            end
            
            field_name = sprintf('sess%d_%s', sess, type_label);
            session_run_data.(field_name) = struct();
            
            mouse_food_runs = cell(length(mice_indices), 1);
            max_runs = 0;
            for mi = 1:length(mice_indices)
                m = mice_indices(mi);
                mouse_runs = run_data(m).runs;
                run_indices = [];
                for r = 1:length(mouse_runs)
                    if any(strcmp(mouse_runs(r).type, valid_types))
                        run_indices = [run_indices, r]; %#ok<AGROW>
                    end
                end
                mouse_food_runs{mi} = run_indices;
                max_runs = max(max_runs, length(run_indices));
            end
            
            if max_runs == 0
                warning('No %s runs found for session %d%s', type_label, sess, category_str);
                session_run_data.(field_name).valid_runs = [];
                session_run_data.(field_name).mouse_food_runs = {};
                continue;
            end
            
            mice_per_run = zeros(max_runs, 1);
            for mi = 1:length(mice_indices)
                n_r = min(length(mouse_food_runs{mi}), max_runs);
                mice_per_run(1:n_r) = mice_per_run(1:n_r) + 1;
            end
            
            valid_runs = find(mice_per_run >= options.min_mice);
            valid_runs = valid_runs(valid_runs <= options.max_runs);
            
            session_run_data.(field_name).valid_runs = valid_runs;
            session_run_data.(field_name).mouse_food_runs = mouse_food_runs;
            session_run_data.(field_name).mice_indices = mice_indices;
            session_run_data.(field_name).max_runs = max_runs;
        end
    end
    
    all_valid_runs = unique([
        session_run_data.sess0_towards.valid_runs(:);
        session_run_data.sess0_away.valid_runs(:);
        session_run_data.sess1_towards.valid_runs(:);
        session_run_data.sess1_away.valid_runs(:)
    ]);
    
    if isempty(all_valid_runs)
        warning('No runs with >= %d mice found for sessions 0 and 1 within first %d runs%s', options.min_mice, options.max_runs, category_str);
        return;
    end
    
    all_valid_runs = sort(all_valid_runs);
    all_valid_runs = all_valid_runs(all_valid_runs <= options.max_runs);
    num_valid_runs = length(all_valid_runs);
    
    fprintf('Plotting %d runs (runs %s)%s\n', num_valid_runs, mat2str(all_valid_runs), category_str);
    
    % Global dF/F limits
    [min_dff, max_dff] = calculate_global_dff_limits(run_data, session_run_data, session_nums, ...
                                                      all_valid_runs, num_valid_runs, options);
    if ~isempty(options.ylim_dff)
        min_dff = options.ylim_dff(1);
        max_dff = options.ylim_dff(2);
    end
    
    % Per-mouse slopes for every session / direction / run
    [slope_results, slope_lim] = calculate_all_slopes(run_data, session_run_data, session_nums, ...
                                                      all_valid_runs, options);
    results.slopes = slope_results;
    results.runs = all_valid_runs;
    
    % ---------------- Main figure: traces + slope bars ----------------
    fig = figure('Name', sprintf('dF/F vs Distance - Sessions Comparison for %s (First %d Runs)%s', ...
                                 options.food_source, options.max_runs, category_str), ...
                 'Position', [100, 100, max(260*num_valid_runs, 700), 850], ...
                 'Color', 'w');
    
    layout = compute_layout(num_valid_runs);
    
    [first_subplot_handle, trace_data] = plot_all_sessions_and_runs(run_data, session_run_data, session_nums, ...
                                                       all_valid_runs, num_valid_runs, ...
                                                       min_dff, max_dff, options, ...
                                                       slope_results, slope_lim, layout);
    
    add_scale_bars(first_subplot_handle, options, min_dff, max_dff);
    
    figure_title = sprintf('dF/F vs Distance Comparison - %s (First %d Runs)%s', ...
                           options.food_source, options.max_runs, category_str);
    if options.use_zscore
        figure_title = [figure_title, ' (Z-scored)'];
    end
    sgtitle(fig, figure_title, 'FontSize', 14);
    
    add_session_labels(fig, layout, session_nums);
    
    % ---------------- Slopes across runs: one figure per session x direction ----------------
    if options.plot_run_anova
        results.anova = plot_slopes_across_runs(slope_results, all_valid_runs, session_nums, ...
                                                options, category_str);
    end
    
    % ---------------- Collect (and optionally save) everything plotted ----------------
    plotted_data = build_plotted_data(run_data, trace_data, slope_results, all_valid_runs, ...
                                      session_nums, results, options, category_str, ...
                                      [min_dff, max_dff], slope_lim);
    results.plotted_data = plotted_data;
    
    if options.save_data
        save_plotted_data(plotted_data, options);
    end
end

% =========================================================================
% Layout
% =========================================================================

function layout = compute_layout(n_cols)
    left = 0.12; right = 0.02; top = 0.10; bottom = 0.05;
    block_gap = 0.075;
    inner_gap = 0.035;
    col_gap = 0.015;
    
    block_h = (1 - top - bottom - block_gap) / 2;
    trace_h = (block_h - inner_gap) * 0.65;
    bar_h   = (block_h - inner_gap) * 0.35;
    col_w   = (1 - left - right - col_gap * (n_cols - 1)) / n_cols;
    
    layout.trace_pos = cell(2, n_cols);
    layout.bar_pos = cell(2, n_cols);
    layout.block_top = zeros(1, 2);
    layout.block_bottom = zeros(1, 2);
    
    block_top = 1 - top;
    for s = 1:2
        trace_bottom = block_top - trace_h;
        bar_bottom = trace_bottom - inner_gap - bar_h;
        for c = 1:n_cols
            x = left + (c - 1) * (col_w + col_gap);
            layout.trace_pos{s, c} = [x, trace_bottom, col_w, trace_h];
            layout.bar_pos{s, c}   = [x, bar_bottom,   col_w, bar_h];
        end
        layout.block_top(s) = block_top;
        layout.block_bottom(s) = bar_bottom;
        block_top = bar_bottom - block_gap;
    end
    
    layout.left = left;
    layout.right = right;
    layout.divider_y = (layout.block_bottom(1) + layout.block_top(2)) / 2;
end

% =========================================================================
% Limits
% =========================================================================

function [min_dff, max_dff] = calculate_global_dff_limits(run_data, session_run_data, session_nums, ...
                                                           all_valid_runs, num_valid_runs, options)
    min_dff = Inf;
    max_dff = -Inf;
    type_labels = {'towards', 'away'};
    
    for sess = session_nums
        for run_idx = 1:num_valid_runs
            run_num = all_valid_runs(run_idx);
            for d = 1:2
                type_label = type_labels{d};
                field_name = sprintf('sess%d_%s', sess, type_label);
                if ~isfield(session_run_data, field_name) || ...
                   isempty(session_run_data.(field_name).valid_runs) || ...
                   ~ismember(run_num, session_run_data.(field_name).valid_runs)
                    continue;
                end
                avg_dff = process_run_data(run_data, session_run_data.(field_name), run_num, ...
                                          type_label, options);
                valid_avg = avg_dff(~isnan(avg_dff));
                if ~isempty(valid_avg)
                    min_dff = min(min_dff, min(valid_avg));
                    max_dff = max(max_dff, max(valid_avg));
                end
            end
        end
    end
    
    if ~isfinite(min_dff) || ~isfinite(max_dff)
        min_dff = -1; max_dff = 1;
    end
    padding = (max_dff - min_dff) * 0.1;
    min_dff = min_dff - padding;
    max_dff = max_dff + padding;
end

% =========================================================================
% Slopes
% =========================================================================

function [slope_results, slope_lim] = calculate_all_slopes(run_data, session_run_data, session_nums, ...
                                                           all_valid_runs, options)
    slope_results = struct();
    all_vals = [];
    type_labels = {'towards', 'away'};
    
    fprintf('\n===== Per-run slopes (within %d cm of food, bin %.1f cm) =====\n', ...
            options.slope_range, options.bin_width);
    
    for sess = session_nums
        for d = 1:2
            type_label = type_labels{d};
            field_name = sprintf('sess%d_%s', sess, type_label);
            slope_results.(field_name) = cell(1, numel(all_valid_runs));
            
            if ~isfield(session_run_data, field_name) || ...
               isempty(session_run_data.(field_name).valid_runs)
                continue;
            end
            
            for k = 1:numel(all_valid_runs)
                run_num = all_valid_runs(k);
                if ~ismember(run_num, session_run_data.(field_name).valid_runs)
                    continue;
                end
                [slopes, mouse_ids] = compute_run_slopes(run_data, session_run_data.(field_name), ...
                                                         run_num, type_label, options);
                p_vs0 = NaN;
                if numel(slopes) >= 2
                    [~, p_vs0] = ttest(slopes, 0);
                end
                slope_results.(field_name){k} = struct('slopes', slopes, 'mice', mouse_ids, ...
                                                       'p_vs_zero', p_vs0);
                all_vals = [all_vals; slopes(:)]; %#ok<AGROW>
                
                if ~isempty(slopes)
                    fprintf('Session %d, run %2d, %-7s: slope = %8.4f +/- %.4f (n=%d), p vs 0 = %.4f%s\n', ...
                            sess, run_num, type_label, mean(slopes), ...
                            std(slopes) / sqrt(numel(slopes)), numel(slopes), p_vs0, ...
                            sig_stars(p_vs0, options.slope_alpha, true));
                end
            end
        end
    end
    
    if ~isempty(options.ylim_slope)
        slope_lim = options.ylim_slope;
    elseif isempty(all_vals)
        slope_lim = [-1, 1];
    else
        lo = min([all_vals; 0]);
        hi = max([all_vals; 0]);
        rng = hi - lo;
        if rng == 0, rng = 1; end
        slope_lim = [lo - 0.1 * rng, hi + 0.2 * rng]; % extra room on top for stars
    end
end

function [slopes, mouse_ids] = compute_run_slopes(run_data, field_data, run_num, type_label, options)
    % Slope per mouse for one run, same method as analyze_single_arm_slopes:
    % bin dF/F by distance from food within slope_range, fit line to bin means.
    slopes = [];
    mouse_ids = [];
    
    edges = 0:options.bin_width:options.slope_range;
    if edges(end) < options.slope_range
        edges(end+1) = options.slope_range;
    end
    centers = edges(1:end-1) + diff(edges) / 2;
    if strcmp(type_label, 'towards')
        x_centers = -centers;   % approaching food: x goes from -range up to 0
    else
        x_centers = centers;    % leaving food: x goes from 0 up to +range
    end
    n_bins = numel(centers);
    
    for mi = 1:numel(field_data.mice_indices)
        m = field_data.mice_indices(mi);
        mouse_run_indices = field_data.mouse_food_runs{mi};
        if run_num > numel(mouse_run_indices)
            continue;
        end
        run_idx = mouse_run_indices(run_num);
        
        dff = get_dff(run_data, m, run_idx, options);
        if isempty(dff)
            continue;
        end
        dist = run_data(m).runs(run_idx).distance(:);
        dff = dff(:);
        
        b = discretize(dist, edges);
        ok = ~isnan(b) & ~isnan(dff);
        if ~any(ok)
            continue;
        end
        bin_means = accumarray(b(ok), dff(ok), [n_bins, 1], @mean, NaN);
        
        filled = ~isnan(bin_means);
        if sum(filled) < options.slope_min_bins
            continue;
        end
        
        X = [ones(sum(filled), 1), x_centers(filled)'];
        beta = X \ bin_means(filled);
        slopes(end+1, 1) = beta(2) * options.slope_per_distance; %#ok<AGROW>
        mouse_ids(end+1, 1) = m; %#ok<AGROW>
    end
end

function str = slope_units_str(options)
    if options.use_zscore
        unit_str = 'Z';
    else
        unit_str = 'dF/F';
    end
    if options.slope_per_distance == 1
        str = sprintf('Slope (%s/cm)', unit_str);
    else
        str = sprintf('Slope (%s/%g cm)', unit_str, options.slope_per_distance);
    end
end

function s = sig_stars(p, alpha, with_space)
    if nargin < 3, with_space = false; end
    if isnan(p) || p >= alpha
        s = '';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    else
        s = '*';
    end
    if with_space && ~isempty(s)
        s = [' ', s];
    end
end

% =========================================================================
% Main figure plotting
% =========================================================================

function [first_subplot_handle, trace_data] = plot_all_sessions_and_runs(run_data, session_run_data, session_nums, ...
                                                            all_valid_runs, num_valid_runs, ...
                                                            min_dff, max_dff, options, ...
                                                            slope_results, slope_lim, layout)
    first_subplot_handle = [];
    
    % trace_data.(sessX_dir){k} holds exactly what is drawn for run all_valid_runs(k)
    trace_data = struct();
    for sess = session_nums
        trace_data.(sprintf('sess%d_towards', sess)) = cell(1, num_valid_runs);
        trace_data.(sprintf('sess%d_away', sess)) = cell(1, num_valid_runs);
    end
    xlim_extended = [-options.max_distance * 1.6, options.max_distance];
    
    for sess_idx = 1:2
        sess = session_nums(sess_idx);
        
        for run_idx = 1:num_valid_runs
            run_num = all_valid_runs(run_idx);
            
            % ---------- Trace axes ----------
            ax = axes('Position', layout.trace_pos{sess_idx, run_idx});
            if sess_idx == 1 && run_idx == 1
                first_subplot_handle = ax;
            end
            hold(ax, 'on');
            box(ax, 'off');
            set(ax, 'XColor', 'none', 'YColor', 'none');
            
            for direction_idx = 1:2
                if direction_idx == 1
                    type_label = 'towards';
                    color = [0, 0, 0.8];
                else
                    type_label = 'away';
                    color = [0.8, 0, 0];
                end
                
                field_name = sprintf('sess%d_%s', sess, type_label);
                if ~isfield(session_run_data, field_name) || ...
                   isempty(session_run_data.(field_name).valid_runs) || ...
                   ~ismember(run_num, session_run_data.(field_name).valid_runs)
                    continue;
                end
                
                [distance_grid, avg_dff, sem_dff, mice_with_run, interp_dff, mice_used] = ...
                    get_processed_run_data(run_data, session_run_data.(field_name), run_num, ...
                                          type_label, options);
                if isempty(distance_grid)
                    continue;
                end
                
                trace_data.(field_name){run_idx} = struct( ...
                    'run', run_num, ...
                    'distance', distance_grid(:), ...
                    'mean', avg_dff(:), ...
                    'sem', sem_dff(:), ...
                    'n_mice', mice_with_run, ...
                    'mice_index', mice_used(:), ...
                    'per_mouse', interp_dff');   % points x mice, unsmoothed
                
                shadedErrorBar(distance_grid, avg_dff, sem_dff, {'Color', color, 'LineWidth', 1.5}, 0.2);
                
                if strcmp(type_label, 'towards')
                    text_x = -options.max_distance * 0.95;
                else
                    text_x = options.max_distance * 0.05;
                end
                text_y = max_dff - 0.1 * (max_dff - min_dff);
                text(ax, text_x, text_y, sprintf('n=%d', mice_with_run), 'Color', color, ...
                     'FontSize', 8, 'VerticalAlignment', 'top');
            end
            
            xlim(ax, xlim_extended);
            ylim(ax, [min_dff, max_dff]);
            plot(ax, [-options.max_distance, options.max_distance], [0, 0], 'k--', 'LineWidth', 0.5);
            title(ax, sprintf('Run %d', run_num), 'FontSize', 10);
            grid(ax, 'off');
            
            % ---------- Slope bar axes ----------
            ax_bar = axes('Position', layout.bar_pos{sess_idx, run_idx});
            plot_slope_bars(ax_bar, slope_results, sess, run_idx, slope_lim, options, run_idx == 1);
        end
    end
end

function plot_slope_bars(ax, slope_results, sess, k, slope_lim, options, show_yaxis)
    % Towards vs away slope bars for one session/run
    hold(ax, 'on');
    colors = {[0, 0, 0.8], [0.8, 0, 0]};
    type_labels = {'towards', 'away'};
    data = cell(1, 2);
    star_y = slope_lim(2) - 0.08 * diff(slope_lim);
    
    for d = 1:2
        field_name = sprintf('sess%d_%s', sess, type_labels{d});
        if isfield(slope_results, field_name) && k <= numel(slope_results.(field_name))
            data{d} = slope_results.(field_name){k};
        end
        if isempty(data{d}) || isempty(data{d}.slopes)
            continue;
        end
        s = data{d}.slopes;
        mu = mean(s);
        se = 0;
        if numel(s) > 1
            se = std(s) / sqrt(numel(s));
        end
        bar(ax, d, mu, 0.6, 'FaceColor', colors{d} * 0.35 + 0.65, ...
            'EdgeColor', colors{d}, 'LineWidth', 1);
        errorbar(ax, d, mu, se, 'Color', 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4);
        
        stars = sig_stars(data{d}.p_vs_zero, options.slope_alpha);
        if ~isempty(stars)
            text(ax, d, star_y, stars, 'HorizontalAlignment', 'center', ...
                 'FontSize', 10, 'FontWeight', 'bold', 'Color', colors{d});
        end
    end
    
    % Paired lines for mice with both towards and away slopes
    if ~isempty(data{1}) && ~isempty(data{2}) && ~isempty(data{1}.slopes) && ~isempty(data{2}.slopes)
        [common, i1, i2] = intersect(data{1}.mice, data{2}.mice);
        for j = 1:numel(common)
            plot(ax, [1, 2], [data{1}.slopes(i1(j)), data{2}.slopes(i2(j))], '-', ...
                 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
        end
    end
    
    for d = 1:2
        if isempty(data{d}) || isempty(data{d}.slopes)
            continue;
        end
        scatter(ax, d * ones(size(data{d}.slopes)), data{d}.slopes, 14, colors{d}, ...
                'filled', 'MarkerEdgeColor', 'w');
    end
    
    plot(ax, [0.4, 2.6], [0, 0], 'k-', 'LineWidth', 0.5);
    xlim(ax, [0.4, 2.6]);
    ylim(ax, slope_lim);
    set(ax, 'XTick', [1, 2], 'XTickLabel', {'To', 'Away'}, 'TickLength', [0, 0], ...
        'FontSize', 8, 'Box', 'off', 'XColor', [0.3, 0.3, 0.3]);
    
    if show_yaxis
        ylabel(ax, slope_units_str(options), 'FontSize', 8);
    else
        set(ax, 'YColor', 'none');
    end
end

% =========================================================================
% Slopes across runs (chronological) + ANOVA
% =========================================================================

function anova_results = plot_slopes_across_runs(slope_results, all_valid_runs, session_nums, ...
                                                 options, category_str)
    % One figure per session x direction: slope vs run number,
    % ANOVA across runs and a linear trend test over run number.
    type_labels = {'towards', 'away'};
    nice_labels = {'Towards', 'Away'};
    colors = {[0, 0, 0.8], [0.8, 0, 0]};
    anova_results = struct();
    
    fprintf('\n===== Slopes across runs (chronological) =====\n');
    
    for s = 1:numel(session_nums)
        sess = session_nums(s);
        for d = 1:2
            field_name = sprintf('sess%d_%s', sess, type_labels{d});
            color = colors{d};
            
            % Gather data
            y = []; g = []; mice = [];
            per_run = struct('run', {}, 'slopes', {}, 'mice', {}, 'p_vs_zero', {});
            if isfield(slope_results, field_name)
                for k = 1:numel(all_valid_runs)
                    e = slope_results.(field_name){k};
                    if isempty(e) || isempty(e.slopes)
                        continue;
                    end
                    y = [y; e.slopes(:)]; %#ok<AGROW>
                    g = [g; repmat(all_valid_runs(k), numel(e.slopes), 1)]; %#ok<AGROW>
                    mice = [mice; e.mice(:)]; %#ok<AGROW>
                    per_run(end+1) = struct('run', all_valid_runs(k), 'slopes', e.slopes(:), ...
                                            'mice', e.mice(:), 'p_vs_zero', e.p_vs_zero); %#ok<AGROW>
                end
            end
            
            if isempty(per_run)
                fprintf('Session %d %s: no slopes available.\n', sess, type_labels{d});
                continue;
            end
            
            res = run_slope_anova(y, g, mice);
            res.per_run = per_run;
            anova_results.(field_name) = res;
            
            fprintf('Session %d %-7s: %s: F(%g,%g) = %.2f, p = %.4f%s | trend: %.4f per run, p = %.4f%s (%s)\n', ...
                    sess, type_labels{d}, res.anova_type, res.df(1), res.df(2), res.F, res.p, ...
                    sig_stars(res.p, options.slope_alpha, true), ...
                    res.trend_slope, res.trend_p, sig_stars(res.trend_p, options.slope_alpha, true), ...
                    res.trend_method);
            
            % ---------- Figure ----------
            run_nums = [per_run.run];
            n_runs = numel(run_nums);
            figure('Name', sprintf('Session %d - %s slopes across runs%s', sess, nice_labels{d}, category_str), ...
                   'Color', 'w', 'Position', [150 + 40*(2*(s-1)+d), 150, max(70*n_runs + 250, 450), 420]);
            ax = axes('Position', [0.15, 0.13, 0.8, 0.68]);
            hold(ax, 'on');
            
            % Mouse trajectories across runs (chronological)
            umice = unique(mice);
            for mi = 1:numel(umice)
                idx = mice == umice(mi);
                [rx, order] = sort(g(idx));
                ry = y(idx);
                ry = ry(order);
                if numel(rx) > 1
                    plot(ax, rx, ry, '-', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.5);
                end
            end
            
            % Bars + SEM + dots
            all_y = y;
            for k = 1:n_runs
                sl = per_run(k).slopes;
                mu = mean(sl);
                se = 0;
                if numel(sl) > 1
                    se = std(sl) / sqrt(numel(sl));
                end
                bar(ax, run_nums(k), mu, 0.6, 'FaceColor', color * 0.35 + 0.65, ...
                    'EdgeColor', color, 'LineWidth', 1);
                errorbar(ax, run_nums(k), mu, se, 'Color', 'k', 'LineStyle', 'none', ...
                         'LineWidth', 1, 'CapSize', 5);
                scatter(ax, run_nums(k) * ones(size(sl)), sl, 18, color, 'filled', ...
                        'MarkerEdgeColor', 'w');
                all_y = [all_y; mu + se; mu - se]; %#ok<AGROW>
            end
            
            % Linear trend line
            if ~isnan(res.trend_slope)
                xr = [min(run_nums) - 0.4, max(run_nums) + 0.4];
                plot(ax, xr, res.trend_intercept + res.trend_slope * xr, '--', ...
                     'Color', color * 0.6, 'LineWidth', 1.5);
            end
            
            % Limits and per-run stars (vs 0)
            lo = min([all_y; 0]); hi = max([all_y; 0]);
            rng = hi - lo; if rng == 0, rng = 1; end
            yl = [lo - 0.1 * rng, hi + 0.2 * rng];
            for k = 1:n_runs
                stars = sig_stars(per_run(k).p_vs_zero, options.slope_alpha);
                if ~isempty(stars)
                    text(ax, run_nums(k), yl(2) - 0.07 * diff(yl), stars, ...
                         'HorizontalAlignment', 'center', 'FontSize', 11, ...
                         'FontWeight', 'bold', 'Color', color);
                end
                text(ax, run_nums(k), yl(1) + 0.03 * diff(yl), sprintf('n=%d', numel(per_run(k).slopes)), ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                     'FontSize', 7, 'Color', [0.4, 0.4, 0.4]);
            end
            
            plot(ax, [min(run_nums) - 0.6, max(run_nums) + 0.6], [0, 0], 'k-', 'LineWidth', 0.5);
            xlim(ax, [min(run_nums) - 0.6, max(run_nums) + 0.6]);
            ylim(ax, yl);
            set(ax, 'XTick', run_nums, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
            xlabel(ax, 'Run #');
            ylabel(ax, slope_units_str(options));
            
            title(ax, {sprintf('Session %d - %s runs: slope across runs', sess, nice_labels{d}), ...
                       sprintf('ANOVA: F(%g,%g) = %.2f, p = %.3g%s   |   Trend: %.3g/run, p = %.3g%s', ...
                               res.df(1), res.df(2), res.F, res.p, sig_stars(res.p, options.slope_alpha, true), ...
                               res.trend_slope, res.trend_p, sig_stars(res.trend_p, options.slope_alpha, true))}, ...
                  'FontSize', 10, 'FontWeight', 'normal');
            
            annotation('textbox', [0.15, 0.005, 0.8, 0.05], 'String', ...
                       sprintf('%s; stars above bars = t-test vs 0; grey lines = individual mice', res.anova_type), ...
                       'EdgeColor', 'none', 'FontSize', 7, 'Color', [0.4, 0.4, 0.4], ...
                       'HorizontalAlignment', 'center');
        end
    end
end

function res = run_slope_anova(y, g, mice)
    % ANOVA of slopes across runs.
    % Preferred: Run (fixed) x Mouse (random) via anovan, which accounts for
    % repeated measures and tolerates mice missing some runs.
    % Fallback: one-way ANOVA across runs (anova1).
    % Trend: linear mixed model slope ~ run + (1|mouse); fallback: regression.
    res = struct('anova_type', 'not enough data', 'F', NaN, 'df', [NaN, NaN], 'p', NaN, ...
                 'table', {{}}, 'trend_slope', NaN, 'trend_intercept', NaN, 'trend_p', NaN, ...
                 'trend_method', 'none', 'n_obs', numel(y), 'n_mice', numel(unique(mice)), ...
                 'n_runs', numel(unique(g)));
    
    if numel(unique(g)) < 2 || numel(y) < 3
        return;
    end
    
    col = @(tbl, name) find(strcmp(tbl(1, :), name), 1);
    
    % ---- ANOVA ----
    done = false;
    if numel(unique(mice)) >= 2 && numel(y) > numel(unique(mice)) + numel(unique(g)) - 1
        try
            [p, tbl] = anovan(y, {g, mice}, 'random', 2, 'varnames', {'Run', 'Mouse'}, ...
                              'display', 'off');
            F = tbl{2, col(tbl, 'F')};
            df1 = tbl{2, col(tbl, 'd.f.')};
            c_den = col(tbl, 'd.f. denom');
            if isempty(c_den)
                df2 = tbl{end-1, col(tbl, 'd.f.')};
            else
                df2 = tbl{2, c_den};
            end
            if ~isnan(p(1)) && ~isempty(F) && ~isnan(F)
                res.anova_type = 'Two-way ANOVA (Run fixed, Mouse random)';
                res.F = F; res.df = [df1, df2]; res.p = p(1); res.table = tbl;
                done = true;
            end
        catch
        end
    end
    if ~done
        [p, tbl] = anova1(y, g, 'off');
        res.anova_type = 'One-way ANOVA across runs';
        res.F = tbl{2, col(tbl, 'F')};
        res.df = [tbl{2, col(tbl, 'df')}, tbl{3, col(tbl, 'df')}];
        res.p = p;
        res.table = tbl;
    end
    
    % ---- Linear trend across runs ----
    try
        if numel(unique(mice)) < 2
            error('need >= 2 mice for mixed model');
        end
        T = table(y, g, categorical(mice), 'VariableNames', {'slope', 'run', 'mouse'});
        lme = fitlme(T, 'slope ~ run + (1|mouse)');
        ci = strcmp(lme.Coefficients.Name, 'run');
        ii = strcmp(lme.Coefficients.Name, '(Intercept)');
        res.trend_slope = lme.Coefficients.Estimate(ci);
        res.trend_intercept = lme.Coefficients.Estimate(ii);
        res.trend_p = lme.Coefficients.pValue(ci);
        res.trend_method = 'LMM: slope ~ run + (1|mouse)';
    catch
        P = polyfit(g, y, 1);
        [R, Pv] = corrcoef(g, y);
        res.trend_slope = P(1);
        res.trend_intercept = P(2);
        res.trend_p = Pv(1, 2);
        res.trend_method = sprintf('linear regression (r = %.2f)', R(1, 2));
    end
end

% =========================================================================
% Trace processing
% =========================================================================

function [distance_grid, avg_dff, sem_dff, mice_with_run, interp_dff, mice_used] = get_processed_run_data(run_data, field_data, run_num, type_label, options)
    if strcmp(type_label, 'towards')
        distance_grid = linspace(-options.max_distance, -options.gap/2, options.distance_points);
    else
        distance_grid = linspace(options.gap/2, options.max_distance, options.distance_points);
    end
    
    [interp_dff, mice_used] = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options);
    
    if isempty(interp_dff) || all(isnan(interp_dff(:)))
        distance_grid = [];
        avg_dff = [];
        sem_dff = [];
        mice_with_run = 0;
        interp_dff = [];
        mice_used = [];
        return;
    end
    
    mice_with_run = size(interp_dff, 1);
    avg_dff = mean(interp_dff, 1, 'omitnan');
    
    if mice_with_run > 1
        valid_count = sum(~isnan(interp_dff), 1);
        valid_count(valid_count == 0) = 1;
        sem_dff = std(interp_dff, 0, 1, 'omitnan') ./ sqrt(valid_count);
    else
        sem_dff = zeros(size(avg_dff));
    end
    
    valid_dff = ~isnan(avg_dff);
    if sum(valid_dff) > options.smoothing
        avg_dff(valid_dff) = movmean(avg_dff(valid_dff), options.smoothing);
        sem_dff(valid_dff) = movmean(sem_dff(valid_dff), options.smoothing);
    end
end

function dff = get_dff(run_data, m, run_idx, options)
    % dF/F of one run, z-scored against all runs of that mouse if requested
    dff = run_data(m).runs(run_idx).dff;
    if options.use_zscore
        all_dff = [];
        for r = 1:length(run_data(m).runs)
            all_dff = [all_dff; run_data(m).runs(r).dff(:)]; %#ok<AGROW>
        end
        dff_mean = mean(all_dff, 'omitnan');
        dff_std = std(all_dff, 'omitnan');
        if dff_std < 1e-10
            dff = [];
            return;
        end
        dff = (dff - dff_mean) / dff_std;
    end
end

function [norm_distance, dff_data] = get_mouse_run_trace(run_data, m, run_idx, type_label, options)
    norm_distance = [];
    dff_data = [];
    dff = get_dff(run_data, m, run_idx, options);
    if isempty(dff)
        return;
    end
    run = run_data(m).runs(run_idx);
    run.dff = dff;
    [norm_distance, dff_data] = normalize_distance_data(run, type_label, options);
end

function [interp_dff, mice_used] = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options)
    mouse_food_runs = field_data.mouse_food_runs;
    mice_indices = field_data.mice_indices;
    interp_dff = [];
    mice_used = [];
    
    for mi = 1:length(mice_indices)
        m = mice_indices(mi);
        mouse_run_indices = mouse_food_runs{mi};
        if run_num > length(mouse_run_indices)
            continue;
        end
        
        [norm_distance, dff_data] = get_mouse_run_trace(run_data, m, mouse_run_indices(run_num), ...
                                                        type_label, options);
        if length(norm_distance) < 5
            continue;
        end
        
        interp_dff_mouse = interp1(norm_distance, dff_data, distance_grid, 'linear', NaN);
        interp_dff = [interp_dff; interp_dff_mouse(:)']; %#ok<AGROW>
        mice_used(end+1, 1) = m; %#ok<AGROW>
    end
end

function [norm_distance, dff_data] = normalize_distance_data(run, type_label, options)
    distance_data = run.distance;
    dff_data = run.dff;
    min_dist = min(distance_data);
    
    if strcmp(type_label, 'towards')
        norm_distance = -(distance_data - min_dist) - options.gap/2;
        if min(norm_distance) < -options.max_distance
            valid_indices = norm_distance >= -options.max_distance;
            norm_distance = norm_distance(valid_indices);
            dff_data = dff_data(valid_indices);
        end
        [norm_distance, sort_idx] = sort(norm_distance, 'descend');
    else
        norm_distance = (distance_data - min_dist) + options.gap/2;
        if max(norm_distance) > options.max_distance
            valid_indices = norm_distance <= options.max_distance;
            norm_distance = norm_distance(valid_indices);
            dff_data = dff_data(valid_indices);
        end
        [norm_distance, sort_idx] = sort(norm_distance, 'ascend');
    end
    
    dff_data = dff_data(sort_idx);
    [norm_distance, unique_idx] = unique(norm_distance, 'stable');
    dff_data = dff_data(unique_idx);
end

function avg_dff = process_run_data(run_data, field_data, run_num, type_label, options)
    if strcmp(type_label, 'towards')
        distance_grid = linspace(-options.max_distance, -options.gap/2, options.distance_points);
    else
        distance_grid = linspace(options.gap/2, options.max_distance, options.distance_points);
    end
    
    interp_dff = collect_mouse_data(run_data, field_data, run_num, type_label, distance_grid, options);
    
    if isempty(interp_dff) || all(isnan(interp_dff(:)))
        avg_dff = [];
    else
        avg_dff = mean(interp_dff, 1, 'omitnan');
    end
end

% =========================================================================
% Plotted data: collect + save
% =========================================================================

function plotted_data = build_plotted_data(run_data, trace_data, slope_results, all_valid_runs, ...
                                           session_nums, results, options, category_str, ...
                                           dff_lim, slope_lim)
    % Everything that is drawn, organised as
    %   plotted_data.traces.sessX_towards/away(k)  - main figure traces
    %   plotted_data.bars.sessX_towards/away(k)    - main figure slope bars
    %   plotted_data.runs_anova.sessX_towards/away - slopes-across-runs figures
    type_labels = {'towards', 'away'};
    
    plotted_data = struct();
    plotted_data.info = struct( ...
        'created', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        'food_source', options.food_source, ...
        'run_category', options.run_category, ...
        'category_str', strtrim(category_str), ...
        'zscored', options.use_zscore, ...
        'slope_units', slope_units_str(options), ...
        'distance_note', 'Trace distance is the plotted axis: towards = negative (food at -gap/2), away = positive (food at +gap/2)', ...
        'slope_note', sprintf('Slope of binned dF/F vs signed distance within %g cm of food (bin %g cm)', ...
                              options.slope_range, options.bin_width), ...
        'per_mouse_note', 'per_mouse traces are unsmoothed; mean/sem are smoothed exactly as plotted', ...
        'dff_ylim', dff_lim, ...
        'slope_ylim', slope_lim);
    plotted_data.options = options;
    plotted_data.runs = all_valid_runs(:)';
    
    empty_trace = struct('run', {}, 'distance', {}, 'mean', {}, 'sem', {}, 'n_mice', {}, ...
                         'mice_index', {}, 'mouse_id', {}, 'per_mouse', {});
    empty_bar = struct('run', {}, 'slopes', {}, 'mice_index', {}, 'mouse_id', {}, ...
                       'n', {}, 'mean', {}, 'sem', {}, 'p_vs_zero', {});
    
    for sess = session_nums
        for d = 1:2
            field_name = sprintf('sess%d_%s', sess, type_labels{d});
            
            % ---- Traces ----
            tr = empty_trace;
            if isfield(trace_data, field_name)
                for k = 1:numel(trace_data.(field_name))
                    t = trace_data.(field_name){k};
                    if isempty(t), continue; end
                    t.mouse_id = mouse_ids_from_index(run_data, t.mice_index);
                    tr(end+1) = orderfields(t, empty_trace); %#ok<AGROW>
                end
            end
            plotted_data.traces.(field_name) = tr;
            
            % ---- Bars ----
            br = empty_bar;
            if isfield(slope_results, field_name)
                for k = 1:numel(slope_results.(field_name))
                    e = slope_results.(field_name){k};
                    if isempty(e) || isempty(e.slopes), continue; end
                    sl = e.slopes(:);
                    se = 0;
                    if numel(sl) > 1, se = std(sl) / sqrt(numel(sl)); end
                    br(end+1) = struct('run', all_valid_runs(k), 'slopes', sl, ...
                                       'mice_index', e.mice(:), ...
                                       'mouse_id', {mouse_ids_from_index(run_data, e.mice)}, ...
                                       'n', numel(sl), 'mean', mean(sl), 'sem', se, ...
                                       'p_vs_zero', e.p_vs_zero); %#ok<AGROW>
                end
            end
            plotted_data.bars.(field_name) = br;
        end
    end
    
    if isfield(results, 'anova')
        plotted_data.runs_anova = results.anova;
    end
end

function ids = mouse_ids_from_index(run_data, idx)
    % Cell array of mouse ID strings for run_data indices
    ids = cell(numel(idx), 1);
    for i = 1:numel(idx)
        if isfield(run_data, 'mouse_id')
            v = run_data(idx(i)).mouse_id;
            if isnumeric(v), v = num2str(v); end
            ids{i} = char(string(v));
        else
            ids{i} = sprintf('idx%d', idx(i));
        end
    end
end

function save_plotted_data(plotted_data, options)
    base = resolve_save_base(options);
    
    % Which formats
    fmt = options.save_format;
    if ischar(fmt) || isstring(fmt)
        switch lower(char(fmt))
            case 'both', fmt = {'mat', 'xlsx'};
            case 'all',  fmt = {'mat', 'xlsx', 'csv'};
            otherwise,   fmt = {lower(char(fmt))};
        end
    end
    fmt = lower(fmt);
    unknown = setdiff(fmt, {'mat', 'xlsx', 'csv'});
    if ~isempty(unknown)
        warning('Unknown save_format(s): %s. Ignored.', strjoin(unknown, ', '));
    end
    
    fprintf('\n');
    if ismember('mat', fmt)
        mat_file = [base, '.mat'];
        save(mat_file, 'plotted_data', '-v7.3');
        fprintf('Saved %s\n', mat_file);
    end
    if ismember('xlsx', fmt)
        write_plotted_xlsx(plotted_data, [base, '.xlsx']);
    end
    if ismember('csv', fmt)
        write_plotted_csvs(plotted_data, base);
    end
end

function base = resolve_save_base(options)
    % save_path is a FOLDER unless it ends in a file extension we write.
    default_name = sprintf('dff_vs_distance_%s_%s_%s', options.food_source, ...
                           options.run_category, datestr(now, 'yyyymmdd_HHMMSS'));
    if ~isempty(options.save_name)
        default_name = char(options.save_name);
    end
    
    p = char(options.save_path);
    if isempty(p)
        folder = pwd;
        name = default_name;
    else
        [f, n, ext] = fileparts(p);
        if any(strcmpi(ext, {'.mat', '.xlsx', '.xls', '.csv'}))
            folder = f;
            name = n;
            if isempty(folder), folder = pwd; end
        else
            folder = p;          % treat the whole path as the target folder
            name = default_name;
        end
    end
    
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
    base = fullfile(folder, name);
end

function write_plotted_xlsx(plotted_data, xlsx_file)
    % Excel workbook, one tab per figure block (see header for layout)
    if exist(xlsx_file, 'file')
        delete(xlsx_file);   % avoid leftover tabs from an older file
    end
    
    opts = plotted_data.options;
    runs = plotted_data.runs;
    nR = numel(runs);
    dirs = {'towards', 'away'};
    dir_names = {'Towards', 'Away'};
    grids = {linspace(-opts.max_distance, -opts.gap/2, opts.distance_points), ...
             linspace(opts.gap/2, opts.max_distance, opts.distance_points)};
    
    tr_fields = fieldnames(plotted_data.traces);
    sessions = unique(cellfun(@(f) split_field(f), tr_fields))';
    
    % ---------------- Trace tabs ----------------
    for sess = sessions
        ncol = 2 + 2 * nR;
        h1 = [{'Direction', 'Distance (cm)'}, repmat({''}, 1, 2 * nR)];
        h2 = [{'', ''}, repmat({'Mean', 'SEM'}, 1, nR)];
        n_rows = {[{'n mice', 'Towards'}, cell(1, 2 * nR)], ...
                  [{'n mice', 'Away'},    cell(1, 2 * nR)]};
        for k = 1:nR
            h1{2 + 2*k - 1} = sprintf('Run %d', runs(k));
        end
        
        data = cell(0, ncol);
        for d = 1:2
            g = grids{d}(:);
            np = numel(g);
            block = cell(np, ncol);
            block(:, 1) = {dir_names{d}};
            block(:, 2) = num2cell(g);
            
            tr = plotted_data.traces.(sprintf('sess%d_%s', sess, dirs{d}));
            for k = 1:nR
                idx = find([tr.run] == runs(k), 1);
                if isempty(idx) || numel(tr(idx).mean) ~= np
                    continue;
                end
                block(:, 2 + 2*k - 1) = nan2empty(tr(idx).mean);
                block(:, 2 + 2*k)     = nan2empty(tr(idx).sem);
                n_rows{d}{2 + 2*k - 1} = tr(idx).n_mice;
            end
            data = [data; block]; %#ok<AGROW>
        end
        
        C = [h1; h2; n_rows{1}; n_rows{2}; data];
        writecell(C, xlsx_file, 'Sheet', sprintf('Sess%d traces', sess));
    end
    
    % ---------------- Slope tabs ----------------
    for sess = sessions
        ncol = 1 + 4 * nR;
        h1 = [{''}, repmat({''}, 1, 4 * nR)];
        h2 = [{'Direction'}, repmat({'Mean', 'SEM', 'n', 'p vs 0'}, 1, nR)];
        for k = 1:nR
            h1{1 + 4*(k-1) + 1} = sprintf('Run %d', runs(k));
        end
        
        summary = cell(2, ncol);
        per_mouse = {};
        for d = 1:2
            summary{d, 1} = dir_names{d};
            br = plotted_data.bars.(sprintf('sess%d_%s', sess, dirs{d}));
            for k = 1:nR
                idx = find([br.run] == runs(k), 1);
                if isempty(idx), continue; end
                c0 = 1 + 4*(k-1);
                summary(d, c0 + (1:4)) = nan2empty([br(idx).mean, br(idx).sem, ...
                                                    br(idx).n, br(idx).p_vs_zero])';
            end
            
            % Per-mouse slopes: rows = mice, columns = runs
            all_ids = {};
            for k = 1:numel(br)
                all_ids = [all_ids; br(k).mouse_id(:)]; %#ok<AGROW>
            end
            ids = unique(all_ids);
            block = cell(numel(ids), 2 + nR);
            block(:, 1) = {dir_names{d}};
            block(:, 2) = ids;
            for k = 1:nR
                idx = find([br.run] == runs(k), 1);
                if isempty(idx), continue; end
                [tf, loc] = ismember(br(idx).mouse_id, ids);
                vals = br(idx).slopes;
                for j = find(tf(:))'
                    block{loc(j), 2 + k} = vals(j);
                end
            end
            per_mouse = [per_mouse; block]; %#ok<AGROW>
        end
        
        pm_header = [{'Direction', 'Mouse ID'}, arrayfun(@(r) sprintf('Run %d', r), runs, ...
                                                        'UniformOutput', false)];
        w = max(ncol, 2 + nR);
        C = [pad_cols([h1; h2; summary], w); ...
             cell(1, w); ...
             pad_cols({sprintf('Per-mouse slopes (%s)', plotted_data.info.slope_units)}, w); ...
             pad_cols(pm_header, w); ...
             pad_cols(per_mouse, w)];
        writecell(C, xlsx_file, 'Sheet', sprintf('Sess%d slopes', sess));
    end
    
    % ---------------- ANOVA tab ----------------
    if isfield(plotted_data, 'runs_anova') && ~isempty(fieldnames(plotted_data.runs_anova))
        an_fields = fieldnames(plotted_data.runs_anova);
        C = {'Session', 'Direction', 'Test', 'F', 'df1', 'df2', 'p', ...
             'Trend slope per run', 'Trend p', 'Trend method', 'n obs', 'n mice', 'n runs'};
        for f = 1:numel(an_fields)
            [sess, dirn] = split_field(an_fields{f});
            r = plotted_data.runs_anova.(an_fields{f});
            row = [{sess, dirn, r.anova_type}, ...
                   nan2empty([r.F, r.df(1), r.df(2), r.p, r.trend_slope, r.trend_p])', ...
                   {r.trend_method, r.n_obs, r.n_mice, r.n_runs}];
            C = [C; row]; %#ok<AGROW>
        end
        writecell(C, xlsx_file, 'Sheet', 'Runs ANOVA');
    end
    
    fprintf('Saved %s\n', xlsx_file);
end

function c = nan2empty(v)
    % Numeric vector -> column cell array with NaNs as empty cells
    c = num2cell(v(:));
    c(isnan(v(:))) = {[]};
end

function C = pad_cols(C, w)
    if size(C, 2) < w
        C = [C, cell(size(C, 1), w - size(C, 2))];
    end
end

function write_plotted_csvs(plotted_data, base)
    % Long-format CSVs: one row per point / mouse / test
    tr_fields = fieldnames(plotted_data.traces);
    
    % ---- Trace means (as plotted) and per-mouse traces ----
    M = {}; P = {};
    for f = 1:numel(tr_fields)
        [sess, dirn] = split_field(tr_fields{f});
        tr = plotted_data.traces.(tr_fields{f});
        for k = 1:numel(tr)
            n = numel(tr(k).distance);
            M{end+1} = table(repmat(sess, n, 1), repmat({dirn}, n, 1), repmat(tr(k).run, n, 1), ...
                             tr(k).distance, tr(k).mean, tr(k).sem, repmat(tr(k).n_mice, n, 1), ...
                             'VariableNames', {'session', 'direction', 'run', 'distance', ...
                                               'mean_dff', 'sem_dff', 'n_mice'}); %#ok<AGROW>
            for j = 1:numel(tr(k).mice_index)
                P{end+1} = table(repmat(sess, n, 1), repmat({dirn}, n, 1), repmat(tr(k).run, n, 1), ...
                                 repmat(tr(k).mouse_id(j), n, 1), repmat(tr(k).mice_index(j), n, 1), ...
                                 tr(k).distance, tr(k).per_mouse(:, j), ...
                                 'VariableNames', {'session', 'direction', 'run', 'mouse_id', ...
                                                   'run_data_index', 'distance', 'dff'}); %#ok<AGROW>
            end
        end
    end
    write_if_any(M, [base, '_traces_mean.csv']);
    write_if_any(P, [base, '_traces_per_mouse.csv']);
    
    % ---- Bars: per-mouse slopes and bar summaries ----
    S = {}; B = {};
    br_fields = fieldnames(plotted_data.bars);
    for f = 1:numel(br_fields)
        [sess, dirn] = split_field(br_fields{f});
        br = plotted_data.bars.(br_fields{f});
        for k = 1:numel(br)
            n = br(k).n;
            S{end+1} = table(repmat(sess, n, 1), repmat({dirn}, n, 1), repmat(br(k).run, n, 1), ...
                             br(k).mouse_id(:), br(k).mice_index(:), br(k).slopes(:), ...
                             'VariableNames', {'session', 'direction', 'run', 'mouse_id', ...
                                               'run_data_index', 'slope'}); %#ok<AGROW>
            B{end+1} = table(sess, {dirn}, br(k).run, n, br(k).mean, br(k).sem, br(k).p_vs_zero, ...
                             'VariableNames', {'session', 'direction', 'run', 'n_mice', ...
                                               'mean_slope', 'sem_slope', 'p_vs_zero'}); %#ok<AGROW>
        end
    end
    write_if_any(S, [base, '_slopes_per_mouse.csv']);
    write_if_any(B, [base, '_slopes_summary.csv']);
    
    % ---- ANOVA / trend across runs ----
    if isfield(plotted_data, 'runs_anova')
        A = {};
        an_fields = fieldnames(plotted_data.runs_anova);
        for f = 1:numel(an_fields)
            [sess, dirn] = split_field(an_fields{f});
            r = plotted_data.runs_anova.(an_fields{f});
            A{end+1} = table(sess, {dirn}, {r.anova_type}, r.F, r.df(1), r.df(2), r.p, ...
                             r.trend_slope, r.trend_p, {r.trend_method}, r.n_obs, r.n_mice, r.n_runs, ...
                             'VariableNames', {'session', 'direction', 'anova_type', 'F', 'df1', ...
                                               'df2', 'p', 'trend_slope_per_run', 'trend_p', ...
                                               'trend_method', 'n_obs', 'n_mice', 'n_runs'}); %#ok<AGROW>
        end
        write_if_any(A, [base, '_runs_anova.csv']);
    end
end

function write_if_any(tables, file)
    if isempty(tables)
        return;
    end
    writetable(vertcat(tables{:}), file);
    fprintf('Saved %s\n', file);
end

function [sess, dirn] = split_field(field_name)
    % 'sess0_towards' -> 0, 'towards'
    tok = regexp(field_name, '^sess(\d+)_(\w+)$', 'tokens', 'once');
    sess = str2double(tok{1});
    dirn = tok{2};
end

% =========================================================================
% Decorations
% =========================================================================

function add_scale_bars(ax, options, min_dff, max_dff)
    if isempty(ax)
        return;
    end
    hold(ax, 'on');
    
    scalebar_x = -options.max_distance * 1.5;
    scalebar_y = min_dff + (max_dff - min_dff) * 0.001;
    
    plot(ax, [scalebar_x, scalebar_x + options.scalebar_distance], [scalebar_y, scalebar_y], 'k-', 'LineWidth', 2);
    plot(ax, [scalebar_x, scalebar_x], [scalebar_y, scalebar_y + options.scalebar_dff], 'k-', 'LineWidth', 2);
    
    if options.use_zscore
        lbl = sprintf('%.1f ZdF/F', options.scalebar_dff);
    else
        lbl = sprintf('%.1f dF/F', options.scalebar_dff);
    end
    text(ax, scalebar_x, scalebar_y + options.scalebar_dff * 1.05, lbl, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    text(ax, scalebar_x + options.scalebar_distance, scalebar_y, sprintf('%d cm', options.scalebar_distance), ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 7);
end

function add_session_labels(fig, layout, session_nums)
    overlay = axes('Parent', fig, 'Position', [0, 0, 1, 1], 'Visible', 'off', ...
                   'XLim', [0, 1], 'YLim', [0, 1], 'HitTest', 'off');
    hold(overlay, 'on');
    set(overlay, 'PickableParts', 'none');
    
    lc = [0.15, 0.15, 0.15];
    for s = 1:2
        y_top = layout.block_top(s);
        y_bot = layout.block_bottom(s);
        y_mid = (y_top + y_bot) / 2;
        
        bx = 0.042;
        plot(overlay, [bx, bx], [y_bot, y_top], '-', 'Color', lc, 'LineWidth', 1.5);
        plot(overlay, [bx, bx + 0.008], [y_top, y_top], '-', 'Color', lc, 'LineWidth', 1.5);
        plot(overlay, [bx, bx + 0.008], [y_bot, y_bot], '-', 'Color', lc, 'LineWidth', 1.5);
        
        text(overlay, 0.022, y_mid, sprintf('Session %d', session_nums(s)), ...
             'Rotation', 90, 'FontSize', 14, 'FontWeight', 'bold', 'Color', lc, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    plot(overlay, [0.01, 1 - layout.right], [layout.divider_y, layout.divider_y], '--', ...
         'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.75);
    
    text(overlay, 1 - layout.right, 0.955, 'Towards', 'Color', [0, 0, 0.8], ...
         'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'right');
    text(overlay, 1 - layout.right, 0.93, 'Away', 'Color', [0.8, 0, 0], ...
         'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'right');
end

function h = shadedErrorBar(x, y, errBar, lineProps, patchSaturation)
    if nargin < 4 || isempty(lineProps)
        lineProps = {'b-', 'LineWidth', 1.5};
    end
    if nargin < 5 || isempty(patchSaturation)
        patchSaturation = 0.5;
    end
    
    if isstruct(lineProps)
        lineColor = lineProps.Color;
    elseif iscell(lineProps)
        colorIdx = find(strcmpi('Color', lineProps));
        if ~isempty(colorIdx) && colorIdx < length(lineProps)
            lineColor = lineProps{colorIdx+1};
        else
            lineColor = [0, 0, 1];
        end
    elseif ischar(lineProps)
        if lineProps(1) == 'b', lineColor = [0, 0, 1];
        elseif lineProps(1) == 'r', lineColor = [1, 0, 0];
        elseif lineProps(1) == 'g', lineColor = [0, 1, 0];
        else, lineColor = [0, 0, 1];
        end
    else
        lineColor = [0, 0, 1];
    end
    
    valid = ~isnan(y) & ~isnan(errBar);
    if sum(valid) < 2
        h = [];
        return;
    end
    
    x_valid = x(valid);
    y_valid = y(valid);
    errBar_valid = errBar(valid);
    
    patchColor = lineColor * patchSaturation + (1 - patchSaturation);
    patch([x_valid fliplr(x_valid)], [y_valid+errBar_valid fliplr(y_valid-errBar_valid)], ...
          patchColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    if iscell(lineProps)
        h = plot(x_valid, y_valid, lineProps{:});
    elseif isstruct(lineProps)
        fields = fieldnames(lineProps);
        propPairs = cell(1, 2*length(fields));
        for i = 1:length(fields)
            propPairs{2*i-1} = fields{i};
            propPairs{2*i} = lineProps.(fields{i});
        end
        h = plot(x_valid, y_valid, propPairs{:});
    else
        h = plot(x_valid, y_valid, lineProps, 'LineWidth', 1.5);
    end
end