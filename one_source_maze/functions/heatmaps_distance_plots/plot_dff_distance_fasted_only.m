function [trace_data] = plot_dff_distance_fasted_only(runs_fasted_food, runs_fasted_gel, options)
    % PLOT_DFF_DISTANCE_FASTED_ONLY Creates a figure comparing dF/F vs distance
    % across sessions for fasted conditions only (food vs gel)
    %
    % OUTPUTS
    %   trace_data - Struct with one field per condition (fasted_food / fasted_gel).
    %                Each field is itself a struct with one field per session
    %                (sess0, sess1, …).  Each session field contains:
    %                  .towards.bin_centers  — distance bin centres (1 x B)
    %                  .towards.mean         — z-scored dF/F mean    (1 x B)
    %                  .towards.sem          — SEM                   (1 x B)
    %                  .towards.n            — number of values/bin  (1 x B)
    %                  .away   (same fields)
    %
    % All other inputs/behaviour are unchanged from the original function.

    % Default options
    if nargin < 3
        options = struct();
    end
    
    % Set default run_category
    if ~isfield(options, 'run_category')
        options.run_category = 'food';
    end
    
    % Validate run_category
    valid_categories = {'food', 'not_food', 'all'};
    if ~ismember(options.run_category, valid_categories)
        error('run_category must be one of: ''food'', ''not_food'', or ''all''');
    end
    
    switch options.run_category
        case 'food'
            type_names = {'towards', 'away'};
            category_label = 'Food Runs';
        case 'not_food'
            type_names = {'not_food_towards', 'not_food_away'};
            category_label = 'Non-Food Runs';
        case 'all'
            type_names = {'towards', 'away', 'not_food_towards', 'not_food_away'};
            category_label = 'All Runs';
    end
    
    if ~isfield(options, 'sessions')
        all_sessions = [];
        for data_set = {runs_fasted_food, runs_fasted_gel}
            run_data = data_set{1};
            if ~isempty(run_data)
                for m = 1:length(run_data)
                    if ~ismember(run_data(m).session, all_sessions)
                        all_sessions = [all_sessions, run_data(m).session];
                    end
                end
            end
        end
        sessions_to_plot_nums = sort(all_sessions);
        options.sessions = cell(1, length(sessions_to_plot_nums));
        for i = 1:length(sessions_to_plot_nums)
            options.sessions{i} = ['sess' num2str(sessions_to_plot_nums(i))];
        end
    end
    
    if ~isfield(options, 'conditions')
        options.conditions = {'fasted_food', 'fasted_gel'};
    end

    if ~isfield(options, 'ylim'),          options.ylim          = [-1.5 1.5]; end
    if ~isfield(options, 'xlim_towards') && ~isfield(options, 'xlim_away')
        if ~isfield(options, 'xlim'),      options.xlim          = [0 210];    end
        options.xlim_towards = options.xlim;
        options.xlim_away    = options.xlim;
    else
        if ~isfield(options, 'xlim_towards')
            options.xlim_towards = getfield_default(options, 'xlim', [0 210]);
        end
        if ~isfield(options, 'xlim_away')
            options.xlim_away = getfield_default(options, 'xlim', [0 210]);
        end
    end
    if ~isfield(options, 'smoothing'),     options.smoothing     = 5;          end
    if ~isfield(options, 'figure_size'),   options.figure_size   = [1200 600]; end
    if ~isfield(options, 'title')
        options.title = sprintf('dF/F vs Distance - Fasted Conditions (%s)', category_label);
    end
    if ~isfield(options, 'plot_sem'),      options.plot_sem      = true;       end

    session_map     = containers.Map();
    session_numbers = [];
    for i = 1:length(options.sessions)
        session_str = options.sessions{i};
        session_num = str2double(session_str(5:end));
        session_map(session_str) = i;
        session_numbers = [session_numbers, session_num];
    end
    [session_numbers, sort_idx] = sort(session_numbers);
    options.sessions = options.sessions(sort_idx);
    
    fprintf('Comparing dF/F across sessions: ');
    for i = 1:length(options.sessions), fprintf('%s ', options.sessions{i}); end
    fprintf('\nFor fasted conditions: ');
    for i = 1:length(options.conditions), fprintf('%s ', options.conditions{i}); end
    fprintf('\nRun category: %s\n', options.run_category);
    fprintf('X-axis limits: Towards [%.0f %.0f], Away [%.0f %.0f]\n', ...
            options.xlim_towards(1), options.xlim_towards(2), ...
            options.xlim_away(1), options.xlim_away(2));
    
    session_labels = cell(1, length(session_numbers));
    for i = 1:length(session_numbers)
        session_labels{i} = sprintf('Session %d', session_numbers(i));
    end
    
    figure('Name', sprintf('dF/F vs Distance - Fasted Conditions (%s)', category_label), ...
           'Position', [100, 100, options.figure_size(1), options.figure_size(2)]);
    
    if strcmp(options.run_category, 'all')
        direction_titles = {'Towards (All)', 'Away (All)'};
    else
        direction_titles = {'Towards', 'Away'};
    end
    
    default_colors = [
        0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980;
        0.9290, 0.6940, 0.1250;
        0.4940, 0.1840, 0.5560;
        0.4660, 0.6740, 0.1880;
        0.3010, 0.7450, 0.9330;
        0.6350, 0.0780, 0.1840;
    ];
    
    if ~isfield(options, 'session_colors')
        options.session_colors = cell(1, length(options.sessions));
        for i = 1:length(options.sessions)
            color_idx = mod(i-1, size(default_colors, 1)) + 1;
            options.session_colors{i} = default_colors(color_idx, :);
        end
    end
    
    condition_data_map = containers.Map(...
        {'fasted_food', 'fasted_gel'}, ...
        {runs_fasted_food, runs_fasted_gel});
    
    max_distance = 0;
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        if ~isKey(condition_data_map, cond) || isempty(condition_data_map(cond)), continue; end
        run_data = condition_data_map(cond);
        for m = 1:length(run_data)
            sess_str = ['sess' num2str(run_data(m).session)];
            if ~ismember(sess_str, options.sessions), continue; end
            for r = 1:length(run_data(m).runs)
                if ~ismember(run_data(m).runs(r).type, type_names), continue; end
                max_distance = max(max_distance, max(run_data(m).runs(r).distance));
            end
        end
    end
    
    bin_width   = 1;
    dist_bins   = 0:bin_width:ceil(max_distance + 5);
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    
    % Initialise binned-data accumulator
    all_data = struct();
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        for s = 1:length(session_numbers)
            sess = session_numbers(s);
            if strcmp(options.run_category, 'all')
                for direction = {'towards', 'away'}
                    key = sprintf('%s_sess%d_%s', cond, sess, direction{1});
                    all_data.(key) = cell(length(dist_bins)-1, 1);
                    for i = 1:length(dist_bins)-1, all_data.(key){i} = []; end
                end
            else
                for t = 1:length(type_names)
                    key = sprintf('%s_sess%d_%s', cond, sess, type_names{t});
                    all_data.(key) = cell(length(dist_bins)-1, 1);
                    for i = 1:length(dist_bins)-1, all_data.(key){i} = []; end
                end
            end
        end
    end
    
    % Fill binned data (z-scored per mouse)
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        if ~isKey(condition_data_map, cond) || isempty(condition_data_map(cond))
            warning('No data for condition %s', cond); continue;
        end
        run_data = condition_data_map(cond);
        for m = 1:length(run_data)
            mouse_id   = run_data(m).mouse_id;
            mouse_sess = run_data(m).session;
            sess_str   = ['sess' num2str(mouse_sess)];
            if ~ismember(sess_str, options.sessions), continue; end
            runs = run_data(m).runs;
            all_dff_values = [];
            for r = 1:length(runs)
                if ismember(runs(r).type, type_names)
                    all_dff_values = [all_dff_values; runs(r).dff];
                end
            end
            if isempty(all_dff_values), continue; end
            dff_mean = mean(all_dff_values);
            dff_std  = std(all_dff_values);
            if dff_std < 1e-10
                warning('Mouse %s has near-zero dF/F std, skipping.', mouse_id); continue;
            end
            for r = 1:length(runs)
                run = runs(r);
                if ~ismember(run.type, type_names), continue; end
                z_scored_dff = (run.dff - dff_mean) / dff_std;
                if strcmp(options.run_category, 'all')
                    if contains(run.type, 'towards')
                        data_key = sprintf('%s_sess%d_towards', cond, mouse_sess);
                    else
                        data_key = sprintf('%s_sess%d_away', cond, mouse_sess);
                    end
                else
                    data_key = sprintf('%s_sess%d_%s', cond, mouse_sess, run.type);
                end
                for i = 1:length(dist_bins)-1
                    indices = run.distance >= dist_bins(i) & run.distance < dist_bins(i+1);
                    if any(indices)
                        all_data.(data_key){i} = [all_data.(data_key){i}; z_scored_dff(indices)];
                    end
                end
            end
        end
    end

    % ------------------------------------------------------------------
    % Build trace_data output — mean, SEM, n per bin for every
    % condition x session x direction combination
    % ------------------------------------------------------------------
    trace_data = struct();
    directions = {'towards', 'away'};

    for cond_idx = 1:length(options.conditions)
        cond     = options.conditions{cond_idx};
        cond_key = strrep(cond, ' ', '_');   % safe field name

        for s = 1:length(session_numbers)
            sess     = session_numbers(s);
            sess_key = sprintf('sess%d', sess);

            for di = 1:2
                dir = directions{di};

                % Determine the all_data key for this combination
                if strcmp(options.run_category, 'food') || strcmp(options.run_category, 'all')
                    akey = sprintf('%s_sess%d_%s', cond, sess, dir);
                else  % not_food
                    akey = sprintf('%s_sess%d_not_food_%s', cond, sess, dir);
                end

                n_bins = length(bin_centers);
                mu  = nan(1, n_bins);
                sem = nan(1, n_bins);
                nn  = zeros(1, n_bins);

                if isfield(all_data, akey)
                    for i = 1:n_bins
                        vals = all_data.(akey){i};
                        if ~isempty(vals)
                            mu(i)  = mean(vals);
                            sem(i) = std(vals) / sqrt(length(vals));
                            nn(i)  = length(vals);
                        end
                    end
                end

                % Keep raw (unsmoothed) values before applying smoothing
                mu_raw  = mu;
                sem_raw = sem;

                % Apply same smoothing as the plot
                valid = ~isnan(mu);
                if sum(valid) > options.smoothing
                    mu_sm  = nan(1, n_bins);
                    sem_sm = nan(1, n_bins);
                    mu_sm(valid)  = movmean(mu(valid),  options.smoothing);
                    sem_sm(valid) = movmean(sem(valid), options.smoothing);
                    mu  = mu_sm;
                    sem = sem_sm;
                end

                trace_data.(cond_key).(sess_key).(dir).bin_centers = bin_centers;
                trace_data.(cond_key).(sess_key).(dir).mean        = mu;
                trace_data.(cond_key).(sess_key).(dir).sem         = sem;
                trace_data.(cond_key).(sess_key).(dir).mean_raw    = mu_raw;
                trace_data.(cond_key).(sess_key).(dir).sem_raw     = sem_raw;
                trace_data.(cond_key).(sess_key).(dir).n           = nn;
            end
        end
    end

    % ------------------------------------------------------------------
    % Plotting (unchanged from original)
    % ------------------------------------------------------------------
    num_rows = length(options.conditions);
    
    for cond_idx = 1:length(options.conditions)
        cond = options.conditions{cond_idx};
        
        for direction_idx = 1:2
            if direction_idx == 1
                plot_type   = 'towards';
                type_suffix = 'towards';
            else
                plot_type   = 'away';
                type_suffix = 'away';
            end
            
            sp_idx = (cond_idx-1)*2 + direction_idx;
            subplot(num_rows, 2, sp_idx);
            hold on;
            
            for s = 1:length(session_numbers)
                sess = session_numbers(s);
                
                if strcmp(options.run_category, 'food') || strcmp(options.run_category, 'all')
                    key = sprintf('%s_sess%d_%s', cond, sess, type_suffix);
                else
                    key = sprintf('%s_sess%d_not_food_%s', cond, sess, type_suffix);
                end
                
                if ~isfield(all_data, key), continue; end
                
                means_vec = nan(length(bin_centers), 1);
                sems_vec  = nan(length(bin_centers), 1);
                
                for i = 1:length(dist_bins)-1
                    bin_data = all_data.(key){i};
                    if ~isempty(bin_data)
                        means_vec(i) = mean(bin_data);
                        sems_vec(i)  = std(bin_data) / sqrt(length(bin_data));
                    end
                end
                
                valid = ~isnan(means_vec);
                if sum(valid) > options.smoothing
                    x_valid   = bin_centers(valid);
                    y_valid   = means_vec(valid);
                    sem_valid = sems_vec(valid);
                    y_smoothed   = movmean(y_valid,   options.smoothing);
                    sem_smoothed = movmean(sem_valid,  options.smoothing);
                    sess_label   = session_labels{s};
                    
                    if options.plot_sem
                        for i = 1:length(x_valid)-1
                            x_patch = [x_valid(i), x_valid(i+1), x_valid(i+1), x_valid(i)];
                            y_patch = [y_smoothed(i)   - sem_smoothed(i), ...
                                       y_smoothed(i+1) - sem_smoothed(i+1), ...
                                       y_smoothed(i+1) + sem_smoothed(i+1), ...
                                       y_smoothed(i)   + sem_smoothed(i)];
                            patch(x_patch, y_patch, options.session_colors{s}, ...
                                  'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
                                  'HandleVisibility', 'off');
                        end
                    end
                    plot(x_valid, y_smoothed, 'Color', options.session_colors{s}, ...
                         'LineWidth', 2, 'DisplayName', sess_label);
                end
            end
            
            cond_parts   = strsplit(cond, '_');
            cond_display = [upper(cond_parts{1}(1)) cond_parts{1}(2:end), ' + ', ...
                            upper(cond_parts{2}(1)) cond_parts{2}(2:end)];
            
            title(sprintf('%s - %s', cond_display, direction_titles{direction_idx}));
            xlabel('Distance from Food', 'FontSize', 14);
            ylabel('Z-scored dF/F',      'FontSize', 14);
            ylim(options.ylim);
            yticks(-10:0.5:10);
            
            if strcmp(plot_type, 'towards')
                xlim(options.xlim_towards);
                set(gca, 'XDir', 'reverse');
                xticks_current = get(gca, 'XTick');
                set(gca, 'XTickLabel', arrayfun(@(x) sprintf('-%g', x), ...
                    xticks_current, 'UniformOutput', false));
            else
                xlim(options.xlim_away);
            end
            
            line(get(gca, 'XLim'), [0 0], 'Color', 'k', 'LineStyle', '--', ...
                 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            add_session_comparison_stats(all_data, cond, plot_type, ...
                options.run_category, session_numbers, bin_centers, options);
            
            legend('show', 'Location', 'best');
            legend('boxoff');
            grid off;
            box off;
        end
    end
    
    sgtitle(options.title, 'FontSize', 14);
    
    filename = sprintf('dff_distance_fasted_%s.pdf', options.run_category);
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('Saved plot to: %s\n', filename);

    % ------------------------------------------------------------------
    % Export trace_data to Excel — one tab per subplot
    % Tab name format: <condition>_<direction>  e.g. fasted_food_towards
    % Columns: distance | mean_sess0 | sem_sess0 | n_sess0 | mean_sess1 | …
    % ------------------------------------------------------------------
    xlsx_name = sprintf('dff_distance_fasted_%s.xlsx', options.run_category);

    % Delete existing file so we start clean (writetable appends by default)
    if exist(xlsx_name, 'file'), delete(xlsx_name); end

    for cond_idx = 1:length(options.conditions)
        cond     = options.conditions{cond_idx};
        cond_key = strrep(cond, ' ', '_');

        for di = 1:2
            dir = directions{di};

            % Build sheet name — Excel tab limit is 31 chars
            sheet_name = sprintf('%s_%s', cond, dir);
            sheet_name = sheet_name(1:min(end, 31));

            % Start with bin_centers as the first column
            T = table(bin_centers(:), 'VariableNames', {'distance'});

            for s = 1:length(session_numbers)
                sess     = session_numbers(s);
                sess_key = sprintf('sess%d', sess);

                td = trace_data.(cond_key).(sess_key).(dir);

                T.(sprintf('mean_sess%d',         sess)) = td.mean_raw(:);
                T.(sprintf('sem_sess%d',          sess)) = td.sem_raw(:);
                T.(sprintf('mean_smoothed_sess%d', sess)) = td.mean(:);
                T.(sprintf('sem_smoothed_sess%d',  sess)) = td.sem(:);
                T.(sprintf('n_sess%d',            sess)) = td.n(:);
            end

            writetable(T, xlsx_name, 'Sheet', sheet_name);
            fprintf('  Wrote sheet "%s"  (%d rows, %d sessions)\n', ...
                    sheet_name, height(T), length(session_numbers));
        end
    end

    fprintf('Excel saved to: %s\n', xlsx_name);
end

% --------------------------------------------------------------------------
function v = getfield_default(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end

% --------------------------------------------------------------------------
function add_session_comparison_stats(all_data, cond, plot_type, run_category, session_numbers, bin_centers, options)
    if length(session_numbers) ~= 2, return; end
    sess1 = session_numbers(1);
    sess2 = session_numbers(2);
    if strcmp(run_category, 'not_food')
        key1 = sprintf('%s_sess%d_not_food_%s', cond, sess1, plot_type);
        key2 = sprintf('%s_sess%d_not_food_%s', cond, sess2, plot_type);
    else
        key1 = sprintf('%s_sess%d_%s', cond, sess1, plot_type);
        key2 = sprintf('%s_sess%d_%s', cond, sess2, plot_type);
    end
    if ~isfield(all_data, key1) || ~isfield(all_data, key2), return; end
    
    p_values  = [];
    valid_bins = [];
    for i = 1:length(bin_centers)
        data1 = all_data.(key1){i};
        data2 = all_data.(key2){i};
        if length(data1) >= 3 && length(data2) >= 3
            [~, p] = ttest2(data1, data2);
            if ~isnan(p)
                p_values  = [p_values,  p];
                valid_bins = [valid_bins, i];
            end
        end
    end
    
    if ~isempty(p_values)
        [~, ~, ~, adj_p] = fdr_bh(p_values, 0.05, 'pdep');
        sig_bins = valid_bins(adj_p < 0.05);
        if ~isempty(sig_bins)
            ylims  = ylim;
            y_stat = ylims(2) - 0.1 * (ylims(2) - ylims(1));
            bin_groups = group_consecutive_bins(sig_bins);
            for g = 1:length(bin_groups)
                group   = bin_groups{g};
                x_start = bin_centers(group(1))   - 0.5;
                x_end   = bin_centers(group(end))  + 0.5;
                line([x_start, x_end], [y_stat, y_stat], 'Color', 'black', ...
                     'LineWidth', 3, 'HandleVisibility', 'off');
            end
        end
    end
end

% --------------------------------------------------------------------------
function bin_groups = group_consecutive_bins(sig_bins)
    bin_groups = {};
    if isempty(sig_bins), return; end
    current_group = sig_bins(1);
    for i = 2:length(sig_bins)
        if sig_bins(i) == sig_bins(i-1) + 1
            current_group = [current_group, sig_bins(i)];
        else
            bin_groups{end+1} = current_group;
            current_group = sig_bins(i);
        end
    end
    bin_groups{end+1} = current_group;
end

% --------------------------------------------------------------------------
function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)
    if nargin < 4, report = 'no'; end
    if nargin < 3, method = 'pdep'; end
    if nargin < 2, q = 0.05; end
    pvals = pvals(:);
    s = length(pvals);
    if s == 0
        h = []; crit_p = []; adj_ci_cvrg = []; adj_p = []; return;
    end
    [pvals_sorted, sort_ids] = sort(pvals);
    if strcmp(method, 'pdep')
        crit_vals = (1:s) * q / s;
    else
        crit_vals = (1:s) * q / (s * sum(1./(1:s)));
    end
    h = pvals_sorted <= crit_vals';
    max_id = find(h, 1, 'last');
    if isempty(max_id)
        crit_p = 0;
        h      = false(size(pvals));
        adj_p  = ones(size(pvals));
    else
        crit_p = pvals_sorted(max_id);
        h      = pvals <= crit_p;
        adj_p_sorted = zeros(size(pvals_sorted));
        adj_p_sorted(s) = pvals_sorted(s);
        for i = (s-1):-1:1
            adj_p_sorted(i) = min(adj_p_sorted(i+1), pvals_sorted(i) * s / i);
        end
        adj_p = zeros(size(pvals));
        adj_p(sort_ids) = adj_p_sorted;
    end
    adj_ci_cvrg = 1 - q;
end