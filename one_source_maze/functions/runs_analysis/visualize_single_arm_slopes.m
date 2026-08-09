function visualize_single_arm_slopes(slope_results, options)
    % visualize_single_arm_slopes - Creates comprehensive visualizations of slope analysis for single-arm maze
    %
    % Inputs:
    %   slope_results - Output from analyze_single_arm_slopes function
    %   options - Visualization options (optional):
    %       .plot_types - Cell array of plot types to create (default: {'heatmap', 'scatter'})
    %                    {'heatmap', 'barplot', 'scatter', 'summary'}
    %       .separate_plots - Boolean, if true creates separate bar plots for towards/away (default: false)
    %                        Note: Only applies to 'barplot' type
    %       .save_plots - Whether to save plots (default: false)
    %                     When true, saves each barplot as SVG and exports per-mouse
    %                     slope data to Excel (one sheet per direction).
    %       .figure_path - Path to save figures (default: current directory)
    %       .session_colors - Cell array of colors for sessions {No Access, Food Available}
    %                        Default: {[0.2, 0.4, 0.8], [0.8, 0.2, 0.2]}

    if nargin < 2
        options = struct();
    end

    if ~isfield(options, 'plot_types')
        options.plot_types = {'heatmap','scatter'};
    end
    if ~isfield(options, 'separate_plots'),  options.separate_plots = false; end
    if ~isfield(options, 'save_plots'),      options.save_plots     = false; end
    if ~isfield(options, 'figure_path'),     options.figure_path    = './';  end

    if ~isfield(options, 'session_colors') || isempty(options.session_colors)
        options.session_colors = {[0.2, 0.4, 0.8], [0.8, 0.2, 0.2]};
    end

    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        fprintf('No slope analysis results to visualize.\n');
        return;
    end

    [slope_data, sessions, directions] = parse_single_arm_slope_results(slope_results);

    for i = 1:length(options.plot_types)
        plot_type = options.plot_types{i};
        switch plot_type
            case 'heatmap'
                create_single_arm_slope_heatmap(slope_data, sessions, directions, options);
            case 'barplot'
                if options.separate_plots
                    create_single_arm_slope_barplot_separate(slope_data, sessions, directions, slope_results, options);
                else
                    create_single_arm_slope_barplot(slope_data, sessions, directions, slope_results, options);
                end
            case 'scatter'
                create_single_arm_slope_scatter(slope_data, sessions, directions, options);
            case 'summary'
                create_single_arm_summary_plot(slope_results, options);
            otherwise
                warning('Unknown plot type: %s', plot_type);
        end
    end
end

% ==========================================================================

function [slope_data, sessions, directions] = parse_single_arm_slope_results(slope_results)
    detailed_fields = fieldnames(slope_results.detailed);
    sessions   = [];
    directions = {};

    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        parts = split(field, '_');
        sess_num  = str2double(parts{1}(5:end));
        direction = parts{2};
        if ~ismember(sess_num, sessions),   sessions(end+1)   = sess_num;   end
        if ~ismember(direction, directions), directions{end+1} = direction; end
    end

    sessions   = sort(sessions);
    directions = {'towards', 'away'};

    slope_data = struct();
    slope_data.slopes          = nan(length(sessions), length(directions));
    slope_data.p_values        = nan(length(sessions), length(directions));
    slope_data.r_squared       = nan(length(sessions), length(directions));
    slope_data.n_mice          = nan(length(sessions), length(directions));
    slope_data.is_significant  = false(length(sessions), length(directions));
    slope_data.correct_direction = false(length(sessions), length(directions));

    for i = 1:length(detailed_fields)
        field = detailed_fields{i};
        stats = slope_results.detailed.(field);
        parts = split(field, '_');
        sess_num  = str2double(parts{1}(5:end));
        direction = parts{2};
        s_idx = find(sessions == sess_num);
        d_idx = find(strcmp(directions, direction));
        if ~isempty(s_idx) && ~isempty(d_idx)
            slope_data.slopes(s_idx, d_idx)           = stats.slope;
            slope_data.p_values(s_idx, d_idx)         = stats.p_value;
            slope_data.r_squared(s_idx, d_idx)        = stats.r_squared;
            slope_data.n_mice(s_idx, d_idx)           = stats.n_mice;
            slope_data.is_significant(s_idx, d_idx)   = stats.is_significant;
            slope_data.correct_direction(s_idx, d_idx) = stats.correct_direction;
        end
    end
end

% ==========================================================================

function create_single_arm_slope_barplot_separate(slope_data, sessions, directions, slope_results, options)
    % Separate figures for towards and away; saves SVG + Excel per direction.

    direction_labels     = {'Towards', 'Away'};
    towards_colors       = {[0.6, 0.8, 1.0], [0.2, 0.4, 0.8]};
    away_colors          = {[1.0, 0.6, 0.6], [0.8, 0.2, 0.2]};
    direction_color_sets = {towards_colors, away_colors};

    % ---- Global y-limits ------------------------------------------------
    detailed_fields      = fieldnames(slope_results.detailed);
    global_individual    = [];
    global_means         = [];
    global_sems          = [];

    for d = 1:length(directions)
        for s = 1:length(sessions)
            field_name = sprintf('sess%d_%s', sessions(s), directions{d});
            if any(strcmp(detailed_fields, field_name))
                st = slope_results.detailed.(field_name);
                if isfield(st, 'mouse_slopes') && ~isempty(st.mouse_slopes)
                    global_individual = [global_individual; st.mouse_slopes];
                    global_means      = [global_means;  mean(st.mouse_slopes)];
                    global_sems       = [global_sems;   std(st.mouse_slopes)/sqrt(length(st.mouse_slopes))];
                else
                    global_means = [global_means; st.slope];
                    global_sems  = [global_sems;  st.se_slope];
                end
            end
        end
    end

    if ~isempty(global_individual)
        global_max = max([global_means + global_sems; global_individual]);
        global_min = min([global_means - global_sems; global_individual]);
    else
        global_max = max(global_means + global_sems, [], 'omitnan');
        global_min = min(global_means - global_sems, [], 'omitnan');
    end
    global_y_range = global_max - global_min;
    if global_y_range == 0 || isnan(global_y_range), global_y_range = 0.1; end
    global_ylim = [global_min - global_y_range*0.15, global_max + global_y_range*0.35];

    % ---- Session labels -------------------------------------------------
    session_labels = make_session_labels(sessions);

    % ---- One figure per direction ---------------------------------------
    for d = 1:length(directions)
        direction = directions{d};
        colors    = direction_color_sets{d};

        figure('Name', sprintf('Single-Arm Slope Analysis - %s Direction', direction_labels{d}), ...
               'Position', [100 + (d-1)*700, 100, 300, 450]);

        individual_data = cell(1, length(sessions));
        mouse_ids       = cell(1, length(sessions));
        mean_values     = zeros(1, length(sessions));
        sem_values      = zeros(1, length(sessions));

        for s = 1:length(sessions)
            field_name = sprintf('sess%d_%s', sessions(s), direction);
            if any(strcmp(detailed_fields, field_name))
                st = slope_results.detailed.(field_name);
                if isfield(st, 'mouse_slopes') && ~isempty(st.mouse_slopes)
                    individual_data{s} = st.mouse_slopes;
                    mean_values(s)     = mean(st.mouse_slopes);
                    sem_values(s)      = std(st.mouse_slopes)/sqrt(length(st.mouse_slopes));
                    if isfield(st, 'mouse_ids') && ~isempty(st.mouse_ids)
                        mouse_ids{s} = st.mouse_ids;
                    else
                        mouse_ids{s} = arrayfun(@(x) sprintf('Mouse_%d',x), ...
                            1:length(st.mouse_slopes), 'UniformOutput', false);
                    end
                else
                    mean_values(s) = st.slope;
                    sem_values(s)  = st.se_slope;
                    individual_data{s} = [];
                    mouse_ids{s}       = {};
                end
            else
                mean_values(s) = NaN; sem_values(s) = NaN;
                individual_data{s} = []; mouse_ids{s} = {};
            end
        end

        y_range = global_y_range;

        bar_width  = 0.7;
        x_positions = 1:length(sessions);
        hold on;

        mouse_positions = cell(1, length(sessions));

        for s = 1:length(sessions)
            bar(x_positions(s), mean_values(s), bar_width, ...
                'FaceColor', colors{s}, 'EdgeColor', colors{s}*0.8, ...
                'LineWidth', 1.5, 'FaceAlpha', 0.8);
            errorbar(x_positions(s), mean_values(s), sem_values(s), ...
                     'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');

            if ~isempty(individual_data{s})
                n_mice  = length(individual_data{s});
                x_jit   = repmat(x_positions(s), n_mice, 1);
                mouse_positions{s} = struct('x_pos', x_jit, ...
                                            'y_pos', individual_data{s}, ...
                                            'names', {mouse_ids{s}});
                scatter(x_jit, individual_data{s}, 50, colors{s}, 'o', ...
                        'filled', 'MarkerFaceAlpha', 0.7, ...
                        'MarkerEdgeColor', colors{s}*0.7, 'LineWidth', 1.5);
            end
        end

        % Paired connections + t-test
        if length(sessions)==2 && ~isempty(mouse_positions{1}) && ~isempty(mouse_positions{2})
            [paired0, paired1] = match_mice(mouse_positions{1}, mouse_positions{2});
            if ~isempty(paired0)
                for pi = 1:length(paired0.y)
                    plot([paired0.x(pi), paired1.x(pi)], [paired0.y(pi), paired1.y(pi)], ...
                         '-', 'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 1.0);
                end
                if length(paired0.y) >= 3
                    [~, pv, ~, tst] = ttest(paired1.y, paired0.y);
                    add_significance_bracket(pv, tst, x_positions, mean_values, sem_values, ...
                                             individual_data, y_range);
                    print_ttest_result(upper(direction), tst, pv, length(paired0.y), ...
                                       mean(paired1.y - paired0.y));
                end
            end
        end

        % Per-session significance stars
        for s = 1:length(sessions)
            if slope_data.is_significant(s, d)
                text(x_positions(s), mean_values(s)+sem_values(s)+y_range*0.3, ...
                     '★', 'FontSize', 12, 'HorizontalAlignment', 'center', ...
                     'Color', 'k', 'FontName', 'Arial Unicode MS');
            end
        end

        xlabel('Sessions', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Slope (z-scored dF/F per distance unit)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('Neural Activity Slopes: %s Direction', direction_labels{d}), ...
              'FontSize', 18, 'FontWeight', 'bold');
        set(gca, 'XTick', x_positions, 'XTickLabel', session_labels, 'FontSize', 14);
        line([0.5, length(sessions)+0.5], [0,0], 'Color', [0.3 0.3 0.3], ...
             'LineStyle', '--', 'LineWidth', 1.5);
        xlim([0.5, length(sessions)+0.5]);
        ylim(global_ylim);

        for s = 1:length(sessions)
            if ~isempty(individual_data{s})
                text(x_positions(s), global_min - global_y_range*0.08, ...
                     sprintf('n=%d', length(individual_data{s})), ...
                     'FontSize', 10, 'HorizontalAlignment', 'center', ...
                     'Color', colors{s}, 'FontWeight', 'bold');
            end
        end

        set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
        hold off;

        % ---- Save SVG ---------------------------------------------------
        if options.save_plots
            svg_name = fullfile(options.figure_path, ...
                sprintf('single_arm_slope_barplot_%s.svg', lower(direction)));
            saveas(gcf, svg_name, 'svg');
            fprintf('Saved SVG: %s\n', svg_name);
        end

        % ---- Export Excel for this direction ----------------------------
        if options.save_plots
            xlsx_name = fullfile(options.figure_path, ...
                sprintf('single_arm_slopes_%s.xlsx', lower(direction)));
            export_slopes_to_excel(xlsx_name, direction, sessions, ...
                                   session_labels, individual_data, mouse_ids);
        end
    end
end

% ==========================================================================

function create_single_arm_slope_barplot(slope_data, sessions, directions, slope_results, options)
    % Single combined figure; saves one SVG and one Excel with two sheets.

    figure('Name', 'Single-Arm Slope Analysis: Enhanced Bar Plot', ...
           'Position', [100, 100, 600, 700]);

    session_labels   = make_session_labels(sessions);
    direction_labels = {'Towards', 'Away'};

    towards_colors = {[0.6, 0.8, 1.0], [0.2, 0.4, 0.8]};
    away_colors    = {[1.0, 0.6, 0.6], [0.8, 0.2, 0.2]};

    direction_session_colors       = cell(length(directions), length(sessions));
    direction_session_colors_light = cell(length(directions), length(sessions));
    for d = 1:length(directions)
        if strcmp(directions{d}, 'towards')
            base = towards_colors;
        else
            base = away_colors;
        end
        for s = 1:length(sessions)
            direction_session_colors{d,s}       = base{s};
            direction_session_colors_light{d,s} = min(base{s}*0.8 + 0.2, 1);
        end
    end

    individual_data = cell(length(sessions), length(directions));
    mouse_ids       = cell(length(sessions), length(directions));
    mean_values     = zeros(length(sessions), length(directions));
    sem_values      = zeros(length(sessions), length(directions));
    detailed_fields = fieldnames(slope_results.detailed);

    for s = 1:length(sessions)
        for d = 1:length(directions)
            field_name = sprintf('sess%d_%s', sessions(s), directions{d});
            if any(strcmp(detailed_fields, field_name))
                st = slope_results.detailed.(field_name);
                if isfield(st, 'mouse_slopes') && ~isempty(st.mouse_slopes)
                    individual_data{s,d} = st.mouse_slopes;
                    mean_values(s,d)     = mean(st.mouse_slopes);
                    sem_values(s,d)      = std(st.mouse_slopes)/sqrt(length(st.mouse_slopes));
                    if isfield(st, 'mouse_ids') && ~isempty(st.mouse_ids)
                        mouse_ids{s,d} = st.mouse_ids;
                    else
                        mouse_ids{s,d} = arrayfun(@(x) sprintf('Mouse_%d',x), ...
                            1:length(st.mouse_slopes), 'UniformOutput', false);
                    end
                else
                    mean_values(s,d) = st.slope; sem_values(s,d) = st.se_slope;
                    individual_data{s,d} = []; mouse_ids{s,d} = {};
                end
            else
                mean_values(s,d) = NaN; sem_values(s,d) = NaN;
                individual_data{s,d} = []; mouse_ids{s,d} = {};
            end
        end
    end

    % y-range
    all_ind = [];
    for s = 1:length(sessions)
        for d = 1:length(directions)
            if ~isempty(individual_data{s,d})
                all_ind = [all_ind; individual_data{s,d}];
            end
        end
    end
    g_max = max(mean_values(:)+sem_values(:), [], 'omitnan');
    g_min = min(mean_values(:)-sem_values(:), [], 'omitnan');
    if ~isempty(all_ind)
        max_val = max(g_max, max(all_ind));
        min_val = min(g_min, min(all_ind));
    else
        max_val = g_max; min_val = g_min;
    end
    y_range = max_val - min_val;
    if y_range == 0 || isnan(y_range), y_range = 0.1; end

    bar_width   = 0.3;
    x_positions = 1:length(directions);
    x_offset    = [-0.2, 0.2];
    hold on;

    mouse_positions = cell(length(directions), length(sessions));
    bar_handles     = gobjects(1, length(sessions));

    for s = 1:length(sessions)
        for d = 1:length(directions)
            x_pos = x_positions(d) + x_offset(s);
            h = bar(x_pos, mean_values(s,d), bar_width, ...
                    'FaceColor', direction_session_colors_light{d,s}, ...
                    'EdgeColor', direction_session_colors{d,s}, ...
                    'LineWidth', 1.5, 'FaceAlpha', 0.8);
            if d == 1, bar_handles(s) = h; end
            errorbar(x_pos, mean_values(s,d), sem_values(s,d), ...
                     'k', 'LineWidth', 2, 'CapSize', 8, 'LineStyle', 'none');

            if ~isempty(individual_data{s,d})
                slopes = individual_data{s,d};
                n      = length(slopes);
                x_jit  = x_pos + (rand(n,1)-0.5)*bar_width*0;
                mouse_positions{d,s} = struct('x_pos', x_jit, 'y_pos', slopes, ...
                                              'names', {mouse_ids{s,d}});
                scatter(x_jit, slopes, 50, direction_session_colors{d,s}, 'o', ...
                        'filled', 'MarkerFaceAlpha', 0.7, ...
                        'MarkerEdgeColor', direction_session_colors{d,s}*0.7, 'LineWidth', 1.5);
            end
        end
    end

    % Paired connections + t-test per direction
    session_comparison_results = struct();
    for d = 1:length(directions)
        if length(sessions)==2 && ~isempty(mouse_positions{d,1}) && ~isempty(mouse_positions{d,2})
            [p0, p1] = match_mice(mouse_positions{d,1}, mouse_positions{d,2});
            if ~isempty(p0)
                for pi = 1:length(p0.y)
                    plot([p0.x(pi), p1.x(pi)], [p0.y(pi), p1.y(pi)], ...
                         '-', 'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 1.0);
                end
                if length(p0.y) >= 3
                    [~, pv, ~, tst] = ttest(p1.y, p0.y);

                    % Store for printing
                    dn = directions{d};
                    session_comparison_results.(dn).p_value       = pv;
                    session_comparison_results.(dn).t_stat        = tst.tstat;
                    session_comparison_results.(dn).n_pairs       = length(p0.y);
                    session_comparison_results.(dn).is_significant = pv < 0.05;
                    session_comparison_results.(dn).mean_diff     = mean(p1.y - p0.y);

                    if pv < 0.05
                        dir_max = max([mean_values(:,d)+sem_values(:,d); ...
                            (cellfun(@(c) max([c;-Inf]), individual_data(:,d)))]);
                        star_x    = x_positions(d);
                        bracket_y = dir_max + y_range*0.06;
                        bx1 = x_positions(d)+x_offset(1);
                        bx2 = x_positions(d)+x_offset(2);
                        plot([bx1 bx1 bx2 bx2], ...
                             [bracket_y-y_range*0.01, bracket_y, bracket_y, bracket_y-y_range*0.01], ...
                             'k-', 'LineWidth', 1.5);
                        star_str = pv_to_stars(pv);
                        text(star_x, dir_max+y_range*0.08, star_str, 'FontSize', 14, ...
                             'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
                        text(star_x, dir_max+y_range*0.04, pv_to_str(pv), 'FontSize', 9, ...
                             'HorizontalAlignment', 'center', 'Color', 'k');
                    end
                end
            end
        end
    end

    % Per-bar significance stars
    for d = 1:length(directions)
        for s = 1:length(sessions)
            if slope_data.is_significant(s,d)
                x_pos = x_positions(d)+x_offset(s);
                text(x_pos, mean_values(s,d)+sem_values(s,d)+y_range*0.3, ...
                     '★', 'FontSize', 12, 'HorizontalAlignment', 'center', ...
                     'Color', 'k', 'FontName', 'Arial Unicode MS');
            end
        end
    end

    xlabel('Movement Direction', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Slope (z-scored dF/F per distance unit)', 'FontSize', 16, 'FontWeight', 'bold');
    set(gca, 'XTick', x_positions, 'XTickLabel', direction_labels, 'FontSize', 14);
    line([0.5, length(directions)+0.5], [0,0], 'Color', [0.3 0.3 0.3], ...
         'LineStyle', '--', 'LineWidth', 1.5);

    % Legend with neutral scatter points
    leg_h = gobjects(1, length(sessions));
    for s = 1:length(sessions)
        c = (s==1)*[0.7 0.7 0.7] + (s==2)*[0.3 0.3 0.3];
        leg_h(s) = scatter(NaN, NaN, 100, c, 'o', 'filled', ...
                           'MarkerEdgeColor', c*0.7, 'LineWidth', 1.5);
    end
    legend(leg_h, session_labels, 'Location', 'best', 'FontSize', 14, ...
           'Box', 'off', 'EdgeColor', 'none');

    xlim([0.5, length(directions)+0.5]);
    ylim([min_val - y_range*0.15, max_val + y_range*0.30]);

    for d = 1:length(directions)
        for s = 1:length(sessions)
            if ~isempty(individual_data{s,d})
                text(x_positions(d)+x_offset(s), min_val - y_range*0.08, ...
                     sprintf('n=%d', length(individual_data{s,d})), ...
                     'FontSize', 10, 'HorizontalAlignment', 'center', ...
                     'Color', direction_session_colors{d,s}, 'FontWeight', 'bold');
            end
        end
    end

    set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
    hold off;

    % Print t-test results
    if ~isempty(fieldnames(session_comparison_results))
        fprintf('\n===== SESSION COMPARISON RESULTS (T-TESTS) =====\n');
        for dn = fieldnames(session_comparison_results)'
            r = session_comparison_results.(dn{1});
            fprintf('%s: t(%d)=%.3f, p=%.4f, mean diff=%.4f%s\n', upper(dn{1}), ...
                r.n_pairs-1, r.t_stat, r.p_value, r.mean_diff, ...
                ternary(r.is_significant, [' ' pv_to_stars(r.p_value)], ''));
        end
        fprintf('===============================================\n');
    end

    % ---- Save SVG -------------------------------------------------------
    if options.save_plots
        svg_name = fullfile(options.figure_path, 'single_arm_slope_barplot_combined.svg');
        saveas(gcf, svg_name, 'svg');
        fprintf('Saved SVG: %s\n', svg_name);

        % ---- Export Excel (one sheet per direction) --------------------
        xlsx_name = fullfile(options.figure_path, 'single_arm_slopes_combined.xlsx');
        if exist(xlsx_name, 'file'), delete(xlsx_name); end
        for d = 1:length(directions)
            export_slopes_to_excel(xlsx_name, directions{d}, sessions, session_labels, ...
                                   individual_data(:,d)', mouse_ids(:,d)');
        end
        fprintf('Saved Excel: %s\n', xlsx_name);
    end
end

% ==========================================================================
%  SHARED HELPERS
% ==========================================================================

function export_slopes_to_excel(xlsx_name, direction, sessions, session_labels, ...
                                 individual_data, mouse_ids)
% Write one Excel sheet named <direction> with columns:
%   mouse_id | slope_<session_label_0> | slope_<session_label_1> | ...
%
% Rows = union of all mouse IDs seen across sessions.
% Cells are NaN when a mouse has no data for that session.

    % Collect the full union of mouse IDs
    all_mice = {};
    for s = 1:length(sessions)
        if ~isempty(mouse_ids{s})
            all_mice = union(all_mice, mouse_ids{s});
        end
    end

    if isempty(all_mice)
        fprintf('  No individual mouse data for direction "%s" — sheet skipped.\n', direction);
        return;
    end

    % Build table: rows = mice, cols = sessions
    n_mice = length(all_mice);
    slope_matrix = nan(n_mice, length(sessions));

    for s = 1:length(sessions)
        if isempty(individual_data{s}), continue; end
        for m = 1:length(mouse_ids{s})
            row = find(strcmp(all_mice, mouse_ids{s}{m}));
            if ~isempty(row)
                slope_matrix(row, s) = individual_data{s}(m);
            end
        end
    end

    % Column names: slope_NoAccess, slope_FoodAvailable, etc.
    col_names = {'mouse_id'};
    for s = 1:length(sessions)
        col_names{end+1} = ['slope_' strrep(session_labels{s}, ' ', '_')];
    end

    T = table(all_mice(:), 'VariableNames', {'mouse_id'});
    for s = 1:length(sessions)
        T.(col_names{s+1}) = slope_matrix(:, s);
    end

    sheet_name = direction(1:min(end,31));  % Excel tab limit
    writetable(T, xlsx_name, 'Sheet', sheet_name);
    fprintf('  Wrote sheet "%s": %d mice x %d sessions\n', ...
            sheet_name, n_mice, length(sessions));
end

% --------------------------------------------------------------------------
function [p0, p1] = match_mice(data0, data1)
% Return matched (paired) x/y values for mice present in both sessions.
    p0 = struct('x', [], 'y', [], 'names', {{}});
    p1 = struct('x', [], 'y', [], 'names', {{}});
    for i = 1:length(data0.names)
        idx1 = find(strcmp(data1.names, data0.names{i}));
        if ~isempty(idx1)
            p0.x(end+1) = data0.x_pos(i);
            p0.y(end+1) = data0.y_pos(i);
            p0.names{end+1} = data0.names{i};
            p1.x(end+1) = data1.x_pos(idx1);
            p1.y(end+1) = data1.y_pos(idx1);
            p1.names{end+1} = data1.names{idx1};
        end
    end
end

% --------------------------------------------------------------------------
function add_significance_bracket(pv, tst, x_positions, mean_values, ...
                                   sem_values, individual_data, y_range)
    if pv >= 0.05, return; end
    all_vals = mean_values + sem_values;
    for s = 1:length(individual_data)
        if ~isempty(individual_data{s})
            all_vals = [all_vals, max(individual_data{s})];
        end
    end
    max_y    = max(all_vals);
    star_x   = mean(x_positions);
    bracket_y = max_y + y_range*0.06;
    plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], ...
         [bracket_y-y_range*0.01, bracket_y, bracket_y, bracket_y-y_range*0.01], ...
         'k-', 'LineWidth', 1.5);
    text(star_x, max_y+y_range*0.08, pv_to_stars(pv), 'FontSize', 14, ...
         'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
    text(star_x, max_y+y_range*0.04, pv_to_str(pv), 'FontSize', 9, ...
         'HorizontalAlignment', 'center', 'Color', 'k');
end

% --------------------------------------------------------------------------
function print_ttest_result(label, tst, pv, n, mean_diff)
    fprintf('\n===== %s DIRECTION: SESSION COMPARISON =====\n', label);
    fprintf('Paired t-test: t(%d)=%.3f, p=%.4f\n', n-1, tst.tstat, pv);
    fprintf('Mean difference: %.4f\n', mean_diff);
    fprintf('==========================================\n');
end

% --------------------------------------------------------------------------
function labels = make_session_labels(sessions)
    labels = cell(1, length(sessions));
    for i = 1:length(sessions)
        switch sessions(i)
            case 0, labels{i} = 'No Access';
            case 1, labels{i} = 'Food Available';
            otherwise, labels{i} = sprintf('Session %d', sessions(i));
        end
    end
end

% --------------------------------------------------------------------------
function s = pv_to_stars(pv)
    if pv < 0.001, s = '***';
    elseif pv < 0.01, s = '**';
    else, s = '*';
    end
end

function s = pv_to_str(pv)
    if pv < 0.001, s = 'p<0.001';
    elseif pv < 0.01, s = 'p<0.01';
    else, s = sprintf('p=%.3f', pv);
    end
end

function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end