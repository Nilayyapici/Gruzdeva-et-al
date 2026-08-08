function plot_slope_bars_by_session(slope_results, options)
% plot_slope_bars_by_session  Bar plots of average slope: Food vs Non-Food per session
%
% Produces 6 subplots in a 2-row x 3-col grid:
%
%              Before      Learning      Test
%  Towards   [Food|NF]   [Food|NF]   [Food|NF]
%  Away      [Food|NF]   [Food|NF]   [Food|NF]
%
% Each subplot:
%   - Two bars: Food arm slope and (averaged) Non-Food arm slope
%   - Error bar: SEM across mice
%   - Dots: individual per-mouse slopes (when available from mixed-effects path)
%   - Stars above each bar: one-sample t-test vs 0 (from analyze function)
%   - Bracket + stars between bars: unpaired t-test Food vs Non-Food
%
% Input:
%   slope_results  - output of analyze_combined_nonfood_slopes
%   options        - optional struct:
%       .alpha             significance level (default 0.05)
%       .food_color        RGB for food bar       (default [0.20 0.50 0.90])
%       .nonfood_color     RGB for non-food bar   (default [0.85 0.33 0.10])
%       .dot_size          scatter marker size     (default 40)
%       .bar_width         bar width               (default 0.50)
%       .show_ns_bar       show 'ns' above each bar vs-zero test (default true)
%       .show_ns_bracket   show 'ns' on Food-vs-NF bracket       (default false)
%       .fig_width         figure width  px        (default 700)
%       .fig_height        figure height px        (default 480)
%       .save_excel        save dot data to Excel  (default true)
%       .excel_file        output filename         (default 'slope_dot_data.xlsx')

    if nargin < 2, options = struct(); end

    if ~isfield(options, 'alpha'),           options.alpha           = 0.05;              end
    if ~isfield(options, 'food_color'),      options.food_color      = [0.20 0.50 0.90];  end
    if ~isfield(options, 'nonfood_color'),   options.nonfood_color   = [0.85 0.33 0.10];  end
    if ~isfield(options, 'dot_size'),        options.dot_size        = 40;                end
    if ~isfield(options, 'bar_width'),       options.bar_width       = 0.50;              end
    if ~isfield(options, 'show_ns_bar'),     options.show_ns_bar     = true;              end
    if ~isfield(options, 'show_ns_bracket'), options.show_ns_bracket = false;             end
    if ~isfield(options, 'fig_width'),       options.fig_width       = 700;               end
    if ~isfield(options, 'fig_height'),      options.fig_height      = 480;               end
    if ~isfield(options, 'save_excel'),      options.save_excel      = true;              end
    if ~isfield(options, 'excel_file'),      options.excel_file      = 'slope_dot_data.xlsx'; end

    % ------------------------------------------------------------------
    % Determine sessions present in results
    % ------------------------------------------------------------------
    detailed_fields = fieldnames(slope_results.detailed);
    if isempty(detailed_fields)
        error('slope_results.detailed is empty — run analyze_combined_nonfood_slopes first.');
    end

    sess_nums = [];
    for i = 1:numel(detailed_fields)
        tok = regexp(detailed_fields{i}, '^sess(\d+)_', 'tokens');
        if ~isempty(tok)
            sess_nums(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
        end
    end
    sess_nums   = sort(unique(sess_nums));
    sess_labels = arrayfun(@(n) sessLabel(n), sess_nums, 'UniformOutput', false);

    directions = {'towards', 'away'};
    arms       = {'food',    'nonfood'};
    arm_labels = {'Food',    'Non-Food'};
    arm_colors = {options.food_color, options.nonfood_color};

    n_sess = numel(sess_nums);
    n_dir  = numel(directions);   % always 2

    % ------------------------------------------------------------------
    % Figure
    % ------------------------------------------------------------------
    figure('Name',     'Slope Bars: Food vs Non-Food', ...
           'Position', [80 80 options.fig_width options.fig_height], ...
           'Color',    'w');

    plot_idx = 1;
    for d = 1:n_dir
        dir_str = directions{d};

        for s = 1:n_sess
            snum = sess_nums(s);
            slbl = sess_labels{s};

            ax = subplot(n_dir, n_sess, plot_idx);
            hold(ax, 'on');

            % --------------------------------------------------------------
            % Collect data for this subplot
            % --------------------------------------------------------------
            bar_h   = nan(1,2);   % bar height  = stats.slope  (authoritative)
            bar_sem = nan(1,2);   % error bar   = stats.se_slope (authoritative)
            bar_p   = nan(1,2);   % p vs 0      = stats.p_value
            bar_sig = false(1,2); % stats.is_significant
            dots    = {[], []};   % per-mouse slopes (for dots + Food-vs-NF test)

            for a = 1:2
                key = sprintf('sess%d_%s_%s', snum, dir_str, arms{a});
                if ~isfield(slope_results.detailed, key), continue; end
                st = slope_results.detailed.(key);

                % These come directly from analyze_combined_nonfood_slopes —
                % do NOT recompute them here.
                bar_h(a)   = st.slope;
                bar_sem(a) = st.se_slope;
                bar_p(a)   = st.p_value;
                bar_sig(a) = st.is_significant;

                % Per-mouse slopes only exist in the mixed-effects path.
                % Used only for: (1) drawing dots, (2) Food-vs-NF bracket test.
                if isfield(st, 'mouse_slopes') && numel(st.mouse_slopes) >= 2
                    dots{a} = st.mouse_slopes(:);
                end
            end

            % --------------------------------------------------------------
            % Food vs Non-Food bracket test
            % Use unpaired t-test: mouse_slopes order is NOT guaranteed to
            % match across arms (mice with <3 pts are dropped independently
            % per arm inside performMixedEffectsRegression).
            % --------------------------------------------------------------
            cmp_p   = NaN;
            cmp_sig = false;
            if numel(dots{1}) >= 2 && numel(dots{2}) >= 2
                [~, cmp_p] = ttest2(dots{1}, dots{2});
                cmp_sig    = cmp_p < options.alpha;
            end

            % --------------------------------------------------------------
            % Draw bars + error bars + dots + per-bar significance
            % --------------------------------------------------------------
            for a = 1:2
                if isnan(bar_h(a)), continue; end

                % Bar
                bar(ax, a, bar_h(a), options.bar_width, ...
                    'FaceColor', arm_colors{a}, ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.75);

                % SEM error bar
                if bar_sem(a) > 0
                    errorbar(ax, a, bar_h(a), bar_sem(a), ...
                        'k', 'LineWidth', 1.2, 'CapSize', 6, 'LineStyle', 'none');
                end

                % Per-mouse dots (jittered)
                if ~isempty(dots{a})
                    jitter = (rand(size(dots{a})) - 0.5) * 0.18;
                    scatter(ax, a + jitter, dots{a}, options.dot_size, ...
                        'MarkerFaceColor', arm_colors{a} * 0.6, ...
                        'MarkerEdgeColor', 'w', ...
                        'MarkerFaceAlpha', 0.9, ...
                        'LineWidth', 0.5);
                end

                % Stars / ns above bar (one-sample t-test vs 0, from analyze function)
                y_annot = bar_h(a) + bar_sem(a);
                if ~isempty(dots{a})
                    y_annot = max(y_annot, max(dots{a}));
                end
                str = sigStars(bar_p(a), bar_sig(a), options.show_ns_bar);
                if ~isempty(str)
                    text(ax, a, y_annot, str, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment',   'bottom', ...
                        'FontSize', 11, 'FontWeight', 'bold');
                end
            end

            % Zero reference line
            yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);

            % --------------------------------------------------------------
            % Comparison bracket between the two bars
            % --------------------------------------------------------------
            str = sigStars(cmp_p, cmp_sig, options.show_ns_bracket);
            if ~isempty(str) && ~any(isnan(bar_h))
                all_dots = [dots{1}; dots{2}];
                y_top = max(bar_h + bar_sem);
                if ~isempty(all_dots)
                    y_top = max(y_top, max(all_dots));
                end
                gap    = (max(abs(bar_h)) + max(bar_sem)) * 0.15 + 0.01;
                bkt_y  = y_top + gap;
                tick_h = gap * 0.35;

                plot(ax, [1 2], [bkt_y bkt_y], 'k-', 'LineWidth', 1.0);
                plot(ax, [1 1], [bkt_y-tick_h bkt_y], 'k-', 'LineWidth', 1.0);
                plot(ax, [2 2], [bkt_y-tick_h bkt_y], 'k-', 'LineWidth', 1.0);
                text(ax, 1.5, bkt_y, str, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'bottom', ...
                    'FontSize', 11, 'FontWeight', 'bold');
            end

            % --------------------------------------------------------------
            % Axes cosmetics
            % --------------------------------------------------------------
            ax.XTick      = [1 2];
            ax.XTickLabel = arm_labels;
            ax.XLim       = [0.4 2.6];
            ax.Box        = 'off';
            ax.TickDir    = 'out';
            ax.FontSize   = 10;

            title(ax, slbl, 'FontSize', 11, 'FontWeight', 'bold');

            if s == 1
                dir_cap = [upper(dir_str(1)) dir_str(2:end)];
                ylabel(ax, sprintf('%s\nSlope (z-score / dist)', dir_cap), 'FontSize', 10);
            end

            hold(ax, 'off');
            plot_idx = plot_idx + 1;
        end
    end

    sgtitle('Average Slope: Food vs Non-Food Arms', 'FontSize', 13, 'FontWeight', 'bold');

    % ------------------------------------------------------------------
    % Save per-mouse slope dots to Excel
    % ------------------------------------------------------------------
    if options.save_excel
        rows = {};
        for d = 1:n_dir
            dir_str = directions{d};
            for s = 1:n_sess
                snum = sess_nums(s);
                slbl = sess_labels{s};
                for a = 1:2
                    key = sprintf('sess%d_%s_%s', snum, dir_str, arms{a});
                    if ~isfield(slope_results.detailed, key), continue; end
                    st = slope_results.detailed.(key);
                    if ~isfield(st, 'mouse_slopes') || isempty(st.mouse_slopes), continue; end
                    ms = st.mouse_slopes(:);
                    for j = 1:numel(ms)
                        rows(end+1, :) = {slbl, dir_str, arm_labels{a}, ms(j)}; %#ok<AGROW>
                    end
                end
            end
        end
        if ~isempty(rows)
            T = cell2table(rows, 'VariableNames', {'Session','Direction','Arm','MouseSlope'});
            writetable(T, options.excel_file);
            fprintf('Slope dot data saved to: %s\n', options.excel_file);
        else
            fprintf('No per-mouse slopes available to save (all conditions used simple regression fallback).\n');
        end
    end
end


% =========================================================================
% Local helpers
% =========================================================================

function lbl = sessLabel(n)
    switch n
        case 0,    lbl = 'Before';
        case 1,    lbl = 'Learning';
        case 2,    lbl = 'Test';
        otherwise, lbl = sprintf('Session %d', n);
    end
end

function str = sigStars(p, is_sig, show_ns)
    if ~is_sig
        if show_ns, str = 'ns'; else, str = ''; end
        return
    end
    if     p < 0.001, str = '***';
    elseif p < 0.01,  str = '**';
    elseif p < 0.05,  str = '*';
    else,             str = '';
    end
end