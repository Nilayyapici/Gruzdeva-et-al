function plot_traces_by_session(run_data, options)
% plot_traces_by_session - Plots dF/F traces (Food vs Non-Food) separately per session
%
% Layout: 2 figures (Towards / Away), each with 1 x n_sessions subplots
%
%   Figure 1 - Towards:   [Before]   [Learning]   [Test]
%   Figure 2 - Away:      [Before]   [Learning]   [Test]
%
% Each subplot shows two traces: Food arm and (combined) Non-Food arm,
% with SEM shading, for that session only.
%
% Inputs:
%   run_data  - structure from analyze_mouse_runs
%   options   - optional struct:
%       .sessions            cell array e.g. {'sess0','sess1','sess2'} (default: all)
%       .food_color          RGB for food trace        (default [0.10 0.40 0.80])
%       .nonfood_color       RGB for non-food trace    (default [0.80 0.30 0.10])
%       .ylim                y-axis limits             (default [-2 2])
%       .smoothing           moving-average window     (default 5)
%       .max_distance        max distance to plot      (default 200)
%       .plot_sem            shade SEM                 (default true)
%       .line_width          trace line width          (default 2.5)
%       .axis_label_font_size                          (default 12)
%       .tick_font_size                                (default 10)
%       .subplot_title_font_size                       (default 13)
%       .x_tick_spacing      x-axis tick spacing       (default 20)
%       .y_tick_spacing      y-axis tick spacing       (default auto)
%       .fig_width           width per figure px       (default 900)
%       .fig_height          height per figure px      (default 320)
%       .save_excel          save mean+SEM to Excel    (default true)
%       .excel_file          output filename           (default 'traces_by_session.xlsx')
%       .save_pdf            save figures as vector PDF (default true)
%       .pdf_prefix          prefix for PDF filenames  (default 'traces_by_session')

    if nargin < 2, options = struct(); end

    % ---- defaults --------------------------------------------------------
    if ~isfield(options, 'sessions')
        all_sess = [];
        for m = 1:numel(run_data)
            if ~ismember(run_data(m).session, all_sess)
                all_sess(end+1) = run_data(m).session; %#ok<AGROW>
            end
        end
        all_sess = sort(all_sess);
        options.sessions = arrayfun(@(n) sprintf('sess%d',n), all_sess, 'UniformOutput', false);
    end

    if ~isfield(options, 'food_color'),             options.food_color            = [0.10 0.40 0.80]; end
    if ~isfield(options, 'nonfood_color'),          options.nonfood_color         = [0.80 0.30 0.10]; end
    if ~isfield(options, 'ylim'),                   options.ylim                  = [-2 2];            end
    if ~isfield(options, 'smoothing'),              options.smoothing             = 5;                 end
    if ~isfield(options, 'max_distance'),           options.max_distance          = 200;               end
    if ~isfield(options, 'plot_sem'),               options.plot_sem              = true;              end
    if ~isfield(options, 'line_width'),             options.line_width            = 2.5;               end
    if ~isfield(options, 'axis_label_font_size'),   options.axis_label_font_size  = 12;                end
    if ~isfield(options, 'tick_font_size'),         options.tick_font_size        = 10;                end
    if ~isfield(options, 'subplot_title_font_size'),options.subplot_title_font_size = 13;              end
    if ~isfield(options, 'x_tick_spacing'),         options.x_tick_spacing        = 20;               end
    if ~isfield(options, 'y_tick_spacing'),         options.y_tick_spacing        = [];               end
    if ~isfield(options, 'fig_width'),              options.fig_width             = 900;               end
    if ~isfield(options, 'fig_height'),             options.fig_height            = 320;               end
    if ~isfield(options, 'save_excel'),             options.save_excel            = true;              end
    if ~isfield(options, 'excel_file'),             options.excel_file            = 'traces_by_session.xlsx'; end
    if ~isfield(options, 'save_pdf'),               options.save_pdf              = true;              end
    if ~isfield(options, 'pdf_prefix'),             options.pdf_prefix            = 'traces_by_session'; end

    % ---- session metadata ------------------------------------------------
    session_numbers = [];
    for i = 1:numel(options.sessions)
        session_numbers(end+1) = str2double(options.sessions{i}(5:end)); %#ok<AGROW>
    end
    [session_numbers, si] = sort(session_numbers);
    options.sessions = options.sessions(si);

    sess_labels = arrayfun(@(n) sessLabel(n), session_numbers, 'UniformOutput', false);

    n_sess      = numel(session_numbers);
    directions  = {'towards', 'away'};
    dir_titles  = {'Towards Food', 'Away from Food'};
    arms        = {'food', 'nonfood'};
    arm_labels  = {'Food', 'Non-Food'};
    arm_colors  = {options.food_color, options.nonfood_color};

    % ---- binning ---------------------------------------------------------
    bin_width   = 1;
    dist_bins   = 0 : bin_width : options.max_distance;
    bin_centers = dist_bins(1:end-1) + bin_width/2;
    n_bins      = numel(bin_centers);

    % ---- data structures -------------------------------------------------
    % mouse_bin_data.(key){bin} = containers.Map(mouse_id -> mean_z)
    mouse_bin_data = struct();
    for s = 1:n_sess
        for d = 1:2
            for a = 1:2
                key = makeKey(session_numbers(s), directions{d}, arms{a});
                mouse_bin_data.(key) = cell(n_bins, 1);
                for i = 1:n_bins
                    mouse_bin_data.(key){i} = containers.Map();
                end
            end
        end
    end

    % ---- process mice ----------------------------------------------------
    fprintf('Processing mice...\n');
    for m = 1:numel(run_data)
        mouse_id   = run_data(m).mouse_id;
        mouse_sess = run_data(m).session;

        sess_str = sprintf('sess%d', mouse_sess);
        if ~ismember(sess_str, options.sessions), continue; end

        runs = run_data(m).runs;

        % z-score parameters from all runs of this mouse
        all_dff = [];
        for r = 1:numel(runs), all_dff = [all_dff; runs(r).dff]; end %#ok<AGROW>
        dff_mean = mean(all_dff);
        dff_std  = std(all_dff);
        if dff_std < 1e-10
            warning('Mouse %s: near-zero std, skipping.', mouse_id);
            continue;
        end

        % temporary per-mouse bin accumulator
        tmp = struct();
        for d = 1:2
            for a = 1:2
                key = makeKey(mouse_sess, directions{d}, arms{a});
                tmp.(key) = cell(n_bins, 1);
                for i = 1:n_bins, tmp.(key){i} = []; end
            end
        end

        for r = 1:numel(runs)
            run = runs(r);
            if strcmp(run.arm, 'food')
                carm = 'food';
            elseif ismember(run.arm, {'nonfood1','nonfood2'})
                carm = 'nonfood';
            else
                continue;
            end
            if ~ismember(run.type, directions), continue; end

            valid_idx = run.distance <= options.max_distance;
            if ~any(valid_idx), continue; end

            dist_f = run.distance(valid_idx);
            dff_z  = (run.dff(valid_idx) - dff_mean) / dff_std;

            key = makeKey(mouse_sess, run.type, carm);
            for i = 1:n_bins
                idx = dist_f >= dist_bins(i) & dist_f < dist_bins(i+1);
                if any(idx)
                    tmp.(key){i} = [tmp.(key){i}; dff_z(idx)];
                end
            end
        end

        % store per-mouse bin averages
        for d = 1:2
            for a = 1:2
                key = makeKey(mouse_sess, directions{d}, arms{a});
                for i = 1:n_bins
                    if ~isempty(tmp.(key){i})
                        bin_avg = mean(tmp.(key){i});
                        mp = mouse_bin_data.(key){i};
                        if mp.isKey(mouse_id)
                            mp(mouse_id) = (mp(mouse_id) + bin_avg) / 2;
                        else
                            mp(mouse_id) = bin_avg;
                        end
                        mouse_bin_data.(key){i} = mp;
                    end
                end
            end
        end
    end

    % ---- compute mean + SEM traces ---------------------------------------
    % trace_data.(key).x / .mean / .sem  (already smoothed, x flipped for towards)
    trace_data = struct();
    for s = 1:n_sess
        for d = 1:2
            for a = 1:2
                key = makeKey(session_numbers(s), directions{d}, arms{a});

                raw_mean = nan(n_bins, 1);
                raw_sem  = nan(n_bins, 1);

                for i = 1:n_bins
                    mp = mouse_bin_data.(key){i};
                    if mp.Count > 0
                        vals = cell2mat(values(mp));
                        raw_mean(i) = mean(vals);
                        raw_sem(i)  = std(vals) / sqrt(numel(vals));
                    end
                end

                valid = ~isnan(raw_mean);
                if sum(valid) > options.smoothing
                    x_v    = bin_centers(valid);
                    y_v    = raw_mean(valid);
                    sem_v  = raw_sem(valid);
                    y_sm   = movmean(y_v,   options.smoothing);
                    sem_sm = movmean(sem_v, options.smoothing);
                else
                    x_v    = bin_centers(valid);
                    y_sm   = raw_mean(valid);
                    sem_sm = raw_sem(valid);
                end

                if strcmp(directions{d}, 'towards')
                    x_v = -x_v;
                end

                trace_data.(key).x    = x_v(:);
                trace_data.(key).mean = y_sm(:);
                trace_data.(key).sem  = sem_sm(:);
            end
        end
    end

    % ---- plot ------------------------------------------------------------
    for d = 1:2
        dir_str = directions{d};

        figure('Name',     sprintf('dF/F Traces - %s', dir_titles{d}), ...
               'Position', [80 80 options.fig_width options.fig_height], ...
               'Color',    'w');

        for s = 1:n_sess
            ax = subplot(1, n_sess, s);
            hold(ax, 'on');

            for a = 1:2
                key   = makeKey(session_numbers(s), dir_str, arms{a});
                if ~isfield(trace_data, key), continue; end
                td    = trace_data.(key);
                x     = td.x;
                ymean = td.mean;
                ysem  = td.sem;

                if isempty(x), continue; end

                col = arm_colors{a};

                % SEM shading
                if options.plot_sem
                    x_patch = [x(:)', fliplr(x(:)')];
                    y_patch = [(ymean-ysem)', fliplr((ymean+ysem)')];
                    fill(ax, x_patch, y_patch, col, ...
                        'EdgeColor', 'none', 'FaceAlpha', 0.20, ...
                        'HandleVisibility', 'off');
                end

                % Mean trace
                plot(ax, x, ymean, 'Color', col, ...
                    'LineWidth', options.line_width, ...
                    'DisplayName', arm_labels{a});
            end

            % Zero line
            if strcmp(dir_str, 'towards')
                xl = [-options.max_distance 0];
            else
                xl = [0 options.max_distance];
            end
            plot(ax, xl, [0 0], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

            % Axes formatting
            xlim(ax, xl);
            ylim(ax, options.ylim);

            x_ticks = xl(1) : options.x_tick_spacing : xl(2);
            xticks(ax, x_ticks);
            if ~isempty(options.y_tick_spacing)
                yticks(ax, options.ylim(1) : options.y_tick_spacing : options.ylim(2));
            end

            ax.Box    = 'off';
            ax.TickDir = 'out';
            ax.FontSize = options.tick_font_size;

            title(ax, sess_labels{s}, 'FontSize', options.subplot_title_font_size, 'FontWeight', 'bold');

            if strcmp(dir_str, 'towards')
                xlabel(ax, 'Distance to Food (cm)', 'FontSize', options.axis_label_font_size);
            else
                xlabel(ax, 'Distance from Food (cm)', 'FontSize', options.axis_label_font_size);
            end
            if s == 1
                ylabel(ax, 'Z-scored dF/F', 'FontSize', options.axis_label_font_size);
            end

            lg = legend(ax, 'show', 'Location', 'best');
            lg.Box = 'off';
            lg.FontSize = options.tick_font_size;
            lg.ItemTokenSize = [15 8];

            hold(ax, 'off');
        end

        sgtitle(dir_titles{d}, 'FontSize', options.subplot_title_font_size + 1, 'FontWeight', 'bold');

        % Save figure as vector PDF
        if options.save_pdf
            dir_cap  = [upper(dir_str(1)) dir_str(2:end)];
            pdf_name = sprintf('%s_%s.pdf', options.pdf_prefix, dir_cap);
            exportgraphics(gcf, pdf_name, 'ContentType', 'vector');
            fprintf('Figure saved to: %s\n', pdf_name);
        end
    end

    % ---- save Excel ------------------------------------------------------
    if options.save_excel
        fprintf('Saving trace data to %s...\n', options.excel_file);

        % Delete existing file so we write fresh tabs
        if exist(options.excel_file, 'file')
            delete(options.excel_file);
        end

        for d = 1:2
            dir_str = directions{d};
            for s = 1:n_sess
                slbl = sess_labels{s};

                % Sheet name: e.g. "Towards_Before"
                dir_cap   = [upper(dir_str(1)) dir_str(2:end)];
                sheet_name = sprintf('%s_%s', dir_cap, slbl);

                % Build table: Distance | Food_Mean | Food_SEM | NonFood_Mean | NonFood_SEM
                key_food = makeKey(session_numbers(s), dir_str, 'food');
                key_nf   = makeKey(session_numbers(s), dir_str, 'nonfood');

                td_food = trace_data.(key_food);
                td_nf   = trace_data.(key_nf);

                % Use the food x-axis as reference (both should be same bins)
                x_ref = td_food.x;
                if isempty(x_ref), x_ref = td_nf.x; end
                if isempty(x_ref), continue; end

                % Interpolate non-food onto same x grid if needed
                food_mean = interp1(td_food.x, td_food.mean, x_ref, 'linear', NaN);
                food_sem  = interp1(td_food.x, td_food.sem,  x_ref, 'linear', NaN);
                nf_mean   = interp1(td_nf.x,   td_nf.mean,  x_ref, 'linear', NaN);
                nf_sem    = interp1(td_nf.x,   td_nf.sem,   x_ref, 'linear', NaN);

                dist_col_name = 'Distance_cm';
                T = table(x_ref, food_mean, food_sem, nf_mean, nf_sem, ...
                    'VariableNames', {dist_col_name, 'Food_Mean', 'Food_SEM', 'NonFood_Mean', 'NonFood_SEM'});

                writetable(T, options.excel_file, 'Sheet', sheet_name);
                fprintf('  Written sheet: %s\n', sheet_name);
            end
        end
        fprintf('Done. Excel saved to: %s\n', options.excel_file);
    end
end


% =========================================================================
% Helpers
% =========================================================================

function key = makeKey(sess_num, direction, arm)
    key = sprintf('sess%d_%s_%s', sess_num, direction, arm);
end

function lbl = sessLabel(n)
    switch n
        case 0,    lbl = 'Before';
        case 1,    lbl = 'Learning';
        case 2,    lbl = 'Test';
        otherwise, lbl = sprintf('Session%d', n);
    end
end