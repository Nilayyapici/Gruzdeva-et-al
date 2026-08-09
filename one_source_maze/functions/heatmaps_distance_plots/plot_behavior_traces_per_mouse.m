function plot_behavior_traces_per_mouse(mice_all, opts)

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'state'), opts.state = 'fasted'; end
if ~isfield(opts,'source'), opts.source = 'food'; end
if ~isfield(opts,'session'), opts.session = 1; end   % 0 or 1
if ~isfield(opts,'food_dist_thresh'), opts.food_dist_thresh = 10; end
if ~isfield(opts,'save_dir'), opts.save_dir = pwd; end

COL_TIME  = 1;
COL_SPEED = 4;
COL_DIST  = 5;
COL_FOOD  = 8;
COL_EAT   = 9;
COL_GROOM = 10;

sess_tag = sprintf('_sess%d', opts.session);

keep = strcmp(mice_all(:,2), opts.state) & ...
       strcmp(mice_all(:,3), opts.source) & ...
       contains(mice_all(:,1), sess_tag);

mice_sub = mice_all(keep,:);

if isempty(mice_sub)
    error('No matching sessions found.');
end

if ~exist(opts.save_dir, 'dir')
    mkdir(opts.save_dir);
end

for i = 1:size(mice_sub,1)

    sname = mice_sub{i,1};
    data  = mice_sub{i,4};
    data(isinf(data)) = NaN;

    time_vec = data(:,COL_TIME);
    t_sec = time_vec - time_vec(1);

    spatial_distance = data(:,COL_DIST);
    speed = data(:,COL_SPEED);

    food_visit = data(:,COL_FOOD) > 0;
    eating     = data(:,COL_EAT) > 0;
    grooming   = data(:,COL_GROOM) > 0;

    % temporal distance = min(time since food, time to food)
    at_food = spatial_distance <= opts.food_dist_thresh;

    time_since_food = NaN(size(time_vec));
    last_t = NaN;
    for t = 1:length(time_vec)
        if at_food(t)
            last_t = time_vec(t);
        end
        if ~isnan(last_t)
            time_since_food(t) = time_vec(t) - last_t;
        end
    end

    time_to_food = NaN(size(time_vec));
    next_t = NaN;
    for t = length(time_vec):-1:1
        if at_food(t)
            next_t = time_vec(t);
        end
        if ~isnan(next_t)
            time_to_food(t) = next_t - time_vec(t);
        end
    end

    tsf = time_since_food;
    ttf = time_to_food;
    tsf(isnan(tsf)) = Inf;
    ttf(isnan(ttf)) = Inf;

    temporal_distance = min(tsf, ttf);
    temporal_distance(isinf(temporal_distance)) = NaN;

    % z-score continuous traces
    spatial_z = zscore_nan(spatial_distance);
    temporal_z = zscore_nan(temporal_distance);
    speed_z = zscore_nan(speed);

    fig = figure('Name', sname, 'Position', [50 50 1300 650]);
    hold on;

    % plot traces
    plot(t_sec, spatial_z + 4, 'LineWidth', 1.2);
    plot(t_sec, temporal_z + 2, 'LineWidth', 1.2);
    plot(t_sec, speed_z,       'LineWidth', 1.0);

    % event patches
    add_event_patches(t_sec, food_visit, [-2.0 -1.4], [0.2 0.7 0.2], 0.35);
    add_event_patches(t_sec, eating,     [-3.0 -2.4], [0.9 0.5 0.1], 0.35);
    add_event_patches(t_sec, grooming,   [-4.0 -3.4], [0.6 0.2 0.8], 0.35);

    yticks([-3.7 -2.7 -1.7 0 2 4]);
    yticklabels({'Grooming','Eating','Visits','Speed','Temporal dist','Spatial dist'});

    xlabel('Time in session (s)');
    ylabel('Behavior / predictor traces');
    title(strrep(sname, '_', '\_'));

    grid on;
    box off;
    xlim([0 max(t_sec)]);

    sgtitle(sprintf('%s | %s | sess%d', opts.state, opts.source, opts.session));

    pdf_name = sprintf('behavior_traces_%s_%s_%s_sess%d.pdf', ...
        clean_name(sname), opts.state, opts.source, opts.session);

    pdf_path = fullfile(opts.save_dir, pdf_name);

    set(fig, 'PaperOrientation', 'landscape', ...
             'PaperUnits', 'normalized', ...
             'PaperPosition', [0 0 1 1]);

    print(fig, pdf_path, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved: %s\n', pdf_path);

    close(fig);
end

end


function z = zscore_nan(x)
    mu = mean(x, 'omitnan');
    sd = std(x, 'omitnan');
    if sd > 0
        z = (x - mu) ./ sd;
    else
        z = x * NaN;
    end
end


function add_event_patches(t_sec, event_vec, y_range, color_rgb, alpha_val)

    event_vec = logical(event_vec(:));
    t_sec = t_sec(:);

    d = diff([false; event_vec; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for k = 1:length(starts)
        x1 = t_sec(starts(k));
        x2 = t_sec(stops(k));

        if x2 <= x1
            x2 = x1 + median(diff(t_sec), 'omitnan');
        end

        patch([x1 x2 x2 x1], ...
              [y_range(1) y_range(1) y_range(2) y_range(2)], ...
              color_rgb, ...
              'FaceAlpha', alpha_val, ...
              'EdgeColor', 'none');
    end
end


function out = clean_name(str)
    out = regexprep(str, '[^a-zA-Z0-9_]', '_');
end