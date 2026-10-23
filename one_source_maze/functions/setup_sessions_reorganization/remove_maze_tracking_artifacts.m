function mice_all = remove_maze_tracking_artifacts(mice_all, corners, opts)
% REMOVE_MAZE_TRACKING_ARTIFACTS  Drop frames whose nose and/or centroid
%   tracking falls outside the physical maze boundary. Run this FIRST,
%   before add_acceleration_column, add_head_direction_column, or
%   glm_dff_model_comparison -- those all trust the position columns to
%   be inside the arena.
%
%   mice_all = REMOVE_MAZE_TRACKING_ARTIFACTS(mice_all, corners)
%   mice_all = REMOVE_MAZE_TRACKING_ARTIFACTS(mice_all, corners, opts)
%
%   DATA COLUMNS USED (mice_all{row,4}):
%     2-x nose   3-y nose   12-x centroid   13-y centroid
%
%   MAZE BOUNDARY
%     corners is a struct with fields x1..xN, y1..yN defining the maze
%     outline as an ordered polygon (the last vertex connects back to
%     the first). N is inferred from however many xK/yK fields are
%     present, e.g.:
%       corners.x1 = 45;  corners.y1 = 10;
%       corners.x2 = 155; corners.y2 = 10;
%       ...
%       corners.x12 = 50; corners.y12 = 508;
%
%     A frame is removed if any checked point (nose and/or centroid,
%     per opts.check_nose / opts.check_centroid) lies outside this
%     polygon, or is NaN. Points exactly on the boundary are kept.
%
%   OPTS (all optional):
%     opts.check_nose     : true/false -- reject a frame if the nose
%                            falls outside the maze. (default: true)
%     opts.check_centroid : true/false -- reject a frame if the
%                            centroid falls outside the maze.
%                            (default: true)
%     opts.margin          : expand (positive) or shrink (negative) the
%                            maze polygon by this many units before
%                            testing, to allow a little slack near the
%                            walls (uses polyshape/polybuffer).
%                            (default: 0)
%     opts.example          : true/false -- plot QC figures for a
%                            handful of example sessions showing the
%                            maze outline with kept points (gray) and
%                            removed points (red) overlaid.
%                            (default: false)
%     opts.n_example         : number of example sessions to plot.
%                            (default: 5)
%     opts.example_rows      : specific row indices into mice_all to
%                            plot instead of a random pick -- pass
%                            1:size(mice_all,1) to review every
%                            session. (default: [] = auto-select)
%     opts.example_seed      : RNG seed for reproducible random example
%                            selection. (default: [])
%
%   OUTPUT
%     mice_all : same cell array, with the offending rows removed from
%                each session's data matrix (mice_all{row,4} shrinks
%                row-wise -- every column of that matrix, including any
%                columns already added by other functions, is trimmed
%                together so everything stays aligned; other columns of
%                mice_all, e.g. name/state/source, are untouched)
%
%   EXAMPLE:
%     corners.x1 = 45;  corners.y1 = 10;
%     corners.x2 = 155; corners.y2 = 10;
%     corners.x3 = 170; corners.y3 = 420;
%     corners.x4 = 210; corners.y4 = 420;
%     corners.x5 = 200; corners.y5 = 10;
%     corners.x6 = 320; corners.y6 = 7;
%     corners.x7 = 330; corners.y7 = 420;
%     corners.x8 = 365; corners.y8 = 420;
%     corners.x9 = 370; corners.y9 = 10;
%     corners.x10 = 490; corners.y10 = 2;
%     corners.x11 = 490; corners.y11 = 510;
%     corners.x12 = 50;  corners.y12 = 508;
%     mice_all = remove_maze_tracking_artifacts(mice_all, corners, ...
%         struct('example', true));

% ------------------------------------------------------------------ %
%  COLUMN INDICES                                                      %
% ------------------------------------------------------------------ %
COL_X_NOSE = 2;
COL_Y_NOSE = 3;
COL_X_CENT = 12;
COL_Y_CENT = 13;

% ------------------------------------------------------------------ %
%  PARSE OPTIONS                                                       %
% ------------------------------------------------------------------ %
if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'check_nose'),     opts.check_nose     = true;  end
if ~isfield(opts,'check_centroid'), opts.check_centroid = true;  end
if ~isfield(opts,'margin'),         opts.margin         = 0;     end
if ~isfield(opts,'example'),        opts.example        = false; end
if ~isfield(opts,'n_example'),      opts.n_example      = 15;     end
if ~isfield(opts,'example_rows'),   opts.example_rows   = [];    end
if ~isfield(opts,'example_seed'),   opts.example_seed   = [];    end

if ~opts.check_nose && ~opts.check_centroid
    error('At least one of opts.check_nose / opts.check_centroid must be true.');
end

% ------------------------------------------------------------------ %
%  BUILD MAZE POLYGON FROM corners STRUCT                             %
% ------------------------------------------------------------------ %
fn   = fieldnames(corners);
is_x = ~cellfun(@isempty, regexp(fn, '^x\d+$', 'once'));
nums = cellfun(@(s) str2double(s(2:end)), fn(is_x));
n_pts = max(nums);

poly_x = NaN(1, n_pts);
poly_y = NaN(1, n_pts);
for k = 1:n_pts
    xf = sprintf('x%d', k);
    yf = sprintf('y%d', k);
    if ~isfield(corners, xf) || ~isfield(corners, yf)
        error('corners is missing field %s or %s (expected fields x1..x%d, y1..y%d).', ...
            xf, yf, n_pts, n_pts);
    end
    poly_x(k) = corners.(xf);
    poly_y(k) = corners.(yf);
end

pgon = polyshape(poly_x, poly_y, 'Simplify', false);
if opts.margin ~= 0
    pgon = polybuffer(pgon, opts.margin);
end

% ------------------------------------------------------------------ %
%  REMOVE OUT-OF-BOUNDS FRAMES, SESSION BY SESSION                    %
% ------------------------------------------------------------------ %
n_rows = size(mice_all,1);
fprintf('Removing tracking artifacts outside the maze boundary (%d-vertex polygon) for %d sessions...\n\n', ...
    n_pts, n_rows);
fprintf('%-28s  %10s  %10s  %8s\n', 'Session', 'n_frames', 'n_removed', 'pct');
fprintf('%s\n', repmat('-',1,62));

if opts.example
    diag_data = cell(n_rows,1);
else
    diag_data = {};
end

n_removed_total = 0;
n_frames_total  = 0;
for r = 1:n_rows
    data = mice_all{r,4};
    name = mice_all{r,1};
    if isempty(data)
        continue;
    end

    x_nose = data(:, COL_X_NOSE);
    y_nose = data(:, COL_Y_NOSE);
    x_cent = data(:, COL_X_CENT);
    y_cent = data(:, COL_Y_CENT);

    bad = false(size(data,1),1);
    if opts.check_nose
        bad = bad | ~isinterior(pgon, x_nose, y_nose) | isnan(x_nose) | isnan(y_nose);
    end
    if opts.check_centroid
        bad = bad | ~isinterior(pgon, x_cent, y_cent) | isnan(x_cent) | isnan(y_cent);
    end

    n_bad = sum(bad);
    n_removed_total = n_removed_total + n_bad;
    n_frames_total  = n_frames_total + size(data,1);
    fprintf('%-28s  %10d  %10d  %7.2f%%\n', name, size(data,1), n_bad, 100*n_bad/size(data,1));

    if opts.example
        diag_data{r} = struct('name', name, ...
            'x_nose', x_nose, 'y_nose', y_nose, ...
            'x_cent', x_cent, 'y_cent', y_cent, ...
            'bad', bad);
    end

    data(bad,:)    = [];
    mice_all{r,4}  = data;
end

fprintf('\nDone. Removed %d/%d frames (%.2f%%) across %d sessions.\n', ...
    n_removed_total, n_frames_total, 100*n_removed_total/max(n_frames_total,1), n_rows);

% ------------------------------------------------------------------ %
%  OPTIONAL: EXAMPLE QC PLOT                                          %
% ------------------------------------------------------------------ %
if opts.example
    plot_example_artifact_removal(diag_data, pgon, opts);
end

end


% ====================================================================== %
%  EXAMPLE QC PLOT: maze outline + kept (gray) vs removed (red) points   %
% ====================================================================== %
function plot_example_artifact_removal(diag_data, pgon, opts)

has_diag   = ~cellfun(@isempty, diag_data);
valid_rows = find(has_diag);

if isempty(valid_rows)
    fprintf('No valid sessions available for artifact-removal example plot.\n');
    return;
end

if ~isempty(opts.example_rows)
    ex_rows = opts.example_rows(:)';
    ex_rows = ex_rows(ismember(ex_rows, valid_rows));
    if isempty(ex_rows)
        fprintf('None of opts.example_rows have diagnostic data -- skipping example plot.\n');
        return;
    end
else
    if ~isempty(opts.example_seed)
        rng(opts.example_seed);
    end
    n_ex    = min(opts.n_example, length(valid_rows));
    perm    = randperm(length(valid_rows), n_ex);
    ex_rows = valid_rows(perm)';
end

n_ex = length(ex_rows);

for i = 1:n_ex
    r = ex_rows(i);
    d = diag_data{r};

    figure('Name', sprintf('Artifact removal QC: %s', strrep(d.name,'_',' ')), ...
        'Position', [50 50 900 750]);
    hold on;

    plot(pgon, 'FaceColor',[0.85 0.92 1], 'FaceAlpha',0.35, ...
        'EdgeColor',[0.2 0.3 0.6], 'LineWidth',1.5, 'DisplayName','Maze boundary');

    good = ~d.bad;
    plot(d.x_cent(good), d.y_cent(good), '.', 'Color',[0.35 0.35 0.35], ...
        'MarkerSize',4, 'DisplayName','Kept (centroid)');
    if any(d.bad)
        plot(d.x_cent(d.bad), d.y_cent(d.bad), '.', 'Color',[0.9 0.15 0.10], ...
            'MarkerSize',9, 'DisplayName','Removed (centroid)');
    end
    plot(d.x_nose(good), d.y_nose(good), '.', 'Color',[0.55 0.75 0.95], ...
        'MarkerSize',2, 'DisplayName','Kept (nose)');
    if any(d.bad)
        plot(d.x_nose(d.bad), d.y_nose(d.bad), '.', 'Color',[1.0 0.55 0.10], ...
            'MarkerSize',6, 'DisplayName','Removed (nose)');
    end

    axis equal; grid on; box off;
    xlabel('x'); ylabel('y');
    n_bad = sum(d.bad);
    n_tot = length(d.bad);
    title(sprintf('%s: %d/%d frames removed (%.1f%%)', ...
        strrep(d.name,'_',' '), n_bad, n_tot, 100*n_bad/max(n_tot,1)), ...
        'FontSize',10, 'Interpreter','none');
    legend('Location','bestoutside', 'FontSize',8);
end

end