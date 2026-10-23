function plotBeforeAfterCorr_water_obj(corr_data, options)
% PLOTBEFOREAFTERCORRELATIONS Creates two separate figures comparing correlations
% before and after discovery for water and object distances
%
% Parameters:
%   corr_data - Structure with correlation data from analyzeDffDistanceCorr_water_object
%   options - Struct with visualization parameters

% Default options
if nargin < 2
    options = struct();
end

% Set default title if not provided
if ~isfield(options, 'title_prefix'), options.title_prefix = 'dF/F with Distance'; end

% Define colors for water (blue tones like in plotSingleGroupComparison)
water_colors = [0.4, 0.7, 0.9; 0.2, 0.4, 0.7];

% Define colors for object (warm tones)
object_colors = [0.9, 0.7, 0.4; 0.7, 0.4, 0.2];

% Create two figures: one for water correlations, one for object correlations
createBeforeAfterFigure('water', corr_data, water_colors, options);
createBeforeAfterFigure('object', corr_data, object_colors, options);
end

function createBeforeAfterFigure(correlation_type, corr_data, colors, options)
% Create a figure comparing before vs after for either water or object correlations

% Get the appropriate correlation data
corr_field = [correlation_type '_corr'];
if strcmp(correlation_type, 'water')
    title_text = [options.title_prefix ' to Water'];
    after_session = 'water';  % Use water session for "after" data
else
    title_text = [options.title_prefix ' to Object'];
    after_session = 'object'; % Use object session for "after" data
end

% Get correlations
before_corrs = corr_data.before.(corr_field);
after_corrs = corr_data.(after_session).(corr_field);

% Apply Fisher z-transformation for averaging
before_z = fisher_z(before_corrs);
after_z = fisher_z(after_corrs);

% Compute means and SEMs in correlation space (like plotSingleGroupComparison)
before_mean = mean(before_corrs, 'omitnan');
after_mean = mean(after_corrs, 'omitnan');
before_sem = std(before_corrs, 'omitnan') / sqrt(sum(~isnan(before_corrs)));
after_sem = std(after_corrs, 'omitnan') / sqrt(sum(~isnan(after_corrs)));

% Calculate means for plotting
means = [before_mean, after_mean];
sems = [before_sem, after_sem];

% Count valid samples
n_before = sum(~isnan(before_corrs));
n_after = sum(~isnan(after_corrs));

% Create figure with same dimensions as plotSingleGroupComparison
figure('Position', [100, 100, 300, 500]);

% Create bar plot with flat face color (same style as plotSingleGroupComparison)
b = bar(1:2, means, 'FaceColor', 'flat');

% Set colors for each bar
for k = 1:2
    b.CData(k,:) = colors(k,:);
end

hold on;

% Add error bars
errorbar(1:2, means, sems, 'k.', 'LineWidth', 1);

% Add individual data points using scatter (same style as plotSingleGroupComparison)
jitter = 0; % No jitter for cleaner look
scatter(ones(size(before_corrs))+randn(size(before_corrs))*jitter, before_corrs, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(2*ones(size(after_corrs))+randn(size(after_corrs))*jitter, after_corrs, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

% Add connecting lines for matched data points (same mouse)
before_mice = corr_data.before.mouse_ids;
after_mice = corr_data.(after_session).mouse_ids;

% Find common mice and create paired data
common_mice = intersect(before_mice, after_mice);
paired_data = [];
paired_z_data = [];

for i = 1:length(common_mice)
    mouse_id = common_mice{i};
    idx_before = find(strcmp(before_mice, mouse_id));
    idx_after = find(strcmp(after_mice, mouse_id));

    if ~isempty(idx_before) && ~isempty(idx_after) && ...
            ~isnan(before_corrs(idx_before)) && ~isnan(after_corrs(idx_after))
        paired_data(end+1, :) = [before_corrs(idx_before), after_corrs(idx_after)];
        paired_z_data(end+1, :) = [before_z(idx_before), after_z(idx_after)];
    end
end

% Plot connecting lines with transparency (same style as plotSingleGroupComparison)
for i = 1:size(paired_data, 1)
    if ~any(isnan(paired_data(i,:)))
        plot([1, 2]+randn(1,2)*jitter*0.5, paired_data(i,:), 'k-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5 0.5]);
    end
end

% Perform statistical test (paired t-test on z-transformed values if possible)
if size(paired_z_data, 1) > 1
    valid_pairs = ~isnan(paired_z_data(:,1)) & ~isnan(paired_z_data(:,2));
    if sum(valid_pairs) > 1
        [~, p_value] = ttest(paired_z_data(valid_pairs,1), paired_z_data(valid_pairs,2));
    else
        % Fallback to unpaired test
        [~, p_value] = ttest2(before_z, after_z, 'Vartype', 'unequal');
    end
else
    % Fallback to unpaired test
    [~, p_value] = ttest2(before_z, after_z, 'Vartype', 'unequal');
end

% Add significance indicator using sigstar (same as plotSingleGroupComparison)
sigstar([1, 2], p_value);

% Format plot (same style as plotSingleGroupComparison)
set(gca, 'XTick', 1:2, 'XTickLabel', {'Before Discovery', 'After Discovery'});
ylabel('Pearson Correlation (r)');
title(title_text, 'FontWeight', 'bold');
ylim([-0.3, 0.6]);
grid off;
box off;

% Print statistics
fprintf('\n----- %s Correlations (Before vs After) -----\n', upper(correlation_type(1)) + correlation_type(2:end));
fprintf('Before: Mean = %.3f, SEM = %.3f, n = %d\n', before_mean, before_sem, n_before);
fprintf('After: Mean = %.3f, SEM = %.3f, n = %d\n', after_mean, after_sem, n_after);
fprintf('Statistical test: p = %.4f\n', p_value);
fprintf('Paired data points: n = %d\n', size(paired_data, 1));

% ── Export individual mouse before/after correlations to Excel ──────────
% Build table for all mice (not just paired), NaN where missing
all_mice_export = union(before_mice, after_mice);
n_exp = length(all_mice_export);

r_before_exp  = NaN(n_exp, 1);
r_after_exp   = NaN(n_exp, 1);
fz_before_exp = NaN(n_exp, 1);
fz_after_exp  = NaN(n_exp, 1);

for i = 1:n_exp
    mouse_id = all_mice_export{i};
    idx_b = find(strcmp(before_mice, mouse_id));
    idx_a = find(strcmp(after_mice,  mouse_id));
    if ~isempty(idx_b)
        r_before_exp(i)  = before_corrs(idx_b);
        fz_before_exp(i) = before_z(idx_b);
    end
    if ~isempty(idx_a)
        r_after_exp(i)  = after_corrs(idx_a);
        fz_after_exp(i) = after_z(idx_a);
    end
end

delta_r  = r_after_exp  - r_before_exp;
delta_fz = fz_after_exp - fz_before_exp;

mouse_ids_exp = reshape(all_mice_export, n_exp, 1);

export_table = table(mouse_ids_exp, r_before_exp, r_after_exp, delta_r, ...
    fz_before_exp, fz_after_exp, delta_fz, ...
    'VariableNames', {'Mouse_ID', 'r_Before', 'r_After', 'Delta_r', ...
    'FisherZ_Before', 'FisherZ_After', 'Delta_FisherZ'});

excel_filename = sprintf('before_after_%s_%s.xlsx', ...
    correlation_type, datestr(now, 'yyyymmdd'));

writetable(export_table, excel_filename);
fprintf('Data exported to: %s\n', excel_filename);
disp(export_table);
end

function z = fisher_z(r)
% Fisher's z-transformation for correlation coefficients
z = NaN(size(r));
valid = ~isnan(r) & abs(r) < 1; % Avoid edge cases where |r| = 1
z(valid) = 0.5 * log((1 + r(valid)) ./ (1 - r(valid)));
end