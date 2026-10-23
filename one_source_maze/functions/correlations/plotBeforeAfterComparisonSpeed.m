function plotBeforeAfterComparisonSpeed(grouped, options)
    % Plot before vs after comparisons for each group
    
    % Create one figure per state-source combination
    for i = 1:length(options.state)
        for j = 1:length(options.source)
            state = options.state{i};
            source = options.source{j};
            field_name = [state, '_', source];
            
            if isfield(grouped, field_name) && ~isempty(grouped.(field_name))
                figure('Position', [100 + (i-1)*400, 100 + (j-1)*400, 450, 600]);
                plotSingleGroupComparisonSpeed(grouped.(field_name), state, source);
            end
        end
    end
end

function plotSingleGroupComparisonSpeed(group_data, state, source)
    % Plot before vs after comparison for a single group
    
    % Extract correlation data
    corr_before = [group_data.corr_before];
    corr_after = [group_data.corr_after];
    
    % Extract just the correlation values (first column)
    r_before = corr_before(1:2:end);
    r_after = corr_after(1:2:end);
    
    % Create paired data for plot
    paired_data = [r_before', r_after'];
    
    % Calculate means and SEMs
    means = [mean(r_before, 'omitnan'), mean(r_after, 'omitnan')];
    sems = [std(r_before, 'omitnan')/sqrt(sum(~isnan(r_before))), ...
            std(r_after, 'omitnan')/sqrt(sum(~isnan(r_after)))];
    
    % Create bar plot
    b = bar(1:2, means, 'FaceColor', 'flat');
    colors = [0.2, 0.6, 0.8; 0.8, 0.4, 0.2];
    for k = 1:2
        b.CData(k,:) = colors(k,:);
    end
    hold on;
    
    % Add error bars
    errorbar(1:2, means, sems, 'k.', 'LineWidth', 1);
    
    % Add individual data points
    jitter = 0;
    scatter(ones(size(r_before))+randn(size(r_before))*jitter, r_before, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2*ones(size(r_after))+randn(size(r_after))*jitter, r_after, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Add connecting lines for each mouse (paired data)
    for i = 1:size(paired_data, 1)
        if ~any(isnan(paired_data(i,:)))
            plot([1, 2]+randn(1,2)*jitter*0.5, paired_data(i,:), 'k-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.5]);
        end
    end
    
    % Perform statistical test (paired t-test on z-transformed values)
    z_before = [group_data.z_before];
    z_after = [group_data.z_after];
    valid_pairs = ~isnan(z_before) & ~isnan(z_after);
    if sum(valid_pairs) > 1
        [~, p_value] = ttest(z_before(valid_pairs), z_after(valid_pairs));
    else
        p_value = NaN;
    end
    
    % Add significance indicator
    sigstar([1, 2], p_value);
    
    % Format plot
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Before Discovery', 'After Discovery'});
    ylabel('Speed-DFF Pearson Correlation (r)');
    title([upper(state(1)) state(2:end), ' + ', upper(source(1)) source(2:end), ...
           ': DFF-Speed Correlation Before vs After Food Discovery'], 'FontWeight', 'bold');
    ylim([-0.5, 0.5]);
    grid on;
    box off;
end