function plotBeforeAfterBothSources(grouped, state)
    % Plot before vs after comparisons for both sources (food and gel) for a single state
    % 
    % Parameters:
    %   grouped: struct with grouped data
    %   state: string, 'fasted' or 'fed'
    %   options: struct with options
    
    % Create figure
    figure('Position', [100, 100, 500, 600]);
    
    % Extract data for each source
    food_data = grouped.([state, '_food']);
    gel_data = grouped.([state, '_gel']);
    
    % Extract correlation data
    % Food data
    food_corr_before = [food_data.corr_before];
    food_corr_after = [food_data.corr_after];
    food_r_before = food_corr_before(1:2:end);
    food_r_after = food_corr_after(1:2:end);
    
    % Gel data
    gel_corr_before = [gel_data.corr_before];
    gel_corr_after = [gel_data.corr_after];
    gel_r_before = gel_corr_before(1:2:end);
    gel_r_after = gel_corr_after(1:2:end);
    
    % Calculate means
    means = [mean(food_r_before, 'omitnan'), mean(food_r_after, 'omitnan'), ...
             mean(gel_r_before, 'omitnan'), mean(gel_r_after, 'omitnan')];
    
    % Calculate SEMs
    sems = [std(food_r_before, 'omitnan')/sqrt(sum(~isnan(food_r_before))), ...
            std(food_r_after, 'omitnan')/sqrt(sum(~isnan(food_r_after))), ...
            std(gel_r_before, 'omitnan')/sqrt(sum(~isnan(gel_r_before))), ...
            std(gel_r_after, 'omitnan')/sqrt(sum(~isnan(gel_r_after)))];
    
    % Create grouped bar plot with gaps
    % Create x positions for bars with gap between groups
    bar_width = 0.8;
    group_gap = 1;
    x_pos = [1, 2, 3+group_gap, 4+group_gap];
    
    % Create bar plot
    b = bar(x_pos, means, bar_width, 'FaceColor', 'flat');
    
    % Set colors for bars
    colors = [0.4, 0.7, 0.9;    % Light blue for Before Food
              0.2, 0.4, 0.7;    % Dark blue for Food After
              0.4, 0.7, 0.9;    % Light blue for Gel Before
              0.2, 0.4, 0.7];   % Dark blue for Gel After
    
    for k = 1:4
        b.CData(k,:) = colors(k,:);
    end
    
    hold on;
    
    % Add error bars
    errorbar(x_pos, means, sems, 'k.', 'LineWidth', 0.8);
    
    % Add individual data points and connecting lines
    jitter = 0;
    
    % Food data points and connections
    scatter(ones(size(food_r_before))+randn(size(food_r_before))*jitter, food_r_before, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(2*ones(size(food_r_after))+randn(size(food_r_after))*jitter, food_r_after, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Connect food before-after points for each mouse
    food_paired_data = [food_r_before', food_r_after'];
    for i = 1:size(food_paired_data, 1)
        if ~any(isnan(food_paired_data(i,:)))
            plot([1, 2]+randn(1,2)*jitter*0.5, food_paired_data(i,:), 'k-', 'LineWidth', 0.75, 'Color', [0.5 0.5 0.5 0.7]);
        end
    end
    
    % Gel data points and connections
    scatter((3+group_gap)*ones(size(gel_r_before))+randn(size(gel_r_before))*jitter, gel_r_before, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter((4+group_gap)*ones(size(gel_r_after))+randn(size(gel_r_after))*jitter, gel_r_after, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Connect gel before-after points for each mouse
    gel_paired_data = [gel_r_before', gel_r_after'];
    for i = 1:size(gel_paired_data, 1)
        if ~any(isnan(gel_paired_data(i,:)))
            plot([3+group_gap, 4+group_gap]+randn(1,2)*jitter*0.5, gel_paired_data(i,:), 'k-', 'LineWidth', 0.75, 'Color', [0.5 0.5 0.5 0.7]);
        end
    end
    
    % Perform statistical tests (paired t-test on z-transformed values) for food
    food_z_before = [food_data.z_before];
    food_z_after = [food_data.z_after];
    food_valid_pairs = ~isnan(food_z_before) & ~isnan(food_z_after);
    if sum(food_valid_pairs) > 1
        [~, food_p_value] = ttest(food_z_before(food_valid_pairs), food_z_after(food_valid_pairs));
    else
        food_p_value = NaN;
    end
    
    % Perform statistical tests for gel
    gel_z_before = [gel_data.z_before];
    gel_z_after = [gel_data.z_after];
    gel_valid_pairs = ~isnan(gel_z_before) & ~isnan(gel_z_after);
    if sum(gel_valid_pairs) > 1
        [~, gel_p_value] = ttest(gel_z_before(gel_valid_pairs), gel_z_after(gel_valid_pairs));
    else
        gel_p_value = NaN;
    end
    
    % Add significance indicators
    % Food comparison
    sigstar([1, 2], food_p_value);
    
    % Gel comparison
    sigstar([3+group_gap, 4+group_gap], gel_p_value);
    
    % Format plot
    set(gca, 'XTick', x_pos, 'XTickLabel', {'Before', 'After', 'Before', 'After'});
    ylabel('Pearson Correlation (r)');
    set(get(gca, 'YLabel'), 'FontSize', 18);

    % Add state to title
    title([upper(state(1)) state(2:end), ': Before vs After Food/Gel Discovery'], 'FontWeight', 'bold', 'FontSize', 16);
    
    % Set y-axis limits based on data range
    all_data = [food_r_before, food_r_after, gel_r_before, gel_r_after];
    y_min = min(all_data) - 0.1;
    y_max = max(all_data) + 0.2; % Extra space for significance markers
    ylim([y_min, y_max]);
    
    % Make plot look nicer
    grid off;
    box off;
    set(gca, 'LineWidth', 1);
    box off
    
    % Create a dummy plot for legend
    h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 10);
    
    % Add legend with only before/after colors
    legend_handle = legend([h1, h2], {'Before', 'After'}, 'Location', 'northeast','FontSize', 14);
    legend_handle.Box = 'off';

     % Add group labels
    text(1.5, y_min - 0.08*(y_max-y_min), 'Food', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);
    text(3.5+group_gap, y_min - 0.08*(y_max-y_min), 'Gel', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14);

end