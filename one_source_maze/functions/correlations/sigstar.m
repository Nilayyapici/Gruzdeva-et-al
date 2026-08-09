function sigstar(group_positions, p_value)
    % Add significance stars to plot
    % Simplified version of the sigstar function
    
    % If p-value is NaN, don't add annotation
    if isnan(p_value)
        return;
    end
    
    % Height of the significance bar
    y_range = ylim;
    bar_height = y_range(2) * 0.98;
    
    % Draw the bar
    plot(group_positions, [bar_height, bar_height], 'k-', 'LineWidth', 1.5);
    
    % Add p-value text
    if p_value < 0.001
        p_text = 'p < 0.001***';
    elseif p_value < 0.01
        p_text = 'p < 0.01**';
    elseif p_value < 0.05
        p_text = 'p < 0.05*';
    else
        p_text = sprintf('p = %.3f', p_value);
    end
    
    text(mean(group_positions), bar_height*1.05, p_text, 'HorizontalAlignment', 'center', 'FontSize', 14);
end