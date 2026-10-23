function plotBeforeAfterComparisons(grouped, options)
    % Plot before vs after comparisons for each group
    
    % Create one figure per state-source combination
    for i = 1:length(options.state)
        for j = 1:length(options.source)
            state = options.state{i};
            source = options.source{j};
            field_name = [state, '_', source];
            
            if isfield(grouped, field_name) && ~isempty(grouped.(field_name))
                % figure('Position', [100 + (i-1)*400, 100 + (j-1)*400, 450, 600]);
                plotSingleGroupComparison(grouped.(field_name), state, source);
            end
        end
    end
end