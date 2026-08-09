function grouped = groupResultsByCondition(results)
    % Group results by condition and stimulus
    grouped = struct();
    
    % Create masks for each group and extract data
    for condition = {'fasted', 'fed'}
        for stimulus = {'food', 'gel'}
            % Create field name
            field_name = [condition{1}, '_', stimulus{1}];
            
            % Find matching mice
            mask = strcmp({results.condition}, condition{1}) & ...
                   strcmp({results.stimulus}, stimulus{1});
            
            % Store matched results
            grouped.(field_name) = results(mask);
        end
    end
end