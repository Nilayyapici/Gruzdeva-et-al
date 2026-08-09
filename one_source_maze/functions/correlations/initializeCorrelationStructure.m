% Initialize the correlation data structure with all required fields
function corr_data = initializeCorrelationStructure(options)
    corr_data = struct();
    
    % Groups, sources, and directions
    groups = {'fasted', 'fed'};
    sources = {'food', 'gel'};
    directions = {'towards', 'away'};
    
    % Create nested structure for all combinations
    for g = 1:length(groups)
        for s = 1:length(sources)
            group_name = groups{g};
            source_name = sources{s};
            field_name = [group_name, '_', source_name];
            
            % Initialize the field for this group-source combination
            corr_data.(field_name) = struct();
            
            for d = 1:length(directions)
                direction = directions{d};
                corr_data.(field_name).(direction) = struct();
                
                % For each session
                for sess = 1:length(options.sessions)
                    session = options.sessions{sess};
                    corr_data.(field_name).(direction).(session) = struct(...
                        'corr', [], ...       % Correlation coefficients
                        'pval', [], ...       % p-values
                        'z_corr', [], ...     % Fisher z-transformed correlations
                        'mouse_ids', {{}});   % Mouse IDs for tracking individual mice
                end
            end
        end
    end
end