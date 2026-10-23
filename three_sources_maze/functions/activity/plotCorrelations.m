%% Function to plot correlation results
function plotCorrelations(corr_data, options)
    % Extract options
    valid_sessions = options.sessions;
    arm_labels = {'Food Arm', 'Non-Food Arm 1', 'Non-Food Arm 2'};
    
    % Check if we have enough sessions for comparison
    if length(valid_sessions) < 1
        error('Need at least one session for plotting.');
    end
    
    % Calculate means and SEMs in z-space for each session
    mean_z = zeros(length(valid_sessions), 3);
    sem_z = zeros(length(valid_sessions), 3);
    mean_corr = zeros(length(valid_sessions), 3);
    sem_corr = zeros(length(valid_sessions), 3);
    
    for s = 1:length(valid_sessions)
        session = valid_sessions{s};
        if ~isempty(corr_data.z_corr.(session))
            % Calculate means and SEMs in z-space
            mean_z(s,:) = nanmean(corr_data.z_corr.(session), 1);
            sem_z(s,:) = nanstd(corr_data.z_corr.(session), [], 1) / sqrt(sum(~isnan(corr_data.z_corr.(session)(:,1))));
            
            % Convert back to correlation values for plotting
            mean_corr(s,:) = inverse_fisher_z(mean_z(s,:));
            sem_corr(s,:) = inverse_fisher_z(mean_z(s,:) + sem_z(s,:)) - mean_corr(s,:);
        end
    end
    
    % Define colors for the sessions
    colors = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0.6350 0.0780 0.1840]};
    
    % Data preparation for grouped bars
    means = mean_corr';
    sems = sem_corr';
    
    % Number of groups (arms) and number of bars (sessions) in each group
    numGroups = size(means, 1);
    numBars = size(means, 2);
    
    % Create figure
    figure('Position', [100, 100, 800, 500]);
    
    % Plot grouped bars
    b = bar(means, 'grouped');
    
    % Set colors for each session
    for s = 1:numBars
        b(s).FaceColor = colors{mod(s-1, length(colors))+1};
    end
    
    hold on;
    
    % Calculate the width of the group of bars
    groupWidth = min(0.8, numBars/(numBars + 1.5));
    
    % Set x-axis labels
    xticks(1:numGroups);
    xticklabels(arm_labels);
    
    % Add error bars
    for i = 1:numBars
        % Calculate x position for error bars
        x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
        
        % Plot error bars
        errorbar(x, means(:, i), sems(:, i), 'k.', 'LineWidth', 1.2);
    end
    
    % Add individual data points and connect them if requested
    if options.connect_points && length(valid_sessions) >= 2
        % Get the individual correlation values for each session
        session_data = cell(1, length(valid_sessions));
        for s = 1:length(valid_sessions)
            session_data{s} = corr_data.corr.(valid_sessions{s});
        end
        
        % Get mouse IDs for matching across sessions
        mouse_ids = cell(1, length(valid_sessions));
        for s = 1:length(valid_sessions)
            mouse_ids{s} = corr_data.mouse_ids.(valid_sessions{s});
        end
        
        % Find common mice across all sessions
        common_mice = mouse_ids{1};
        for s = 2:length(valid_sessions)
            common_mice = intersect(common_mice, mouse_ids{s});
        end
        
        % If we have common mice, plot connected points
        if ~isempty(common_mice)
            % For each arm
            for j = 1:numGroups
                % For each pair of consecutive sessions
                for s = 1:(length(valid_sessions)-1)
                    % Get indices for common mice in these sessions
                    [~, idx1, ~] = intersect(mouse_ids{s}, common_mice);
                    [~, idx2, ~] = intersect(mouse_ids{s+1}, common_mice);
                    
                    % Calculate x positions
                    x1 = j - groupWidth/2 + (2*s-1) * groupWidth / (2*numBars);
                    x2 = j - groupWidth/2 + (2*(s+1)-1) * groupWidth / (2*numBars);
                    
                    % Plot individual points
                    scatter(repmat(x1, length(idx1), 1), session_data{s}(idx1, j), 8, 'ko', 'filled');
                    scatter(repmat(x2, length(idx2), 1), session_data{s+1}(idx2, j), 8, 'ko', 'filled');
                    
                    % Connect points from consecutive sessions
                    for k = 1:length(common_mice)
                        plot([x1, x2], [session_data{s}(idx1(k), j), session_data{s+1}(idx2(k), j)], 'k-', 'LineWidth', 0.5);
                    end
                end
            end
        else
            % If no common mice, just plot individual points without connections
            for s = 1:length(valid_sessions)
                for j = 1:numGroups
                    x = j - groupWidth/2 + (2*s-1) * groupWidth / (2*numBars);
                    scatter(repmat(x, size(session_data{s}, 1), 1), session_data{s}(:, j), 8, 'ko', 'filled');
                end
            end
        end
    else
        % Just plot individual points without connections
        for s = 1:length(valid_sessions)
            session = valid_sessions{s};
            for j = 1:numGroups
                x = j - groupWidth/2 + (2*s-1) * groupWidth / (2*numBars);
                scatter(repmat(x, size(corr_data.corr.(session), 1), 1), corr_data.corr.(session)(:, j), 8, 'ko', 'filled');
            end
        end
    end
    
    % Labels and formatting
    ylabel('Pearson Correlation', 'FontSize', 12);
    if isfield(options, 'ylimit') && ~isempty(options.ylimit)
        ylim(options.ylimit);
    end
    
    title(['Correlation of dF/F with Distance (', upper(options.group), ')'], 'FontSize', 14);
    
    % Create session labels for legend
    session_labels = cell(1, length(valid_sessions));
    for s = 1:length(valid_sessions)
        session_labels{s} = strrep(valid_sessions{s}, 'sess', 'Session ');
    end
    
    legend(b, session_labels, 'Location', 'best');
    legend('boxoff');
    box off;
    
    hold off;
end