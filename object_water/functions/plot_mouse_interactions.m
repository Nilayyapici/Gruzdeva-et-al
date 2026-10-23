function plot_mouse_interactions(data, mouse_id, interaction_type)
    % Create figure with two subplots
    figure;
    
    % Choose which data to use based on interaction_type
    if strcmp(interaction_type, 'water')
        interaction_data = data(:,15);
        distance_data = data(:,9);  % distance to water
        title_text = ['Mouse Position: ' mouse_id ' - Water'];
        colorbar_label_1 = 'Water Interaction';
        colorbar_label_2 = 'Distance to Water (cm)';
    elseif strcmp(interaction_type, 'object')
        interaction_data = data(:,16);
        distance_data = data(:,10);  % distance to object
        title_text = ['Mouse Position: ' mouse_id ' - Object'];
        colorbar_label_1 = 'Object Interaction';
        colorbar_label_2 = 'Distance to Object (cm)';
    end
    
    % Subplot 1: Interaction
    subplot(1, 2, 1);
    scatter(data(:,2), data(:,3), 20, interaction_data, 'filled');
    colormap(gca, [0 0 1; 1 0 0]); % Blue for 0 (no interaction), Red for 1 (interaction)
    c1 = colorbar;
    c1.Ticks = [0.25, 0.75];
    c1.TickLabels = {'No Interaction', colorbar_label_1};
    xlabel('X Position');
    ylabel('Y Position');
    title(['Interaction: ' title_text]);
    axis equal;
    grid on;
    
    % Subplot 2: Distance
    subplot(1, 2, 2);
    scatter(data(:,2), data(:,3), 20, distance_data, 'filled');
    colormap(gca, 'jet'); % Use jet colormap for distances
    c2 = colorbar;
    c2.Label.String = colorbar_label_2;
    xlabel('X Position');
    ylabel('Y Position');
    title(['Distance: ' title_text]);
    axis equal;
    grid on;
    
    % Adjust figure size
    set(gcf, 'Position', [100, 100, 1200, 500]);
end