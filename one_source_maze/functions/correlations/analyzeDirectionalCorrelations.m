function [corr_data] = analyzeDirectionalCorrelations(mice_all, options)
    % ANALYZEDIRECTIONALCORRELATIONS Analyze correlations between dF/F and distance
    % during towards/away movements for food sources across sessions
    %
    % Inputs:
    %   mice_all - Cell array containing preprocessed mouse data
    %     Column 1: Mouse ID (with session suffix)
    %     Column 2: Condition ('fasted' or 'fed')
    %     Column 3: Source type ('food' or 'gel')
    %     Column 4: Data matrix with columns:
    %       1-time, 2-x, 3-y, 4-speed, 5-dist to food, 6-path from last visit,
    %       7-door (0-closed/1-open), 8-food interaction, 9-eating (0/1),
    %       10-grooming (0/1), 11-dff
    %     Column 5: Signal quality
    %     Column 6: Discovery frame
    %
    %   options - Structure with analysis parameters
    %
    % Outputs:
    %   corr_data - Structure with correlation results organized by condition,
    %               source, direction, and session
    
    % Set default options if not provided
    if nargin < 2
        options = struct();
    end
    
    % Default options
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2', 'sess3'}; end
    if ~isfield(options, 'group'), options.group = 'all'; end % 'all', 'fasted', or 'fed'
    if ~isfield(options, 'source'), options.source = 'all'; end % 'all', 'food', or 'gel'
    if ~isfield(options, 'dist_limit'), options.dist_limit = 200; end % For correlation calculations
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'remove_eating'), options.remove_eating = true; end
    if ~isfield(options, 'filter_order'), options.filter_order = 1; end
    if ~isfield(options, 'cutoff_freq'), options.cutoff_freq = 0.08; end
    if ~isfield(options, 'movement_threshold'), options.movement_threshold = 0.4; end % Threshold for towards/away detection
    if ~isfield(options, 'ylimit'), options.ylimit = [-0.5, 0.5]; end % For plots
    if ~isfield(options, 'run_stats'), options.run_stats = true; end
    if ~isfield(options, 'debug'), options.debug = false; end
    if ~isfield(options, 'debug_plot_selection'), options.debug_plot_selection = 'all'; end % 'all', 'first_n', 'random_n', specific mouse IDs
    if ~isfield(options, 'debug_plot_count'), options.debug_plot_count = 3; end % Number of mice to plot when using 'first_n' or 'random_n'
    
    % Relevant column indices
    COL_TIME = 1;
    COL_DIST = 5;
    COL_DOOR = 7;
    COL_EATING = 9;
    COL_GROOMING = 10;
    COL_DFF = 11;
    
    % Initialize correlation data structure with all needed fields
    corr_data = initializeCorrelationStructure(options);
    
    % Process each mouse/session
    if options.debug
        fprintf('Processing %d entries...\n', size(mice_all, 1));
        
        % Prepare debugging plot selection
        mice_to_plot = [];
        if ischar(options.debug_plot_selection)
            if strcmp(options.debug_plot_selection, 'all')
                mice_to_plot = 1:size(mice_all, 1);
            elseif strcmp(options.debug_plot_selection, 'first_n')
                mice_to_plot = 1:min(options.debug_plot_count, size(mice_all, 1));
            elseif strcmp(options.debug_plot_selection, 'random_n')
                % Select random indices without replacement
                n_mice = min(options.debug_plot_count, size(mice_all, 1));
                mice_to_plot = randperm(size(mice_all, 1), n_mice);
            end
        elseif iscell(options.debug_plot_selection)
            % User provided specific mouse IDs
            for i = 1:size(mice_all, 1)
                if any(strcmp(extractMouseID(mice_all{i, 1}), options.debug_plot_selection))
                    mice_to_plot = [mice_to_plot, i];
                end
            end
        end
    end
    
    for i = 1:size(mice_all, 1)
        % Extract information
        session_info = mice_all{i, 1};  % Mouse session identifier
        group = mice_all{i, 2};         % 'fasted' or 'fed'
        source = mice_all{i, 3};        % 'food' or 'gel'
        raw_data = mice_all{i, 4};      % Actual data matrix
        
        % Check if we should include this entry based on group filter
        if ~strcmp(options.group, 'all') && ~strcmp(group, options.group)
            continue;
        end
        
        % Check if we should include this entry based on source filter
        if ~strcmp(options.source, 'all') && ~strcmp(source, options.source)
            continue;
        end
        
        % Determine field name for this combination
        field_name = [group, '_', source];
        
        % Determine which session we're processing
        if contains(session_info, 'sess0')
            session_type = 'sess0';
        elseif contains(session_info, 'sess1')
            session_type = 'sess1';
        elseif contains(session_info, 'sess2')
            session_type = 'sess2';
        elseif contains(session_info, 'sess3')
            session_type = 'sess3';
        else
            continue;  % Skip if session type is unknown
        end
        
        % Check if this session is in our analysis list
        if ~ismember(session_type, options.sessions)
            continue;
        end
        
        % Extract mouse ID from the session info
        mouse_id = session_info;
        if contains(session_info, '_sess')
            parts = strsplit(session_info, '_sess');
            mouse_id = parts{1};
        end
        
        % Filter and analyze the data
        try
            % Create a copy of the data
            data = raw_data;
            
            % Remove grooming and eating if requested
            if options.remove_grooming && size(data, 2) >= COL_GROOMING
                data = data(data(:, COL_GROOMING) == 0, :);
            end
            
            if options.remove_eating && size(data, 2) >= COL_EATING
                data = data(data(:, COL_EATING) == 0, :);
            end
            
            % Skip if not enough data after filtering
            if size(data, 1) < 10
                if options.debug
                    fprintf('Warning: Not enough data points for mouse %s after filtering\n', mouse_id);
                end
                continue;
            end
            
            % Filter the distance signal to smooth it
            try
                % Minimum data length required for filtfilt (3 times filter order)
                min_data_length = 3 * options.filter_order;
                
                if size(data, 1) > min_data_length
                    % Create Butterworth filter
                    [b, a] = butter(options.filter_order, options.cutoff_freq, 'low');
                    
                    % Apply filter to distance data
                    data(:, COL_DIST) = filtfilt(b, a, data(:, COL_DIST));
                end
            catch filter_err
                if options.debug
                    fprintf('Warning: Filtering error for mouse %s: %s\n', mouse_id, filter_err.message);
                end
                % Continue with unfiltered data
            end
            
            % Calculate distance derivative to determine approach/retreat
            dist_derivative = [0; diff(data(:, COL_DIST))]; % Append zero at start
            
            % Identify towards and away movements based on distance derivative
            towards_idx = dist_derivative < -options.movement_threshold & ...
                          data(:, COL_DIST) <= options.dist_limit;
            away_idx = dist_derivative > options.movement_threshold & ...
                       data(:, COL_DIST) <= options.dist_limit;
            
            % Create debug plot if requested
            if options.debug && ismember(i, mice_to_plot)
                plotTowardsAwayDetection(data, dist_derivative, towards_idx, away_idx, ...
                                        mouse_id, group, source, session_type, options);
            end
            
            % Optional debugging output
            if options.debug
                fprintf('Mouse %s (%s, %s, %s): Total=%d, Towards=%d, Away=%d points\n', ...
                        mouse_id, group, source, session_type, ...
                        size(data, 1), sum(towards_idx), sum(away_idx));
            end
            
            % Calculate correlation for towards movement
            if sum(towards_idx) >= 4  % Need at least 4 points for meaningful correlation
                [rho_towards, pval_towards] = corr(data(towards_idx, COL_DFF), ...
                                                   data(towards_idx, COL_DIST), ...
                                                   'Type', 'Pearson');
                
                % Convert to Fisher's z-transform for better statistics
                z_towards = fisher_z(rho_towards);
                
                % Store results
                corr_data.(field_name).towards.(session_type).corr = ...
                    [corr_data.(field_name).towards.(session_type).corr; rho_towards];
                corr_data.(field_name).towards.(session_type).pval = ...
                    [corr_data.(field_name).towards.(session_type).pval; pval_towards];
                corr_data.(field_name).towards.(session_type).z_corr = ...
                    [corr_data.(field_name).towards.(session_type).z_corr; z_towards];
                corr_data.(field_name).towards.(session_type).mouse_ids = ...
                    [corr_data.(field_name).towards.(session_type).mouse_ids; {mouse_id}];
            end
            
            % Calculate correlation for away movement
            if sum(away_idx) >= 4  % Need at least 4 points for meaningful correlation
                [rho_away, pval_away] = corr(data(away_idx, COL_DFF), ...
                                             data(away_idx, COL_DIST), ...
                                             'Type', 'Pearson');
                
                % Convert to Fisher's z-transform for better statistics
                z_away = fisher_z(rho_away);
                
                % Store results
                corr_data.(field_name).away.(session_type).corr = ...
                    [corr_data.(field_name).away.(session_type).corr; rho_away];
                corr_data.(field_name).away.(session_type).pval = ...
                    [corr_data.(field_name).away.(session_type).pval; pval_away];
                corr_data.(field_name).away.(session_type).z_corr = ...
                    [corr_data.(field_name).away.(session_type).z_corr; z_away];
                corr_data.(field_name).away.(session_type).mouse_ids = ...
                    [corr_data.(field_name).away.(session_type).mouse_ids; {mouse_id}];
            end
            
        catch error_msg
            if options.debug
                fprintf('Error processing mouse %s: %s\n', mouse_id, error_msg.message);
            end
        end
    end
    
    % Visualization
    plotDirectionalCorrelations(corr_data, options);
    
    % Statistical analysis
    if options.run_stats
        performCorrelationStats(corr_data, options);
    end
end