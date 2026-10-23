%% Function to analyze correlations between dF/F and distance
function [corr_data] = analyzeDffDistanceCorrelations(mice_all, options)
    % Default options
    if nargin < 2
        options = struct();
    end
    
    % Basic options
    if ~isfield(options, 'sessions'), options.sessions = {'sess0', 'sess1', 'sess2'}; end
    if ~isfield(options, 'group'), options.group = 'all'; end
    if ~isfield(options, 'title'), options.title = 'dF/F and Distance Correlation Analysis'; end
    if ~isfield(options, 'use_limits'), options.use_limits = true; end
    
    % Time and distance limits
    if ~isfield(options, 'time_lim'), options.time_lim = 7; end  % in minutes
    if ~isfield(options, 'dist_lim'), options.dist_lim = 150; end  % for positions as far as this limit
    if ~isfield(options, 'dist_too_close'), options.dist_too_close = 0; end  % for removing too close to food positions
    % Filter parameters
    if ~isfield(options, 'filter_order_dist'), options.filter_order_dist = 1; end
    if ~isfield(options, 'cutoff_freq_dist'), options.cutoff_freq_dist = 0.08; end
    if ~isfield(options, 'filter_order_dff'), options.filter_order_dff = 1; end
    if ~isfield(options, 'cutoff_freq_dff'), options.cutoff_freq_dff = 0.05; end
    
    % Plotting options
    if ~isfield(options, 'connect_points'), options.connect_points = true; end
    if ~isfield(options, 'ylimit'), options.ylimit = [-0.4, 0.4]; end
    if ~isfield(options, 'run_stats'), options.run_stats = true; end
    if ~isfield(options, 'remove_grooming'), options.remove_grooming = true; end
    if ~isfield(options, 'speed_threshold'), options.speed_threshold = 0; end  % 0 means no threshold
    
    % Initialize correlation data structure
    corr_data = struct();
    
    % Initialize session data
    for s = 1:length(options.sessions)
        session = options.sessions{s};
        corr_data.corr.(session) = [];
        corr_data.pval.(session) = [];
        corr_data.z_corr.(session) = [];
        corr_data.groups.(session) = [];
        corr_data.mouse_ids.(session) = {};
    end
    
    % Process each mouse
    for i = 1:size(mice_all, 1)
        session_info = mice_all{i, 1};  % Mouse session identifier
        group = mice_all{i, 2};        % 'control' or 'CNO'
        food_arm = mice_all{i, 3};     % Food arm location
        data = mice_all{i, 4};         % Actual session data
        
        % Check if we should include this entry based on group filter
        if ~strcmp(options.group, 'all') && ~strcmp(group, options.group)
            continue;  % Skip if it doesn't match the filter
        end
        
        % Determine which session we're processing
        if contains(session_info, 'sess0')
            session_type = 'sess0';
        elseif contains(session_info, 'sess1')
            session_type = 'sess1';
        elseif contains(session_info, 'sess2')
            session_type = 'sess2';
        else
            continue;  % Skip if session type is unknown
        end
        
        % Check if this session is in our analysis list
        if ~ismember(session_type, options.sessions)
            continue;
        end
        
        % Extract mouse ID
        mouse_id = extractMouseID(session_info);
        
        % Process session data based on options
        if strcmp(session_type, 'sess2')
            % For session 2, keep only closed door condition
            data = data(data(:,12) == 0, :);
        elseif strcmp(session_type, 'sess1') && ~isempty(data)
            % For session 1, make sure we're only using the door=1 (open) condition
            data = data(data(:,12) == 1, :);
            
            % Fix time for sess1 by making it relative to the start
            if ~isempty(data)
                data(:,1) = data(:,1) - data(1,1);
            end
        elseif strcmp(session_type, 'sess0') && ~isempty(data)
            % For session 0, make sure we're only using the door=0 (closed) condition
            data = data(data(:,12) == 0, :);
        end
        
        % Skip if no data after filtering
        if isempty(data)
            continue;
        end
        
        % Remove grooming periods if requested
        if options.remove_grooming && size(data, 2) >= 13
            data = data(data(:,13) == 0, :);  % Keep only rows where grooming is 0
        end
        
        % Apply speed threshold if requested
        if options.speed_threshold > 0 && size(data, 2) >= 7
            data = data(data(:,7) >= options.speed_threshold, :);
        end
        
        % Skip if no data after additional filtering
        if isempty(data)
            continue;
        end
        
        % Apply time limit if needed
        if options.use_limits
            % Find the index corresponding to the time limit
            limit_idx = findnearest(options.time_lim * 60, data(:,1));
            data_limited = data(1:limit_idx, :);
        else
            data_limited = data;
        end
        
        % Skip if no data after time limiting
        if isempty(data_limited)
            continue;
        end
        
        % Design Butterworth filters for distance and dF/F
        [b_dist, a_dist] = butter(options.filter_order_dist, options.cutoff_freq_dist, 'low');
        [b_dff, a_dff] = butter(options.filter_order_dff, options.cutoff_freq_dff, 'low');
        
        % Create filtered data
        filt_data = data_limited;
        
        % Filter distance measurements
        filt_data(:,9) = filtfilt(b_dist, a_dist, data_limited(:,9));   % Food arm
        filt_data(:,10) = filtfilt(b_dist, a_dist, data_limited(:,10)); % Non-food arm 2
        filt_data(:,11) = filtfilt(b_dist, a_dist, data_limited(:,11)); % Non-food arm 1
        
        % Filter dF/F
        filt_data(:,6) = filtfilt(b_dff, a_dff, data_limited(:,6));
        
        % Extract filtered data
        dist_food = filt_data(:,9);     % Distance to food arm
        dist_nonfood1 = filt_data(:,11); % Distance to non-food arm 1
        dist_nonfood2 = filt_data(:,10); % Distance to non-food arm 2
        dff = filt_data(:,6);          % Filtered dF/F

        % Apply distance limit filter for correlations
        % Apply distance limit filter for correlations - DIFFERENT FOR SESS1
        if strcmp(session_type, 'sess1')
            % For session 1, apply the dist_too_close limit
            valid_indices_food = (dist_food > options.dist_too_close) & (dist_food <= options.dist_lim);
        else
            % For other sessions, use only the upper limit
            valid_indices_food = dist_food <= options.dist_lim;
        end
        valid_indices_nonfood1 = dist_nonfood1 <= options.dist_lim;
        valid_indices_nonfood2 = dist_nonfood2 <= options.dist_lim;


        % Calculate correlations only for distances <= dist_lim
        % If no valid data points, set correlation and p-value to NaN
        if sum(valid_indices_food) > 2
            [rho_food, pval_food] = corr(dist_food(valid_indices_food), dff(valid_indices_food), 'Type', 'Pearson');
        else
            rho_food = NaN;
            pval_food = NaN;
        end
        
        if sum(valid_indices_nonfood1) > 2
            [rho_nonfood1, pval_nonfood1] = corr(dist_nonfood1(valid_indices_nonfood1), dff(valid_indices_nonfood1), 'Type', 'Pearson');
        else
            rho_nonfood1 = NaN;
            pval_nonfood1 = NaN;
        end
        
        if sum(valid_indices_nonfood2) > 2
            [rho_nonfood2, pval_nonfood2] = corr(dist_nonfood2(valid_indices_nonfood2), dff(valid_indices_nonfood2), 'Type', 'Pearson');
        else
            rho_nonfood2 = NaN;
            pval_nonfood2 = NaN;
        end
        
        % Store correlations and p-values
        corr_values = [rho_food, rho_nonfood1, rho_nonfood2];
        pval_values = [pval_food, pval_nonfood1, pval_nonfood2];
        
        % Convert to Fisher's z for averaging
        z_corr_values = fisher_z(corr_values);
        
        % Store results
        corr_data.corr.(session_type) = [corr_data.corr.(session_type); corr_values];
        corr_data.pval.(session_type) = [corr_data.pval.(session_type); pval_values];
        corr_data.z_corr.(session_type) = [corr_data.z_corr.(session_type); z_corr_values];
        corr_data.groups.(session_type) = [corr_data.groups.(session_type); {group}];
        corr_data.mouse_ids.(session_type){end+1} = mouse_id;
    end

    % Plot the results
    plotCorrelations(corr_data, options);

    % Run statistical tests if requested
    if options.run_stats
        performCorrStats(corr_data, options);
    end
end