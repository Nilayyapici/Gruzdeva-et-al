%% Convert photometry files for all the learning trials
% Set the base directory where all your folders are located
% baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard\MDRI17\MDRI17_saline';
% batchProcessDoricFiles(baseDir);

% %% Transfer DLC files to the folders (from hard drive to the cheeseboard folders)
% organizeFilesByMouseAndTrial();

% %% Batch trials 
% 
% mouseName = 'MDRI17';
% % baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard\MDRE8\MDRE8_saline2';
% x_food = 414;
% y_food = 216;
% 
% processBatchTrials(mouseName, baseDir,x_food,y_food);

% %% Plot trajectories to check everything
% plotMouseTracksWithDFF(mouseName, baseDir);   
    
%% Learning analysis all mice together
% Learning curves for all the mice
mouseNames = {'MDRE8', 'F13', 'M21', 'M22', 'FDRE14','FDRE15','FDRE4', 'MDRE27', 'MDRE29', 'MDRI24'};
baseDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';

options.num_trials = 10;
options.use_zscore = 1;           % Use z-scored dF/F (default: false)
options.distance_limit = 30;        % Maximum distance to analyze in cm (default: 20)
options.bin_size = 0.5;             % Bin size for distance binning in cm (default: 0.5)
options.dff_ylim = [-1, 1.2];
options.food_discovery_threshold = 1.5;
options.smoothing_method = 'sgolay';
options.smoothing_window = 10;
options.save_xlsx = 1;
options.speed_threshold = 50;           % cm/s
options.exclude_first_n_frames = 20;     % frames
options.min_speed_threshold = 5; 

options.early_trials = [1:10];      % Trials in early phase
options.middle_trials = [4:6];  % Trials in middle phase  
options.late_trials = [7:10];      % Trials in late phase
options.plot_phases = {'early'};  % Which phases to plot {'early','middle', 'late'};

% Run main analysis and get the data
allMiceData = analyzeAllMice_learning_cheeseboard(baseDir, mouseNames, options);

%% Individual mice

options.num_trials = 10;
options.food_discovery_threshold = 1;
validate_path_lengths(baseDir, mouseNames, options);

%% Validation&Trajectories

options.trial_number = 2;
options.speed_threshold = 50;           % cm/s — default unchanged from before
options.exclude_first_n_frames = 20;
validateTrialData(allMiceData, options);

%% Individual trajectories

opts.min_speed_threshold = 0;
opts.food_discovery_threshold = 1.2;
plotMouseTrajectories_allTrials(allMiceData, opts)