function organizeFilesByMouseAndTrial()
    % Automatically organize files from hard drive to appropriate folders
    % based on mouse name and trial information in filename
    
    % Set source directory (where files are currently located)
    sourceDir = 'D:\temp';  % UPDATE THIS PATH
    
    % Set base destination directory
    baseDestDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';
    
    % Get all files in source directory (adjust pattern as needed)
    filePattern = '*.csv';  % Change to appropriate file extension
    files = dir(fullfile(sourceDir, filePattern));
    
    fprintf('Found %d files to organize...\n', length(files));
    
    % Process each file
    for i = 1:length(files)
        try
            filename = files(i).name;
            fprintf('Processing %d/%d: %s\n', i, length(files), filename);
            
            % Extract mouse name and trial from filename
            [mouseName, trialName] = parseFilename(filename);
            
            if isempty(mouseName) || isempty(trialName)
                fprintf('  - Warning: Could not parse mouse/trial from %s\n', filename);
                continue;
            end
            
            % Find destination folder
            destFolder = findDestinationFolder(baseDestDir, mouseName, trialName);
            
            if isempty(destFolder)
                fprintf('  - Warning: Could not find destination folder for %s/%s\n', mouseName, trialName);
                continue;
            end
            
            % Move/copy file
            sourceFile = fullfile(sourceDir, filename);
            destFile = fullfile(destFolder, filename);
            
            % Create destination directory if it doesn't exist
            if ~exist(destFolder, 'dir')
                mkdir(destFolder);
                fprintf('  - Created directory: %s\n', destFolder);
            end
            
            % Copy file (change to movefile if you want to move instead of copy)
            copyfile(sourceFile, destFile);
            fprintf('  - Moved to: %s\n', destFolder);
            
        catch ME
            fprintf('  - Error processing %s: %s\n', filename, ME.message);
        end
    end
    
    fprintf('File organization complete!\n');
end

function [mouseName, trialName] = parseFilename(filename)
    % Parse mouse name and trial from filename
    % Examples:
    % F13_saline_A12024-09-11T11_39_00DLC_resnet50_cheeseboardJun4shuffle1_200000
    % FDRE14_C22024-08-21T14_02_14DLC_resnet50_cheeseboardJun4shuffle1_200000
    
    mouseName = '';
    trialName = '';
    
    % Remove file extension
    [~, name, ~] = fileparts(filename);
    
    % Method 1: Look for patterns like F13_saline_A1 or FDRE14_C2
    % Split by underscore
    parts = split(name, '_');
    
    if length(parts) >= 2
        % First part should be mouse name (F13, FDRE14, etc.)
        mouseName = parts{1};
        
        % Look for trial pattern in the concatenated string
        % Pattern: Mouse_something_Trial followed by date
        fullString = name;
        
        % Use regex to find trial pattern
        % Look for pattern: letter(s) followed by number(s) before the date
        trialPattern = '([A-Z]+\d+)(\d{4}-\d{2}-\d{2})';
        tokens = regexp(fullString, trialPattern, 'tokens');
        
        if ~isempty(tokens)
            trialName = tokens{1}{1};
        else
            % Alternative: look for last part before date
            datePattern = '\d{4}-\d{2}-\d{2}';
            beforeDate = regexprep(fullString, datePattern, '');
            
            % Extract trial from end of string before date
            trialPattern2 = '([A-Z]\d+)$';
            trialTokens = regexp(beforeDate, trialPattern2, 'tokens');
            if ~isempty(trialTokens)
                trialName = trialTokens{1}{1};
            end
        end
    end
    
    fprintf('  - Parsed: Mouse=%s, Trial=%s\n', mouseName, trialName);
end

function destFolder = findDestinationFolder(baseDir, mouseName, trialName)
    % Find the appropriate destination folder for given mouse and trial
    % Structure: MouseName/MouseName_saline*/TrialName
    
    destFolder = '';
    
    % Look for mouse folder
    mouseFolder = fullfile(baseDir, mouseName);
    if ~exist(mouseFolder, 'dir')
        fprintf('    - Mouse folder not found: %s\n', mouseFolder);
        return;
    end
    
    % Look for saline folder (with potential variations)
    salineFolders = dir(fullfile(mouseFolder, sprintf('%s*saline*', mouseName)));
    
    if isempty(salineFolders)
        % Try without mouse name prefix
        salineFolders = dir(fullfile(mouseFolder, '*saline*'));
    end
    
    if isempty(salineFolders)
        fprintf('    - No saline folder found in: %s\n', mouseFolder);
        return;
    end
    
    % Use first saline folder found
    salineFolder = fullfile(mouseFolder, salineFolders(1).name);
    
    % Look for trial folder
    trialFolder = fullfile(salineFolder, trialName);
    
    if exist(trialFolder, 'dir')
        destFolder = trialFolder;
    else
        fprintf('    - Trial folder not found: %s\n', trialFolder);
        fprintf('    - Available folders in %s:\n', salineFolder);
        trialFolders = dir(salineFolder);
        for j = 1:length(trialFolders)
            if trialFolders(j).isdir && ~startsWith(trialFolders(j).name, '.')
                fprintf('      %s\n', trialFolders(j).name);
            end
        end
    end
end

% Alternative version with manual mapping for complex cases
function organizeFilesWithMapping()
    % Version with manual mapping for complex naming patterns
    
    sourceDir = 'D:\YourFiles';  % UPDATE THIS PATH
    baseDestDir = 'C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard';
    
    % Define manual mapping for special cases
    mappings = containers.Map();
    
    % Add manual mappings as needed
    % mappings('F13_A1') = 'F13\F13_saline_walls\A1';
    % mappings('FDRE14_C2') = 'FDRE14\FDRE14_saline\C2';
    
    files = dir(fullfile(sourceDir, '*.csv'));  % Adjust file pattern
    
    for i = 1:length(files)
        filename = files(i).name;
        [mouseName, trialName] = parseFilename(filename);
        
        key = sprintf('%s_%s', mouseName, trialName);
        
        if isKey(mappings, key)
            % Use manual mapping
            relativePath = mappings(key);
            destFolder = fullfile(baseDestDir, relativePath);
        else
            % Use automatic detection
            destFolder = findDestinationFolder(baseDestDir, mouseName, trialName);
        end
        
        if ~isempty(destFolder)
            sourceFile = fullfile(sourceDir, filename);
            destFile = fullfile(destFolder, filename);
            
            if ~exist(destFolder, 'dir')
                mkdir(destFolder);
            end
            
            copyfile(sourceFile, destFile);
            fprintf('Moved %s to %s\n', filename, destFolder);
        else
            fprintf('Could not find destination for %s\n', filename);
        end
    end
end