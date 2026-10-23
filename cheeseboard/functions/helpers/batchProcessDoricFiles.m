function batchProcessDoricFiles(baseDir)
    % Batch process all .doric files in subdirectories and convert to CSV
    % Then convert to old photometry format
   
    % Find all .doric files recursively
    doricFiles = dir(fullfile(baseDir, '**', '*.doric'));
    
    fprintf('Found %d .doric files to process...\n', length(doricFiles));
    
    % Process each file
    for i = 1:length(doricFiles)
        try
            % Get file info
            filename = fullfile(doricFiles(i).folder, doricFiles(i).name);
            [filepath, name, ~] = fileparts(filename);
            
            fprintf('Processing %d/%d: %s\n', i, length(doricFiles), name);
            
            % Extract data from .doric file
            extractDoricData(filename, filepath, name);
            
            % Convert to old photometry format
            convertToPhotometryFormat(filepath, name);
            
        catch ME
            fprintf('Error processing %s: %s\n', doricFiles(i).name, ME.message);
        end
    end
    
    fprintf('Batch processing complete!\n');
end

function extractDoricData(filename, outputDir, baseName)
    % Extract data from single .doric file and save as CSV
    
    try
        % Read AnalogIn data
        AnalogIn_time = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AnalogIn/Time');
        AnalogIn_values = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AnalogIn/AIN01');
        AnalogIn_data = [AnalogIn_time, AnalogIn_values];
        
        % Read AIN01xAOUT01-LockIn data
        AIN01xAOUT01_time = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time');
        AIN01xAOUT01_values = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values');
        AIN01xAOUT01_data = [AIN01xAOUT01_time, AIN01xAOUT01_values];
        
        % Read AIN01xAOUT02-LockIn data
        AIN01xAOUT02_time = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Time');
        AIN01xAOUT02_values = h5read(filename, '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Values');
        AIN01xAOUT02_data = [AIN01xAOUT02_time, AIN01xAOUT02_values];
        
        % Save as CSV files
        csvwrite(fullfile(outputDir, 'AnalogIn_0000.csv'), AnalogIn_data);
        csvwrite(fullfile(outputDir, 'AIN01xAOUT01_LockIn_0000.csv'), AIN01xAOUT01_data);
        csvwrite(fullfile(outputDir, 'AIN01xAOUT02_LockIn_0000.csv'), AIN01xAOUT02_data);
        
        fprintf('  - Saved CSV files for %s\n', baseName);
        
    catch ME
        fprintf('  - Error extracting data from %s: %s\n', baseName, ME.message);
        rethrow(ME);
    end
end

function convertToPhotometryFormat(folderPath, baseName)
    % Convert CSV files to old photometry format using your existing function
    
    try
        % Define CSV file paths
        ain1File = fullfile(folderPath, 'AIN01xAOUT01_LockIn_0000.csv');
        ain2File = fullfile(folderPath, 'AIN01xAOUT02_LockIn_0000.csv');
        ainFile = fullfile(folderPath, 'AnalogIn_0000.csv');
        
        % Check if files exist
        if ~exist(ain1File, 'file') || ~exist(ain2File, 'file') || ~exist(ainFile, 'file')
            fprintf('  - Warning: Some CSV files missing for %s\n', baseName);
            return;
        end
        
        % Convert using your existing function
        [photometry, ain] = convertToOldPhotometryFormat(ain1File, ain2File, ainFile);
        
        % Save with unique names in the same folder
        save(fullfile(folderPath, sprintf('photometry_%s.mat', baseName)), 'photometry');
        save(fullfile(folderPath, sprintf('ain_%s.mat', baseName)), 'ain');
        
        fprintf('  - Converted to photometry format for %s\n', baseName);
        
    catch ME
        fprintf('  - Error converting %s: %s\n', baseName, ME.message);
    end
end

% Your existing conversion function (modified to work within the script)
function [photometry, ain] = convertToOldPhotometryFormat(ain1File, ain2File, ainFile)
    % convertToOldPhotometryFormat This function converts data from CSV files
    % to a specific old photometry format and saves it as .mat files, and
    % returns the matrices.
    %
    % Inputs:
    % ain1File - CSV file name for AIN1 data
    % ain2File - CSV file name for AIN2 data
    % ainFile - CSV file name for AIN data
    %
    % Outputs:
    % photometry - Combined photometry data
    % ain - AIN data
    
    % Read and convert AIN1 data from the specified file
    ain1_table = readtable(ain1File);
    ain1 = table2array(ain1_table(:, [1, 2]));
    
    % Read and convert AIN2 data from the specified file
    ain2_table = readtable(ain2File);
    ain2 = table2array(ain2_table(:, [1, 2]));
    
    % Read and convert AIN data from the specified file
    ain_table = readtable(ainFile);
    ain = table2array(ain_table(:, [1, 2]));
    
    % Combine AIN1 and AIN2 data into one array for photometry
    photometry = [ain1(:, 1), ain1(:, 2), ain2(:, 2)];
    
    % Note: Individual .mat files are saved in the calling function
    % with unique names to avoid overwriting
end