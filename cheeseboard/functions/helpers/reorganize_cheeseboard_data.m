function mice = reorganize_cheeseboard_data(Cheeseboard)
% REORGANIZE_CHEESEBOARD_DATA - Convert Cheeseboard array to mice array format
%
% Usage:
%   mice = reorganize_cheeseboard_data(Cheeseboard);
%
% Input:
%   Cheeseboard - cell array with structure:
%     Cheeseboard{i,1} - mouse name
%     Cheeseboard{i,2} - [x,y] food position saline
%     Cheeseboard{i,3} - pre-test data saline (Nx10 double)
%     Cheeseboard{i,4} - test data saline (Nx10 double)
%     Cheeseboard{i,5} - [x,y] food position CNO
%     Cheeseboard{i,6} - pre-test data CNO (Nx10 double)
%     Cheeseboard{i,7} - test data CNO (Nx10 double)
%     Cheeseboard{i,8} - [x0,y0] starting position saline
%     Cheeseboard{i,9} - [x0,y0] starting position CNO
%     Cheeseboard{i,10} - latencies saline [pre;test]
%     Cheeseboard{i,11} - latencies CNO [pre;test]
%
% Output:
%   mice - cell array with structure:
%     mice{i,1} - mouse name (with _saline or _CNO suffix)
%     mice{i,2} - [x,y] food position
%     mice{i,3} - pre-test data (Nx10 double)
%     mice{i,4} - test data (Nx10 double)
%     mice{i,5} - [x0,y0] starting position
%     mice{i,6} - latencies [pre;test]
%     mice{i,7} - Time;Percent Pre (empty, to be filled later)
%     mice{i,8} - Time;Percent Test (empty, to be filled later)

%% Input validation
if nargin < 1
    error('Cheeseboard data array is required');
end

if size(Cheeseboard, 2) < 11
    error('Cheeseboard array must have at least 11 columns');
end

%% Initialize
n_cheeseboard_mice = size(Cheeseboard, 1);
mice = {};
mouse_count = 0;

fprintf('Converting Cheeseboard data to mice array format...\n');
fprintf('Processing %d mice from Cheeseboard array\n\n', n_cheeseboard_mice);

%% Process each mouse in Cheeseboard array
for i = 1:n_cheeseboard_mice
    mouse_name = Cheeseboard{i,1};
    fprintf('Processing mouse %d/%d: %s\n', i, n_cheeseboard_mice, mouse_name);
    
    %% Process SALINE session (always present)
    has_saline_data = ~isempty(Cheeseboard{i,2}) && ~isempty(Cheeseboard{i,3}) && ~isempty(Cheeseboard{i,4});
    
    if has_saline_data
        mouse_count = mouse_count + 1;
        
        % Create saline session entry
        mice{mouse_count, 1} = [mouse_name '_saline'];
        mice{mouse_count, 2} = Cheeseboard{i,2};  % food position saline
        mice{mouse_count, 3} = Cheeseboard{i,3};  % pre data saline
        mice{mouse_count, 4} = Cheeseboard{i,4};  % test data saline
        mice{mouse_count, 5} = Cheeseboard{i,8};  % starting position saline
        mice{mouse_count, 6} = Cheeseboard{i,10}; % latencies saline
        mice{mouse_count, 7} = [];  % Time;Percent Pre (empty)
        mice{mouse_count, 8} = [];  % Time;Percent Test (empty)
        
        % Validate saline data
        if isempty(mice{mouse_count, 2})
            warning('Mouse %s: Missing food position for saline session', mouse_name);
        end
        if isempty(mice{mouse_count, 5})
            warning('Mouse %s: Missing starting position for saline session', mouse_name);
        end
        if isempty(mice{mouse_count, 6})
            warning('Mouse %s: Missing latencies for saline session', mouse_name);
        end
        
        fprintf('  → Added saline session: %s\n', mice{mouse_count, 1});
        fprintf('    Food position: [%.1f, %.1f]\n', mice{mouse_count, 2}(1), mice{mouse_count, 2}(2));
        fprintf('    Pre data: %dx%d, Test data: %dx%d\n', ...
                size(mice{mouse_count, 3}, 1), size(mice{mouse_count, 3}, 2), ...
                size(mice{mouse_count, 4}, 1), size(mice{mouse_count, 4}, 2));
        if ~isempty(mice{mouse_count, 5})
            fprintf('    Start position: [%.1f, %.1f]\n', mice{mouse_count, 5}(1), mice{mouse_count, 5}(2));
        end
        if ~isempty(mice{mouse_count, 6})
            fprintf('    Latencies: [%.1f, %.1f] seconds\n', mice{mouse_count, 6}(1), mice{mouse_count, 6}(2));
        end
    else
        fprintf('  → No saline data found, skipping\n');
    end
    
    %% Process CNO session (if present)
    has_cno_data = ~isempty(Cheeseboard{i,5}) && ~isempty(Cheeseboard{i,6}) && ~isempty(Cheeseboard{i,7});
    
    if has_cno_data
        mouse_count = mouse_count + 1;
        
        % Create CNO session entry
        mice{mouse_count, 1} = [mouse_name '_CNO'];
        mice{mouse_count, 2} = Cheeseboard{i,5};  % food position CNO
        mice{mouse_count, 3} = Cheeseboard{i,6};  % pre data CNO
        mice{mouse_count, 4} = Cheeseboard{i,7};  % test data CNO
        mice{mouse_count, 5} = Cheeseboard{i,9};  % starting position CNO
        mice{mouse_count, 6} = Cheeseboard{i,11}; % latencies CNO
        mice{mouse_count, 7} = [];  % Time;Percent Pre (empty)
        mice{mouse_count, 8} = [];  % Time;Percent Test (empty)
        
        % Validate CNO data
        if isempty(mice{mouse_count, 2})
            warning('Mouse %s: Missing food position for CNO session', mouse_name);
        end
        if isempty(mice{mouse_count, 5})
            warning('Mouse %s: Missing starting position for CNO session', mouse_name);
        end
        if isempty(mice{mouse_count, 6})
            warning('Mouse %s: Missing latencies for CNO session', mouse_name);
        end
        
        fprintf('  → Added CNO session: %s\n', mice{mouse_count, 1});
        fprintf('    Food position: [%.1f, %.1f]\n', mice{mouse_count, 2}(1), mice{mouse_count, 2}(2));
        fprintf('    Pre data: %dx%d, Test data: %dx%d\n', ...
                size(mice{mouse_count, 3}, 1), size(mice{mouse_count, 3}, 2), ...
                size(mice{mouse_count, 4}, 1), size(mice{mouse_count, 4}, 2));
        if ~isempty(mice{mouse_count, 5})
            fprintf('    Start position: [%.1f, %.1f]\n', mice{mouse_count, 5}(1), mice{mouse_count, 5}(2));
        end
        if ~isempty(mice{mouse_count, 6})
            fprintf('    Latencies: [%.1f, %.1f] seconds\n', mice{mouse_count, 6}(1), mice{mouse_count, 6}(2));
        end
    else
        fprintf('  → No CNO data found, saline-only mouse\n');
    end
    
    fprintf('\n');
end

%% Summary
fprintf('=== CONVERSION SUMMARY ===\n');
fprintf('Original Cheeseboard mice: %d\n', n_cheeseboard_mice);
fprintf('Total sessions in mice array: %d\n', mouse_count);

% Count session types
saline_sessions = sum(contains(mice(:,1), '_saline'));
cno_sessions = sum(contains(mice(:,1), '_CNO'));

fprintf('  - Saline sessions: %d\n', saline_sessions);
fprintf('  - CNO sessions: %d\n', cno_sessions);
fprintf('  - Mice with both sessions: %d\n', min(saline_sessions, cno_sessions));
fprintf('  - Saline-only mice: %d\n', saline_sessions - min(saline_sessions, cno_sessions));

%% Validation
fprintf('\n=== DATA VALIDATION ===\n');

% Check for any mice with missing critical data
missing_data_count = 0;
for i = 1:size(mice, 1)
    mouse_name = mice{i,1};
    issues = {};
    
    if isempty(mice{i,2})
        issues{end+1} = 'missing food position';
    end
    if isempty(mice{i,3}) || isempty(mice{i,4})
        issues{end+1} = 'missing behavioral data';
    end
    if isempty(mice{i,5})
        issues{end+1} = 'missing start position';
    end
    if isempty(mice{i,6})
        issues{end+1} = 'missing latencies';
    end
    
    if ~isempty(issues)
        fprintf('WARNING - %s: %s\n', mouse_name, strjoin(issues, ', '));
        missing_data_count = missing_data_count + 1;
    end
end

if missing_data_count == 0
    fprintf('✓ All mice have complete data\n');
else
    fprintf('⚠ %d mice have missing data (see warnings above)\n', missing_data_count);
end

%% Display final mice array structure
fprintf('\n=== FINAL MICE ARRAY STRUCTURE ===\n');
fprintf('Size: %dx%d\n', size(mice, 1), size(mice, 2));
fprintf('Columns:\n');
fprintf('  1: Mouse name (with _saline or _CNO suffix)\n');
fprintf('  2: [x,y] food position\n');
fprintf('  3: Pre-test behavioral data (Nx10 double)\n');
fprintf('  4: Test behavioral data (Nx10 double)\n');
fprintf('  5: [x0,y0] starting position\n');
fprintf('  6: Latencies [pre;test]\n');
fprintf('  7: Time;Percent Pre (empty - for future use)\n');
fprintf('  8: Time;Percent Test (empty - for future use)\n');

fprintf('\nFirst few entries:\n');
for i = 1:min(5, size(mice, 1))
    fprintf('  %s: Food [%.1f,%.1f], Data %dx%d & %dx%d\n', ...
            mice{i,1}, mice{i,2}(1), mice{i,2}(2), ...
            size(mice{i,3}, 1), size(mice{i,3}, 2), ...
            size(mice{i,4}, 1), size(mice{i,4}, 2));
end

if size(mice, 1) > 5
    fprintf('  ... and %d more\n', size(mice, 1) - 5);
end

fprintf('\nData reorganization completed!\n');
fprintf('You can now use this mice array with the behavioral analysis functions.\n');

end