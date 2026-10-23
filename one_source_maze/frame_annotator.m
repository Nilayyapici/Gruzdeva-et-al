function frame_annotator(videoFile, outCsv)
% FRAME_ANNOTATOR  Step through a video frame by frame and mark event onsets/offsets.
%
%   frame_annotator                      % opens a file picker
%   frame_annotator('myvideo.avi')       % saves to myvideo_events.csv
%   frame_annotator('myvideo.avi', 'out.csv')
%
% Controls (click on the figure first so it has keyboard focus):
%   Right / Left arrow   : next / previous frame
%   Shift + arrow        : jump 10 frames
%   Ctrl  + arrow        : jump 100 frames
%   Home / End           : first / last frame
%   O                    : mark ONSET at current frame
%   F                    : mark OFFSET at current frame (completes the event)
%   U                    : undo (clears pending onset, or removes last event)
%   G                    : go to a specific frame number
%   S                    : save CSV (also auto-saves after every change)
%   Q                    : save and quit
%   The slider at the bottom can also be dragged to navigate.
%
% Output CSV: two columns, Onset and Offset (1-based frame numbers),
% one row per event.

    % ---------- Inputs ----------
    if nargin < 1 || isempty(videoFile)
        [f, p] = uigetfile({'*.avi;*.mp4;*.mov', 'Video files'}, 'Select video');
        if isequal(f, 0), return; end
        videoFile = fullfile(p, f);
    end
    if nargin < 2 || isempty(outCsv)
        [p, n] = fileparts(videoFile);
        outCsv = fullfile(p, [n '_events.csv']);
    end

    v = VideoReader(videoFile);
    try
        nFrames = v.NumFrames;
    catch
        nFrames = floor(v.Duration * v.FrameRate);
    end

    % ---------- State ----------
    events  = zeros(0, 2);   % [onset offset]
    pending = [];            % onset waiting for its offset
    k       = 1;             % current frame

    if isfile(outCsv)
        choice = questdlg(sprintf('%s already exists.', outCsv), 'Existing file', ...
                          'Load existing', 'Start fresh', 'Load existing');
        if strcmp(choice, 'Load existing')
            T = readtable(outCsv);
            events = T{:, 1:2};
        end
    end

    % ---------- GUI ----------
    fig = figure('Name', ['Frame annotator - ' videoFile], 'NumberTitle', 'off', ...
                 'Color', 'w', 'KeyPressFcn', @onKey, 'CloseRequestFcn', @onClose);
    ax = axes('Parent', fig, 'Position', [0.03 0.14 0.94 0.78]);
    hImg = image(ax, getFrame(1));
    axis(ax, 'image', 'off');

    sld = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.03 0.06 0.94 0.04], 'Min', 1, 'Max', max(nFrames, 2), ...
        'Value', 1, 'SliderStep', [1 10] / max(nFrames - 1, 1), ...
        'Callback', @onSlider, 'KeyPressFcn', @onKey);

    uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.03 0.005 0.94 0.045], 'BackgroundColor', 'w', ...
        'String', ['<-/->: 1 frame | Shift: 10 | Ctrl: 100 | O: onset | F: offset | ' ...
                   'U: undo | G: go to | S: save | Q: quit']);

    show();

    % ---------- Nested functions ----------
    function img = getFrame(idx)
        img = read(v, idx);
        if size(img, 3) == 1
            img = repmat(img, 1, 1, 3);   % grayscale -> RGB for display
        end
    end

    function goTo(idx)
        k = min(max(round(idx), 1), nFrames);
        show();
    end

    function show()
        hImg.CData = getFrame(k);
        sld.Value = k;

        if isempty(pending), pStr = '-'; else, pStr = num2str(pending); end
        tag = '';
        if any(events(:, 1) == k), tag = [tag '  [ONSET]'];  end
        if any(events(:, 2) == k), tag = [tag '  [OFFSET]']; end
        if ~isempty(pending) && pending == k, tag = [tag '  [PENDING ONSET]']; end

        if isempty(tag), col = 'k'; else, col = 'r'; end
        title(ax, sprintf('Frame %d / %d   (t = %.3f s)   |   events: %d   |   pending onset: %s%s', ...
              k, nFrames, (k - 1) / v.FrameRate, size(events, 1), pStr, tag), ...
              'Color', col, 'Interpreter', 'none');
        drawnow limitrate;
    end

    function onSlider(src, ~)
        goTo(src.Value);
    end

    function onKey(~, evt)
        step = 1;
        if ismember('shift', evt.Modifier), step = 10; end
        if any(ismember({'control', 'command'}, evt.Modifier)), step = 100; end

        switch lower(evt.Key)
            case 'rightarrow', goTo(k + step);
            case 'leftarrow',  goTo(k - step);
            case 'home',       goTo(1);
            case 'end',        goTo(nFrames);
            case 'o',          markOnset();
            case 'f',          markOffset();
            case 'u',          undoLast();
            case 'g',          askFrame();
            case 's',          saveCsv(true);
            case 'q',          onClose();
        end
    end

    function markOnset()
        pending = k;
        show();
    end

    function markOffset()
        if isempty(pending)
            beep; warndlg('Mark an onset (O) before marking an offset (F).', 'No onset');
            return;
        end
        if k < pending
            beep; warndlg(sprintf('Offset (%d) is before onset (%d).', k, pending), 'Invalid offset');
            return;
        end
        events(end + 1, :) = [pending k];
        pending = [];
        saveCsv(false);
        show();
    end

    function undoLast()
        if ~isempty(pending)
            pending = [];
        elseif ~isempty(events)
            events(end, :) = [];
            saveCsv(false);
        end
        show();
    end

    function askFrame()
        a = inputdlg(sprintf('Go to frame (1-%d):', nFrames), 'Go to', 1, {num2str(k)});
        if ~isempty(a)
            val = str2double(a{1});
            if ~isnan(val), goTo(val); end
        end
        figure(fig);   % return keyboard focus to the figure
    end

    function saveCsv(verbose)
        T = array2table(events, 'VariableNames', {'Onset', 'Offset'});
        writetable(T, outCsv);
        if verbose
            fprintf('Saved %d events to %s\n', size(events, 1), outCsv);
        end
    end

    function onClose(~, ~)
        if ~isempty(pending)
            fprintf('Note: unpaired onset at frame %d was not saved.\n', pending);
        end
        saveCsv(true);
        delete(fig);
    end
end