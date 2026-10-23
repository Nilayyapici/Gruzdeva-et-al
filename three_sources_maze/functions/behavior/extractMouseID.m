%% Helper function to extract mouse ID from session_info
function mouse_id = extractMouseID(session_info)
    % Extract mouse ID from session_info string
    parts = strsplit(session_info, '_');
    if length(parts) >= 1
        mouse_id = parts{1};
    else
        mouse_id = session_info;
    end
end
