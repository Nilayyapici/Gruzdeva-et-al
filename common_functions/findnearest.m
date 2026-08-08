%% Helper function to find nearest value (replacement for findnearest)
function [idx, val] = findNearest(target, array)
    [~, idx] = min(abs(array - target));
    val = array(idx);
end