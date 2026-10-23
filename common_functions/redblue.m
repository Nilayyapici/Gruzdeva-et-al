function cmap = redblue(n)
    % Create red-blue diverging colormap
    
    if nargin < 1
        n = 256;
    end
    
    % Create diverging colormap: blue -> white -> red
    half = ceil(n/2);
    
    % Blue to white
    r1 = linspace(0, 1, half)';
    g1 = linspace(0, 1, half)';
    b1 = ones(half, 1);
    
    % White to red
    r2 = ones(n-half, 1);
    g2 = linspace(1, 0, n-half)';
    b2 = linspace(1, 0, n-half)';
    
    cmap = [r1, g1, b1; r2, g2, b2];
end