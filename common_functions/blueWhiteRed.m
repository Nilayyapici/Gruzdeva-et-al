function cmap = blueWhiteRed(m)
% BLUEWHITERED creates a blue-to-white-to-red colormap
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a colormap
%   that transitions from dark blue through a small white region to true red

if nargin < 1
   m = 100;
end

% Define the color transitions
colors = [
    0.0, 0.0, 0.5;  % Dark blue
    0.1, 0.3, 0.7;  % Softer blue
    0.4, 0.6, 0.9;  % Medium blue
    0.7, 0.85, 1.0; % Very light blue
    1.0, 1.0, 1.0;  % White (center point)
    1.0, 0.85, 0.85; % Very light red
    0.9, 0.4, 0.4;  % Medium red
    0.8, 0.1, 0.1   % True red
];

n = size(colors, 1);
cmap = zeros(m, 3);

% Create non-uniform spacing to make white region smaller
positions = [0, 0.15, 0.3, 0.45, 0.5, 0.55, 0.7, 1.0];
idx = round(positions * (m-1)) + 1;

% Interpolate colors between key points
for i = 1:n-1
    i1 = idx(i);
    i2 = idx(i+1);
    ni = i2 - i1 + 1;
    
    r = linspace(colors(i,1), colors(i+1,1), ni);
    g = linspace(colors(i,2), colors(i+1,2), ni);
    b = linspace(colors(i,3), colors(i+1,3), ni);
    
    cmap(i1:i2, :) = [r', g', b'];
end