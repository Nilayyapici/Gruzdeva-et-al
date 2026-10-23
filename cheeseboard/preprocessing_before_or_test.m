%% Data input
% Fill in separately for Before (Pre) or Test
cd C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\Cheeseboard\M22\M22_saline_solid_walls\Pre

clear all
% Example
tracking = readtable('M22_pre2024-09-18T15_42_09.csv'); % write a name of the file with time and x,y centroid; folder should be in the matlab path
Analog_table = readtable('M22_pre_AI2024-09-18T15_42_09.csv');
DLC = readtable('M22_pre2024-09-18T15_42_10DLC_resnet50_cheeseboardJun4shuffle1_200000.csv');

load('photometry.mat')
load('ain.mat')

%% Synch TTL from Doric

txy = [];
txy(:,1) = table2array(tracking(:,1)); % time in system time
txy(:,2:5) = table2array(DLC(:,[2,3,5,6])); % coordinates Nose, Centroid
Analog_copy = table2array(Analog_table);
Analog_time = Analog_copy(:,1);
Analog_data = Analog_copy(:,2);
Copy_photometry = find(Analog_data <= 2);
Analog_time(Copy_photometry)=[];
Analog_data(Copy_photometry)=[];
Analog_bonsai= [Analog_time,Analog_data];
start_phot = Analog_bonsai(2,1);
frame_beh = findnearest(start_phot, txy(:,1)); % find fist timepoint when photometry starts
dt_synch = (start_phot - txy(1, 1))/1000;
txy(:, 1) = (txy(:, 1) - txy(1, 1))/1000; % converting time of the system to ms and then ms to sec
Analog_bonsai(:, 1) = (Analog_bonsai(:, 1) - Analog_bonsai(1, 1))/1000;
time_tracking = txy(:,1);
txy(1:frame_beh-1, :) = []; % remove all frames before photometry started
txy(:, 1) = (txy(:, 1) - txy(1, 1)); % start from 0

%% Correction for the photometry cumulative delay

photometry(:,1) = photometry(:,1)*1.000138;
ain(:,1) = ain(:,1)*1.000138;

figure
hold all
plot(Analog_bonsai(:,1), Analog_bonsai(:,2))
plot(photometry(:,1), photometry(:,3)*100)
plot(ain(:,1), ain(:,2)*100)

%% coordinates of the arena

x0 = 336;% coordinates of center
y0 = 262;

x_right_extern = 563;
r_extern = x_right_extern - x0; % external circle

x_food=380;
y_food=320;

r_food = 20;

pix_to_cm =75/400;

%% Track in the arena
figure;
hold on
scatter(txy(1:end, 2), txy(1:end, 3), 5, txy(1:end, 1), 'filled'); 

extern = [x0-r_extern, y0-r_extern, 2*r_extern, 2*r_extern];
rectangle('Position', extern, 'Curvature',[1 1]);

food_loc = [x_food-r_food, y_food-r_food, 2*r_food, 2*r_food];
rectangle('Position', food_loc, 'EdgeColor','r', 'Curvature',[1 1]);

%% speed

speed = sqrt(diff(txy(:,4)).^2 + diff(txy(:,5)).^2) ./ diff(txy(:,1));;
speed_mm = 75/400*speed; % speed in cm/sec
xy_speed = txy;
xy_speed(2:end, 6) = speed_mm;
xy_speed([1,2], :) = [];
figure;
plot(txy(:, 2), txy(:, 3))
hold on
scatter(xy_speed(:, 2), xy_speed(:, 3), 15, xy_speed(:, 6), 'filled');
axis([0 400 0 450]);
% caxis([0 10]);
colormap(jet);
colorbar
figure;
plot(xy_speed(:, 1), xy_speed(:, 6));
xlim([0, 600]);

%% Artefacts of tracking

xy_corrected = xy_speed;
for i = 1:length(xy_speed)
    if xy_speed(i, 6) >= 60 
       xy_corrected(i, :) = NaN;
    end
end

xy_corrected = rmmissing(xy_corrected); % removing all the NaNs

figure;
plot(xy_corrected(:, 4), xy_corrected(:, 5))
hold on;
scatter(xy_corrected(:, 4), xy_corrected(:, 5), 15, xy_corrected(:, 6), 'filled');
axis([0 400 0 450]);
% caxis([0 10]);
colormap(jet);
colorbar
figure;
plot(xy_corrected(:, 1), xy_corrected(:, 6));


% %% manual removing artefacts
% 
% for i = 8410:length(xy_corrected)
%     if xy_corrected(i, 2) > 280
%         xy_corrected(i, 5) = 1;
%     end
% end
% xy_corr = [];
% for i = 1:length(xy_corrected)
%     if xy_corrected(i, 5) ~= 1
%        xy_corr = [xy_corr; xy_corrected(i, :)];
%     end
% end

% %% Heat map
% 
% figure;
% dscatter(xy_speed(:,2), xy_speed(:,3), 'plottype', 'image');
% colormap('jet')

%% Photometry

auto = photometry(:, 3);
gcamp = photometry(:, 2);
time_phot = photometry(:, 1);

% fitting like in LERNER paper
reg = polyfit(auto, gcamp, 1);
%reg_coeff (:,file) = reg;
a = reg(1);
b = reg(2);
controlFit = a.*auto + b; % fiting
dff = (gcamp - controlFit)./controlFit; %this gives deltaF/F % motion and photo
dff = dff * 100; % get %

tautodff = [time_phot, auto, gcamp, dff];

order = 2;
framelen = 51;

tautodff(:,3) = sgolayfilt(tautodff(:, 3),order,framelen);
tautodff(:,2) = sgolayfilt(tautodff(:, 2),order,framelen);
tautodff(:,4) = dff;

figure
subplot(2,1,1);
plot(tautodff(:,1), tautodff(:,2), 'b', tautodff(:,1), tautodff(:,3), 'r');
legend('auto', 'GCaMP')
legend('boxoff');
ylabel('Voltage, V');

subplot(2,1,2);
plot(tautodff(:,1), tautodff(:,4));
ylabel('dff, %');
xlabel('Time, sec');

%% Time step

time_step = mean(diff(tautodff(:,1)));

%% interpolation

txy_phot_speed = [];
interp_coord = interp1(xy_corrected(:, 1), xy_corrected(:, [2,3,6]), tautodff(:,1));
txy_phot_speed = time_phot;
txy_phot_speed(:, [2,3,7]) = interp_coord; %2-x, 3-y, 7- speed
txy_phot_speed(:, 4) = tautodff(:,3); % gcamp
txy_phot_speed(:, 5) = tautodff(:,2); % auto
txy_phot_speed(:, 6) = tautodff(:,4); % dff
interp_coord_centroid = interp1(xy_corrected(:, 1), xy_corrected(:, [4,5]), tautodff(:,1));
txy_phot_speed(:, [8 9]) = interp_coord_centroid; % x,y of centroid

start_index = findnearest(20, txy_phot_speed(:,1)); % after first 20 sec
txy_phot_speed(1:start_index, :) = []; % removing first 20 sec
figure
plot(txy_phot_speed(:,6))

%% binning to n sec

n = 0.05;

time_step = mean(diff(txy_phot_speed(:,1)));

txy_phot_sp_bin = zeros([round(length(txy_phot_speed)/round(n/time_step))-1, size(txy_phot_speed,2)]);
for i = round(n/time_step):round(n/time_step):length(txy_phot_speed)
    txy_phot_sp_bin(i/round(n/time_step),1) = txy_phot_speed(i,1);
    txy_phot_sp_bin(i/round(n/time_step),[2,3]) = txy_phot_speed(i,[2,3]);
    txy_phot_sp_bin(i/round(n/time_step),4:end) = mean(txy_phot_speed(i-round(n/time_step)+1:i,4:end));
end

%% Check the track in the arena
txy_zones = txy_phot_sp_bin;

txy_zones_norm = txy_zones;
txy_zones_norm(:,6) = (txy_zones(:,6)-min(txy_zones(:,6)))...
    /(max(txy_zones(:,6))-min(txy_zones(:,6)))*100;

figure
% plot(txy_phot_sp_bin(:, 2), txy_phot_sp_bin(:, 3));
hold on
scatter(txy_zones_norm(:, 2), txy_zones_norm(:, 3), 7, txy_zones_norm(:, 6), 'filled'); 
% scatter(txy_phot_sp_bin(end/2:end, 2), txy_phot_sp_bin(end/2:end, 3), 5, txy_phot_sp_bin(end/2:end, 6), 'filled'); 
% colormap(multigradient(rgb));
cal = colormap(turbo);
caxis([10,80])
axis off

rectangle('Position', extern, 'Curvature',[1 1]);
rectangle('Position', food_loc, 'EdgeColor','r','LineWidth', 1.2,'Curvature',[1 1]);

%%
% grooming_csv = readtable('Grooming.csv');
grooming_frames = [1597,1811];

groom_start = time_tracking(grooming_frames(:,1)) - dt_synch;
groom_end = time_tracking(grooming_frames(:,2)) - dt_synch;
groom_duration = time_tracking(grooming_frames(:,2))-time_tracking(grooming_frames(:,1));

time_groom_start = [];
time_groom_end = [];

for i = 1:length(groom_start)
    time_groom_start_i = findnearest(groom_start(i), txy_zones(:, 1));
    time_groom_start = [time_groom_start;time_groom_start_i];
    time_groom_end_i = findnearest(groom_end(i), txy_zones(:, 1));
    time_groom_end = [time_groom_end;time_groom_end_i];
    
    % put the episodes in the table with everything (13th column):
    txy_zones(time_groom_start(i):time_groom_end(i), 13) = 1;
end


figure
hold on
plot(txy_zones(:,1),txy_zones(:,6))
for i=1:length(groom_start)
    g_int =rectangle('Position',[groom_start(i) min(txy_zones(:,6)) groom_duration(i) 26], 'FaceColor', [0 0 0],'FaceAlpha', 0.2,'EdgeColor',[0 0 0 0.3],...
    'LineWidth',0.1); % grooming
end

%% Distance too food - 9th column
data = txy_zones;
for i = 1:length(data)
    txy_zones(i,9) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2);
end

figure
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 9), 'filled');

%% Time spent in the zones - 8th column

% x_food=x_food_cno;
% y_food=y_food_cno;

% Calculate the radius of the inner circle (distance from the arena center to the food center)
radius_inner = sqrt((x_food - x0)^2 + (y_food - y0)^2);

% Angle to the food center relative to the circle's center
angle_food = atan2(y_food - y0, x_food - x0);

% Calculate the angles for area 2 and area 3
angle_area2 = angle_food + 2*pi/3;
angle_area3 = angle_food - 2*pi/3;

% Calculate the coordinates for area 2
x_area2 = x0-30 + radius_inner * cos(angle_area2);
y_area2 = y0 + radius_inner * sin(angle_area2);

% Calculate the coordinates for area 3
x_area3 = x0 + radius_inner * cos(angle_area3);
y_area3 = y0-30 + radius_inner * sin(angle_area3);

% Define the radius of the zones (Change)
radius = 84;

% Calculate the distance from the mouse to the center of each zone
distance_to_food = sqrt((txy_zones(:, 2) - x_food).^2 + (txy_zones(:, 3) - y_food).^2);
distance_to_area2 = sqrt((txy_zones(:, 2) - x_area2).^2 + (txy_zones(:, 3) - y_area2).^2);
distance_to_area3 = sqrt((txy_zones(:, 2) - x_area3).^2 + (txy_zones(:, 3) - y_area3).^2);

% Initialize a new column for zone labels
txy_zones(:, 8) = 0; % Assuming column 8 is currently empty or can be overwritten

% Assign labels based on proximity to the zone centers within the specified radius
for i = 1:length(txy_zones)
    if distance_to_food(i) <= radius
        txy_zones(i, 8) = 1; % Food zone
    elseif distance_to_area2(i) <= radius
        txy_zones(i, 8) = 2; % Area 2
    elseif distance_to_area3(i) <= radius
        txy_zones(i, 8) = 3; % Area 3
    else
        txy_zones(i, 8) = 0; % Outside all defined zones
    end
end


% Initialize time counters for each zone
time_in_food_zone = 0;
time_in_area2 = 0;
time_in_area3 = 0;

% Calculate time differences
time_differences = diff(txy_zones(:, 1));

% Append an additional time difference for the last entry to keep the length consistent
time_differences = [time_differences; 0];

% Sum the time spent in each zone
for i = 1:length(txy_zones)
    if txy_zones(i, 8) == 1
        time_in_food_zone = time_in_food_zone + time_differences(i);
    elseif txy_zones(i, 8) == 2
        time_in_area2 = time_in_area2 + time_differences(i);
    elseif txy_zones(i, 8) == 3
        time_in_area3 = time_in_area3 + time_differences(i);
    end
end

figure('Position', [100, 50, 800, 600]); % Set the size of the figure to 800x600 pixels
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 8), 'filled');
axis off
colormap([0.6 0.6 0.6; 1 0 0; 0 1 0; 0 0 1]); % Red for zone 1, Green for zone 2, Blue for zone 3
% colorbar;
title('Mouse Trajectory and Zone Times', 'FontSize', 14);
% Add text annotations
x_limits = xlim;
y_limits = ylim;

% Position the text annotations in the plot
text_x_pos = x_limits(1) - 0.05 * (x_limits(2) - x_limits(1));
text_y_pos = y_limits(2) - 0.02 * (y_limits(2) - y_limits(1));

text(text_x_pos, text_y_pos, sprintf('Time in food zone: %.2f s', time_in_food_zone), 'Color', 'red', 'FontSize', 12);
text(text_x_pos, text_y_pos - 0.05 * (y_limits(2) - y_limits(1)), sprintf('Time in area 2: %.2f s', time_in_area2), 'Color', 'green', 'FontSize', 12);
text(text_x_pos, text_y_pos - 0.10 * (y_limits(2) - y_limits(1)), sprintf('Time in area 3: %.2f s', time_in_area3), 'Color', 'blue', 'FontSize', 12);

% Optional: Save the figure
saveas(gcf, 'empty_trajectory_zones.png');

%%
txy_zones_wo_groom = [];

for i = 1:length(txy_zones)
    if txy_zones(i, 13) == 0 && txy_zones(i, 10) == 0
        txy_zones_wo_groom = [txy_zones_wo_groom; txy_zones(i, :)];
    end
end

figure
hold on
plot(txy_zones_wo_groom(:,1),txy_zones_wo_groom(:,6))

%% Save data

cd 'path'

load("Cheeseboard.mat")
%%
txy_zones(:, [10,11,12]) = [];
Cheeseboard{9,2} = [x_food;y_food];
Cheeseboard{9,3} = txy_zones;
%%
save('Cheeseboard.mat', 'Cheeseboard')
