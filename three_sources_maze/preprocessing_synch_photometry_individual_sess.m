%% Data input
% Fill in individually for session 1 (Before+Leraning) and for session 2 (Test)
clear all
% Example:
tracking = readtable('MDRE8_3sources4tr_top_sess12024-06-11T15_45_14.csv'); % write a name of the file with time and x,y centroid; folder should be in the matlab path
Analog_table = readtable('MDRE8_3sources4tr_AI_sess12024-06-11T15_45_14.csv');

DLC = readtable('MDRE8_3sources4tr_top_sess12024-06-11T15_45_15DLC_resnet50_MmazeMay11shuffle1_480000.csv');
% Food cam only for session 1:
food_cam = readtable('MDRE8_foodcam_time2024-06-11T15_45_14.csv');
time_food = table2array(food_cam(:,1));
time_food = (time_food(:, 1) - time_food(1, 1))/1000; % ms to sec

load('photometry.mat')
load('ain.mat')

%% Synch TTL from Doric

txy = [];
txy(:,1) = table2array(tracking(:,1)); % time in system time
txy(:,2:5) = table2array(DLC(:,[2,3,5,6])); % coordinates Nose, Centroid

Analog_copy = table2array(Analog_table);
Analog_time = Analog_copy(:,1);
Analog_data = Analog_copy(:,2);
Copy_photometry = find(Analog_data <= 3);
Analog_time(Copy_photometry)=[];
Analog_data(Copy_photometry)=[];
Analog_bonsai= [Analog_time,Analog_data];
start_phot = Analog_bonsai(2,1);
frame_beh = findnearest(start_phot, txy(:,1)); % find fist timepoint when photometry starts !!?? DELAY 30 frames?
dt_synch = (start_phot - txy(1, 1))/1000;

txy(:, 1) = (txy(:, 1) - txy(1, 1))/1000; % converting time of the system to ms and then ms to sec

Analog_bonsai(:, 1) = (Analog_bonsai(:, 1) - Analog_bonsai(1, 1))/1000;

time_tracking = txy(:,1);
txy(1:frame_beh-1, :) = []; % remove all frames before photometry started
txy(:, 1) = (txy(:, 1) - txy(1, 1)); % start from 0
%% Correction for the photometry cumulative delay

photometry(:,1) = photometry(:,1)*1.000138;
ain(:,1) = ain(:,1)*1.000138;

%% Plotting

figure
hold all
plot(Analog_bonsai(:,1), Analog_bonsai(:,2))
plot(photometry(:,1), photometry(:,3)*100)
plot(ain(:,1), ain(:,2)*100)

%% coordinates of the arena

x1 = 78; % all the corners starting from left top
y1 = 13;
x2 = 144;
y2 = 16;
x3 = 158;
y3 = 414;
x4 = 245;
y4 = 413;
x5 = 244;
y5 = 13;
x6 = 310;
y6 = 16;
x7 = 310;
y7 = 419;
x8 = 401;
y8 = 416;
x9 = 413;
y9 = 17;
x10 = 481;
y10 = 14;
x11 = 468;
y11 = 495;
x12 = 85;
y12 = 485;

x1_door_right = 408;
y1_door_right = 79;
x2_door_right = 470;
y2_door_right = 175;

x1_door_center = 245;
y1_door_center = 66;
x2_door_center = 310;
y2_door_center = 63;

x1_door_left = 73;
y1_door_left = 79;
x2_door_left = 139;
y2_door_left = 79;

% r_cent = 3*r_of/4; % if walls are 10 cm  
x_food_right = 437;
y_food_right = 64;

x_food_center = 283;
y_food_center = 53;

x_food_left = 127;
y_food_left = 61;

r_food = 14; %radius of food

pix_to_cm =75/400;

%% Track in the arena
figure;
plot(txy(1:end, 2), txy(1:end,3));
hold on

plot([x1 x2], [y1 y2], 'k','LineWidth', 1); % low close left wall
plot([x2 x3], [y2 y3], 'k','LineWidth', 1); % low close left wall
plot([x3 x4], [y3 y4], 'k', 'LineWidth', 1); % low close right wall
plot([x4 x5], [y4 y5], 'k','LineWidth', 1); % low close left wall
plot([x5 x6], [y5 y6], 'k', 'LineWidth', 1); % low close low wall
plot([x6 x7], [y6 y7], 'k','LineWidth', 1); % low close left wall
plot([x7 x8], [y7 y8], 'k','LineWidth', 1); % upper close left wall
plot([x8 x9], [y8 y9], 'k','LineWidth', 1); % low close left wall
plot([x9 x10], [y9 y10], 'k', 'LineWidth', 1); % upper close right wall
plot([x10 x11], [y10 y11], 'k','LineWidth', 1); % low close left wall
plot([x11 x12], [y11 y12], 'k', 'LineWidth', 1); % upper close upper wall
plot([x12 x1], [y12 y1], 'k','LineWidth', 1);

% plot([x1_door x2_door], [y1_door y2_door], 'k','LineWidth', 1); 

%% speed

speed = sqrt(diff(txy(:,4)).^2 + diff(txy(:,5)).^2) ./ diff(txy(:,1));;
speed_mm = 75/400*speed; % speed in cm/sec
xy_speed = txy;
xy_speed(2:end, 6) = speed_mm;
xy_speed([1,2], :) = [];
figure;
plot(txy(:, 2), txy(:, 3))
hold on;
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
    if xy_speed(i, 6) >= 80 
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
    txy_phot_sp_bin(i/round(n/time_step),4:end) = mean(txy_phot_speed(i-round(n/time_step)+1:i,[4:end]));
end

%% Check the track in the arena

x1_door = x1_door_center;
x2_door = x2_door_center;
y1_door = y1_door_center;
y2_door = y2_door_center;

figure
% plot(txy_phot_sp_bin(:, 2), txy_phot_sp_bin(:, 3));
hold on
scatter(txy_phot_sp_bin(1:end/2-300, 2), txy_phot_sp_bin(1:end/2-300, 3), 5, txy_phot_sp_bin(1:end/2-300, 6), 'filled'); 

plot([x1 x2], [y1 y2], 'k','LineWidth', 1); % low close left wall
plot([x2 x3], [y2 y3], 'k','LineWidth', 1); % low close left wall
plot([x3 x4], [y3 y4], 'k', 'LineWidth', 1); % low close right wall
plot([x4 x5], [y4 y5], 'k','LineWidth', 1); % low close left wall
plot([x5 x6], [y5 y6], 'k', 'LineWidth', 1); % low close low wall
plot([x6 x7], [y6 y7], 'k','LineWidth', 1); % low close left wall
plot([x7 x8], [y7 y8], 'k','LineWidth', 1); % upper close left wall
plot([x8 x9], [y8 y9], 'k','LineWidth', 1); % low close left wall
plot([x9 x10], [y9 y10], 'k', 'LineWidth', 1); % upper close right wall
plot([x10 x11], [y10 y11], 'k','LineWidth', 1); % low close left wall
plot([x11 x12], [y11 y12], 'k', 'LineWidth', 1); % upper close upper wall
plot([x12 x1], [y12 y1], 'k','LineWidth', 1);

plot([x1_door x2_door], [y1_door y2_door], 'k','LineWidth', 2); 

%% Defining zones - 8th column
% The food is in the opposite arm (left) - change in the code

coefficients1 = polyfit([x1, x12], [y1, y12], 1);
a1 = coefficients1(1); % left side out
b1 = coefficients1(2);
coefficients2 = polyfit([x2, x3], [y2, y3], 1);
a2 = coefficients2(1); % left inner
b2 = coefficients2(2);
coefficients3 = polyfit([x4, x5], [y4, y5], 1);
a3 = coefficients3(1); % left center 
b3 = coefficients3(2);
coefficients4 = polyfit([x6, x7], [y6, y7], 1);
a4 = coefficients4(1); % right center
b4 = coefficients4(2);
coefficients5 = polyfit([x8, x9], [y8, y9], 1);
a5 = coefficients5(1); % right inner
b5 = coefficients5(2);
coefficients6 = polyfit([x10, x11], [y10, y11], 1);
a6 = coefficients6(1); % right outside 
b6 = coefficients6(2);
coefficients7 = polyfit([x3, x8], [y3, y8], 1);
a7 = coefficients7(1); % horizontal inner
b7 = coefficients7(2);
coefficients8 = polyfit([x11, x12], [y11, y12], 1);
a8 = coefficients8(1); % horizontal outside
b8 = coefficients8(2);
coefficients9 = polyfit([x1_door, x2_door], [y1_door, y2_door], 1);
a9 = coefficients9(1); % food zone
b9 = coefficients9(2);
coefficients10 = polyfit([x1, x2], [y1, y2], 1);
a10 = coefficients10(1); % food zone
b10 = coefficients10(2);

txy_zones = txy_phot_sp_bin;
for i = 1:length(txy_zones)
    if (txy_zones(i, 3) - a1*txy_zones(i, 2)) <= b1 && (txy_zones(i, 3) - a2*txy_zones(i, 2)) ...
            >= b2 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7 && (txy_zones(i, 3) - a10*txy_zones(i, 2)) > b10 % added everything after && (for the opposite arm)
        txy_zones(i, 8) = 1; % left arm
    elseif (txy_zones(i, 3) - a3*txy_zones(i, 2)) <= b3 && (txy_zones(i, 3) - a4*txy_zones(i, 2)) ...
             <= b4 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7 && (txy_zones(i, 3) - a9*txy_zones(i, 2)) > b9 
        txy_zones(i, 8) = 2; % center arm
    elseif (txy_zones(i, 3) - a5*txy_zones(i, 2)) >= b5 && (txy_zones(i, 3) - a6*txy_zones(i, 2)) ...
            < b6 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7 && (txy_zones(i, 3) - a5*txy_zones(i, 2)) > b5% && (txy_zones(i, 3) - a9*txy_zones(i, 2)) > b9 % commented for the food in the opposite arm
        txy_zones(i, 8) = 3; % right arm
    elseif (txy_zones(i, 3) - a7*txy_zones(i, 2)) > b7 && (txy_zones(i, 3) - a8*txy_zones(i, 2)) ...
            <= b8 && (txy_zones(i, 3) - a6*txy_zones(i, 2)) < b6 && (txy_zones(i, 3) - a1*txy_zones(i, 2)) <= b1
        txy_zones(i, 8) = 4;% horizontal corridor
    elseif (txy_zones(i, 3) - a9*txy_zones(i, 2)) <= b9 && (txy_zones(i, 3) - a6*txy_zones(i, 2)) <= b6% ...
            %&& (txy_zones(i, 3) - a10*txy_zones(i, 2)) > b10 && (txy_zones(i, 3) - a3*txy_zones(i, 2)) > b3
        txy_zones(i, 8) = 5; % (food)
    else
        txy_zones(i, 8) = NaN;
    end
end

txy_zones_corr =[];
for i = 1:length(txy_zones)
    if  ~(isnan(txy_zones(i, 8)))
        txy_zones_corr = [txy_zones_corr;txy_zones(i,:)];
    end
end       
        
txy_zones = txy_zones_corr;

figure
hold all
% plot(txy_zones(:, 2), txy_zones(:, 3));
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 8), 'filled');
title('zones')

%% Find distance to food (shortest) - 9th column
x_food = x_food_center;
y_food = y_food_center;

length_left_pix = sqrt((x3-x_food).^2 + (y3-y_food).^2);
length_horizontal_pix = sqrt((x11-x12).^2 + (y11-y12).^2);
length_center_pix = sqrt(((x4+(x7-x4)/2)-x_food).^2 + (y4-y_food).^2);

for i = 1:length(txy_zones)
    if txy_zones(i, 8) == 2 || txy_zones(i, 8) == 5 % center arm with food
        txy_zones(i, 9) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2); % a distance from the coordinat to center of food
    elseif txy_zones(i, 8) == 4 % horizontal arm
        txy_zones(i, 9) = pix_to_cm*(length_center_pix + sqrt((txy_zones(i, 2)-(x4+(x7-x4)/2)).^2 + (txy_zones(i, 3)-y3).^2)); % length of the right arm + distance from the coordinate to the beginning of the horizontal arm
    elseif txy_zones(i, 8) == 3  % right arm
        txy_zones(i, 9) = pix_to_cm*(length_center_pix + length_horizontal_pix/2 + sqrt((txy_zones(i, 2)-(x8+(x11-x8)/2)).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + length of the horizontal arm + distance from the coordinate to the beginning of the right arm
    elseif txy_zones(i, 8) == 1 % left arm
        txy_zones(i, 9) = pix_to_cm*(length_center_pix + length_horizontal_pix/2 + sqrt((txy_zones(i, 2)-(x12+(x3-x12)/2)).^2 + (txy_zones(i, 3)-y3).^2)); % length of the right arm + half of the length of the horizontal arm + distance from the coordinate to the beginning of the central arm
    elseif txy_zones(i, 8) == 0 % center arm
        txy_zones(i, :) = NaN;
    end
end


txy_zones = rmmissing(txy_zones); % removing all the NaNs

figure
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 9), 'filled');

%% Door open(1)/closed(0) - 10th column
%%% CHANGE THE DOOR INDEXIES

door_open1 = 19077;
door_remov1 = time_tracking(door_open1)-dt_synch;
door_remov1_indx = findnearest(door_remov1, txy_zones(:,1));

door_closed1 = 25919;
door_add1 = time_tracking(door_closed1)-dt_synch;
door_add1_indx = findnearest(door_add1, txy_zones(:,1));

door_open2 = 30202;
door_remov2 = time_tracking(door_open2)-dt_synch;
door_remov2_indx = findnearest(door_remov2, txy_zones(:,1));

door_closed2 = 38097;
door_add2 = time_tracking(door_closed2)-dt_synch;
door_add2_indx = findnearest(door_add2, txy_zones(:,1));

door_open3 = 42814;
door_remov3 = time_tracking(door_open3)-dt_synch;
door_remov3_indx = findnearest(door_remov3, txy_zones(:,1));

door_closed3 = 48205;
door_add3 = time_tracking(door_closed3)-dt_synch;
door_add3_indx = findnearest(door_add3, txy_zones(:,1));

door_open4 = 52325;
door_remov4 = time_tracking(door_open4)-dt_synch;
door_remov4_indx = findnearest(door_remov4, txy_zones(:,1));

txy_zones(1:door_remov1_indx,12) = 0;
txy_zones(door_remov1_indx:door_add1_indx,12) = 1;
txy_zones(door_add1_indx:door_remov2_indx,12) = 2;
txy_zones(door_remov2_indx:door_add2_indx,12) = 1;
txy_zones(door_add2_indx:door_remov3_indx,12) = 2;
txy_zones(door_remov3_indx:door_add3_indx,12) = 1;
txy_zones(door_add3_indx:door_remov4_indx,12) = 2;
txy_zones(door_remov4_indx:end,12) = 1;


%Check:
figure
plot(txy_zones(1:end,1), txy_zones(1:end,6), 'LineWidth',1);
hold on
plot(txy_zones(1:end,1), txy_zones(1:end,12), 'LineWidth',1);

plot([1 1]*(door_remov1), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door_remov2), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2);
plot([1 1]*(door_remov3), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2);
plot([1 1]*(door_remov4), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2);
plot([1 1]*(door_add1), ylim, 'Color', 'b','lineWidth',2);
plot([1 1]*(door_add2), ylim, 'Color', 'b','lineWidth',2);
plot([1 1]*(door_add3), ylim, 'Color', 'b','lineWidth',2);

%% Check the correlation in both conditions, open/closed door

txy_open_door = [];
txy_closed_door = [];

for i = 1:length(txy_zones)
        if txy_zones(i, 12) == 1 && txy_zones(i,9) > 0 % added a condition to proximity to food
        txy_open_door = [txy_open_door; txy_zones(i, :)];
    elseif txy_zones(i, 12) == 0
        txy_closed_door = [txy_closed_door; txy_zones(i, :)];
    end
end

fig = figure;
hold all

subplot(1,2,1)
plot(txy_closed_door(:,9), txy_closed_door(:,6), 'ko', 'MarkerSize', 2)
axis([0 230 -6 20])% for path
% axis([0 20 -3 5]) % for speed
[rho, pval1] = corr(txy_closed_door(:,9), txy_closed_door(:,6),'Type','Pearson');
title({'Closed door', sprintf('Pearson correlation = %f', rho)}, 'FontSize', 14)

subplot(1,2,2)
x = txy_open_door(:,9);
y = txy_open_door(:,6);

coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

plot(txy_open_door(:,9), txy_open_door(:,6), 'ko', 'MarkerSize', 2)
hold on
plot(xFit, yFit, 'r-', 'LineWidth', 2);
axis([0 230 -6 20]) % for path
% axis([0 20 -3 5]) % for speed
[rho, pval2] = corr(txy_open_door(:,9), txy_open_door(:,6),'Type','Pearson');
title({'Open door', sprintf('Pearson correlation = %f', rho)}, 'FontSize', 14)


han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Distance to food, cm');
% xlabel(han,'Speed, cm/s');
ylabel(han,'dFF, %');
%%
figure
scatter(txy_open_door(1:end, 2), txy_open_door(1:end, 3), 5, txy_open_door(1:end, 6), 'filled');

%% Left arm

x_food = x_food_left;
y_food = y_food_left;

length_left_pix = sqrt(((x12+(x3-x12)/2)-x_food).^2 + (y3-y_food).^2);
length_horizontal_pix = sqrt((x11-x12).^2 + (y11-y12).^2);

for i = 1:length(txy_zones)
    if txy_zones(i, 8) == 1 % left arm with food
        txy_zones(i, 10) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2); % a distance from the coordinat to center of food
    elseif txy_zones(i, 8) == 4 % horizontal arm
        txy_zones(i, 10) = pix_to_cm*(length_left_pix + sqrt((txy_zones(i, 2)-(x12+(x3-x12)/2)).^2 + (txy_zones(i, 3)-y3).^2)); % length of the right arm + distance from the coordinate to the beginning of the horizontal arm
    elseif txy_zones(i, 8) == 3 % right arm
        txy_zones(i, 10) = pix_to_cm*(length_left_pix + length_horizontal_pix + sqrt((txy_zones(i, 2)-(x8+(x11-x8)/2)).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + length of the horizontal arm + distance from the coordinate to the beginning of the right arm
    elseif txy_zones(i, 8) == 2 || txy_zones(i, 8) == 5 % center arm
        txy_zones(i, 10) = pix_to_cm*(length_left_pix + 1/2*length_horizontal_pix + sqrt((txy_zones(i, 2)-(x4+(x7-x4)/2)).^2 + (txy_zones(i, 3)-y4).^2)); % length of the right arm + half of the length of the horizontal arm + distance from the coordinate to the beginning of the central arm        
    elseif txy_zones(i, 8) == 0 % center arm
        txy_zones(i, :) = NaN;
    end
end

txy_zones = rmmissing(txy_zones); % removing all the NaNs

figure
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 10), 'filled');
% colormap(multigradient(rgb));

%% Find distance to the second arm (right)

x_food = x_food_right;
y_food = y_food_right;

length_right_pix = sqrt((x8-x_food).^2 + (y8-y_food).^2);
length_horizontal_pix = sqrt((x11-x12).^2 + (y11-y12).^2);

for i = 1:length(txy_zones)
    if txy_zones(i, 8) == 3  % right arm with food
        txy_zones(i, 11) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2); % a distance from the coordinat to center of food
    elseif txy_zones(i, 8) == 4 % horizontal arm
        txy_zones(i, 11) = pix_to_cm*(length_right_pix + sqrt((txy_zones(i, 2)-x8).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + distance from the coordinate to the beginning of the horizontal arm
    elseif txy_zones(i, 8) == 1 % left arm
        txy_zones(i, 11) = pix_to_cm*(length_right_pix + length_horizontal_pix + sqrt((txy_zones(i, 2)-x3).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + length of the horizontal arm + distance from the coordinate to the beginning of the right arm
    elseif txy_zones(i, 8) == 2 || txy_zones(i, 8) == 5 % center arm
        txy_zones(i, 11) = pix_to_cm*(length_right_pix + 1/2*length_horizontal_pix + sqrt((txy_zones(i, 2)-x7).^2 + (txy_zones(i, 3)-y7).^2)); % length of the right arm + half of the length of the horizontal arm + distance from the coordinate to the beginning of the central arm        
    elseif txy_zones(i, 8) == 0 % center arm
        txy_zones(i, :) = NaN;
    end
end

txy_zones = rmmissing(txy_zones); % removing all the NaNs

figure
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 11), 'filled');
%%
% grooming_csv = readtable('Grooming.csv');
grooming_frames = [1433,1644;6218,6294;8276,8433;13204,13709;21764,21813;25147,25438;39098,39314;44428,45099;57356,57644]; %table2array(grooming_csv);

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
plot(txy_zones(:,1),txy_zones(:,6))
for i=1:length(groom_start)
    g_int =rectangle('Position',[groom_start(i) min(txy_zones(:,6)) groom_duration(i) 26], 'FaceColor', [0 0 0],'FaceAlpha', 0.2,'EdgeColor',[0 0 0 0.3],...
    'LineWidth',0.1); % grooming
end

%% Exclude grooming

txy_zones_wo_groom_c = [];
txy_zones_wo_groom_o = [];

for i = 1:length(txy_zones)
    if txy_zones(i, 12) == 0 && txy_zones(i, 13) == 0
        txy_zones_wo_groom_c = [txy_zones_wo_groom_c; txy_zones(i, :)];
    elseif txy_zones(i, 12) == 1 && txy_zones(i, 13) == 0 && txy_zones(i, 9) >=5
        txy_zones_wo_groom_o = [txy_zones_wo_groom_o; txy_zones(i, :)];
    end
end

fig = figure;
hold all
column = 11;

subplot(1,2,1)
x = txy_zones_wo_groom_c(:,column);
y = txy_zones_wo_groom_c(:,6);

coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

plot(txy_zones_wo_groom_c(:,column), txy_zones_wo_groom_c(:,6), 'ko', 'MarkerSize', 1)
hold on
plot(xFit, yFit, 'r-', 'LineWidth', 2);
axis([0 230 -8 15]) % for path
% axis([0 20 -3 5]) % for speed
[rho, pval] = corr(txy_zones_wo_groom_c(:,column), txy_zones_wo_groom_c(:,6),'Type','Pearson');
title({'No access to food (no groom)', sprintf('Pearson correlation = %f', rho)}, 'FontSize', 14)

subplot(1,2,2)
x = txy_zones_wo_groom_o(:,column);
y = txy_zones_wo_groom_o(:,6);

coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

plot(txy_zones_wo_groom_o(:,column), txy_zones_wo_groom_o(:,6), 'ko', 'MarkerSize', 1)
hold on
plot(xFit, yFit, 'r-', 'LineWidth', 2);
axis([0 230 -8 15]) % for path
% axis([0 20 -3 5]) % for speed
[rho, pval] = corr(txy_zones_wo_groom_o(:,column), txy_zones_wo_groom_o(:,6),'Type','Pearson');
title({'Access to food (no groom)', sprintf('Pearson correlation = %f', rho)}, 'FontSize', 14)



han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Distance to food, cm');
% xlabel(han,'Speed, cm/s');
ylabel(han,'dFF, %');


% %% Saving
% 
% MDRE8_sess1 = txy_zones;
% 
% cd C:\Users\Anna\Dropbox\PhD\Cornell\Nilay_Antonio\Photometry\AgRP\M_maze\3sources_4tr_memory_task_all_mice
% 
% save('MDRE8_sess1.mat', 'MDRE8_sess1')

% 1-time
% 2-x
% 3-y
% 4-465
% 5-405
% 6-dF/F
% 7-speed
% 8-zones
% 9-distance to food
% 10-distance to food diff1 (left for center, center for side)
% 11-distance to food diff2 (right for center, opposite for side)
% 12-door 0 closed/1 open/2 closed for the second time
% 13-grooming