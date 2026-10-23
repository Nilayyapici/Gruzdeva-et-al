%% Data input
clear all

cd % example data dir

% Example data:
tracking = readtable('.csv'); % file name for the top view camera
Analog_table = readtable('.csv'); % file name for analog input from Doric system for synchronization
DLC = readtable('.csv'); % file name for deeplabcut tracking
food_cam = readtable('.csv'); % file name for food camera

time_food = table2array(food_cam(:,1));
time_food = (time_food(:, 1) - time_food(1, 1))/1000; % ms to sec

load('photometry.mat') % from Doric
load('ain.mat') % from Doric

%for color map
rgb = [0 0 1; 0.2 0.2 1; 0.3 0.3 1; 0.5 0.5 1; 1 0.95 1; 1 0.5 0.5; 1 0.3 0.3; 1 0.2 0.2; 1 0 0];

%% Synch TTL from Doric

txy = [];
txy(:,1) = table2array(tracking(:,1)); % time in system time
txy(:,[2:5]) = table2array(DLC(:,[2,3,5,6])); % coordinates Nose, Centroid

Analog_copy = table2array(Analog_table);

Analog_time = Analog_copy(:,1);
Analog_data = Analog_copy(:,2);
Copy_photometry = find(Analog_data <= 3);
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

%% Plotting

figure
hold all
plot(Analog_bonsai(:,1), Analog_bonsai(:,2))
plot(photometry(:,1), photometry(:,3)*100)
plot(ain(:,1), ain(:,2)*100)

%% arena 

x1 = 82; % all the corners starting from left top
y1 = 19;
x2 = 150;
y2 = 20;
x3 = 169;
y3 = 424;
x4 = 261;
y4 = 424;
x5 = 244;
y5 =14;
x6 = 312;
y6 = 12;
x7 = 326;
y7 = 420;
x8 = 420;
y8 = 418;
x9 = 410;
y9 = 12;
x10 = 481;
y10 = 7;
x11 = 493;
y11 = 488;
x12 = 101;
y12 = 498;

x1_door = 415;
y1_door = 100;
x2_door = 483;
y2_door = 99;

x_food = 440;
y_food = 75;

pix_to_cm =0.17;

%% Track in the arena
figure;
plot(txy(1:end/4, 4), txy(1:end/4, 5));
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

plot([x1_door x2_door], [y1_door y2_door], 'k','LineWidth', 1); 

%% speed

speed = sqrt(diff(txy(:,4)).^2 + diff(txy(:,5)).^2) ./ diff(txy(:,1));;
speed_cm = pix_to_cm*speed; % speed in cm/sec
xy_speed = txy;
xy_speed(2:end, 6) = speed_cm;
xy_speed([1,2], :) = [];
figure;
plot(txy(:, 2), txy(:, 3))
hold on;
scatter(xy_speed(:, 2), xy_speed(:, 3), 15, xy_speed(:, 6), 'filled');
axis([0 400 0 450]);
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

% Remove first 20s of the recording
start_index = findnearest(20, txy_phot_speed(:,1)); % after first 20 sec
txy_phot_speed(1:start_index, :) = []; % removing first 20 sec
figure
plot(txy_phot_speed(:,6))

%% Plotting with manually labeled events

%Fill in for each trial
%events form the notes (from the video):
interaction_food_csv = readtable('.csv'); % food visits file
interaction_food_frames = table2array(interaction_food_csv);

eating_csv = readtable('.csv'); % eating file
eating_frames = table2array(eating_csv);

food_start = time_food(interaction_food_frames(:,1)) - dt_synch;
food_end = time_food(interaction_food_frames(:,2)) - dt_synch;
food_duration = time_food(interaction_food_frames(:,2))-time_food(interaction_food_frames(:,1));

eating_start = time_food(eating_frames(:,1)) - dt_synch;
eating_end = time_food(eating_frames(:,2)) - dt_synch;
eating_duration = time_food(eating_frames(:,2))-time_food(eating_frames(:,1));

% Grooming
grooming_csv = readtable('.csv'); % grooming file
grooming_frames = table2array(grooming_csv);

groom_start = time_tracking(grooming_frames(:,1)) - dt_synch;
groom_end = time_tracking(grooming_frames(:,2)) - dt_synch;
groom_duration = time_tracking(grooming_frames(:,2))-time_tracking(grooming_frames(:,1));

% Fill the frames from video when the door was open/closed
door_remov1 = time_tracking(29114)-dt_synch;
door1 = time_tracking(45542)-dt_synch;
door_remov2 = time_tracking(65288)-dt_synch;
door2 = time_tracking(128468)-dt_synch;

% Plotting

figure

txy_phot_precise = txy_phot_speed;

plot(txy_phot_precise(1:end,1), txy_phot_precise(1:end,6), 'LineWidth',1);
ylabel('dF/F, %', 'Fontsize', 16);
xlabel('Time, s', 'Fontsize', 14);
xticks(0:200:4000);
hold on

plot([1 1]*(door_remov1), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door_remov2), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door1), ylim, 'Color', 'b','lineWidth',2);

for i=1:length(food_start)
    rectangle('Position',[food_start(i) min(txy_phot_precise(:,6)) food_duration(i) 70], 'FaceColor', [0.4 0.7 0.2],'FaceAlpha',0.1,'EdgeColor',[0.4 0.7 0.2 0.3],...
    'LineWidth',0.1); % interaction with food % x y w h
end
% 
for i=1:length(eating_start)
    rectangle('Position',[eating_start(i) min(txy_phot_precise(:,6)) eating_duration(i) 70], 'FaceColor', [1 0 0],'FaceAlpha',0.1,'EdgeColor',[1 0 0 0.3],...
    'LineWidth',0.1); % interaction with food
end

for i=1:length(groom_start)
    g_int =rectangle('Position',[groom_start(i) min(txy_phot_precise(:,6)) groom_duration(i) 70], 'FaceColor', [0 0 0],'FaceAlpha',0.1,'EdgeColor',[0 0 0 0.3],...
    'LineWidth',0.1); % grooming
end

hold off


%% For any onset and offset
% new way to calculate df/f
% 
% The 470 and 405 nm signals were independently processed and normalized to baseline signals to
% determine ∆F/F, where ∆F/F= (F-Fbaseline)/Fbaseline and Fbaseline is the median of pre-stimulus signal. - I use mean!! (Anna) No isosbestic
% normalization was introduced. Data were down-sampled to 1 Hz in MATLAB.

%change time before/after here
time_before_inter = 15;
time_after_inter = 15; % 15 for interaction and grooming 25 for eating

%change for trigger
temp_on = groom_start;
temp_off = groom_end;
temp_duration = groom_duration;
% temp_color = [0.1 0.4 0.2]; % interaction with food
% temp_color = [1 0 0]; % eating
temp_color = [0.4940 0.1840 0.5560]; % grooming

%onset of interaction
index_start = [];
for i = 1:length(temp_on)
    k = findnearest(temp_on(i), txy_phot_precise(:, 1), 0); %index of starting eating
    index_start = [index_start; k];
end

pre_inter = zeros(length(index_start), 1);
pre_inter_405 = zeros(length(index_start), 1);
after_inter = zeros(length(index_start), 1);
after_inter_405 = zeros(length(index_start), 1);
pre_inter_dff = [];
after_inter_dff = [];
pre_inter_dff_405 = [];
after_inter_dff_405 = [];

for i = 1:length(index_start)
    for j = 1:round(time_before_inter/time_step)
        pre_inter(j, i) = txy_phot_precise((index_start(i)-j), 4);
        pre_inter_405(j, i) = txy_phot_precise((index_start(i)-j), 5);
    end
    pre_inter_dff_i = (pre_inter(:,i)-mean(pre_inter(:,i)))/mean(pre_inter(:,i))*100;
    pre_inter_dff = [pre_inter_dff, pre_inter_dff_i];
    pre_inter_dff_i_405 = (pre_inter_405(:,i)-mean(pre_inter_405(:,i)))/mean(pre_inter_405(:,i))*100;
    pre_inter_dff_405 = [pre_inter_dff_405, pre_inter_dff_i_405];
    
    after_inter(1, i) = txy_phot_precise((index_start(i)), 4);
    after_inter_405(1, i) = txy_phot_precise((index_start(i)), 5);
    for jj = 1:round(time_after_inter/time_step)
        if index_start(i)+jj < length(txy_phot_precise)
            after_inter(jj+1, i) = txy_phot_precise((index_start(i)+jj), 4);
            after_inter_405(jj+1, i) = txy_phot_precise((index_start(i)+jj), 5);
        end
    end
    after_inter_dff_i = (after_inter(:,i)-mean(pre_inter(:,i)))/mean(pre_inter(:,i))*100;
    after_inter_dff = [after_inter_dff, after_inter_dff_i];
    after_inter_dff_i_405 = (after_inter_405(:,i)-mean(pre_inter_405(:,i)))/mean(pre_inter_405(:,i))*100;
    after_inter_dff_405 = [after_inter_dff_405, after_inter_dff_i_405];
    
end

around_inter = [flip(pre_inter_dff); after_inter_dff];
around_inter_405 = [flip(pre_inter_dff_405); after_inter_dff_405];
time_around_inter = [time_step:time_step:(length(around_inter)*time_step)]';

%Offset of interaction

index_end = [];
for i = 1:length(temp_off)
    k = findnearest(temp_off(i), txy_phot_precise(:, 1), 0); %index of starting eating
    index_end = [index_end; k];
end

pre_off = zeros(length(index_end), 1);
after_off = zeros(length(index_end), 1);
pre_off_dff = [];
after_off_dff = [];
pre_off_405 = zeros(length(index_end), 1);
after_off_405 = zeros(length(index_end), 1);
pre_off_dff_405 = [];
after_off_dff_405 = [];

for i = 1:length(index_end)
    for j = 1:round(time_before_inter/time_step)
        pre_off(j, i) = txy_phot_precise((index_end(i)-j), 4);
        pre_off_405(j, i) = txy_phot_precise((index_end(i)-j), 5);
    end
    pre_off_dff_i = (pre_off(:,i)-mean(pre_off(:,i)))/mean(pre_off(:,i))*100;
    pre_off_dff = [pre_off_dff, pre_off_dff_i];
    pre_off_dff_i_405 = (pre_off_405(:,i)-mean(pre_off_405(:,i)))/mean(pre_off_405(:,i))*100;
    pre_off_dff_405 = [pre_off_dff_405, pre_off_dff_i_405];
    
    after_off(1, i) = txy_phot_precise((index_end(i)), 4);
    after_off_405(1, i) = txy_phot_precise((index_end(i)), 5);
    
    for jj = 1:round(time_after_inter/time_step)
        if index_end(i)+jj < length(txy_phot_precise)
            after_off(jj+1, i) = txy_phot_precise((index_end(i)+jj), 4);
            after_off_405(jj+1, i) = txy_phot_precise((index_end(i)+jj), 5);
        end
    end
    after_off_dff_i = (after_off(:,i)-mean(pre_off(:,i)))/mean(pre_off(:,i))*100;
    after_off_dff = [after_off_dff, after_off_dff_i];
    after_off_dff_i_405 = (after_off_405(:,i)-mean(pre_off_405(:,i)))/mean(pre_off_405(:,i))*100;
    after_off_dff_405 = [after_off_dff_405, after_off_dff_i_405];
    
end

around_off = [flip(pre_off_dff); after_off_dff];
around_off_405 = [flip(pre_off_dff_405); after_off_dff_405];
time_around_off = [time_step:time_step:(length(around_off)*time_step)]';

figure
subplot(1,2,1)
hold all
plot(time_around_inter-time_before_inter, mean(around_inter, 2),'Color', temp_color,'lineWidth', 2.5);
plot(time_around_inter-time_before_inter, mean(around_inter_405, 2),'Color', [0 0 0],'lineWidth', 1.5);
plot(time_around_inter-time_before_inter, around_inter, 'Color', [0, 0, 0, 0.04], 'lineWidth',1.5)

plot([1 1]*((length(pre_inter)+1)*time_step-time_before_inter), [-10 20], '--r','lineWidth',1.5);

% rectangle('Position',[0 -2 max(food_duration) 10], 'FaceColor', [1 0 0 0.05],'EdgeColor',[1 0 0 0.3],...
%     'LineWidth',0.1) % position - x y w h

rectangle('Position',[0 -10 mean(temp_duration) 100], 'FaceColor', [1 0 0 0.1],'EdgeColor',[0 0 0 0],...
    'LineWidth',0.1) % position - x y w h

curve1 = mean(around_inter, 2) + (std(around_inter,0,2)/sqrt(length(temp_on)));
curve2 = mean(around_inter, 2) - (std(around_inter,0,2)/sqrt(length(temp_on)));
patch([(time_around_inter-time_before_inter)' fliplr((time_around_inter-time_before_inter)')], ...
    [curve1' fliplr(curve2')], temp_color, 'Facealpha', 0.2, 'Edgecolor', 'None');

curve1_405 = mean(around_inter_405, 2) + (std(around_inter_405,0,2)/sqrt(length(temp_on)));
curve2_405 = mean(around_inter_405, 2) - (std(around_inter_405,0,2)/sqrt(length(temp_on)));
patch([(time_around_inter-time_before_inter)' fliplr((time_around_inter-time_before_inter)')], ...
    [curve1_405' fliplr(curve2_405')], [0 0 0], 'Facealpha', 0.2, 'Edgecolor', 'None');

xlabel('Time, sec', 'FontSize', 14);
ylabel('dff, %', 'FontSize', 14);
ylim([-10, 10]);
xlim([(-1)*(time_before_inter+5),time_after_inter+5]);
xticks((-1)*time_before_inter:5:time_after_inter)
leg = legend('465', '405', 'Fontsize', 14);
leg.ItemTokenSize = [12,12];
legend('box', 'off');

title('Onset', 'FontSize', 14)

subplot(1,2,2)
hold all
plot(time_around_off-time_before_inter, mean(around_off, 2),'Color', temp_color,'lineWidth', 2.5);
plot(time_around_off-time_before_inter, mean(around_off_405, 2),'Color', [0 0 0],'lineWidth', 1.5);
plot(time_around_off-time_before_inter, around_off, 'Color', [0, 0, 0, 0.04], 'lineWidth',1.5)

plot([1 1]*((length(pre_off)+1)*time_step-time_before_inter), [-20 10], '--k','lineWidth',1.5);

rectangle('Position',[0-mean(temp_duration) -20 mean(temp_duration) 100], 'FaceColor', [1 0 0 0.1],'EdgeColor',[0 0 0 0],...
    'LineWidth',0.1) % position - x y w h

curve1 = mean(around_off, 2) + (std(around_off,0,2)/sqrt(length(temp_off)));
curve2 = mean(around_off, 2) - (std(around_off,0,2)/sqrt(length(temp_off)));
patch([(time_around_off-time_before_inter)' fliplr((time_around_off-time_before_inter)')], [curve1' fliplr(curve2')], ...
    temp_color, 'Facealpha', 0.2, 'Edgecolor', 'None');
curve1_405 = mean(around_off_405, 2) + (std(around_off_405,0,2)/sqrt(length(temp_off)));
curve2_405 = mean(around_off_405, 2) - (std(around_off_405,0,2)/sqrt(length(temp_off)));
patch([(time_around_off-time_before_inter)' fliplr((time_around_off-time_before_inter)')], [curve1_405' fliplr(curve2_405')], ...
    [0 0 0], 'Facealpha', 0.2, 'Edgecolor', 'None');

xlabel('Time, sec', 'FontSize', 14);
ylabel('dff, %', 'FontSize', 14);
ylim([-10, 10]);
xlim([(-1)*(time_before_inter+5),time_after_inter+5]);
xticks((-1)*time_before_inter:5:time_after_inter)
leg = legend('490', '405');
leg.ItemTokenSize = [12,12];
legend('box', 'off');

title('Offset', 'FontSize', 14)

% sgtitle('Interaction with food')
% sgtitle('Eating')
% sgtitle('Grooming')


%% Reorganization for farther analysis

txy_phot_speed(:,16) = txy_phot_speed(:,8);
txy_phot_speed(:,17) = txy_phot_speed(:,9);

%% binning to n sec

n = 0.2; % 0.2 for one source maze and 0.05 for three sources maze and for cheeseboard

time_step = mean(diff(txy_phot_speed(:,1)));

txy_phot_sp_bin = zeros([round(length(txy_phot_speed)/round(n/time_step))-1, size(txy_phot_speed,2)]);
for i = round(n/time_step):round(n/time_step):length(txy_phot_speed)
    txy_phot_sp_bin(i/round(n/time_step),1) = txy_phot_speed(i,1);
    txy_phot_sp_bin(i/round(n/time_step),[2,3]) = txy_phot_speed(i,[2,3]);
    txy_phot_sp_bin(i/round(n/time_step),[4:end]) = mean(txy_phot_speed(i-round(n/time_step)+1:i,[4:end]));
end

%% Check the track in the arena

figure
hold on
scatter(txy_phot_sp_bin(:, 2), txy_phot_sp_bin(:, 3), 5, txy_phot_sp_bin(:, 6), 'filled'); 

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

txy_zones = txy_phot_sp_bin;
for i = 1:length(txy_zones)
    if (txy_zones(i, 3) - a1*txy_zones(i, 2)) <= b1 && (txy_zones(i, 3) - a2*txy_zones(i, 2)) ...
            >= b2 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7
        txy_zones(i, 8) = 1; % left arm
    elseif (txy_zones(i, 3) - a3*txy_zones(i, 2)) <= b3 && (txy_zones(i, 3) - a4*txy_zones(i, 2)) ...
            >= b4 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7
        txy_zones(i, 8) = 2; % center arm
    elseif (txy_zones(i, 3) - a5*txy_zones(i, 2)) <= b5 && (txy_zones(i, 3) - a6*txy_zones(i, 2)) ...
            > b6 && (txy_zones(i, 3) - a7*txy_zones(i, 2)) <= b7 && (txy_zones(i, 3) - a9*txy_zones(i, 2)) > b9
        txy_zones(i, 8) = 3; % right arm (food)
    elseif (txy_zones(i, 3) - a7*txy_zones(i, 2)) > b7 && (txy_zones(i, 3) - a8*txy_zones(i, 2)) ...
            <= b8 && (txy_zones(i, 3) - a6*txy_zones(i, 2)) > b6 && (txy_zones(i, 3) - a1*txy_zones(i, 2)) <= b1
        txy_zones(i, 8) = 4;% horizontal corridor
    elseif (txy_zones(i, 3) - a9*txy_zones(i, 2)) <= b9 && (txy_zones(i, 3) - a5*txy_zones(i, 2)) < b5
        txy_zones(i, 8) = 5; % right arm (food)
    else
        txy_zones(i, 8) = NaN;
    end
end

figure
hold all
% plot(txy_phot(:, 2), txy_phot(:, 3));
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 8), 'filled');
title('zones')

%% Find distance to food (shortest) - 9th column

length_right_pix = sqrt((x8-x_food).^2 + (y8-y_food).^2);
length_horizontal_pix = sqrt((x11-x12).^2 + (y11-y12).^2);

for i = 1:length(txy_zones)
    if txy_zones(i, 8) == 3 || txy_zones(i, 8) == 5 % right arm with food
        txy_zones(i, 9) = pix_to_cm*sqrt((txy_zones(i, 2)-x_food).^2 + (txy_zones(i, 3)-y_food).^2); % a distance from the coordinat to center of food
    elseif txy_zones(i, 8) == 4 % horizontal arm
        txy_zones(i, 9) = pix_to_cm*(length_right_pix + sqrt((txy_zones(i, 2)-x8).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + distance from the coordinate to the beginning of the horizontal arm
    elseif txy_zones(i, 8) == 1 % left arm
        txy_zones(i, 9) = pix_to_cm*(length_right_pix + length_horizontal_pix + sqrt((txy_zones(i, 2)-x3).^2 + (txy_zones(i, 3)-y8).^2)); % length of the right arm + length of the horizontal arm + distance from the coordinate to the beginning of the right arm
    elseif txy_zones(i, 8) == 2 % center arm
        txy_zones(i, 9) = pix_to_cm*(length_right_pix + 1/2*length_horizontal_pix + sqrt((txy_zones(i, 2)-x7).^2 + (txy_zones(i, 3)-y7).^2)); % length of the right arm + half of the length of the horizontal arm + distance from the coordinate to the beginning of the central arm        
    elseif txy_zones(i, 8) == 0 % center arm
        txy_zones(i, :) = NaN;
    end
end

txy_zones = rmmissing(txy_zones); % removing all the NaNs

figure
scatter(txy_zones(:, 2), txy_zones(:, 3), 5, txy_zones(:, 9), 'filled')

%% Door open(1)/closed(0) - 10th column

door_remov1_indx = findnearest(door_remov1, txy_zones(:,1));
door1_indx = findnearest(door1, txy_zones(:,1));
door_remov2_indx = findnearest(door_remov2, txy_zones(:,1));

txy_zones(1:door_remov1_indx,10) = 0;
txy_zones(door_remov1_indx:door1_indx,10) = 1;
txy_zones(door1_indx:door_remov2_indx,10) = 0;
txy_zones(door_remov2_indx:end,10) = 1;

%% Interaction with food 0/1 - 11th column

%to find time of the episodes in a newly binned array:
time_food_start = [];
time_food_end = [];

for i = 1:length(food_start)
    time_food_start_i = findnearest(food_start(i), txy_zones(:, 1));
    time_food_start = [time_food_start;time_food_start_i];
    time_food_end_i = findnearest(food_end(i), txy_zones(:, 1));
    time_food_end = [time_food_end;time_food_end_i];
    
    % put the episodes in the table with everything (11th column):
    txy_zones(time_food_start(i):time_food_end(i), 11) = 1;
end

%% Check correlation before the door is open for a first time

txy_before_1food_inter = [];

for i = 1:length(txy_zones)
    if txy_zones(i,11) == 0
        txy_before_1food_inter = [txy_before_1food_inter; txy_zones(i, :)];
    elseif txy_zones(i,11) == 1
        break
    end
end

ind_Iinter = i;
txy_after_1food_inter = txy_zones(ind_Iinter:end,:);
txy_after_1food_inter_Iopen = [];
txy_after_1food_inter_IIclosed = [];

for i = 1:length(txy_after_1food_inter)
    if txy_after_1food_inter(i,10) == 1
        txy_after_1food_inter_Iopen = [txy_after_1food_inter_Iopen; txy_after_1food_inter(i, :)];
    elseif txy_after_1food_inter(i,10) == 0
        break
    end
end

ind_close1 = i;
txy_after_1food_inter_closed = txy_zones(ind_Iinter+ind_close1:end,:);

for i = 1:length(txy_after_1food_inter_closed)
    if txy_after_1food_inter_closed(i,10) == 0
        txy_after_1food_inter_IIclosed = [txy_after_1food_inter_IIclosed; txy_after_1food_inter_closed(i, :)];
    elseif txy_after_1food_inter_closed(i,10) == 1
        break
    end
end

ind_open2 = i;
txy_after_1food_inter_IIopen = txy_zones(ind_Iinter+ind_close1+ind_open2:end,:);

%check:

figure
hold on
plot(txy_before_1food_inter(:,1), txy_before_1food_inter(:,6))
plot(txy_after_1food_inter_Iopen(:,1), txy_after_1food_inter_Iopen(:,6))
plot(txy_after_1food_inter_IIclosed(:,1), txy_after_1food_inter_IIclosed(:,6))
plot(txy_after_1food_inter_IIopen(:,1), txy_after_1food_inter_IIopen(:,6))

plot([1 1]*(door_remov1), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door_remov2), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door1), ylim, 'Color', 'b','lineWidth',2);

%% 

fig = figure;
hold all
subplot(1,2,1)
x = txy_before_1food_inter(:,9);
y = txy_before_1food_inter(:,6);

coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

plot(txy_before_1food_inter(:,9), txy_before_1food_inter(:,6), 'ko', 'MarkerSize', 1)
hold on
plot(xFit, yFit, 'r-', 'LineWidth', 2);
[rho, pval] = corr(txy_before_1food_inter(:,9), txy_before_1food_inter(:,6),'Type','Spearman');
title({'Before 1st food visit', sprintf('Spearman correlation = %f', rho)}, 'FontSize', 16)

subplot(1,2,2)
x = txy_after_1food_inter(:,9);
y = txy_after_1food_inter(:,6);

coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
plot(txy_after_1food_inter(:,9), txy_after_1food_inter(:,6), 'ko', 'MarkerSize', 1)
hold on
plot(xFit, yFit, 'r-', 'LineWidth', 2);
[rho, pval] = corr(txy_after_1food_inter(:,9), txy_after_1food_inter(:,6),'Type','Spearman');
title({'After 1st food visit', sprintf('Spearman correlation = %f', rho)}, 'FontSize', 16)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% xlabel(han,'Distance to food, cm');
xlabel(han,'Distance, cm', 'Fontsize', 16);
ylabel(han,'dFF, %', 'Fontsize', 16);

%% Eating 0/1 - 12 column

%to find time of the episodes in a newly binned array:
time_eating_start = [];
time_eating_end = [];

for i = 1:length(eating_start)
    time_eating_start_i = findnearest(eating_start(i), txy_zones(:, 1));
    time_eating_start = [time_eating_start;time_eating_start_i];
    time_eating_end_i = findnearest(eating_end(i), txy_zones(:, 1));
    time_eating_end = [time_eating_end;time_eating_end_i];
    
    % put the episodes in the table with everything (12th column):
    txy_zones(time_eating_start(i):time_eating_end(i), 12) = 1;
end

%% Grooming 0/1 - 13 column

%to find time of the episodes in a newly binned array:
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

%% Find the path from the interaction with food/eating - 15th column

food_eating_start = sort([time_food_start;time_eating_start]);
food_eating_end = sort([time_food_end;time_eating_end]);

txy_zones(:,14) = zeros(length(txy_zones),1);
txy_zones(:,15) = zeros(length(txy_zones),1);

for j = 1:length(food_eating_start)-1
    for i = food_eating_end(j):food_eating_start(j+1)
        txy_zones(i, 14) = (pix_to_cm*sqrt((txy_zones(i, 2)-txy_zones(i-1, 2)).^2 + (txy_zones(i, 3)-txy_zones(i-1, 3)).^2));
        txy_zones(i, 15) = (pix_to_cm*sqrt((txy_zones(i, 2)-txy_zones(i-1, 2)).^2 + (txy_zones(i, 3)-txy_zones(i-1, 3)).^2))+txy_zones(i-1, 15);
    end
end

for i = food_eating_end(end-1):length(txy_zones)
    txy_zones(i, 14) = (pix_to_cm*sqrt((txy_zones(i, 2)-txy_zones(i-1, 2)).^2 + (txy_zones(i, 3)-txy_zones(i-1, 3)).^2));
    txy_zones(i, 15) = (pix_to_cm*sqrt((txy_zones(i, 2)-txy_zones(i-1, 2)).^2 + (txy_zones(i, 3)-txy_zones(i-1, 3)).^2))+txy_zones(i-1, 15);
end

%% Checking everything

figure
hold on
plot(txy_zones(:,1), txy_zones(:,6), 'LineWidth',1);
plot(txy_zones(:,1), txy_zones(:,11), 'g', 'LineWidth',1);
plot(txy_zones(:,1), txy_zones(:,9)/20, 'r', 'LineWidth',1);
plot(txy_zones(:,1), txy_zones(:,13), 'k', 'LineWidth',1);
plot(txy_zones(:,1), txy_zones(:,15)/800, 'c', 'LineWidth',1);


plot([1 1]*(door_remov1), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door_remov2), ylim, 'Color', [0.95 0.7 0.2],'lineWidth',2); 
plot([1 1]*(door1), ylim, 'Color', 'b','lineWidth',2);


for i=1:length(food_start)
    rectangle('Position',[food_start(i) min(txy_zones(:,6)) food_duration(i) 10], 'FaceColor', [0.4 0.7 0.2 0.2],'EdgeColor',[0.4 0.7 0.2 0.3],...
    'LineWidth',0.1); % interaction with food % x y w h
end

for i=1:length(eating_start)
    rectangle('Position',[eating_start(i) min(txy_zones(:,6)) eating_duration(i) 25], 'FaceColor', [1 0 0 0.1],'EdgeColor',[1 0 0 0.3],...
    'LineWidth',0.1); % eating
end


for i=1:length(groom_start)
    g_int =rectangle('Position',[groom_start(i) min(txy_zones(:,6)) groom_duration(i) 25], 'FaceColor', [0 0 0 0.2],'EdgeColor',[0 0 0 0.3],...
    'LineWidth',0.1); % grooming
end

hold off


%% making a matrix for further analysis

% 1-time
% 2-x
% 3-y
% 4-speed
% 5-dist to food (shortest)
% 6-path from the last visit
% 7-door closed-1/open-2
% 8-food interaction
% 9-eating 0/1
% 10-grooming 0/1
% 11-hanging 0/1
% 12-dff
% 13-x centroid
% 14-y centroid

temp = [txy_zones(:,1:3),txy_zones(:,7),txy_zones(:,9),txy_zones(:,15),txy_zones(:,10),...
    txy_zones(:,11:13), txy_zones(:,6), txy_zones(:,16), txy_zones(:,17)];

%% saving
% cd %path
% F9I_GLM_fasted_food = temp;
% save('F9I_GLM_fasted_food.mat', 'F9I_GLM_fasted_food')

% 1-time
% 2-x
% 3-y
% 4-speed
% 5-dist to food (shortest)
% 6-path from the last visit
% 7-door closed-1/open-2
% 8-food interaction
% 9-eating 0/1
% 10-grooming 0/1
% 11-dff
% 12-x-centroid
% 13-y-centroid