clc
close all %closes all the figures
clear all
%you can type close all in command window to get rid of figures
%Typing whos variable_name to check properties of the variables to debug
%syntax stuff
%now we're automating...he wants us to automate the 4th degree fitting
%note: the 4th degree is best bc its a balance between accuracy and
%efficiency, also prevents overfitting for the algorithm. Go through each
%block one by one

%there's 8 files of cell data, each one has 7 cycles :(
%%
cell_id = '1p4'; %p = point; 1.1-1.4; 2.1-2.4
% Specify the full path to the .mat file
filename = strcat('rpt_cycle_raw\cell_', cell_id ,'_rpt_cycle_raw.mat');

% Load the .mat file using the full path
load(filename);

% load('rpt_cycle_raw/cell_1p1_rpt_cycle_raw.mat')
hppc_idx = find(raw.rpt_type == "hppc");
cycle = raw.super_cycle_count(hppc_idx);

time_hppc = datetime(raw.RealTime_sec_(hppc_idx), 'ConvertFrom', 'posixtime');
voltage_hppc = raw.SampleVoltage_mV_(hppc_idx);
current_hppc = raw.SampleCurrent_mA_(hppc_idx);

figure()
plot(time_hppc,voltage_hppc*(1e-3))
ylabel('voltage')
xlabel('time')
title(strcat('Voltage in all RPT cycles cell ', cell_id));

grid on
figure()
plot(time_hppc,current_hppc*(1e-3))
xlabel('time')
ylabel('current')
grid on
title(strcat('Current in all cycles cell  ', cell_id)); %note in the powerpoint how many cycles there are


%% hppc discharge
cycle_num = 1; %there are 7 cycles total, for some there are only 6, check the graph
cycle = raw.super_cycle_count(hppc_idx);
discharge = raw.StepNo_(hppc_idx);
cycle_cnt = raw.CycleCount(hppc_idx);
% discharge
h_idx = find(cycle == cycle_num & cycle_cnt < 14 & cycle_cnt > 1 & discharge > 4);

t_cycle = time_hppc(h_idx);
v_cycle = voltage_hppc(h_idx)*(1e-3);
i_cycle = current_hppc(h_idx)*(1e-3);
figure()
plot(t_cycle, v_cycle)
xlabel('time')
ylabel('voltage')
grid on
title(strcat('Voltage in cycle ', num2str(cycle_num), ' cell  ', cell_id));
figure()
plot(t_cycle, i_cycle)
xlabel('time')
ylabel('voltage')
grid on
title(strcat('Current in cycle ', num2str( cycle_num), ' cell  ', cell_id));
%% save normalized_Ah
% save('normalized_Ah.mat', 'normalized_Ah')
%%
% save('hppc_v_i_test.mat', 'v_cycle','i_cycle')
%% identify the voltage drops - soc levels

%ok the code works. Depending on the cycle and the cell tho, ur gonna need
%to tune it.

%so locs1 should be the indices and peaks is the actual value. 
%The issue is definitely findpeaks, there's too many peaks being stored...
%The negative in v_cycle makes sense since findpeaks finds maximas but in
%our case we want minimas.
significant_dip_threshold = 0.11; %this isn't a part of the findpeaks function
%i will add fancy parameters to findpeaks i think
%I need to convert t_cycle into 0,1,2,...59084 lol. I can proly use a
%function that's like zeroes but is ascending instead lol.
t_cycledum = [1:size(t_cycle,1)];

%[peaks1, locs1] = findpeaks(-v_cycle); %ok even tho this makes no sense, it works so don't fix it lmao.
%nvm it must be removed lol

%the code doesn't use locs to graph the mins, it only uses it to find the
%points b4 each dip lol
%%[peaks1, locs1, dummy1, dummy2] = findpeaks(-v_cycle, t_cycledum,'MinPeakProminence',0.09, 'Threshold', 0.00002, 'MinPeakWidth', 100);
%ok so for finding the dips, this works vv well. 100% accuracy i would say.
%Adjust the threshold value to get the right points i think
%%peaks1 = abs(peaks1);

%I'm going to make a tweak. Voltage is too complicated, let's use the
%current curve :))))))) hella big brain :))))))
%the advantage is now i can use minheight insted of prominence :)
[peaks1, locs1, dummy1, dummy2] = findpeaks(-i_cycle, t_cycledum,'MinPeakHeight', 30,'MinPeakProminence',20);
locs1 = sort(locs1);
baddd = [];
for i = 2:length(locs1)
    if abs(locs1(i) -  locs1(i-1)) < 20 && locs1(i-1)
        baddd(end+1) = i-1;
    end
end
baddd = unique(baddd);
j = 0;
for i = 1:length(baddd)
    locs1(baddd(i) - j) = [];
    j = j + 1;
end

peaks1 = v_cycle(locs1);

%so we want to make sure peaks and locs are sorted correctly so the
%algorithms later on work well. We also want to make sure no values get
%repeated but im lazy so only add once an error happens lol


%hrmm ok i thought of an algorithm. We know that the dip index is pretty
%close to SOC level index. and we know what the difference is from tuning,
%its 0.11. Let's use those two facts. This is something I'll have to tune
%as well. For now let's say if its within 0.1 of 0.11. Should prolly
%decrease. 
% Initialize array to hold the SOC levels
v_soc_10 = [];
v_idxb4dip = [];

% figure()
% plot(t_cycle, v_cycle,t_cycle(locs1),v_cycle(locs1),'r*') %unnecessary for presenting
% title(['Peaks locations']);
% ylabel('Voltage');
% grid on
%oh ok, if there's ever an error around here, it means the threshold at the
%while loop is too high, you need to lower it
for i = 1:length(locs1)
    %refreshing myself on loops, i = 1 corresponds to the 1st value found
    %in the array lol.
    j = 1;
    while abs(i_cycle(locs1(i) - j) - i_cycle(locs1(i))) <= 20
        %addition bc peaks stored negative values while v is positive
        j = j + 1;
    end
    v_idxb4dip(end+1) = locs1(i) - j;
    v_soc_10(end+1) = v_cycle(locs1(i) - j);
end
%to tune the min distance, choose a reasonable value, then look at the data
%stored and adjust from there


%the code uses local_min to graph
%both locs and local_min are incorrect lol
%I think locs is more useful, so lets keep local_min in the backwater for
%now lol.
local_min = islocalmin(v_cycle);

% Loop over each peak to check if it is followed by a significant dip
%This chunk of code does not work well, I need to change it.
%what if i use findpeaks, and use the identify flat edges function to my
%advantage? wait no ive tried that, i cant even identify the sharp peaks :(


% for i = 1:length(peaks1)
%     % Check if the peak is followed by a significant dip
%     if i < length(peaks1) && abs(peaks1(i)) > abs(peaks1(i+1)) && (abs(peaks1(i)) - abs(peaks1(i+1))) > significant_dip_threshold %(peaks1(i) - voltage_cut(locs1(i)+1) 
%         v_soc_10(end+1) = abs(peaks1(i)); % Get the voltage right before the dip
%         %Yeah i've deemed that v_soc_10 is useless. We only care about the
%         %idx lol.
%         idx(end+1) = locs1(i);
% 
%     end
% end

%what exactly is peaks1 and locs1?
%what exactly is v_soc_10 and idx? 
%Once i figure those things out, debugging should be no issue.

% Output the SOC levels
disp('The voltage of SOC levels before each significant dip are:');
disp(v_idxb4dip);
figure()
subplot(2,1,1)
plot(t_cycle, v_cycle,t_cycle(v_idxb4dip),v_cycle(v_idxb4dip),'o')
title(['Voltage of SOC Levels before Significant Dip']); %I feel like this title is inaccurate lol.
ylabel('Voltage');
grid on
subplot(2,1,2)
plot(t_cycle, i_cycle,t_cycle(v_idxb4dip),i_cycle(v_idxb4dip),'o')
ylabel('Current');
grid on
figure()
plot(t_cycle(v_idxb4dip), v_cycle(v_idxb4dip), '*-')
title(['V levels']);
ylabel('Voltage');
xlabel('Time');
grid on

figure()
plot(t_cycle, v_cycle,t_cycle(locs1),v_cycle(locs1),'r*') %unnecessary for presenting
title(['Peaks locations']);
ylabel('Voltage');
grid on

figure()
plot(t_cycle, i_cycle,t_cycle(locs1),i_cycle(locs1),'r*') %unnecessary for presenting
title(['Peaks locations']);
ylabel('Current');
grid on




%% Calculate the SoC at each point in the cycle using the current

% Calculate the charge (Q) discharged during the cycle
Qcap_cut = trapz(i_cycle); % Charge capacity for the cycle
%trapz is matlabs builtin integration thing

% Initialize an array to store the SoC at each point in the cycle
soc_cut = [];
Q_cut = [];
for i = 1:length(v_idxb4dip)
    Q_cut(i) = trapz(i_cycle(1:v_idxb4dip(i)));
    soc_cut(i) = abs(abs(Qcap_cut) - abs(Q_cut(i))) / abs(Qcap_cut);
end
disp(soc_cut)

% Perform min-max normalization
min_value = min(soc_cut);
max_value = max(soc_cut);
normalized_soc = (soc_cut - min_value) / (max_value - min_value);

% To ensure the maximum value is exactly 1, we can apply a small fix due to precision issues:
normalized_soc = (normalized_soc - min(normalized_soc)) / (max(normalized_soc) - min(normalized_soc));
disp(normalized_soc)
[soc,idxx] = sort(normalized_soc);
ocv = sort(v_soc_10);

[t_soc,i_s] = sort(normalized_soc);
t_ocv = sort(v_soc_10);

% Get coefficients of a line fit through the data.
coefficients = polyfit(t_soc, t_ocv, 4);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(t_soc), max(t_soc), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything

figure()
plot(t_soc,t_ocv, 'r*', 'DisplayName', 'V-SOC point')
hold on;
%{' '} ruins your code...
plot(xFit, yFit, 'b-', 'LineWidth', 1, 'DisplayName', sprintf('Fit: y = %.1fx^4 + %.1fx^3 + %.1fx^2 + %.1fx + %.1f', coefficients(1), coefficients(2), coefficients(3), coefficients(4), coefficients(5))); % Plot fitted line.
title(strcat('V SoC curve of cycle ', num2str( cycle_num), ' cell  ', cell_id));
%for this one, we need to manually add a line of best fit
%if we wanted to, we can run fitline through a function to remove plus
%signs and replace with the negative signs since its been entirely
%transformed into a string
hold off;
legend('show', 'Location', 'southeast')
ylabel('Voltage');
xlabel('SoC');
grid on

%% save capacity Qcapt(t) where t is cycles - incorrect - need to log the cap time
Qcap_1_As = abs(Qcap_cut);
Qcap_1_Ah = Qcap_1_As/3600;
disp('Capacity in Ah: ')
disp(Qcap_1_Ah)   %we need to present this value, figure out a way to automate saving this value
time_start_cycle1 = t_cycle(1);

% % % for cycle two, we would need an algorithm to identify the starting
% index...actually, since h_index is being plugged into time_hppc,
% t_cycle(1) is always going to be the start of the cycle regardless of
% which cycle it is
%fileName = 'Qcap_t.mat';

% % since there's only one Qcap_t.mat, i think we want to have two array
% for each cell. One time array and one Qcap array for the 6-7 cycles.
% Overwrite the file with an empty variable

%save('Qcap_t.mat', strcat('Qcap_cycle', num2str(cycle_num), '_cell  ', cell_id, '_Ah'), strcat('time_start_cycle', num2str(cycle_num), ' cell  ', cell_id),'-append')
%capacity_file = sprintf('cell_%s_capacity.mat', cell_id);
%cell 1 => cycle 1
%cell 7 => cycle 7

% emptyData = [];
% save(capacity_file, 'emptyData');

% Check if the file already exists
% if isfile(capacity_file)
%     % If it exists, load the existing capacities
%     load(capacity_file, 'Qcap_all_cycles');
%     load(capacity_file, 'Time_all_cycles');
% 
% else
%     % If the file doesn't exist, initialize an empty array
%     Qcap_all_cycles = [];
%     Time_all_cycles = datetime.empty;
% end
% % Add the new capacity value for the current cycle
% 
% Qcap_all_cycles(end + 1) = Qcap_1_Ah;
% Time_all_cycles(end + 1) = time_start_cycle1;
% 
% if isfile(capacity_file)
% % Save the updated capacities back to the file
%     save(capacity_file, 'Qcap_all_cycles', 'Time_all_cycles', '-append');
% else
%     save(capacity_file, 'Qcap_all_cycles', 'Time_all_cycles');
% end

%% for clearing specific elements
% Qcap_all_cycles = load(capacity_file, 'Qcap_all_cycles'); 
% Time_all_cycles = load(capacity_file, 'Time_all_cycles'); 
% 
% Qcap_all_cycles = Qcap_all_cycles.Qcap_all_cycles;
% Qcap_all_cycles(4) = [];
% Time_all_cycles = Time_all_cycles.Time_all_cycles;
% Time_all_cycles(4) = [];
% save(capacity_file, 'Qcap_all_cycles', 'Time_all_cycles');

%%
% % Initialize array to hold the SOC levels
% i need to change all of this code, its not robust enough
%good practice, don't change any of nhat's variables unless they're useless
%or you need one he hasn't created

%ok lol i_peaks does absolutely nothing lmao. What Nhat's code does is it
%selects the less prominent or shorter peaks for v_cycle. Then it stores the index
%right after
%i_peaks = i_cycle(local_min);
significant_r_dip = 0.01;
R_soc = [];
idx_r = [];
v_peak_tau_2 = [];
idx_vpk_calc_time_2 = [];


%I need to tune this now after fixing the previous bug
%actually, i need to replace this with my new code that focuses on the
%current curve
[v_peak_tau, idx_vpk_calc_time, dummy1, dummy2] = findpeaks(-v_cycle, t_cycledum,'MinPeakProminence',0.030, 'MinPeakDistance',600);
% 
% idx_vpk_calc_time(2) = []; %there's no shame in hardcoding lmao
% idx_vpk_calc_time(4) = []; %there's no shame in hardcoding lmao
% idx_vpk_calc_time(6) = []; %there's no shame in hardcoding lmao
% idx_vpk_calc_time(12) = []; %there's no shame in hardcoding lmao
% idx_vpk_calc_time(16) = []; %there's no shame in hardcoding lmao
% 

% [v_peak_tau, idx_vpk_calc_time, dummy1, dummy2] = findpeaks(-i_cycle, t_cycledum,'MinPeakHeight', 30,'MinPeakProminence',20);
% v_peak_tau = v_cycle(idx_vpk_calc_time);
%so now i have both the smol dips and the big dips. I already have the big
%dips stored in locs1 so I can easily remove all the big dips v_peak_tau
%via a simple comparison function.
% % 
% figure()
% plot(t_cycledum, v_cycle,t_cycledum(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o',t_cycledum(idx_vpk_calc_time_2),v_cycle(idx_vpk_calc_time_2),'r*' )
% title(['Peaks locations']);
% ylabel('Voltage');
% grid on

%b4 we use set diff, we need to do some rounding for v peak and idx. The
%margin of error for indeces i think is 5 but let's say 10 to be safe. When
%we round we always want to round to the number that's closer to the edge
%ie the greater value.

%ok so we need to deal with repeat values...
%this algorithm below, only use it for unique cases such as cycle 4 for
%cell 1p2
% wait i have an ingenius idea: screw using v_peak_tau at all, let's just
% idx, and stick it into v_cycle and call that v_tau omfg so much easier
% lol
baddd = [];
for i = 2:length(idx_vpk_calc_time)
    if abs(idx_vpk_calc_time(i) -  idx_vpk_calc_time(i-1)) < 20 && idx_vpk_calc_time(i-1)
        baddd(end+1) = i-1;
    end
end
baddd = unique(baddd);
j = 0;
for i = 1:length(baddd)
    idx_vpk_calc_time(baddd(i) - j) = [];
    j = j + 1;
end

% figure()
% plot(t_cycledum, v_cycle,t_cycledum(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o' )
% title(['Peaks locations']);
% ylabel('Voltage');
% grid on


%Actually we know its just going to be every other point do to nature of
%the curve/peaks so we don't need to use find i think
%check each cycle to see whether the big or small dip comes first in the
%graph. Based on the order, you may need to adjust this for loop in terms
%of which index we need to first check ie first or second and the
%subsequent ones.
j = 0;
for i = 1:length(idx_vpk_calc_time)
    if mod(i, 2) == 1
        j = j + 1;
        if j > length(locs1)
            break;
        end
        if idx_vpk_calc_time(i)  ~= locs1(j) %this line of code is unnecessary lol
            if idx_vpk_calc_time(i)  < locs1(j)
                idx_vpk_calc_time(i) = locs1(j);
            elseif idx_vpk_calc_time(i)  > locs1(j)
                locs1(j) = idx_vpk_calc_time(i);
            end
        end
    end
end

idx_vpk_calc_time = setdiff(idx_vpk_calc_time, locs1);

% figure()
% plot(t_cycledum, v_cycle,t_cycledum(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o' )
% title(['Peaks locations']);
% ylabel('Voltage');
% grid on


% j = 0;
% for i = 1:length(v_peak_tau)
%     if mod(i, 2) == 1
%         j = j + 1;
%         if j > length(peaks1)
%             break;
%         end
%         if v_peak_tau(i)  ~= peaks1(j) %this line of code is unnecessary lol
%             if v_peak_tau(i)  < peaks1(j)
%                 v_peak_tau(i) = peaks1(j);
%             elseif v_peak_tau(i)  > peaks1(j)
%                 peaks1(j) = v_peak_tau(i);
%             end
%         end
%     end
% end
v_peak_tau = [];
v_peak_tau = v_cycle(idx_vpk_calc_time);

%these two lines of code indicate we want it to be in acsending order? 
v_peak_tau = v_peak_tau.';
% v_peak_tau = flip(v_peak_tau);

%The double comment below is used for debugging

% figure()
% plot(t_cycle, v_cycle,t_cycle(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o' )
% title(['Peaks locations']);
% ylabel('Voltage');
% grid on

for i = 1:length(idx_vpk_calc_time)
    %refreshing myself on loops, i = 1 corresponds to the 1st value found
    %in the array lol.
    j = 1;
    if idx_vpk_calc_time(i) + j <= length(i_cycle)
        while abs(i_cycle(idx_vpk_calc_time(i) + j) - i_cycle(idx_vpk_calc_time(i))) < 7 || abs(i_cycle(idx_vpk_calc_time(i) + j) - i_cycle(idx_vpk_calc_time(i))) > 11
            %addition bc peaks stored negative values while v is positive
            j = j + 1;
            %i need another if in herebc j is going to increment still the
            %sum will go past the max
            %ok maybe the value that gets stored when the loop breaks out
            %is wrong, but that's okay bc later on we splice the array from
            %1 to 9. 
            if idx_vpk_calc_time(i) + j >= length(i_cycle)
                break;
            end

        end
    end

    v_peak_tau_2(end+1) = v_cycle(idx_vpk_calc_time(i)+j);
    idx_vpk_calc_time_2(end+1) = t_cycledum(idx_vpk_calc_time(i)+j);
end

idx_r = idx_vpk_calc_time;

%delta_v0 = abs(peaks1(i+1))- abs(peaks1(i));
delta_v0 = abs(v_peak_tau_2 - v_peak_tau);

%delta_i = abs(i_cycle(locs1(i+1)))- abs(i_cycle(locs1(i)));
delta_i = abs(i_cycle(idx_vpk_calc_time_2)) - abs(i_cycle(idx_vpk_calc_time));
delta_i = delta_i.';

R_soc = abs(delta_v0 ./ delta_i);

% Loop over each peak to check if it is followed by a significant dip
% for i = 1:length(peaks1)
%     % Check if the peak is followed by a significant dip
%     if i < length(peaks1) && abs(peaks1(i)) < abs(peaks1(i+1)) && (abs(peaks1(i+1)) - abs(peaks1(i))) > significant_r_dip && (abs(peaks1(i+1)) - abs(peaks1(i))) < significant_r_dip*10 && abs(peaks1(i+3)) > abs(peaks1(i+1)) 
%         % && (abs(i_peaks(i+1)) - abs(i_peaks(i))) < 10
%         delta_v0 = abs(peaks1(i+1))- abs(peaks1(i));
%         delta_i = abs(i_cycle(locs1(i+1)))- abs(i_cycle(locs1(i)));
%         R_soc(end+1) = abs(delta_v0 / delta_i);
%         v_peak_tau(end+1) = abs(peaks1(i));
%         v_peak_tau_2(end+1) = abs(peaks1(i+1));
%         idx_vpk_calc_time(end+1) = locs1(i);
%         idx_vpk_calc_time_2(end+1) = locs1(i+1);
%         idx_r(end+1) = locs1(i);
%     end
% end

disp(R_soc)
%0.0034    0.0036    0.0035    0.0035    0.0036    0.0034    0.0032    0.0034    0.0036

disp('1')
disp(v_peak_tau_2)
 % remove below 3V ripples
disp('2')
disp(v_peak_tau)
% disp(idx_vpk_calc_time(1:11))

figure()
plot(t_cycle, v_cycle,t_cycle(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o',t_cycle(idx_vpk_calc_time_2),v_cycle(idx_vpk_calc_time_2),'r*' )
title(['R0 Markers']);
ylabel('Voltage');
grid on

figure()
plot(t_cycle, i_cycle,t_cycle(idx_vpk_calc_time),i_cycle(idx_vpk_calc_time),'o',t_cycle(idx_vpk_calc_time_2),i_cycle(idx_vpk_calc_time_2),'r*' )
title(['R0 Markers']);
ylabel('Current');
grid on


%% calculate the RC discharge
% use soc V value and R soc value

% Initialize array to hold the values
% tau_RC = [];
% idx_tau = [];
v_soc_tau = v_soc_10(2:10); % only care about 90 % to 10 % SoC
idx_tau = v_idxb4dip(2:10);

v_R_soc_tau =  v_peak_tau;
idx_R0 = idx_vpk_calc_time;
disp('1')
disp(v_soc_tau)
 % remove below 3V ripples
disp('2')
disp(v_peak_tau)

% visualize RC discharge
figure()
plot(t_cycle, v_cycle,t_cycle(idx_vpk_calc_time),v_cycle(idx_vpk_calc_time),'o',t_cycle(idx_tau),v_cycle(idx_tau),'b*',t_cycle(idx_vpk_calc_time_2),v_cycle(idx_vpk_calc_time_2),'r*'  )
title(['ECM Markers']);
ylabel('Voltage');
grid on

figure()
plot(t_cycle, i_cycle,t_cycle(idx_vpk_calc_time),i_cycle(idx_vpk_calc_time),'o',t_cycle(idx_tau),i_cycle(idx_tau),'b*',t_cycle(idx_vpk_calc_time_2),i_cycle(idx_vpk_calc_time_2),'r*'  )
title(['ECM Markers']);
ylabel('Voltage');
grid on

%% find R0
% calculate and plot R0
subst = i_cycle(idx_vpk_calc_time_2) - i_cycle(idx_vpk_calc_time);
delta_current_R0 = subst(1:9);
subst2 = v_peak_tau_2(1:9) - v_peak_tau(1:9);
delta_voltage_R0 = reshape(subst2, [9,1]);
R_0 = delta_voltage_R0./delta_current_R0;
disp(R_0)
figure
plot(R_0)
title(strcat('R0 in cycle ', num2str(cycle_num), ' cell  ', cell_id));
% title('R0 cycle 1')
ylabel('R (Ohm)');
xlabel('SoC');
grid on

%% calculate R1

subst3 = v_peak_tau(1:9);
delta_v_inf = reshape((v_soc_tau - subst3),[9,1]);
subst4 = idx_R0(1:9);
delta_i_inf = i_cycle(idx_tau) - i_cycle(subst4);

R_1 = abs(delta_v_inf./delta_i_inf) - R_0;

disp(R_1)
figure
plot(R_1)
% title('R1 cycle 1')
title(strcat('R1 in cycle ', num2str(cycle_num), ' cell  ', cell_id));
ylabel('R1 (Ohm)');
xlabel('SoC');
grid on


%% calculate C1

% find time
delta_t = seconds(t_cycle(idx_tau) - t_cycle(subst4));
disp(delta_t)

C_1 = delta_t./(4*R_1);
disp(C_1)
figure
plot(C_1)
% title('C1 cycle 1')
title(strcat('C1 in cycle ', num2str(cycle_num), ' cell  ', cell_id));
ylabel('C1 (Farad)');
xlabel('SoC');
grid on









%% ignore this part
%% simulink data calculate rmse
% voltage_cut 
% simout_cut = interp1(voltage_cut, ScopeData(:,2));


% %%
% % hc = load('hppc_v_i_exp_data/hppc_v_i_cycle1.mat', 'voltage_cut','current_cut');
% % voltage_cut = hc.voltage_cut;
% % v_cut_int = interp1(1:size(simout,1)/size(voltage_cut,1):size(simout,1), voltage_cut, 1:size(simout,1), 'linear', 'extrap')';
% % time_cut_int = interp1(1:size(simout,1)/size(time_cut,1):size(simout,1), time_cut, 1:size(simout,1), 'linear', 'extrap')';
% v_cut_int = interp1(1:size(simout,1)/size(voltage_cut,1):size(simout,1), voltage_cut, 1:size(simout,1), 'pchip')';
% time_cut_int = interp1(1:size(simout,1)/size(time_cut,1):size(simout,1), time_cut, 1:size(simout,1), 'pchip')';
% rms_error = rmse(simout,v_cut_int);
% map_error = mape(simout,v_cut_int);
% disp('The Root Mean Square Error:');
% disp(rms_error)
% disp('The Mean Absolute Percentage Error:');
% disp(map_error)
% 
% 
% % %%
% % v_cut_int = interp1(1:size(voltage_cut,1)/size(simout,1):size(voltage_cut,1), simout, 1:size(voltage_cut,1), 'linear', 'extrap')';
% % time_cut_int = interp1(1:size(voltage_cut,1)/size(simout,1):size(voltage_cut,1), 1:length(simout), 1:size(voltage_cut,1), 'linear', 'extrap')';
% % rms_error = rmse(simout,v_cut_int); % replace simout with voltage_cut
% %%
% 
% % Original size of simout
% original_length = length(simout);
% 
% % Target size - the length of voltage_cut
% target_length = length(voltage_cut);
% 
% % Generate a vector of points at which to interpolate simout
% % linspace creates a linearly spaced vector between 1 and the length of simout
% new_indices = linspace(1, original_length, target_length);
% 
% % Perform the interpolation
% simout_resized = interp1(1:original_length, simout, new_indices);
% % simout_fake = voltage_cut*0.99;
% 
% figure;
% plot(time_cut,voltage_cut, linewidth = 3.5);hold on;
% plot(time_cut, simout_resized, '--', linewidth = 3.5);
% 
% legend('Experiment','ECM Model')
% % xlabel('Time')
% ylabel('Voltage (V)')
% grid on
% title('c1p1 cycle 1')
% txt = ['RMSE: ' num2str(round(rms_error,4))];
% txt2 = ['MAPE: ' num2str(round(map_error,4))];
% text(time_cut_int(1), v_cut_int(1)-0.965, txt)
% text(time_cut_int(1), v_cut_int(1)-0.765, txt2)
% xlim([datetime('31-Jan-2022 05:34:00'), datetime('31-Jan-2022 21:52:05')]);
% ylim([2.9, 4.3]); % Adjust the values based on your expected range
% % figure;
% % plot(time_cut,voltage_cut, time_cut, simout_fake)
% %%
% % simout_diff = 
% 
% figure()
% plot(time_cut_int, v_cut_int, '-.',  linewidth = 2);hold on; % real
% plot(time_cut_int,simout,'--', linewidth = 2);
% legend('Experiment','ECM Model')
% % xlabel('Time')
% ylabel('Voltage (V)')
% grid on
% title('c1p1 cycle 1')
% txt = ['RMSE: ' num2str(round(rms_error,4))];
% txt2 = ['MAPE: ' num2str(round(map_error,4))];
% text(time_cut_int(1), v_cut_int(1)-0.965, txt)
% text(time_cut_int(1), v_cut_int(1)-0.765, txt2)
% xlim([datetime('31-Jan-2022 05:34:08'), datetime('31-Jan-2022 21:52:02')]);
% ylim([3, 4.1]); % Adjust the values based on your expected range

%% d 
%dis chunk of code is LEGENDARY
% FolderName = "C:\Users\tas87\OneDrive\Documents\Cui's Research Stuff\UG_hppc_extract_rev2\UG_hppc_extract\figures\cell2p4\Cell2p4\Cell2p4\2p4cycle4";
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = strcat('Cycle ', num2str( cycle_num), ' cell  ', cell_id, ' Figure ', num2str(get(FigHandle, 'Number')));
%   set(0, 'CurrentFigure', FigHandle);
%   saveas(FigHandle,fullfile(FolderName, [FigName '.fig'])); %Specify format for the figure
%   saveas(FigHandle,fullfile(FolderName, [FigName '.png'])); %having both figs and pngs is the play
% end

%% 

% 
% Rs_file = sprintf('cell_%s_cycle__%i_Rs_R1_C1.mat', cell_id, cycle_num);
% if isfile(Rs_file)
%     % If it exists, load the existing capacities
%     load(Rs_file, 'Rs', 'R1', 'C1');
% else
%     % If the file doesn't exist, initialize an empty array
%     Rs = [];
% end
% % Add the new capacity value for the current cycle
% % 
% Rs = mean(R_0);
% R1 = mean(R_1);
% C1 = mean(C_1);
% % 
% if isfile(Rs_file)
%     save(Rs_file, 'Rs', 'R1', 'C1', '-append');
% else
%     save(Rs_file, 'Rs', 'R1', 'C1');
% end
% 
% coeff_file = sprintf('cell_%s_cycle_%i_coeff.mat', cell_id, cycle_num);
% if isfile(coeff_file)
%     % If it exists, load the existing capacities
%     load(coeff_file, 'coeff');
% else
%     % If the file doesn't exist, initialize an empty array
%     coeff = [];
% end
% % Add the new capacity value for the current cycle
% % 
% coeff = coefficients;
% % 
% if isfile(coeff_file)
%     save(coeff_file, 'coeff', '-append');
% else
%     save(coeff_file, 'coeff');
% end
