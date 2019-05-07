clear all
close all
clc

% Author - Adrian S Roman: asroman@ucdavis.edu
% Simulation of asynchronies between partners with 
% matching and mismatching SPRs jointly tapping a melody 
% over four consecutive time periods of equal length [3]. 
% Our model replicates the behavioral data, showing greater 
% asynchrony for the “mismatching” condition. 


fs = 10000;
T = 1/fs;
dur = 52.915;
t = 0:T:(dur-T);
ntime = length(t);

% z_a - conditions
a = 1;
b = -1;

z = 0.5*exp(1i*2*pi)*ones(size(t));
y = 0.5*exp(1i*2*pi)*ones(size(t));


%%%%%%%%%%%%%%%%%%%% Group Mismatch %%%%%%%%%%%%%%%%%%%%%%%
% SPR diff > 110 ms

freqs_miss1 = [180 210 215 220 280 310 330 350 340 420];
freqs_miss2 = [310 340 330 380 410 460 462 470 500 520];
%freqs_miss2 = 320:25:550;        % Group missmatch frequencies

%%%%%%%%%%%%%%%%%%%% Group Match %%%%%%%%%%%%%%%%%%%%%%%
% SPR diff < 10 ms
freqs_match1 = [269 280 350 343 373 376 398 420 433 451];
freqs_match2 = [278 289 359 352 280 384 407 429 440 460];


%%%%%%%%%%%%%%%%%% HEBBIAN LEARNING PARAMETERS %%%%%%%%%%%%%%

lambda_1 = 2.7;          % learning parameter
lambda_2 = 0.76;         % flexibility parameter

f_1 = zeros(size(t));    % Adaptive frequency Osc1 (Hebbian)
f_2 = zeros(size(t));    % Adaptive frequency Osc2 (Hebbian)

f_s = 1000/400; 
F = exp(1i*2*pi*t*(f_s));  

mean_SPR_miss_pairs = zeros(length(freqs_miss1), 4);
mean_SPR_match_pairs = zeros(length(freqs_miss1), 4);


% iterate over SPRs of musicians or non-musicians
tr_time = 0.400 * 4; % in seconds
tr_samps = tr_time * 1000;


for i = 1:length(freqs_miss1)
    
            f_1(1) = 1000/freqs_miss1(i);
            f_2(1) = 1000/freqs_miss2(i); 
            

                % forward euler implementation
                for j = 2:ntime
                    if t(j - 1) <= tr_time  % train during four metronome cycles:
                        z(j) = (z(j-1) + T*f_2(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + F(j-1)));
                        f_2(j) = f_2(j-1) + T*(1/(2*pi))*f_2(j-1)*(-lambda_1*real(F(j-1))*sin(angle(z(j-1))) - lambda_2*(f_2(j-1) - f_2(1))/f_2(1));
                        %phase_shft1(j) = angle(z(j-1)) - angle(F(j-1));

                        y(j) = (y(j-1) + T*f_1(j-1)*(y(j-1)*(a + 1i*2*pi + b*(abs(y(j-1)).^2)) + F(j-1))); 
                        f_1(j) = f_1(j-1) + T*(1/(2*pi))*f_1(j-1)*(-lambda_1*real(F(j-1))*sin(angle(y(j-1))) - lambda_2*(f_1(j-1) - f_1(1))/f_1(1));  
                     else % interaction between two oscillators while now learning from each other:
                        z(j) = (z(j-1) + T*f_2(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + y(j-1)));
                        f_2(j) = f_2(j-1) + T*(1/(2*pi))*f_2(j-1)*(-lambda_1*real(y(j-1))*sin(angle(z(j-1))) - lambda_2*(f_2(j-1) - f_2(1))/f_2(1));

                        y(j) = (y(j-1) + T*f_1(j-1)*(y(j-1)*(a + 1i*2*pi + b*(abs(y(j-1)).^2)) + z(j-1))); 
                        f_1(j) = f_1(j-1) + T*(1/(2*pi))*f_1(j-1)*(-lambda_1*real(z(j-1))*sin(angle(y(j-1))) - lambda_2*(f_1(j-1) - f_1(1))/f_1(1));
                        %phase_shft2(j) = angle(y(j-1)) - angle(F(j-1));    
                    end
                end
                
                [pks_z,locs_z] = findpeaks(real(z));
                [pks_y,locs_y] = findpeaks(real(y));
                
                find_leader = [abs(f_s - real(f_2(tr_samps))) abs(f_s - real(f_1(tr_samps)))];
                [~, which_min] = min(find_leader);
                
                % Find which oscillator is more similar to stimulus
                if which_min == 1
                    locs_lead = locs_z;
                    locs_follow = locs_y;
                else
                    locs_lead = locs_y;
                    locs_follow = locs_z;
                end           
                    
                new_followlocs = zeros(1, length(locs_lead));
                for iloc = 1:length(locs_lead)
                  locsy_diff = abs(locs_follow-locs_lead(iloc));
                  [~,nypeak_index] = min(locsy_diff);
                  new_followlocs(iloc) = locs_follow(nypeak_index);
                end
                % find the index after training  
                tr_samps_locsz_diff = abs(tr_samps-locs_lead);
                [~,nzpeak_index] = min(tr_samps_locsz_diff);
                mid_nzpeak = locs_lead(nzpeak_index);
                % eliminate the training part of the simulation for z
                locs_lead = locs_lead(nzpeak_index:end);
                % eliminate the training part of the simulation for y
                new_followlocs = new_followlocs(nzpeak_index:end);

                % calculate number of peaks divisible by 4
                mod_four = mod(length(locs_lead(1:end)), 4);
                % eliminate extra peaks
                locs_lead = locs_lead(1:end-mod_four);
                new_followlocs = new_followlocs(1:end-mod_four);
                
                % Recover Variable names for computation
                if which_min == 1
                    locs_z = locs_lead;
                    locs_y = new_followlocs;
                else
                    locs_y = locs_lead;
                    locs_z = new_followlocs;
                end   
                  
                  
                % Mean Asynchrony   
                z_locsFourCycles = reshape(locs_z, [], 4);
                y_locsFourCycles = reshape(locs_y, [], 4);
                mean_SPR_miss_pairs(i,:) = mean(abs((t(z_locsFourCycles) - t(y_locsFourCycles))));

end
mean_abs_mismch_asyn = 1000*mean(mean_SPR_miss_pairs);  % recaling to ms
f_1 = zeros(size(t));    % Adaptive frequency Osc1 (Hebbian)
f_2 = zeros(size(t));    % Adaptive frequency Osc2 (Hebbian)

for i = 1:length(freqs_miss1)
    
            f_1(1) = 1000/freqs_match1(i);
            f_2(1) = 1000/freqs_match2(i); 
            

                % forward euler integration
                for j = 2:ntime
                    if t(j - 1) <= tr_time  % train during four metronome cycles:
                        z(j) = (z(j-1) + T*f_2(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + F(j-1)));
                        f_2(j) = f_2(j-1) + T*(1/(2*pi))*f_2(j-1)*(-lambda_1*real(F(j-1))*sin(angle(z(j-1))) - lambda_2*(f_2(j-1) - f_2(1))/f_2(1));
                        %phase_shft1(j) = angle(z(j-1)) - angle(F(j-1));

                        y(j) = (y(j-1) + T*f_1(j-1)*(y(j-1)*(a + 1i*2*pi + b*(abs(y(j-1)).^2)) + F(j-1))); 
                        f_1(j) = f_1(j-1) + T*(1/(2*pi))*f_1(j-1)*(-lambda_1*real(F(j-1))*sin(angle(y(j-1))) - lambda_2*(f_1(j-1) - f_1(1))/f_1(1));  
                     else % interaction between two oscillators while now learning from each other:
                        z(j) = (z(j-1) + T*f_2(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + y(j-1)));
                        f_2(j) = f_2(j-1) + T*(1/(2*pi))*f_2(j-1)*(-lambda_1*real(y(j-1))*sin(angle(z(j-1))) - lambda_2*(f_2(j-1) - f_2(1))/f_2(1));

                        y(j) = (y(j-1) + T*f_1(j-1)*(y(j-1)*(a + 1i*2*pi + b*(abs(y(j-1)).^2)) + z(j-1))); 
                        f_1(j) = f_1(j-1) + T*(1/(2*pi))*f_1(j-1)*(-lambda_1*real(z(j-1))*sin(angle(y(j-1))) - lambda_2*(f_1(j-1) - f_1(1))/f_1(1));
                        %phase_shft2(j) = angle(y(j-1)) - angle(F(j-1));    
                    end
                end
                
                [pks_z,locs_z] = findpeaks(real(z));
                [pks_y,locs_y] = findpeaks(real(y));
                
                % Find which oscillator is more similar to stimulus
                find_leader = [abs(f_s - real(f_2(tr_samps))) abs(f_s - real(f_1(tr_samps)))];
                [~, which_min] = min(find_leader);
             
                if which_min == 1
                    locs_lead = locs_z;
                    locs_follow = locs_y;
                else
                    locs_lead = locs_y;
                    locs_follow = locs_z;
                end           
                    
                new_followlocs = zeros(1, length(locs_lead));
                for iloc = 1:length(locs_lead)
                  locsy_diff = abs(locs_follow-locs_lead(iloc));
                  [~,nypeak_index] = min(locsy_diff);
                  new_followlocs(iloc) = locs_follow(nypeak_index);
                end
                % find the index after training  
                tr_samps_locsz_diff = abs(tr_samps-locs_lead);
                [~,nzpeak_index] = min(tr_samps_locsz_diff);
                mid_nzpeak = locs_lead(nzpeak_index);
                % eliminate the training part of the simulation for z
                locs_lead = locs_lead(nzpeak_index:end);
                % eliminate the training part of the simulation for y
                new_followlocs = new_followlocs(nzpeak_index:end);

                % calculate number of peaks divisible by 4
                mod_four = mod(length(locs_lead(1:end)), 4);
                % eliminate extra peaks
                locs_lead = locs_lead(1:end-mod_four);
                new_followlocs = new_followlocs(1:end-mod_four);
                
                % Recover Variable names for computation
                if which_min == 1
                    locs_z = locs_lead;
                    locs_y = new_followlocs;
                else
                    locs_y = locs_lead;
                    locs_z = new_followlocs;
                end   
                  
                  
                % Mean Asynchrony   
                z_locsFourCycles = reshape(locs_z, [], 4);
                y_locsFourCycles = reshape(locs_y, [], 4);
                mean_SPR_match_pairs(i,:) = mean(abs((t(z_locsFourCycles) - t(y_locsFourCycles))));
end
mean_abs_match_asyn = 1000*mean(mean_SPR_match_pairs);  % recaling to ms

data = [mean_abs_match_asyn' mean_abs_mismch_asyn'];
% 
bar_handle = bar(data, 'grouped')
set(bar_handle(1),'FaceColor','k');
set(bar_handle(2),'FaceColor','w');
set(bar_handle, 'LineWidth', 2);
names = {'Rate Match', 'Rate Mismatch'};
ylabel('Mean Absolute Duet Asynchrony (ms)')
ylim([0 35])
xlabel('Cycle of Duet Performance')
leg = legend(names);
title(leg,'SPR-Group')
set(gca,'FontSize', 34)

