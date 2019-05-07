clear all
close all
clc

% Author - Adrian S Roman: asroman@ucdavis.edu
% Simulation of human asynchronies when tapping with four different
% metronome rates deviating around an individual's SPR. Due to faster or
% slower SCR compared to the stimulus frequency, our model shows the
% anticipating and lagging bevaior oberved in humans.

fs = 10000;
T = 1/fs;
dur = 250.915;
t = 0:T:(dur-T);
ntime = length(t);
halfsamps = floor(ntime/2);

% z - conditions
a = 1;
b = -1;

z = 0.5*exp(1i*2*pi)*ones(size(t));

%%%%%%%%%%%%%%%%%%%% Group Musicians %%%%%%%%%%%%%%%%%%%%%%%
% Mean Group SPR - (404ms)
musicians = [250 260 300 310 325 340 345 350 380 400 410 430 440 450 460 465 475 480 600 650];

%%%%%%%%%%%%%%%%%% Group Non-Musicians %%%%%%%%%%%%%%%%%%%%%%
% Mean Group SPR - (306ms)
%musicians = [200 250 255 260 265 275 275 285 300 305 310 315 315 315 320 340 345 360 380 450];

%%%%%%%%%%%%%%%%%% HEBBIAN LEARNING PARAMETERS %%%%%%%%%%%%%%
lambda_1 = 2.7;%9.6316;      % learning parameter
lambda_2 = 0.76;%2.2368;       % flexibility parameter

mean_indiv = zeros(length(musicians), 4); 

% iterate over SPRs of musicians or non-musicians
for ispr = 1:length(musicians)
            %%%%%%%%%%%%%%%%%%%%%% STIMULUS %%%%%%%%%%%%%%%%%%%%%%%%%%
            % generate period lengths ms 30% faster 15% ... to a given SPR
            f_s = (0.7:0.15:1.30)*musicians(ispr);   % Metronome's period length in miliseconds
            f_s(3) = [];                             % remove central element
            mean_asyn_freqs = zeros(1, length(f_s));
            
            % iterate over stimulus frequencies
            for i = 1:length(f_s)
                f = zeros(size(t));    % Adaptive frequency (Hebbian)
                f(1) = 1000/(musicians(ispr));
                F = exp(1i*2*pi*t*(1000/f_s(i)));  % Stimulus "Metronome"
                
                % forward euler implementation
                for j = 2:ntime
                    z(j) = (z(j-1) + T*f(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + F(j-1)));
                    f(j) = f(j-1) + T*(1/(2*pi))*f(j-1)*(-lambda_1*real(F(j-1))*sin(angle(z(j-1))) - lambda_2*(f(j-1) - f(1))/f(1));
                end
                
                %%%%%%%%%%% Finding peaks %%%%%%%%%
                % Peaks for oscillator and stimilus
                [pks_z,locs_z] = findpeaks(real(z));
                [pks_F,locs_F] = findpeaks(real(F));
                locs_F = [1 locs_F];
                
               
                try
                    % which z peak is closest to the midpoint of the
                    % simulation?
                    halfsamps_locsz_diff = abs(halfsamps-locs_z);
                    [~,mid_nzpeak_index] = min(halfsamps_locsz_diff);
                    mid_nzpeak = locs_z(mid_nzpeak_index);
                    
                    % eliminate the first half of the simulation for z
                    locs_z = locs_z(mid_nzpeak_index:end);
                    
                    % which F peak is closest to mid_nzpeak?
                    mid_nzpeak_locs_F_diff = abs(locs_F - mid_nzpeak);
                    [~,mid_F_peaks_index] = min(mid_nzpeak_locs_F_diff);
                    
                    % which z peak is the penultimate one?
                    pen_nzpeak = locs_z(end-1);
                    % which F peak is closest to the penultimate z peak?
                    pen_nzpeak_locs_F_diff = abs(locs_F - pen_nzpeak);
                    [~,pen_F_peaks_index] = min(pen_nzpeak_locs_F_diff);
                    
                    % compute the mean asynchrony
                    mean_asynchrony = locs_z(1:end-1) - locs_F(mid_F_peaks_index:pen_F_peaks_index);
                    
                    mean_asyn_freqs(i) = 1000*mean(mean_asynchrony)/fs;
                    
                catch
                end
            end

            mean_indiv(ispr,:) = mean_asyn_freqs;

    
end
mean_asyncrhonies = mean(mean_indiv);
bar(mean_asyncrhonies,'FaceColor','k');
ylabel('Mean Adjusted Asynchrony (ms)')
ylim([-30,30])
xlabel('Rate Condition')
name = {'F30', 'F15', 'S15', 'S30'};
set(gca,'xticklabel',name)
set(gca,'FontSize', 34)


