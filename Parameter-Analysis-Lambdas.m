clear all
%close all
clc

fs = 10000;
T = 1/fs;
dur = 100.915;
t = 0:T:(dur-T);
ntime = length(t);
halfsamps = floor(ntime/2);

% z_a - conditions
a = 1;
b = -1;

z = 0.5*exp(1i*2*pi)*ones(size(t));

musicians = 350;%[200 250 255 260 265 275 275 285 300 305 310 315 315 315 320 340 345 360 380 450];

%%%%%%%%%%%%%%%%%%%%%%%%%%% STIMULUS %%%%%%%%%%%%%%%%%%%%%%%%
f_s = musicians * 1.3;  % Metronome's period length in miliseconds

% CHANGE TO 0.7 or 1.3 to get matrices 

mean_indiv = nan(length(musicians), length(f_s));
std_indiv = nan(length(musicians), length(f_s));

% iterate over SPRs of musicians or non-musicians
for ispr = 1:length(musicians)
    lambda_1 = linspace(1,3,5);%22;%17.2; % 3.6%3.5;       % learning parameter
    lambda_2 = linspace(0.1,6,10);%1.45;%1.9; % 2.4%2.5;       % flexibility parameter
    figure;
    X = nan(length(lambda_1),length(lambda_2));
    X_std = nan(length(lambda_1),length(lambda_2));
    X_hebb = nan(length(lambda_1),length(lambda_2));
    
    for ilamb1 = 1:length(lambda_1)
        
        for ilamb2 = 1:length(lambda_2)
            
            % mean_asyn_freqs = zeros(1, length(f_s));
            %std_asyn_freqs = zeros(1, length(f_s));
            
            % iterate over stimulus frequencies
            for i = 1:length(f_s)
                f = zeros(size(t));    % Adaptive frequency (Hebbian)
                f(1) = 1000/(musicians(ispr));%1000/musicians(k);
                F = exp(1i*2*pi*t*(1000/f_s));  % Stinlus Metronome
                
                % forward euler integration
                for j = 2:ntime
                    z(j) = (z(j-1) + T*f(j-1)*(z(j-1)*(a + 1i*2*pi + b*(abs(z(j-1)).^2)) + F(j-1)));
                    f(j) = f(j-1) + T*(1/(2*pi))*f(j-1)*(-lambda_1(ilamb1)*real(F(j-1))*sin(angle(z(j-1))) - lambda_2(ilamb2)*(f(j-1) - f(1))/f(1));
                end
                
                
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
                    
                    X(ilamb1,ilamb2) = 1000*mean(mean_asynchrony)/fs;
                    X_std(ilamb1,ilamb2) = 1000*std(mean_asynchrony)/fs;
                    X_hebb(ilamb1,ilamb2) = mean(f(halfsamps:end));
                    
                catch
                end
            end
            
        end
        
    end
    
    lambda_1 = round(lambda_1, 2, 'significant');
    lambda_2 = round(lambda_2, 2, 'significant');

    figure;
    imagesc(X)
    xticklabels(lambda_2)
    yticklabels(lambda_1)
    c = colorbar;
    caxis([-120 55]);
    xlabel('\lambda_1', 'FontSize', 36, 'FontWeight', 'bold');
    ylabel('\lambda_2', 'FontSize', 36, 'FontWeight', 'bold');
    set(get(gca,'ylabel'),'Rotation', 90);
    c.Label.String = 'Mean Asynchrony (ms)';
    set(gca,'FontSize', 34); 
    
    
end


