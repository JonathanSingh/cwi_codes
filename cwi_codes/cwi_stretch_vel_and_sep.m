function [epsilon,variance] = cwi_stretch_vel_and_sep(sig1,sig2,dt,win_start,win_end)

% Author: Jonathan Singh - University of Edinburgh, School of GeoSciences
% Email: jonathan.singh@ed.ac.uk

% Description:
% Function for the joint estimation of a change in velocity of a medium and
% the displacement of a source or receiver's location, when both 
% perturbation occur simultaneously - using the coda wave interferometry. 
% This function is accompanies Singh et al. 2018:
% Coda Wave Interferometry for Velocity Monitoring and Acoustic Source 
% Location in Experimental Rock Physics and Rock Mechanics Applications

% Inputs:
% sig1 = unperturbed signal 
% sig2 = peturbed signal
% dt = sampling rate 
% win_start = window start (index - NOT TIME)
% wind_end = window end (index - NOT TIME)

% NOTE: both signals are required to have the same sampling rate

% Outputs: 
% epsilon - stretching factor which maximises correlation (equal to -dv/v) 
% variance = variance of travel time pertrubations. Proportional to the
%    separation and the velocity of the medium (see Singh et al., 2018)


% Arbirtary time scale (independent of sampling rate)
t0 = [1:length(sig1)];

% FIRST SEARCH
% Coarse steps from 10% to -10%
% Create array of stretching factors to be tested:
s_factors = [-0.1:0.001:0.1];

% Loop through stretching factors to interpolate sig2 by different amounts
for s =1:length(s_factors); 
    
    % sig2new is the stretched version of sig2
    sig2new = interp1(t0,sig2,t0/(1+s_factors(s)));
    
    % convert any nans to 0 
    sig2new(isnan(sig2new))=0;
    
    % maximum correlation between reference trace and stretched sig2
    cor_e(s)=max(xcorr(sig1,sig2new,'coeff'));
end

% Find the stretching factor that resulted in maximum correlation 
S1 = s_factors(cor_e==max(cor_e));  

% SECOND SEARCH
% Using S1 as a reference, fine sampling between -1% and 1% of S1 
s2_factors = [S1-0.01: 0.00005: S1 + 0.01];

% Loop through stretching factors to interpolate sig2 by different amounts
for s2 =1:length(s2_factors); 
    
    % sig2new is the stretched version of sig2
    sig2new = interp1(t0,sig2,t0/(1+s2_factors(s2)));
    
    % convert any nans to 0 
    sig2new(isnan(sig2new))=0;
    
    % maximum correlation between reference trace and stretched sig2
    cor_e2(s2)=max(xcorr(sig1,sig2new,'coeff'));
end

% find best stetching factor which maximises correlation
epsilon = s2_factors(cor_e2==max(cor_e2)); % epsilon = -dV/V


% Interpolate time windowed signals by epsilon
sig2new = interp1(t0,sig2,t0/(1+epsilon));
    
% convert any nans to 0
sig2new(isnan(sig2new))=0;

% windowed signals
sig_unp = sig1(win_start:win_end);
sig_per = sig2new(win_start:win_end);

% new time axis
t = [win_start:win_end]* dt;

% Max correlation:
R = max(xcorr(sig_unp,sig_per,'coeff'));

% Decorrelation:
decor = 1-(max(R));

% To calculate the variance of travel time pertrubations, need the dominant
% mean squared frequency - need to differentiate input signals

% Time axis for differentials 
ti = t(1:end-1)+0.5*dt;

% Differentiate and interpolate onto time axis
sigu_1= interp1(ti,(diff(sig_unp)/dt),t);
sigu_2= interp1(ti,(diff(sig_per)/dt),t);

% remove unwanted nans
sigu_1(isnan(sigu_1))=0;
sigu_2(isnan(sigu_2))=0;

% calculate mean squared frequency as a mean of the two signals
ang_freq = ((sum(sigu_1.^2)/sum(sig_unp.^2))+(sum( sigu_2.^2)/sum(sig_per.^2)))/2;

% varaince calculated from decorrelation and mean squared frequency
variance = (2*decor)/ang_freq;

end