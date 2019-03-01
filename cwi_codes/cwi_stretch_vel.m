function [epsilon] = cwi_stretch_vel(sig1,sig2)

% Author: Jonathan Singh - University of Edinburgh, School of GeoSciences
% Email: jonathan.singh@ed.ac.uk

% Description:
% Function for estimating the change in velocity using the coda wave 
% interferometry stretching method. This function is accompanies 
% Singh et al. 2018:
% Coda Wave Interferometry for Velocity Monitoring and Acoustic Source 
% Location in Experimental Rock Physics and Rock Mechanics Applications


% Inputs:
% sig1 = unperturbed signal 
% sig2 = peturbed signal
% NOTE: both signals are required to have the same sampling rate

% Outputs: 
% epsilon - stretching factor which maximises correlation (equal to -dv/v) 


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
epsilon = -s2_factors(cor_e2==max(cor_e2)); % epsilon = -dV/V

end