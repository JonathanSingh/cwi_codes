function [variance] = cwi_sep(sig1,sig2,dt,win_start,win_end)

% Author: Jonathan Singh - University of Edinburgh, School of GeoSciences
% Email: jonathan.singh@ed.ac.uk

% Description:
% Function for estimating the separation between a pair of sources or a 
% pair of receivers using coda wave interferometry. 
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
% variance = variance of travel time perturbations. Proportional to the
%    separation and the velocity of the medium (see Singh et al., 2018)


% windowed signals
sig_unp = sig1(win_start:win_end);
sig_per = sig2(win_start:win_end);

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
