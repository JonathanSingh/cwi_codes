function [ dv_mrt ] = mov_ref_trace( epsilon_matrix,k )

% Author: Jonathan Singh - University of Edinburgh, School of GeoSciences
% Email: jonathan.singh@ed.ac.uk

% Description:
% Function for combining estimates of velocity change from coda wave
% interferometry (epsilon from cwi_stretch_vel.m), employing a moving
% reference trace - with a user selected step size
% This function is accompanies Singh et al. 2018:
% Coda Wave Interferometry for Velocity Monitoring and Acoustic Source 
% Location in Experimental Rock Physics and Rock Mechanics Applications

% Inputs:
% epsilon_matrix = 
% k = step size for moving reference trace

% Output:
% dv_mrt = 1xN vector of cumulative velocity change

% Number of signals represented in eplison matrix
num_eps = size(epsilon_matrix,1);


% For loop through 1:N, outputing a cumulative dv/v  
for i = 1:num_eps
    
   % If trace index is less than the step size k, use the first trace as
   % the referece, i.e., epsilon_matrix(1,:)
    if i <= k 
        dv_mrt(i) = epsilon_matrix(1,i);
    else
        
        % find new reference trace index s 
        s = k*floor((i-0.5)/k);
        
        % write cumulative velocity change to dv_mrt
        dv_mrt(i) = dv_mrt(s) + epsilon_matrix(s,i);
        
    end
end
end

