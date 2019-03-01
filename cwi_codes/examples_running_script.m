%% examples_running_script.m
% Author: Jonathan Singh - University of Edinburgh, School of GeoSciences
% Email: jonathan.singh@ed.ac.uk

% Description:
% Example scripts for the use of coda wave interferometry coda package.
% Scipts and examples are as follows:
%
%   1) Estimating a change in velocity - using a synthetic dataset
%   generated from perturbing the fluid velocity in a digital rock based on
%   a microtomography scan of a Tivoli Traverine.
%
%   2) Estimating a source location perturbation - using a synthetic
%   dataset of a fracture plane occuring in a Tivoli Traverine
%
%   3) Estimating simultaneous source location and velocity perturbations
%   using synthetic data generated in a Berea Sandstone digital rock
%
%   4) Combining estimates of velocity change using moving reference trace
%   using laboratory data collected during deformation of a laminated
%   carbonate core. 


% This script accompanies Singh et al. 2018:
% Coda Wave Interferometry for Velocity Monitoring and Acoustic Source 
% Location in Experimental Rock Physics and Rock Mechanics Applications

clear all
close all

%% 1) Estimating a change in velocity

% Load example data: Fluid change in Tivoli Travertine digital rock
load('example_data/fluid_change_data.mat')

% Sampling interval 
dt = 5e-5;

% Create time vector
time = [1:size(fluid_change_data,1)].*dt; 

% Fluid P wave velocity used  in finite difference modelling 
fluid_vals = [1500 1501.5 1507.5 1515 1530 1575 1620 1650];

%%%%%%%%%%%%%%%%%%%%    Plot example data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% Full traces comparing unperturbed and perturbed signals
subplot(2,2,1:2)
plot(time,fluid_change_data(:,1))
hold on
plot(time,fluid_change_data(:,2))
legend('Unpertrubed','Perturbed')
title('a) Velocity Change: Full Signal')
xlabel('Time (s)')

% Plot comparison of signals for first arriving waves
subplot(2,2,3)
plot(time(13000:24000),fluid_change_data(13000:24000,1))
hold on
plot(time(13000:24000),fluid_change_data(13000:24000,2))
xlim([13000*dt 24000*dt])
title('b) First Arrival')
xlabel('Time (s)')

% Plot comparison of signals for coda waves
subplot(2,2,4)
plot(time(140000:end),fluid_change_data(140000:end,1))
hold on
plot(time(140000:end),fluid_change_data(140000:end,2))
xlim([140000*dt 150000*dt])
title('c) Coda Waves')
xlabel('Time (s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform CWI

% Unperturbed signal
sig1 = fluid_change_data(:,1);

% Loop through all signals - This may take a few minutes
for i = 1:size(fluid_change_data,2) 
    
    % Perturbed signal
    sig2 = fluid_change_data(:,i);
    
    % CWI stretching method for velocity change
    epsilon=cwi_stretch_vel(sig1,sig2);
    
    % Velocity change dV/V = - epsilon
    dV(i) = -epsilon;
        
end


%%%%%%%%%%%%%%%%%%%%%%    Plot CWI Output    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(fluid_vals,dV*100,'-x')
ylabel('\Delta V/V (%)')
xlabel('Fluid Velocity (m/s)')
title('Change in Velocity \Delta V/V from CWI')

%% 2) Estimating a change in source location

clear all

% Load example data
load('example_data/source_change_data.mat')

% Sampling interval 
dt = 5e-5;

% Create time vector
time = [1:size(source_change_data,2)].*dt; 

%%%%%%%%%%%%%%%%%%%%    Plot example data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot true source locations as a fraction of dominant wavelength
figure
scatter(true_src_locs(1,:),true_src_locs(2,:))
xlabel('X location/ \lambda')
ylabel('Y location/ \lambda')
title('True relative source locations') 

% Compare signals for sources at different locations
figure
subplot(2,2,1:2)
plot(time,source_change_data(1,:))
hold on
plot(time,source_change_data(20,:))
legend('Unperturbed','Perturbed')
title('a) Source Location Change: Full Signal')
xlabel('Time (s)')

% Plot comparison of signals for first arriving waves
subplot(2,2,3)
plot(time(7000:18000),source_change_data(1,7000:18000))
hold on

% Take 20th signal as an example
plot(time(2000:12000),source_change_data(20,2000:12000)) 
xlim([2000*dt 12000*dt])
title('b) First Arrival')
xlabel('Time (s)')

% Plot comparison of signals for coda waves
subplot(2,2,4)
plot(time(140000:end),source_change_data(1,140000:end))
hold on
plot(time(140000:end),source_change_data(20,140000:end))
xlim([140000*dt 150000*dt])
title('c) Coda Waves')
xlabel('Time (s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform CWI

% Number of signals
N = size(source_change_data,1);

% Index for start time window
win_start = 100000;

% Index for end of time window
win_end = 150000;

% Sampling interval 
dt = 5e-5;

% velocity of the medium
vel = 2200;

% Assign reference signal 
sig1 = source_change_data(1,:);    

% Loop through all sources 
for i = 1:N
    
    % 'Perturbed' signal for varying source location  
    sig2 = source_change_data(i,:);
    
   
    % Perform CWI to estimate variance of travel time perturbations
    [var] = cwi_sep(sig1,sig2,dt,win_start,win_end);
    
    % Calculate inter-source separation using relationship between
    % velocity and variance of travel time perturbations.
    % (See Singh et al. 2018) 
    sep_cwi(i) = sqrt(2*vel^2.*var)/67; % Normalised by dominant wavelength
    
    % Calculate true distance between sources
    % X component displacement
    dx = (true_src_locs(1,i)-true_src_locs(1,1))^2;
    
    % Y Component displacement
    dy = (true_src_locs(2,i)-true_src_locs(2,1))^2;
    
    % True displacement vector
    sep_true(i)= sqrt(dx + dy);
    
end


%%%%%%%%%%%%%%%%%%%%%%    Plot CWI Output    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot CWI separations as a function of true estimates
figure
scatter(sep_true,sep_cwi,'o')
hold on
% Plot true estimates line for comparison
plot(sep_true,sep_true)
xlim([0 2.2])
ylim([0 1])
xlabel('True Inter Source Distance /\lambda')
ylabel('Estimate CWI Inter Source Distance / \lambda')
legend('CWI Estimates','True Separation')
title('Estimated source separation from CWI')

%% 3) Estimating simultaneous changes in velocity and source location
clear all

% Load example data
load('example_data/source_and_vel_data.mat')
% Signals for a range of velocity changes and simultaneous source location
% perturbations. Velocity changes: 0 - 1%, source changes: 1-100/lambda


% Sampling interval 
dt = 5e-5;

% Create time vector
time = [1:size(source_and_vel_data,1)].*dt; 


%%%%%%%%%%%%%%%%%%%%    Plot example data    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare unperturbed signal with 0.2% velocity change and ~0.05*lambda
% source displacement
figure
subplot(2,2,1:2)
plot(time,source_and_vel_data(:,1,1))
hold on
plot(time,source_and_vel_data(:,3,3))
legend('Unperturbed','Perturbed')
title('a) Simultaneous Source Location and Velocity Change: Full Signal')
xlabel('Time (s)')

% Plot comparison of signals for first arriving waves
subplot(2,2,3)
plot(time(6000:16000),source_and_vel_data(6000:16000,1,1))
hold on
plot(time(6000:16000),source_and_vel_data(6000:16000,3,3))
xlim([6000*dt 16000*dt])
title('b) First Arrival')
xlabel('Time (s)')

% Plot comparison of signals for coda waves
subplot(2,2,4)
plot(time(90000:end),source_and_vel_data(90000:end,1,1))
hold on
plot(time(90000:end),source_and_vel_data(90000:end,3,3))
xlim([90000*dt 100000*dt])
title('c) Coda Waves')
xlabel('Time (s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform CWI


% Index for start time window
win_start = 500;

% Index for end of time window
win_end = 100000;

% Sampling interval 
dt = 5e-5;

% velocity of the medium
vel = 3500;

% Assign reference signal
sig1 = source_and_vel_data(:,1,1);

% Loop through all all velocity perturbations and source location changes
% NOTE: This will take a few minutes
%       Consider parallelising for large datasets
for i = 1:length(vel_change)
  
    for j = 1:length(source_pos)

        % Assign perturbed signal
        sig2 = source_and_vel_data(:,i,j);
        
        % CWI for simultaneous changes in velocity and source location
         [eps,var]=cwi_stretch_vel_and_sep(sig1,sig2,dt,win_start,win_end); 
         
        % Velocity change = -epsilon:
        vel_change_cwi(i,j) = eps;
        
        % Calculate inter-source separation using relationship between
        % velocity and variance of travel time perturbations:
        sep_cwi(i,j) = sqrt(2*vel^2.*var)/67; 
       
        % Calculate true change in velocity
        true_vel_change(i,j) =  vel_change(i)-vel_change(1);
        
        % Calculate true inter source distance
        true_source_change(i,j) = source_pos(j)-source_pos(1);
           
    end
end

%%%%%%%%%%%%%%%%%%%%%%    Plot CWI Output    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
% Compare true and estimated inter source distances
subplot(2,1,1)
scatter(true_source_change(:),sep_cwi(:))
hold on
plot(true_source_change(:),true_source_change(:))
xlabel('True Inter Source Distance /\lambda')
ylabel('Estimate CWI Inter Source Distance / \lambda')
legend('CWI Estimates','True Separation')
title('a) Estimated source separation')

% Compare true and estimated velocity changes
subplot(2,1,2)
scatter(true_vel_change(:),vel_change_cwi(:))
hold on
plot(true_vel_change(:),true_vel_change(:))
xlabel('True velocity change \Delta V/V')
ylabel('Estimate CWI velocity change \Delta V/V')
legend('CWI Estimates','True velocity change')
title('b) Estimated velocity change')

%% 4) Combine estimates with moving refereance trace (MRT)

clear all

% Combine estimates with moving reference trace function

% Function requires velocity change data to be in an NxN matrix.
% Take the output of cwi_stretch_vel in a nested for loop

% signals: Nxl matrix, where N is the number of signals and l is the signal length

% for i = 1:N
%     for j = i:N
%         
%         epsilon=cwi_stretch_vel(signals(i,:),signals(j,:));
%         epsilon_matrix(i,j)=epsilon;
%         
%     end
% end

% Load epsilon matrix: Experimental deformation of a laminated carbonate
load('example_data/mrt_example_data.mat')

% Step size for moving reference trace
k = 32;

% Use mov_ref_trace.m function
[dv_mrt] = mov_ref_trace(epsilon_matrix,k);


%%%%%%%%%%%%%%%%%%%%%%    Plot MRT Output    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

 % Compare differential stress and dv/v WITHOUT moving reference trace
subplot(2,1,1)
% Plot differential stress
yyaxis left
plot(time_elapsed,differential_stress)
ylabel('Differential Stress (MPa)')

% Plot velocity change
yyaxis right
plot(time_elapsed,epsilon_matrix(1,5:end))
ylabel('Velocity Change \Delta V/V')
xlabel('Time Elapsed (HH:MM:SS)')
legend('Differential Stress','Velocity Change')
title('a) Without moving reference trace')

 % Compare differential stress and dv/v WITH moving reference trace
subplot(2,1,2)
% Plot differential stress
yyaxis left
plot(time_elapsed,differential_stress)

ylabel('Differential Stress (MPa)')

% Plot velocity change
yyaxis right
plot(time_elapsed,dv_mrt(5:end))
legend('Differential Stress','Velocity Change')
ylim([0 0.025])
ylabel('Velocity Change \Delta V/V')
xlabel('Time Elapsed (HH:MM:SS)')
title('b) With moving reference trace')