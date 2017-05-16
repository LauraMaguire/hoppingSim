function [meanMSD, dtime] = runHoppingSimulation(params)
%% Run hopping simulation many times

disp('Beginning simulations.');
tic % begin measuring how long this section takes to run

%varied_params.binding_energy=[ 0 logspace(-2,2)];

% Set plot_flag to 0 so hundreds of plots don't pop up.
plot_flag = 0;

% Initialize x-array:
%   Dimension 1 indexes the run number.
%   Dimension 2 indexes the timestep.
%   Dimension 3 gives (1) the position and (2) the tether location, if
%   bound to a tether.  If unbound, (2) is zero.
clear x
x = zeros(params.runs,params.timesteps+1,2);

%size(varied_params.binding_energy,2)

% Loop over all runs.
for i=1:params.runs%:size(varied_params.binding_energy,2);
    %params.binding_energy = varied_params.binding_energy(i);
    
    % Run hopping simulation and store results.
    % tether_locs is an array giving the tether location for each tether.
    [ x(i,:,:), tether_locs ] = NumericalHoppingTether( params, plot_flag );
end
    
%varied_parameter = varied_params.binding_energy;

toc % finish time measurement

%% Process the results

% Re-format x-array so that Mike's MSD calculator can use it.
xx=zeros(1,params.runs,params.timesteps+1);
for i=1:params.runs
    xx(1,i,:) = x(i,:,1);
end

disp('Beginning MSD calculations.');
tic % Begin time measurement.

% Initialize the msd-array.
%   Dimension 1 is the MSD.
%   Dimension 2 is the standard deviation of the MSD.
%   Dimension 3 is the number of intervals used in the calculation.
msd = zeros(params.runs, params.timesteps,3);

% Call the MSD computer.
for i=1:params.runs
    [msd(i,:,:),dtime] = computeMSD(xx(1,i,:), min(1e5,params.timesteps), 0, 1);
end
toc % End time measurement.

% Take the mean MSD over all runs.
meanMSD = mean(squeeze(msd(:,:,1)),1);
close all
figure
plot(dtime,meanMSD);
title ('Mean MSD vs time');
figure
plot(dtime,meanMSD.'./dtime);
title('Diffusion coefficient vs time');
end
