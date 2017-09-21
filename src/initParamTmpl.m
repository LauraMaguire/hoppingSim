% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. This should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

clear param

% Variable nondimensional parameters
param.kHop = [0]; % non-dim hopping parameter
param.Ef = [4]; % Free energy in units of kT (non-dim)
%param.hop_probability = [0]; % Hopping attempt rate (non-dim)

% Variable dimensional parameters
param.lc = [100]; % Nup contour length (nm)
param.koff = [1e-2]; % macroscopic off-rate (us^-1)
param.deltaT = [0.005 0.01 0.02 0.04]; % timestep (units of us)

% Constant nondimensional parameters
param.lp = 1; % Nup persistence length (nm)
param.c = 0.1; % number of tethers per nm

% Constant dimensional parameters
param.lp = 1; % Nup persistence length (nm)
param.c = 0.1; % number of tethers per nm
param.L = 500; % length of simulated medium in nm
param.D = 1; % diffusion coefficient (nm^2/us)

% Parameters relating to the simulation logistics
param.timesteps = 10^5; % Total number of timesteps.
param.recsteps = 10; % number of steps before recording
param.maxComputeMsdPnts = min(1e5,param.timesteps);
param.runs = 10; % number of times to run the simulation
param.right_probability = 0.5; % in case I want to include drift
param.numTr = 1; % number of trials
param.trID = 1; % trial ID

% calculated time parameters
param.numrec = round( param.timesteps/param.recsteps) + 1 ;
% Fix timesteps to be integer number of recsteps
param.timesteps = round( param.timesteps/param.recsteps) * param.recsteps;
% Save it
save('Params', 'param');
