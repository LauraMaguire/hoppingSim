% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. This should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

clear param

% Variable parameters
param.lc = [100 500]; % Nup contour length in nm
param.konSite = [0.1 0.5]; % kon per site (nondimensional)
param.koff = 1; % off-rate (in units of deltaT?)
param.Ef = [0.01 0.1 1]; % Free energy in units of kT (non-dim)
param.hop_probability = [0 0.1 0.9]; % Hopping attempt rate (non-dim)

% Parameters that should stay constant
param.lp = 1; % Nup persistence length (in nm)
param.c = 0.1; % number of tethers per nm

% Parameters relating to the simulation logistics
param.L = 500; % length of simulated medium in nm
param.D = 1; % diffusion coefficient
param.deltaT = 0.01; % timestep (in what units?)
param.timesteps = 10^5; % Total number of timesteps.
param.runs = 10; % number of times to run the simulation
param.right_probability = 0.5; % in case I want to include drift
param.numTr = 1; % number of trials
param.trID = 1; % trial ID

% Save it
save('Params', 'param');
