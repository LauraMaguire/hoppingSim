% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. This should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

clear param

% Variable parameters
param.lc = [100 500]; % Nup contour length in nm
param.r0 = [0.1 0.5]; % Binding attempt rate (non-dim)
param.Ef = [0.01 0.1 1]; % Free energy in units of kT (non-dim)
param.hop_probability = [0 0.1 0.9]; % Hopping attempt rate (non-dim)

% Parameters that should stay constant
param.lp = 1; % Nup persistence length (in nm)
param.c = 0.1; % fraction of lattice sites with tether attached

% Parameters relating to the simulation logistics
param.N = 1e5; % total number of lattice points
param.timesteps = 10^5; % Total number of timesteps.
param.runs = 10; % number of times to run the simulation
param.right_probability = 0.5; % in case I want to include drift
param.numTr = 1; % number of trials
param.trID = 1; % trial ID

% Save it
save('Params', 'param');
