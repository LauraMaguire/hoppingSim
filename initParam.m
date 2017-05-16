% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. This should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

clear param

% set params
param.numTr = 10; % number of trials
param.trID = 1; % trial ID

% Added by Laura to make this into the hopping code
% USER-SET PARAMETERS - physical
%param.koff = logspace(-5,-2); % off-rate (in us^-1)
param.koff = 0.01;
param.kon = 1e-3; % on-rate (in us^-1 uM^-1)
param.Nt = 1e3; % total Nup concentration (in uM)
param.lc = [100 500]; % Nup contour length (in nm)
param.lp = 1; % Nup persistence length (in nm)
param.a = 1; % lattice spacing (in nm)
param.tau = 1; % timestep (in us)
%p.hop_probability = [0 0.1 0.3 0.9]; % percent of time you attempt a hop
param.hop_probability = 0;
param.right_probability = 0.5; % in case I want to include drift

% USER_SET PARAMETERS - simulation
param.N = 1e5; % total number of lattice points
param.timesteps = 10^5; % Total number of timesteps.
param.runs = 2; % number of times to run the simulation

% CALCULATED PARAMETERS (dimensional)
% p.Kd = p.koff./p.kon; % in uM
% p.Df = p.a^2./p.tau; % in nm^2/s
% 
% CALCULATED PARAMETERS (dimensionless)
%c = (Nt*1e-6)*a^3/1.66; % fraction of lattice sites with tether attachment
%point - wrong way to calculate!?
% p.c = p.a.*(p.Nt.*1e-6./1.66).^(1/3); % fraction of lattice sites with tether attachment point
% p.k = (3*p.a.^2)./(2.*p.lc.*p.lp); % n.d. spring constant
% p.nu = sqrt(pi./(2.*p.k)).*erf((1./(2.*p.c)).*sqrt(p.k./2)); % handy constant
%p.Ef = -log((2.*p.c.*p.Kd./p.Nt).*p.nu); % n.d. energy of a free particle (divided by thermal energy)
%p.Eb = (1./(2.*p.c.*p.nu)).*(p.c.*p.nu-exp(-p.k./(8.*p.c.^2))./2); % n.d. avg. energy of a bound particle
%p.Zf = p.N.*exp(-p.Ef); % free partition function
%p.Zb = 2.*p.N.*p.c.*p.nu; % bound partition function
%p.Z = p.Zf+p.Zb; % total partition function
%p.binding_energy = p.Ef;
%p.binding_rate = p.koff.*p.tau.*exp(p.Ef-p.Eb); % binding/unbinding attempt rate (should always be 1?)

% Save it
save('Params', 'param');
