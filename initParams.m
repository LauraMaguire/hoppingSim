% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. this should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

%trial master
trialmaster.tind = 1; % trial indicator
trialmaster.runstrtind = 1; % run indicator
trialmaster.nt = 1; % number of trials
trialmaster.seedShift = 1; % seed shifter resolve cluster issues
trialmaster.verbose = 0; % print things

%key parameters and constants
params.bind_energy_vec = [0]; % binding_energy
params.ffrac_obst_vec= [ 0.1 0.2 ]; %filling fraction of obstacles
params.size_obst = [3]; % size of obst. prgm forces it ot be odd
params.num_tracer = 100; %filling fraction of tracers
params.tr_unbnd_diff = 1; % unbound diffusion
params.tr_bnd_diff = 0; %bound diffusion

%grid stuff
const.dim = 2; % system dimension
const.n_trials = trialmaster.nt; % number of trials
const.n_gridpoints = 40; % number of grid points, same in x and y
const.ntimesteps = round( 10 ^ ( 3.0 ) ); % number of timesteps Note 1e5 gives errors on my laptop.
const.rec_interval = round( 10 ^ ( 1.0 ) ); % time elapsed before recording
const.write_interval = round( 10 ^ ( 2.0 ) ); % time elapsed before writing to file
const.twait = 1; % time waited before recording
const.trPosRecModFlag = 0; % Flag to record tracers or not. modulated
const.trPosRecNoModFlag = 1; % Flag to record tracer state. not modulated
const.trStateRecFlag = 0; % Flag to record tracer state
const.trackOcc = 1; % Flag to track occupancy
const.obsPosRecModFlag = 0; % Flag to record tracers or not. modulated
const.obsPosRecNoModFlag = 0; % Flag to record tracer state. not modulated
const.size_tracer = 1; %tracer linear dimension, MUST BE odd integer

%msd stuff
const.calcQuad = 0; % Flag for calculating quad
const.maxpts_msd = 100; % Flag for calculating quad
const.useStart = 1; % Using t=1 to start windows instead of t=end

%model stuff
modelopt.animate=0;          %1 to show animation, 0 for no animation
modelopt.tpause=0.0;         %pause time in animation, 0.1 s is fast, 1 s is slow
modelopt.movie=0;           %1 to record movie
modelopt.obst_excl=0;       %1 if obstacles sterically exclude each other, 0 if not
modelopt.tracer_excl=0;     %MUST BE 0 so tracers don't interact (ghosts)
modelopt.obst_trace_excl=0;  %1 if obstacles and tracers mutually exclude
modelopt.edges_place=0;   %1 if place tracers on obstacle edges

% save something to const and modelopt 
modelopt.dimension=const.dim; %system dimension
const.obst_excl = modelopt.obst_excl; %system dimension
const.obst_trace_excl = modelopt.obst_trace_excl; %system dimension

% Dont place on edges if obstacles can overlap
if modelopt.tracer_excl == 0; modelopt.edges_place == 0; end;
% Fix time stuff and add some calculated things
if const.twait < 1; const.twait = 1; end;
if const.write_interval > const.ntimesteps; const.write_interval = const.ntimesteps; end;
if const.rec_interval > const.write_interval; const.rec_interval = const.write_interval; end;

const.write_interval  = round( const.write_interval / const.rec_interval ) * const.rec_interval; % fix write_interval
const.NrecChunk  = const.write_interval / const.rec_interval; % Number of recorded points / chunk
const.ntimesteps = floor( const.ntimesteps / const.write_interval ) * const.write_interval; % fix ntimesteps
const.Nreclost   = ceil( const.twait / const.rec_interval) - 1; % Number of recorded points we skip due to waiting
const.Nchunklost = ceil( const.twait / const.write_interval) - 1; % Number of chunks we skip due to waiting
const.twait      = (const.Nreclost + 1) * const.rec_interval; % fix wait time. Add 1 to start recording at twait
const.NrecTot    = const.ntimesteps / const.rec_interval - const.Nreclost; % Number of recorded points
const.NchunkTot  = const.ntimesteps / const.write_interval - const.Nchunklost; % Number of times time to file

% Store verbose
const.verbose = trialmaster.verbose;

% Fix size issues
params.size_obst( ~mod(params.size_obst,2) ) = ...
  params.size_obst(~mod(params.size_obst,2) ) - 1; 
params.size_obst = unique(params.size_obst);
params.size_obst = params.size_obst( params.size_obst <= ...
  round( max(params.ffrac_obst_vec) ^ (1/const.dim) * const.n_gridpoints ) );
if isempty(params.size_obst); params.size_obst = 1; end;

const.size_tracer( ~mod(const.size_tracer,2) ) = ...
  const.size_tracer(~mod(const.size_tracer,2) ) - 1; 
const.size_tracer = unique(const.size_tracer);
if isempty(const.size_tracer); const.size_tracer = 1; end;


% Save it
save('Params', 'const','params','trialmaster','modelopt');
