% runHoppingSimulation()
% Description: executeable that calls main body of model diffusion_model
% Program calls loads parameter file or calls initial parameter file if one
% doesn't exist yet, sets up parallelization, and moves outputs
%

% Current problems:
% The script to find horizontal asymptotes doesn't work.
% Effective diffusion coefficient is wrong.
% Average bound lifetime is wrong (80 when it should be 100, 580 when it
% should be 1000).
% Actual simulated bound particle energy doesn't match the predicted value,
% and the mean distance from the tether's center is too large.  Hopefully
% fixing this problem will fix the two above problems as well.

% Things I've fixed/tried:
% Changed random number generation to fix problem with noise.  Difficult to
% check - will look out for issue in the future.
% Fixed wrapdistance error and tested.
% Fixed tether location error and tested.
% Checked that Deff = 1 for no-hopping, no-binding case.
% Fixed error leading to superdiffusive behavior - now randomly picks
% between two equidistant tethers.
% Average particle energy is very close when calculated two different ways.
% Differences are probably due to discrete vs continuous calculations.
% Changed limits of integration on stat mech calculation.
% Wrote a function to calculate Zb and Eb numerically.

% Things to try: 
% Run through conversion scheme carefully.
% Think physically about why Eb is coming up as 1/2 all the time.
% Check that eCurrent calculator is working properly.
% Redo entire analytical calculation in Matlab to eliminate algebra errors.


function runHoppingSimulation()
try
  addpath('./src');
  StartTime = datestr(now);
  currentdir=pwd;
  fprintf('In dir %s\n',currentdir);
  fprintf('Start time, %s\n', StartTime);
  
  % Allocate params
  param = struct();
  
  %make output directories if they don't exist
  if exist('./output','dir') == 0; mkdir('./output') ;end
  
  %load params. check if it exists, if not, run it, then delete it
  %initparams on tracked, so make it if it's not there
  if exist('Params.mat','file') == 0
    if exist('initParam.m','file') == 0
      cpmatparams
    end
    initParam
  end
  load Params.mat;
  
  %display everything
  fprintf('parameters read in\n');
  disp(param);
  
  %build a parameter matrix - I think these are the ones that get varied
  param_mat = combvec( param.koff, param.lc );
  [~,nparams] = size(param_mat);
  
  % For some reason, param_mat gets "sliced". Create vectors to get arround
  param_koff = param_mat(1,:);
  param_lc = param_mat(2,:);
  
  % print some stuff
  fprintf('Starting paramloop \n')
  fprintf('nparams = %d\n', nparams)
  RunTimeID = tic;
  
  % loop over parameters
  for ii=1:nparams
    % scramble rng in parfor! It's rng is indepedent on ML's current state
    pause(ii); % pause for ii seconds
    rng('shuffle');
    fprintf('for ii = %d Rand num = %f \n', ii, rand() );
    
    % assign temp variables
    
    paramTemp = param;
    paramTemp.koff = param_koff(ii);
    paramTemp.lc = param_lc(ii);
    
    paramTemp.Kd = paramTemp.koff/paramTemp.kon; % in uM
    paramTemp.Df = paramTemp.a^2/paramTemp.tau; % in nm^2/s
    paramTemp.pf = (1+paramTemp.Nt/paramTemp.Kd)^(-1); % free probability
    paramTemp.Db_theo = paramTemp.Df*paramTemp.koff*paramTemp.lc*paramTemp.lp/...
        (paramTemp.koff*paramTemp.lc*paramTemp.lp + 3*paramTemp.Df); % theoretical bound diffusion coefficient
    paramTemp.Deff_theo = paramTemp.pf*paramTemp.Df + (1-paramTemp.pf)*paramTemp.Db_theo; % theoretical effective diffusion
    
    paramTemp.c = paramTemp.a*(paramTemp.Nt*1e-6/1.66)^(1/3); % fraction of lattice sites with tether attachment point
    paramTemp.k = (3*paramTemp.a^2)/(2*paramTemp.lc*paramTemp.lp); % n.d. spring constant
    %paramTemp.nu = sqrt(pi/(2*paramTemp.k))*erf((1/(2*paramTemp.c))*sqrt(paramTemp.k/2)); % handy constant
    paramTemp.nu = sqrt(2*pi/(paramTemp.k)); % handy constant
    paramTemp.Ef = -log((2*paramTemp.c*paramTemp.Kd/paramTemp.Nt)*paramTemp.nu); % n.d. energy of a free particle (divided by thermal energy)
    %paramTemp.Eb = (1./(2.*paramTemp.c.*paramTemp.nu)).*(paramTemp.c.*paramTemp.nu-exp(-paramTemp.k./(8.*paramTemp.c.^2))./2); % n.d. avg. energy of a bound particle
    paramTemp.Eb = 0.5;
    paramTemp.Zf = paramTemp.N.*exp(-paramTemp.Ef); % free partition function
    paramTemp.Zb = paramTemp.N.*paramTemp.c.*paramTemp.nu; % bound partition function
    
    paramTemp.Z = paramTemp.Zf+paramTemp.Zb; % total partition function
    
    paramTemp.binding_energy = paramTemp.Ef;
    paramTemp.binding_rate = paramTemp.koff.*paramTemp.tau.*exp(paramTemp.Ef-paramTemp.Eb); % binding/unbinding attempt rate (should always be 1?)
    
    plot_flag = 0;
    
    filestring=['Kd',num2str(paramTemp.Kd,'%.3f'),...
      '_lc',num2str(paramTemp.lc,'%.0f'),...
      '_Nt',num2str(paramTemp.Nt,'%.0f'),...
      '_hopProb',num2str(paramTemp.hop_probability,'%.2f') ];
    filename=['data_',filestring,'.mat'];
    fprintf('%s\n',filename);

    %   Initialize x-array:
    %   Dimension 1 indexes the run number.
    %   Dimension 2 indexes the timestep.
    %   Dimension 3 gives (1) the position and (2) the tether location, if
    %   bound to a tether.  If unbound, (2) is zero.
    all_x_output = zeros(paramTemp.runs,paramTemp.timesteps+1,2);
    lifetimeList = zeros(1,paramTemp.runs);
    eCurrent = zeros(paramTemp.runs, paramTemp.timesteps);
    distList = zeros(paramTemp.runs, paramTemp.timesteps);
    % Loop over all runs.
    for i=1:paramTemp.runs
        pause(i/100); % pause for i/100 seconds
        rng('shuffle');
        %fprintf('for i = %d Rand num = %f \n', ii, rand() );
        % Run hopping simulation and store results.
        % tether_locs is an array giving the tether location for each tether.
        [ x, ~,eCurrent(i,:),distList(i,:)] = NumericalHoppingTether( paramTemp, plot_flag );
        all_x_output(i,:,:) = x;
        [~,~,lifetimeList(i)] = LifetimeCalculator(x);
    end
    % Process the results.
    % Re-format x-array so that Mike's MSD calculator can use it.
    xx=zeros(1,paramTemp.runs,paramTemp.timesteps+1);
    lr = zeros(1,paramTemp.runs);
    for i=1:paramTemp.runs
        xx(1,i,:) = all_x_output(i,:,1);
        lr(i) = sign(xx(1,i,1)-xx(1,i,end));
    end

    % Initialize the msd-array.
    %   Dimension 1 is the MSD.
    %   Dimension 2 is the standard deviation of the MSD.
    %   Dimension 3 is the number of intervals used in the calculation.
    msd = zeros(paramTemp.runs, paramTemp.timesteps,3);
    timesteps = paramTemp.timesteps;
    % Call the MSD computer.
    parfor i=1:paramTemp.runs
        [msd(i,:,:),~] = computeMSD(xx(1,i,:), min(1e5,timesteps), 0, 1);
    end
    % Take the mean MSD over all runs.
    meanMSD = mean(squeeze(msd(:,:,1)),1);
    meanErr = std(squeeze(msd(:,:,2)),1);
    
    % Save the important results in a .mat file in output directory.
    fileObj = matfile(filename,'Writable',true);
    fileObj.meanMSD = meanMSD;
    fileObj.meanErr = meanErr;
    fileObj.dtime = 1:timesteps;
    fileObj.t = 1:timesteps/2;
    fileObj.Deff = meanMSD(1:end/2)./fileObj.t;
    fileObj.Derr = meanErr(1:end/2)./fileObj.t;
    fileObj.param = param;
    fileObj.paramTemp = paramTemp;
    fileObj.lr = lr;
    fileObj.lifetimeList = lifetimeList;
    fileObj.eCurrent = eCurrent;
    fileObj.distList = distList;
    movefile(filename,'./output'); %why is this giving an error?
  end
  runTime = toc(RunTimeID);
  runHr = floor( runTime / 3600); runTime = runTime - runHr*3600;
  runMin = floor( runTime / 60);  runTime = runTime - runMin*60;
  runSec = floor(runTime);
  fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
  EndTime = datestr(now);
  fprintf('Completed run: %s\n',EndTime);
  movefile('Params.mat','ParamsFinished.mat')
catch err
  fprintf('%s',err.getReport('extended') );
  keyboard
end
