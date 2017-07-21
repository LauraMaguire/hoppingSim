function runHoppingSimulation()
try
    
% Things that need doing:
% Write a better wrapped displacement finder.
% Figure out why plotted Deff is smaller than Dcalc.
% See if we can get to diffusion-limited kon.


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
  param_mat = combvec( param.lc, param.konSite, param.Ef, param.hop_probability );
  [~,nparams] = size(param_mat);
  
  % For some reason, param_mat gets "sliced". Create vectors to get arround
  param_lc = param_mat(1,:);
  param_konSite = param_mat(2,:);
  param_Ef = param_mat(3,:);
  param_hop_probability= param_mat(4,:);
  
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
    paramTemp.lc = param_lc(ii);
    paramTemp.konSite = param_konSite(ii);
    paramTemp.Ef = param_Ef(ii);
    paramTemp.hop_probability = param_hop_probability(ii);
    
    % calculate remaining parameters
    lp = param.lp;
    lc = paramTemp.lc;
    D = paramTemp.D;
    k = 3/(2*lc*lp);
    Ef = paramTemp.Ef;
    koff = paramTemp.koff;
    runs = paramTemp.runs;
    timesteps = paramTemp.timesteps;
    %Nt = 1.66e6*c^3; % tether concentration in uM
    Nt = 1e3;
    c = (Nt/1.66e6)^(1/3);
    paramTemp.c=c;
    
    % make sure all parameters are recorded
    paramTemp.Nt = Nt;
    paramTemp.k = k;
    
    % set plot_flag to zero so lots of plots don't pop up.
    plot_flag = 0;
    
    filestring=['Ef',num2str(paramTemp.Ef,'%.2f'),...
      '_lc',num2str(paramTemp.lc,'%.0f'),...
      '_konSite',num2str(paramTemp.konSite,'%.2f'),...
      '_hopProb',num2str(paramTemp.hop_probability,'%.2f'),...
      '_TrID', num2str(paramTemp.trID)];
    filename=['data_',filestring,'.mat'];
    fprintf('%s\n',filename);
    
    % Set timescale (override deltaT input)
    % Timescale 1: diffusion between adjacent wells.
    t1 = 1/(D*c^2);
    % Timescale 2: diffusion from one side of well to another
    t2 = 8*Ef/(k*D);
    % Timescale 3: bound lifetime (1/koff)
    t3 = 1/koff;
    % Timescale 4: because the on probability keeps being larger than one
    t4 = 1/(koff*exp(Ef));

    % deltaT must be much smaller than each timescale.  If our other
    % assumptions are being met properly, t1 should be much larger than t2.
    %  I'm not sure if I'm calculating t2 correctly, though.
    
    deltaT = min([t1,t2,t3,t4])/100;
    paramTemp.deltaT = deltaT;

    %   Initialize x-array:
    %   Dimension 1 indexes the run number.
    %   Dimension 2 indexes the timestep.
    %   Dimension 3 gives (1) the position and (2) the tether location, if
    %   bound to a tether.  If unbound, (2) is zero.
    all_x_output = zeros(runs,timesteps+1,2);
    boundRecord = zeros(runs, timesteps+1);
    unboundList = zeros(runs, timesteps+1);
    unboundRecord = zeros(runs, timesteps+1);

    % Loop over all runs.
    parfor i=1:runs
        pause(i/100); % pause for i/100 seconds
        rng('shuffle');
        %fprintf('for i = %d Rand num = %f \n', i, rand() );
        % Run hopping simulation and store results.
        [ x, ~,~,~] = NumericalHoppingTether( paramTemp, plot_flag );
        all_x_output(i,:,:) = x;
        [~,br] = listBoundEvents(x);
        br(timesteps+1) = 0;
        boundRecord(i,:) = br;
        [unbound,ur] = listUnboundEvents(x);
        unboundList(i,:) = unbound;
        ur(timesteps+1) = 0;
        unboundRecord(i,:) = ur;
    end
    % Process the results.
    
    % Calculate koff (units of us^-1)
    boundRecord = nonzeros(boundRecord);
    if length(boundRecord) > 50
        [fit, ~] = ExpFit(boundRecord,deltaT);
        koff = -fit.b;
    else
        koff = 0;
    end
    
    % Calculate kon (units of us^-1 uM^-1)
    unboundRecord = nonzeros(unboundRecord);
    if length(unboundRecord) > 50
        [fit, ~] = ExpFit(unboundRecord,deltaT);
        kon = -fit.b/Nt;
    else
        kon = 0;
    end
    close all
    
    % Calculate pf, fraction of time spent free
    pf = sum(sum(unboundList))/(runs*timesteps);
    
    % Reformat the position results to calculate msd.
    xx=zeros(1,paramTemp.runs,paramTemp.timesteps+1);
    for i=1:paramTemp.runs
        xx(1,i,:) = all_x_output(i,:,1);
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
    dtime = param.deltaT*(1:timesteps);
    %Deff = findHorztlAsymp(dtime(1:end/2),meanMSD(1:end/2),meanErr(1:end/2));
    
    % Make results structure
    results = struct();
    results.meanMSD = meanMSD;
    results.meanErr = meanErr;
    results.dtime = dtime;
    
    t = param.deltaT*(1:timesteps/2);
    results.Deff = meanMSD(1:end/2)./(2*t);
    results.Derr = meanErr(1:end/2)./(2*t);

    results.boundRecord = nonzeros(boundRecord);
    results.unboundRecord = nonzeros(unboundRecord);
    results.koffCalc = koff;
    results.konCalc = kon;
    results.pfCalc = pf;

    % Save the important results in a .mat file in output directory.
    fileObj = matfile(filename,'Writable',true);
    fileObj.DeffCalc = (pf*D)+(1-pf)*(koff*lc*lp*D)/(3*D+koff*lc*lp);
    fileObj.paramIn = param;
    fileObj.paramOut = paramTemp;
    fileObj.results = results;
    movefile(filename,'./output');
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
