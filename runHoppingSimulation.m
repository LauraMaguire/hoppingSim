function runHoppingSimulation()
% This is the main function for the hopping simulation.  It uses parameters
% set in initParam and executes several runs of NumericalHoppingTether.  It
% saves an object in "Output" that contains the input and output parameters
% as well as a results object.
try

  % Initialize the run.  Print current directory and start time.

  addpath('./src');
  StartTime = datestr(now);
  currentdir=pwd;
  fprintf('In dir %s\n',currentdir);
  fprintf('Start time, %s\n', StartTime);
  
  % Allocate params.
  param = struct();
  
  % Make output directories if they don't exist.
  if exist('./output','dir') == 0; mkdir('./output') ;end
  
  % Check that Params exists.  If not, make it.
  % initparams is untracked, so make it if it's not there.
  if exist('Params.mat','file') == 0
    if exist('initParam.m','file') == 0
      cpmatparams
    end
    initParam
  end
  
  % Load Params.mat and move it to ParamsFinished.mat
  load Params.mat;
  if exist('Params.mat','file')==2
    movefile('Params.mat','ParamsFinished.mat')
  end
  
  % Display everything.
  fprintf('parameters read in\n');
  disp(param);
  
  % Build a matrix of parameter vectors.
  param_mat = combvec( param.lc, param.kHop, param.koff, param.deltaT);
  [~,nparams] = size(param_mat);
  
  % For some reason, param_mat gets "sliced". Create vectors to get arround
  param_lc = param_mat(1,:);
  param_kHop = param_mat(2,:);
  param_koff= param_mat(3,:);
  param_deltaT = param_mat(4,:);
  
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
    fprintf(['Time is now ' datestr(now) '\n']);
    
    % assign temp variables
    paramTemp = param;
    paramTemp.lc = param_lc(ii);
    paramTemp.kHop = param_kHop(ii);
    paramTemp.koff = param_koff(ii);
    paramTemp.deltaT = param_deltaT(ii);
    
    % rename imported parameters for convenience
    lp = paramTemp.lp;
    lc = paramTemp.lc;
    D = paramTemp.D;
    koff = paramTemp.koff;
    Nt = paramTemp.Nt;
    runs = paramTemp.runs;
    timesteps = paramTemp.timesteps;
    deltaT = paramTemp.deltaT;
    numrec = paramTemp.numrec;

    % Calculate and save remaining parameters
    k = 3/(2*lc*lp);
    paramTemp.k = k;
    Ef = 1.88*koff^(-0.168);
    paramTemp.Ef = Ef;
    c = (Nt/1.66e6)^(1/3);
    paramTemp.c=c;
       
    % set plot_flag to zero so lots of plots don't pop up.
    plot_flag = 0;
    
    % Set name that results will be saved under
    filestring=['TrID', num2str(paramTemp.trID),...
      '_koff',num2str(paramTemp.koff,'%.2e'),...
      '_kHop',num2str(paramTemp.kHop,'%.3f'),...
      '_dT_',num2str(paramTemp.deltaT, '%.3f'),...
      '_lc',num2str(paramTemp.lc,'%.0f')];
    filename=[filestring,'.mat'];
    fprintf('%s\n',filename);
    
    % Set timescale (override deltaT input)
    % Timescale 1: diffusion between adjacent wells.
    %     t1 = 1/(D*c^2);
    %     % Timescale 2: diffusion from one side of well to another
    %     t2 = 8*Ef/(k*D);
    %     % Timescale 3: bound lifetime (1/koff)
    %     t3 = 1/koff;
    %     % Timescale 4: because the on probability keeps being larger than one
    %     t4 = 1/(koff*exp(Ef)*(20*sqrt(2*Ef/k)));
    %
    %     % deltaT must be much smaller than each timescale.  If our other
    %     % assumptions are being met properly, t1 should be much larger than t2.
    %     %  I'm not sure if I'm calculating t2 correctly, though.
    %
    %     %deltaT = min([t1,t2,t3,10*t4])/10;
    %     disp([t1,t2,t3,t4]);
    %     disp(num2str(deltaT));
    %     paramTemp.deltaT = deltaT;
    
    %     % reset number of steps based on timestep
    %     %timesteps = round(5000/deltaT);
    %     paramTemp.timesteps = timesteps;
    %     disp(num2str(timesteps));
    
    %   Initialize x-array:
    %   Dimension 1 indexes the run number.
    %   Dimension 2 indexes the timestep.
    %   Dimension 3 gives (1) the position and (2) the tether location, if
    %   bound to a tether.  If unbound, (2) is zero.
    all_x_output = zeros(runs,numrec,2);
    boundRecord = [];
    unboundRecord = [];
    hopCount = zeros(1,runs);
    hopOverageCount = zeros(1,runs);
    onOverageCount = zeros(1,runs);
    
    parfor i=1:runs
      pause(i/100); % pause for i/100 seconds
      rng('shuffle');
      fprintf(['Trial ' num2str(ii) ', Run ' num2str(i) '. Time is now ' datestr(now) '\n']);
      %fprintf('for i = %d Rand num = %f \n', i, rand() );
      % Run hopping simulation and store results.
      [ x, ~,br,hc,hoc,oo] = NumericalHoppingTether( paramTemp, plot_flag );
      all_x_output(i,:,:) = x;
      % Create list of binding events.
      [bl,ul] = listBindingEvents(br);
      boundRecord = vertcat(boundRecord,bl);
      unboundRecord = vertcat(unboundRecord,ul);
      % Add to hopping and overage counts.
      hopCount(i) = hc;
      hopOverageCount(i) = hoc;
      onOverageCount(i) = oo;
    end
    
    % Process the results.
    % Calculate koff (units of us^-1)
    if length(boundRecord) > 50
      [fit, ~] = ExpFit(boundRecord,deltaT);
      koff = -fit.b;
    else
      koff = NaN;
    end
    
    % Calculate kon (units of us^-1 uM^-1)
    if length(unboundRecord) > 50
      [fit, ~] = ExpFit(unboundRecord,deltaT);
      kon = -fit.b/Nt;
    else
      kon = NaN;
    end
    close all
    
    % Calculate pf, fraction of time spent free
    pf = sum(sum(unboundRecord))/sum(sum(vertcat(boundRecord, unboundRecord)));
    
    % Reformat the position results to calculate msd.
    %numLags = 3e4-1;
    xx=zeros(paramTemp.runs, 1,numrec);
    for i=1:paramTemp.runs
      xx(i,1,:) = reshape( all_x_output(i,:,1), [1 1 numrec] ) ;
    end
    
    % Initialize the msd-array.
    %   Dimension 1 is the MSD.
    %   Dimension 2 is the standard deviation of the MSD.
    %   Dimension 3 is the number of intervals used in the calculation.
%     msd = zeros(paramTemp.runs, numrec-1,3);
%     % Call the MSD computer.
%     for i=1:paramTemp.runs
%       [msdTemp,dtime] = computeMSD(xx(i,:,:), paramTemp.maxComputeMsdPnts, 0, 2);
%       msd(i,:,:) = msdTemp;
%     end
    
    % take average over all msd
    [msdAll,dtime] = computeMSD(xx, paramTemp.maxComputeMsdPnts, 0, 2);
    % Take the mean MSD over all runs.
%     if param.runs>1
%       meanMSD = mean(squeeze(msd(:,:,1)),1);
%       meanErr = (1/sqrt(runs))*mean(msd(:,:,2)./sqrt(msd(:,:,3)),1);
%     else
%       meanMSD = squeeze(msd(:,:,1));
%       meanErr = zeros(1,numrec);
%     end
    % Time stuff
    dtime = deltaT * paramTemp.recsteps * dtime';
    
    % Make results structure
    results = struct();
    results.meanMSD = msdAll(:,1)';
    results.meanErr = msdAll(:,2)'./sqrt(msdAll(:,3)');
    results.dtime = dtime;
%     results.msdAll = msdAll(:,1)';
%     results.msdSigAll = msdAll(:,2)';
%     results.msdNumPtnsAll = msdAll(:,3)';
    if paramTemp.storePos
      xx = reshape( xx, [paramTemp.runs, numrec] );
      results.xx = xx;
    end
    
    % calculate D
    %Deff = findHorztlAsymp(dtime(1:end/2),meanMSD(1:end/2),meanErr(1:end/2));
    results.Deff = results.meanMSD ./ ( 2*dtime );
    results.Derr = results.meanErr ./ ( 2*dtime );
    results.boundRecord = boundRecord;
    results.unboundRecord = unboundRecord;
    results.koffCalc = koff;
    results.konCalc = kon;
    results.pfCalc = pf;
    
    % Calcuate frequency of hopping per bound timestep
    results.hopFreq = sum(hopCount)/sum(boundRecord);
    % Calculate frequency of hop prob too high per hop
    if results.hopFreq ~= 0
      results.hopOverageFreq = mean(hopOverageCount/hopCount);
    else
      results.hopOverageFreq = NaN;
    end
    
    % shows number of timesteps that binding happened with onProb > 1
    results.onOverageCount = onOverageCount;
    
    % Give some warnings about the time scales
    if results.Deff(1) > 1.5
      disp('Check time scale: Deff > 1.5 at t=0.');
    elseif results.Deff(1) < 0.1
      disp('Check time scale: Deff < 0.1 at t=0.');
    end
    
    % Save the important results in a .mat file in output directory.
    fileObj = matfile(filename,'Writable',true);
    fileObj.DeffCalc = (pf*D)+(1-pf)*(koff*lc*lp*D)/(3*D+koff*lc*lp);
    fileObj.paramIn = param;
    fileObj.paramOut = paramTemp;
    fileObj.results = results;
    movefile(filename,'./output');
  end % end of loop over parameters
  
  % Finish up and display runtime.
  runTime = toc(RunTimeID);
  runHr = floor( runTime / 3600); runTime = runTime - runHr*3600;
  runMin = floor( runTime / 60);  runTime = runTime - runMin*60;
  runSec = floor(runTime);
  fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
  EndTime = datestr(now);
  fprintf('Completed run: %s\n',EndTime);
catch err
  fprintf('%s',err.getReport('extended') );
  keyboard
end
