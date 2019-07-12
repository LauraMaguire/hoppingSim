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
  
  %% Build a matrix of parameter vectors.
  param_mat = combvec( param.lc, param.rHop, param.koff, param.deltaT);
  [~,nparams] = size(param_mat);
  
  % For some reason, param_mat gets "sliced". Create vectors to get around.
  param_lc      = param_mat(1,:);
  param_rHop    = param_mat(2,:);
  param_koff    = param_mat(3,:);
  param_deltaT  = param_mat(4,:);
  
  % print some stuff
  fprintf('Starting paramloop \n')
  fprintf('nparams = %d\n', nparams)
  RunTimeID = tic;
  
  %% Loop over parameters
  for ii=1:nparams
    % scramble rng in parfor! It's rng is indepedent on ML's current state
    pause(ii); % pause for ii seconds
    rng('shuffle');
    fprintf('for ii = %d Rand num = %f \n', ii, rand() );
    fprintf(['Time is now ' datestr(now) '\n']);
    
    % assign temp variables
    paramTemp        = param;
    paramTemp.lc     = param_lc(ii);
    paramTemp.rHop   = param_rHop(ii);
    paramTemp.koff   = param_koff(ii);
    paramTemp.deltaT = param_deltaT(ii);
    
    % rename imported parameters for convenience
    lp            = paramTemp.lp;
    lc            = paramTemp.lc;
    D             = paramTemp.D;
    koff          = paramTemp.koff;
    Nt            = paramTemp.Nt;
    runs          = paramTemp.runs;
    timesteps     = paramTemp.timesteps;
    deltaT        = paramTemp.deltaT;
    numrec        = paramTemp.numrec;
    unbindFlag    = paramTemp.unbindFlag;
    bindFlag      = paramTemp.bindFlag;
    springForceFlag = paramTemp.springForceFlag;

    % Calculate and save remaining parameters
    k             = 3/(2*lc*lp);
    paramTemp.k   = k;
    Ef            = 1.88*koff^(-0.168);
    paramTemp.Ef  = Ef;
    c             = (Nt/1.66e6)^(1/3);
    paramTemp.c   = c;
       
    % set plot_flag to zero so lots of plots don't pop up.
    plot_flag = 0;
    
    % Set name that results will be saved under
    filestring=['TrID',  num2str(paramTemp.trID          ),...
                '_koff', num2str(paramTemp.koff,   '%.2e'),...
                '_rHop', num2str(paramTemp.rHop,   '%.2e'),...
                '_dT_',  num2str(paramTemp.deltaT, '%.3f'),...
                '_lc',   num2str(paramTemp.lc,     '%.0f')];
    filename=[filestring,'.mat'];
    fprintf('%s\n',filename);
 
    % Initialize x-array:
    % Dimension 1 indexes the run number.
    % Dimension 2 indexes the timestep.
    % Dimension 3 gives (1) the position and (2) the tether location, if
    % bound to a tether.  If unbound, (2) is zero.
    all_x_output    = zeros(runs,numrec,2);
    
    % Initialize other arrays:
    hopCount        = zeros(1,runs); % count number of hops per run
    hopOverageCount = zeros(1,runs); % count overage hops per run (hopProb > 1)
    tetherLocations = cell(1,runs);  % record tether locations for each run
    dist            = zeros(runs,numrec); % record particle's distance from
    % current tether at each timestep (NaN if not bound)
    
    % if binding AND unbinding are enabled, intialize additional arrays
    onOverageCount  = []; % count on-overage per run (onProb > 1)
    boundRecord     = []; % record binding events (will determine kon)
    unboundRecord   = []; % record unbinding events (will determine koff)
    
    %% Loop over all runs in parallel.
    parfor i=1:runs
        pause(i/100); % pause for i/100 seconds
        rng('shuffle');
        fprintf(['Trial ' num2str(ii) ', Run ' num2str(i) ...
            '. Time is now ' datestr(now) '\n']);
        %fprintf('for i = %d Rand num = %f \n', i, rand() );
      
        % Run hopping simulation and store results.
        [ x, tl,br,hc,hoc,oo] = NumericalHoppingTether(paramTemp, plot_flag);
        all_x_output(i,:,:) = x;
        % If binding AND unbinding enabled:
         if (bindFlag && unbindFlag)
            % Create list of binding and unbinding events
            [bl,ul] = listBindingEvents(br);
            boundRecord = vertcat(boundRecord,bl);
            unboundRecord = vertcat(unboundRecord,ul);
            onOverageCount = vertcat(onOverageCount,oo);
         end
        
        % Add to hopping and overage counts.
        hopCount(i) = hc;
        hopOverageCount(i) = hoc;
      
        % Record tether locations and distances.
        tetherLocations{i} = tl;
        dist(i,:) = distToTether(x,tl);
    end % End of loop over runs.
    
    %% Process the results.
    
    % Additional processing if binding AND unbinding enabled.
    if (bindFlag && unbindFlag)
        % Calculate koff (units of us^-1)
        if length(boundRecord) > 50
            [fit, ~] = ExpFit(boundRecord,deltaT); koff = -fit.b;
        else, koff = NaN;
        end
    
        % Calculate kon (units of us^-1 uM^-1)
        if length(unboundRecord) > 50
            [fit, ~] = ExpFit(unboundRecord,deltaT); kon = -fit.b/Nt;
        else, kon = NaN;
        end
        
        close all
    
        % Calculate pf, fraction of time spent free
        pf = sum(sum(unboundRecord))...
            /sum(sum(vertcat(boundRecord, unboundRecord)));
    end
    
    %% Reformat the position results to calculate msd.
    xx=zeros(paramTemp.runs, 1,numrec);
    for i=1:paramTemp.runs
      xx(i,1,:) = reshape( all_x_output(i,:,1), [1 1 numrec] ) ;
    end
    
    % Initialize the msd-array.
    %   Dimension 1 is the MSD.
    %   Dimension 2 is the standard deviation of the MSD.
    %   Dimension 3 is the number of intervals used in the calculation.

    fprintf('Running msd for all runs, i.e., particles\n')
    % Call the MSD computer.
    distFromTethers = dist;
    if (bindFlag && unbindFlag)
        [msdAll,dtime] = computeMSDFixedTimeOrigin(xx(:,:,:), ...
        paramTemp.maxComputeMsdPnts,0);
    else
        [msdAll,dtime] = computeMSDTethers(xx(:,:,:), ...
        distFromTethers, param.thresholdDistance, ...
        paramTemp.maxComputeMsdPnts);
    end
    fprintf('Finished msd for all runs, i.e., particles\n')
    dtime = deltaT * paramTemp.recsteps * dtime';
    
    %% Make results structure
    results = struct();
    
    % Populate structure
    results.dist                = dist;
    results.tetherLocations     = tetherLocations;
    results.x                   = all_x_output;
    results.meanMSD             = msdAll(:,1)';
    results.meanErr             = msdAll(:,2)'./sqrt(msdAll(:,3)');
    results.dtime               = dtime;
    
    % Calcuate frequency of hopping per bound timestep
    results.hopFreq             = sum(hopCount)/timesteps;
    
    % Calculate frequency of hop prob too high per hop
    if results.hopFreq ~= 0
        results.hopOverageFreq  = mean(hopOverageCount/hopCount);
    else
        results.hopOverageFreq  = NaN;
    end
    
    if (bindFlag && unbindFlag)
    results.boundRecord         = boundRecord;
    results.unboundRecord       = unboundRecord;
    results.koffCalc            = koff;
    results.konCalc             = kon;
    results.pfCalc              = pf;
    results.onOverageCount      = onOverageCount;
    end
    
    %% Save the important results in a .mat file in output directory.
    fileObj             = matfile(filename,'Writable',true);
    fileObj.paramIn     = param;
    fileObj.paramOut    = paramTemp;
    fileObj.results     = results;
    
    movefile(filename,'./output');
  end % end of loop over parameters
  
  % Finish up and display runtime.
  runTime   = toc(RunTimeID);
  runHr     = floor( runTime / 3600); runTime = runTime - runHr*3600;
  runMin    = floor( runTime / 60);  runTime = runTime - runMin*60;
  runSec    = floor(runTime);
  fprintf('RunTime: %.2d:%.2d:%.2d (hr:min:sec)\n', runHr, runMin,runSec);
  EndTime   = datestr(now);
  fprintf('Completed run: %s\n',EndTime);
catch err
  fprintf('%s',err.getReport('extended') );
end
