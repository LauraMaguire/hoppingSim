% runHoppingSimulation()
% Description: executeable that calls main body of model diffusion_model
% Program calls loads parameter file or calls initial parameter file if one
% doesn't exist yet, sets up parallelization, and moves outputs
%

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
    pause(ii);
    rng('shuffle');
    fprintf('for ii = %d Rand num = %f \n', ii, rand() );
    
    % assign temp variables
    
    paramTemp = param;
    paramTemp.koff = param_koff(ii);
    paramTemp.lc = param_lc(ii);
    
    paramTemp.Kd = paramTemp.koff/paramTemp.kon; % in uM
    paramTemp.Df = paramTemp.a^2/paramTemp.tau; % in nm^2/s
    
    paramTemp.c = paramTemp.a*(paramTemp.Nt*1e-6/1.66)^(1/3); % fraction of lattice sites with tether attachment point
    paramTemp.k = (3*paramTemp.a^2)/(2*paramTemp.lc*paramTemp.lp); % n.d. spring constant
    paramTemp.nu = sqrt(pi/(2*paramTemp.k))*erf((1/(2*paramTemp.c))*sqrt(paramTemp.k/2)); % handy constant
    paramTemp.Ef = -log((2*paramTemp.c*paramTemp.Kd/paramTemp.Nt)*paramTemp.nu); % n.d. energy of a free particle (divided by thermal energy)
    paramTemp.Eb = (1./(2.*paramTemp.c.*paramTemp.nu)).*(paramTemp.c.*paramTemp.nu-exp(-paramTemp.k./(8.*paramTemp.c.^2))./2); % n.d. avg. energy of a bound particle
    paramTemp.Zf = paramTemp.N.*exp(-paramTemp.Ef); % free partition function
    paramTemp.Zb = 2.*paramTemp.N.*paramTemp.c.*paramTemp.nu; % bound partition function
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

    % Loop over all runs.
    parfor i=1:paramTemp.runs
        % Run hopping simulation and store results.
        % tether_locs is an array giving the tether location for each tether.
        [ all_x_output(i,:,:), ~] = NumericalHoppingTether( paramTemp, plot_flag );
    end
    % Process the results.
    % Re-format x-array so that Mike's MSD calculator can use it.
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
    meanErr = mean(squeeze(msd(:,:,2)),1);
    
    % Save the important results in a .mat file in output directory.
    fileObj = matfile(filename,'Writable',true);
    fileObj.meanMSD = meanMSD;
    fileObj.meanErr = meanErr;
    fileObj. dtime = 1:timesteps;
    fileObj.param = param;
    fileObj.paramTemp = paramTemp;
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


% function [meanMSD, dtime] = runHoppingSimulation(params)
% %% Run hopping simulation many times
% 
% disp('Beginning simulations.');
% tic % begin measuring how long this section takes to run
% 
% %varied_params.binding_energy=[ 0 logspace(-2,2)];
% 
% % Set plot_flag to 0 so hundreds of plots don't pop up.
% plot_flag = 0;
% 
% % Initialize x-array:
% %   Dimension 1 indexes the run number.
% %   Dimension 2 indexes the timestep.
% %   Dimension 3 gives (1) the position and (2) the tether location, if
% %   bound to a tether.  If unbound, (2) is zero.
% clear x
% x = zeros(params.runs,params.timesteps+1,2);
% 
% %size(varied_params.binding_energy,2)
% 
% % Loop over all runs.
% for i=1:params.runs%:size(varied_params.binding_energy,2);
%     %params.binding_energy = varied_params.binding_energy(i);
%     
%     % Run hopping simulation and store results.
%     % tether_locs is an array giving the tether location for each tether.
%     [ x(i,:,:), tether_locs ] = NumericalHoppingTether( params, plot_flag );
% end
%     
% %varied_parameter = varied_params.binding_energy;
% 
% toc % finish time measurement
% 
% %% Process the results
% 
% % Re-format x-array so that Mike's MSD calculator can use it.
% xx=zeros(1,params.runs,params.timesteps+1);
% for i=1:params.runs
%     xx(1,i,:) = x(i,:,1);
% end
% 
% disp('Beginning MSD calculations.');
% tic % Begin time measurement.
% 
% % Initialize the msd-array.
% %   Dimension 1 is the MSD.
% %   Dimension 2 is the standard deviation of the MSD.
% %   Dimension 3 is the number of intervals used in the calculation.
% msd = zeros(params.runs, params.timesteps,3);
% 
% % Call the MSD computer.
% for i=1:params.runs
%     [msd(i,:,:),dtime] = computeMSD(xx(1,i,:), min(1e5,params.timesteps), 0, 1);
% end
% toc % End time measurement.
% 
% % Take the mean MSD over all runs.
% meanMSD = mean(squeeze(msd(:,:,1)),1);
% close all
% figure
% plot(dtime,meanMSD);
% title ('Mean MSD vs time');
% figure
% plot(dtime,meanMSD.'./dtime);
% title('Diffusion coefficient vs time');
% end
