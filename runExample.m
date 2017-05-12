% runExample()
% Description: executeable that calls main body of model diffusion_model
% Program calls loads parameter file or calls initial parameter file if one
% doesn't exist yet, sets up parallelization, and moves outputs
%

function runExample()
try
  addpath('./src');
  StartTime = datestr(now);
  currentdir=pwd;
  fprintf('In dir %s\n',currentdir);
  fprintf('In run_bindobs, %s\n', StartTime);
  
  % Allocate params
  param = struct();
  
  %make output directories if they don't exist
  if exist('./output','dir') == 0; mkdir('./output') ;end
  
  %load params. check if it exists, if not, run it, then delete it
  %initparams on tracked, so make it if it's not there
  if exist('Params.mat','file') == 0
    if exist('initParamExample.m','file') == 0
      cpmatparams
    end
    initParamExample
  end
  load Params.mat;
  
  %display everything
  fprintf('parameters read in\n');
  disp(param);
  
  %build a parameter matrix
  param_mat = combvec( param.a, param.b );
  [~,nparams] = size(param_mat);
  
  % For some reason, param_mat gets "sliced". Create vectors to get arround
  param_a = param_mat(1,:);
  param_b = param_mat(2,:);
  
  % print some stuff
  fprintf('Starting paramloop \n')
  fprintf('nparams = %d\n', nparams)
  RunTimeID = tic;
  
  % eliminate broadcast warning for fixed parameters
  numTr = param.numTr;
  c = param.c;
  d = param.d;
  
  % loop over parameters
  for ii=1:nparams
    % scramble rng in parfor! It's rng is indepedent on ML's current state
    pause(ii);
    rng('shuffle');
    fprintf('for ii = %d Rand num = %f \n', ii, rand() );
    
    % assign temp variables
    aTemp = param_a(ii);
    bTemp = param_b(ii);
    
    % varying parameters
    pvec=[aTemp bTemp]; %parameter vector
    
    filestring=['a',num2str(aTemp,'%.2f'),...
      '_b',num2str(bTemp,'%.2f'),...
      '_c',num2str(c,'%.2f'),...
      '_d',num2str(d,'%.2f') ];
    filename=['data_',filestring,'.mat'];
    fprintf('%s\n',filename);
    % run for trial
    out = zeros(1, numTr);
    parfor jj = 1:numTr
      [out(jj)] = model(pvec, param);
    end
    % do average
    average = mean( out );
    fileObj = matfile(filename,'Writable',true);
    fileObj.average = average;
    fileObj.a = aTemp;
    fileObj.b = bTemp;
    fileObj.param = param;
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

