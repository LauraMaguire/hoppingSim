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
  params = struct();
  trialmaster = struct();
  const = struct();
  modelopt = struct();
  

