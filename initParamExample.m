% Initialize parameters here and save it to a Params.mat
% This is the tracked copy of the params. this should not be
% edited unless you are adding a new parameter. The parameter
% file that is called, initparams_bindobs, should be a copy of this.
% initparams_bindobs should not be tracked.

% set params
param.a = [ 1 2 3 ]; % parameter 1 that changes
param.b = [ 4 5 ];  % parameter 2 that c
param.c = 4; % parameter that does not change

% Fix things just in case
if param.c > 3 
  param.c = 2;
  fprintf('param.c too large! Setting to it to 2\n');
end
% calculated parameter
param.d = param.c / 2;

% Save it
% save( filename, variables );
save('Params', 'param');
