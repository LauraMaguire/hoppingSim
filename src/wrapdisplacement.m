function disp = wrapdisplacement(x,xt,N)
% returns the signed displacement between x and xt, taking periodic
% boundary conditions into account
% assumes that x and xt are between 0 and N
test = x-xt;
if test > N/2
    disp = test - N;
elseif test < -N/2
    disp = N + test;
else
    disp = test;
end
end