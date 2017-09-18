function disp = wrapdisplacement(x,xt,N)
% returns the signed displacement between x and xt, taking periodic
% boundary conditions into account
% assumes that x and xt are between 0 and N
% test = x-xt;
% if test > N/2
%     disp = test - N;
% elseif test < -N/2
%     disp = N + test;
% else
%     disp = test;
% end

d1 = abs(mod(x,N)-mod(xt,N)); % distance if no wrapping needed
d2 = N-abs(mod(xt,N)-mod(x,N)); % distance if wrapping needed

if d1 <= d2 % no wrapping needed
    s = sign(mod(x,N)-mod(xt,N)); % finding the sign is straightforward
    disp = s*d1; % use the smaller distance
elseif d2 < d1 % wrapping needed
    s = -sign(mod(x,N)-mod(xt,N)); % reverse the sign
    disp = s*d2; % use the smaller distance
end


end