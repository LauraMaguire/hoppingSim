function [ xout ] = wrap( xin,N )
%Enforces periodic boundary condidions
% inputs  = curent position, x,
%           total number of gridpoints = N
% output = new x

xtemp = mod(xin,N);
if xtemp==0
    xout=N;
else 
    xout=xtemp;
end

end





