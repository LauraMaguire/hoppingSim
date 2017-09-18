function [ xout ] = wrap( xin,N )
%Enforces periodic boundary condidions
% inputs  = curent position, x,
%           total number of gridpoints = N
% output = new x
% this function works for floating-point values of xin and N

for j=1:length(xin)
    xtemp = mod(xin(j),N);
    if xtemp==0
        xout=N;
    else 
        xout(j)=xtemp;
    end
end

end





