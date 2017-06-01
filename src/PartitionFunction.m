function [Z,Zsum, E] = PartitionFunction(N,k,c)
% Calculate the bound partition function and average bound energy. Input
% the number N of lattice site, the dimensionless spring constant k, and
% the fraction c of sites with a tether attachement point.  Outputs: Z, the
% partition function (a number); and E, the dimensionless average bound
% energy (in units of kT).

% First create a copy of the tether location vector.
spacing = round(1/c);
% Slightly adjust the total number of lattice sites to make the tether
% spacing fit neatly
remainder = mod(N,spacing);
Nadj = N-remainder;
M = Nadj/spacing;
% Make the tether locations vector.
tether_locations = zeros(1,M);
for i=1:M
    tether_locations(i) = spacing*(i-1)+1;
end
% Reset the total number of lattice sites.
N = Nadj;

% Calculate Zb by looping over all possible bound states
Zsum = zeros(N,length(tether_locations));
% Loop over particle positions
tic
for i=1:N
    % Loop over all tether locations
    for j=tether_locations
        dist = wrapdistance(i,j,N);
        Zsum(i,j) = exp(-(k/2)*dist^2);
    end
end
Z = sum(sum(Zsum));
toc

% Calculate E in a similar manner: E = (1/Z)sum(Ei*exp(-Ei))
ENumerator = zeros(N,length(tether_locations));
% Loop over particle positions
tic
for i=1:N
    % Loop over all tether locations (predictable)
    for j=tether_locations
        ENumerator(i,j) = (k/2)*wrapdistance(i,j,N)^2*exp(-(k/2)*wrapdistance(i,j,N)^2);
    end
end
E = sum(sum(ENumerator))/Z;
toc

end