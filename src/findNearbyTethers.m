function [nearbyTethers]=findNearbyTethers(x,k,Ef,L, tetherLocations, threshold)
% This function takes a position, simulation parameters, and a list of
% tether locations and returns a list of "reasonably nearby" tethers that
% might participate in binding from that location.

% Inputs: (1) x, a position; (2) k, the spring constant; (3) Ef, the free
% energy; (4) L, the length of the box (5) tetherLocations, the list of
% tether locations created earlier in NumericalHoppingTether, and (6) a
% threshold value for multiplying by the width of the well.  For instance,
% a threshold value of 10 will lead to a list of tethers within 10
% well-widths of the particle.

% Output: a subset of indexes from the tetherLocations list corresponding
% to tethers within (threshold*well-width) of the particle.

% Calculate the width of a well using SHO turning points:
w = 2*sqrt(2*Ef/k);
% Initialize nearbyTethers list:
nearbyTethers = [];
% Loop over all tethers:
for i=1:length(tetherLocations)
    % Calculate distance from tether to particle, taking into account
    % periodic boundary conditions:
    dist = wrapdistance(x,tetherLocations(i),L);
    % Add tether index to list if dist is less than (threshold*width):
    if dist <= threshold*w
        nearbyTethers = cat(1,nearbyTethers,i);
    end
end

end