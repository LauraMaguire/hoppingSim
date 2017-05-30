function [ y, tether_locations , energy_for_average, averageBoundE] = BoundEnergyCalculator( params )
% This script calculates the average energy of a bound particle.  It can be
% used as a check against the predicted average bound energy Eb determined
% using the macroscopic parameters.  Input a paramTemp structure from a
% previous simulation run.  The outputs are (1) y, a confusing array, (2)
% tether_locations, an array indexing each tether in the simulation, (3)
% energy_for_average, an array with the bound state energy at each lattice
% site, and (4) the average bound state energy, calculated using the
% bound-state partition function.

% set parameters from input structure
N = params.N;
k = params.k;
c = params.c;   

% Make a tether vector with evenly spaced tethers
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

% For sanity, calculate Sum E exp(-E) explicitly

% begin the y-array, which will be expanded during the for-loop
% dimension 1 will run from 1 to N, indexing location within the simulation
% dimension 2 will be index of nearest tether, in tether_locations array
y(1,:) = [1,1];
% intialize the energy array - will store bound energy for each lattice
% site
energy_for_average = zeros(1,N);
% loop over all lattice sites
for site=1:N
    % find the i+1 tether site, wrapping around for periodic BCs if
    % necessary
    sitePlusOne=wrap(y(site,2)+1,M);
    % find the i-1 tether site
    siteMinusOne = wrap(y(site,2)-1,M);
    
    % find the distances to the possible nearest tethers
    
    dist_plus1 = wrapdistance(tether_locations(sitePlusOne),y(site,1), N);
    dist_minus1 = wrapdistance(tether_locations(siteMinusOne),y(site,1), N);
    dist = wrapdistance(tether_locations(wrap(y(site,2),M)),y(site,1), N);
    test_tethers = [dist_plus1, wrap(y(site,2)+1,M); ...
                    dist_minus1, wrap(y(site,2)-1,M);...
                    dist, wrap(y(site,2),M)];
    mindistance = min(test_tethers(:,1));
    testIndex = find(test_tethers(:,1)==mindistance);
    index = datasample(testIndex,1);
    energy_nearest = 0.5*k*(mindistance)^2;
    energy_for_average(site) = energy_nearest;
    if site<N
        y(site+1,:) = [site+1,test_tethers(index,2)];
    end
denominator = (1/N)*sum(exp(-energy_for_average));
averageBoundE = (1/N)*sum(energy_for_average.*exp(-energy_for_average))/denominator;
end

