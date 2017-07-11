function [ x, tether_locations, Ecurrent,distances] = NumericalHoppingTether( params, plot_flag )
try
 
% This code runs the simulation with a continuous model, taking in only
% non-dimensional parameters.

L = params.L;
D = params.D;
deltaT = params.deltaT;
timesteps = params.timesteps;
k = params.k;
c = params.c;
koff = params.koff;
hop_probability = params.hop_probability; 
konSite = params.konSite;  
Ef = params.Ef; 
right_probability = params.right_probability; 

% Make a tether vector with randomly-spaced tethers from a continuous
% uniform distribution.

% Set number of tethers:
M = round(L*c);
% Make a sorted list of random tethers:
tether_locations = L*sort(rand(M,1));

% x(i,1) = position, x(i,2) is well number (0 if
% unbound.
x = zeros(timesteps,2);
x(1,:) = [L/2 0]; % start at the center.

% intializing arrays for some diagnostic variables
Ecurrent = zeros(timesteps,1); % can remove when lifetime problem is solved
distances = zeros(timesteps,1);

for i=1:timesteps
    randomNumber = rand;    
    % Check for binding/unbinding state changes.
    if x(i,2) ~= 0 % enter this loop is the particle is currently bound
        probOff = koff*deltaT;
        if randomNumber < probOff;
            x(i+1,2) = 0;  % move to unbound state
        else
            x(i+1,2) = x(i,2); % otherwise, stay bound (no hopping yet)
        end
        
    elseif x(i,2)==0 % particle is unbound.  This loop attempts binding.
        % First, pick a random nearby tether to attempt binding to.
        nearbyIndices = findNearbyTethers(x(i,1),k,Ef,L,tether_locations,10);
        tetherIndex = datasample(nearbyIndices,1);
        
        % Calculate Delta G for this position and tether.
        DeltaG = -Ef+0.5*k*wrapdistance(x(i,1),tether_locations(tetherIndex),L)^2;
    
        % Sum up Boltzmann factors for nearby sites.
        denominator = sumNearbyBFs(x(i,1),k,L,nearbyIndices,tether_locations);
    
        % Calculate probability of binding to nearest tether:
        konCurrent = konSite*exp(-DeltaG)/denominator;
        onProb = konCurrent*deltaT;
    
        if randomNumber < onProb % accept binding to new tether
            x(i+1,2) = tetherIndex;
        else % rejected the move, so stay where you are
            x(i+1,2) = x(i,2);
        end
    end
    
    
    % Done with state changes - now deal with diffusion
    % gamma from Robert's Science Advances paper is just drag coefficient
    % for a sphere; gamma = kBT/D; force-dependent term works out to be
    % D*k*deltaT*deltaX (where k = spring constant/kBT)
    
    % Define Gaussian distribution from which to pick step sizes.
    % sigma corresponds to Gaussian solution to diffusion equation.
    sigma = sqrt(2*D*deltaT);
    step = normrnd(0,sigma);
    if x(i+1,2) == 0 % unbound, accept move
        x(i+1,1) = x(i,1) + step;
    else % bound, move in a force-dependent way
        % find displacement from center of well
        dispFromCenter = wrapdisplacement(x(i,1),tether_locations(x(i+1,2)),L);
        % incorporate a term based on spring force.Num
        if rand < right_probability % attempt move to the right.
            x(i+1,1) = x(i,1)-D*k*dispFromCenter*deltaT + step;
        else % or to the left.
            x(i+1,1) = x(i,1)-D*k*dispFromCenter*deltaT-step;
        end
    end
    % if particle has moved off the end of the map, put it back on the
    % other side
    x(i+1,1) = wrap(x(i+1,1),L);
    
end


if plot_flag
    close all
    subplot(2,2,1)
    plot_time = timesteps;% min(10^6, timesteps); %only plot of subset for long runs.
    histogram(nonzeros(x(1:plot_time,2)))
    title('tether locations')
    subplot(2,2,2)
    plot(x(1:plot_time,1))
    title('position vs time')
    subplot(2,2,4)
    plot(nonzeros(x(1:plot_time,2)))
    title('tether locations vs time')
    subplot(2,2,3)
    histfit(x(1:plot_time,1), 100)
    title('histogram of locations')
end
catch err
  fprintf('%s',err.getReport('extended') );
  keyboard
end



