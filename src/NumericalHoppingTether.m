function [ x, tether_locations,hopCount,hopOverageCount, onOverage] = NumericalHoppingTether( params, plot_flag )
try
% currently hard-coded for koff=0!
% This code runs the simulation with a continuous model, taking in only
% non-dimensional parameters.

L = params.L;
D = params.D;
deltaT = params.deltaT;
timesteps = params.timesteps;
numrecpoints = par

k = params.k;
c = params.c;
%koff = params.koff;
koff = 0;
%hop_probability = params.hop_probability; 
kHop = params.kHop;  
Ef = params.Ef; 
%right_probability = params.right_probability;

% Make a tether vector with randomly-spaced tethers from a continuous
% uniform distribution.

% Set number of tethers:
M = round(L*c);
% Make a sorted list of random tethers:
tether_locations = L*sort(rand(M,1));

% x(i,1) = position, x(i,2) is well number (0 if
% unbound.
x = zeros(params.numrec,2);
x(1,:) = [L/2 0]; % start at the center.
jrec = 2; % record index

% intializing arrays for some diagnostic variables
%Ecurrent = zeros(timesteps,1); % can remove when lifetime problem is solved
%distances = zeros(timesteps,1);
hopCount = 0;
hopOverageCount = 0;
onOverage = 0;

% get current position and binding
currPos = x(1,1);
currBind = x(1,2);
nextPos = x(1,1);
nextBind = x(1,2);
for i=0:timesteps
    randomNumber = rand;    
    % update positions and binding
    currPos = nextPos;
    currBind = nextBind;

    % Check for binding/unbinding state changes.
    if currBind ~= 0 % enter this loop if the particle is currently bound
        probOff = koff*deltaT;
        if randomNumber < probOff
            nextBind = 0;  % move to unbound state
        else
            % First, pick a random nearby tether to attempt hopping to.
            nearbyIndices = findNearbyTethers(x(i,1),k,Ef,L,tether_locations,10);
            tetherIndex = datasample(nearbyIndices,1);
            
            % Calculate Delta G to hop between current tether and new
            % tether.
            DeltaG = (0.5*k*(wrapdistance(x(i,1),tether_locations(tetherIndex),L)^2 ...
                -wrapdistance(x(i,1),tether_locations(x(i,2)),L)^2));
            
            % Calculate the probability of a hop.
            kHopCurrent = kHop*length(nearbyIndices)*exp(-DeltaG/2);
%             disp(['M = ' num2str(length(nearbyIndices))]);
%             disp(['exp(-DG/2) = ' num2str(exp(-DeltaG/2))]);
            probHop = kHopCurrent*deltaT;
            
            % Hop if needed, otherwise stay bound to original tether.
            if randomNumber < probOff +probHop
                nextBind = tetherIndex;
%                 disp('Hopping');
                hopCount = hopCount+1;
                if (probOff+probHop) > 1
%                   disp('probOff + probHop > 1');
%                   disp(['probOff = ' num2str(probOff)]);
%                   disp(['probHop = ' num2str(probHop)]);
                    hopOverageCount = hopOverageCount+1;
                end
            else
                nextBind = currBind;
            end
            
        end
        
    elseif currBind==0 % particle is unbound.  This loop attempts binding.
        % First, pick a random nearby tether to attempt binding to.
        nearbyIndices = findNearbyTethers(x(i,1),k,Ef,L,tether_locations,10);
        tetherIndex = datasample(nearbyIndices,1);
        
        % Calculate Delta G for this position and tether.
        DeltaG = -Ef+0.5*k*wrapdistance(x(i,1),tether_locations(tetherIndex),L)^2;
    
        % Sum up Boltzmann factors for nearby sites.
        % denominator = sumNearbyBFs(x(i,1),k,L,Ef,nearbyIndices,tether_locations);
    
        % Calculate probability of binding to nearest tether:
        %konCurrent = konSite*exp(-DeltaG)/denominator;
        konCurrent = koff*length(nearbyIndices)*exp(-DeltaG);
        onProb = 1;%konCurrent*deltaT;
        if onProb > 1
            disp('onProb greater than one');
            disp(['konCurrent = ' num2str(konCurrent)]);
            disp(['deltaT = ' num2str(deltaT)]);
        end
    
        if randomNumber < onProb % accept binding to new tether
            nextBind = tetherIndex;
            if onProb > 1
                onOverage = onOverage+1;
            end
        else % rejected the move, so stay where you are
            nextBind = currBind;
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
    %distances(i) = step; %remove after debugging
    if nextBind == 0 % unbound, accept move
        nextPos = currPos + step;
    else % bound, move in a force-dependent way
        % find displacement from center of well
        dispFromCenter = wrapdisplacement(x(i,1),tether_locations(x(i+1,2)),L);
        %dispFromCenter = 0; % remove after debugging!
        % incorporate a term based on spring force.
        nextPos = currPos-D*k*dispFromCenter*deltaT + step;
    end

    % recording
    if mod(i, particle.recsteps) == 0 
      x(jrec,1) = currPos;
      x(jrec,2) = currBind;
      jrec = jrec+1;
    end
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



