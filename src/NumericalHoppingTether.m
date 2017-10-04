function [ x, tether_locations,binding_record,hopCount,hopOverageCount, onOverage] = NumericalHoppingTether( params, plot_flag )
try
% This code runs the simulation with a continuous model, taking in only
% non-dimensional parameters.

% Import parameters.
L = params.L;
D = params.D;
deltaT = params.deltaT;
timesteps = params.timesteps;

k = params.k;
c = params.c;
koff = params.koff;
kHop = params.kHop;  
Ef = params.Ef; 

% Make a tether vector with randomly-spaced tethers from a continuous
% uniform distribution.

% Set number of tethers:
M = round(L*c);
% Make a sorted list of random tethers:
tether_locations = L*sort(rand(M,1));

% x(i,1) = position, x(i,2) is well number (0 if unbound).
x = zeros(params.numrec,2);
x(1,:) = [L/2 0]; % start at the center.
jrec = 1; % record index (first step gets recorded if t = 0)
recsteps = params.recsteps;

% initialize logical array for binding record
binding_record = false(params.timesteps,1);

% intializing arrays for some diagnostic variables
hopCount = 0;
hopOverageCount = 0;
onOverage = 0;

% get current position and binding for first step (t = 0)
nextPos = x(1,1);
nextBind = x(1,2);
%fprintf('Starting time loop\n')
for i=0:timesteps-1
    randomNumber = rand;    
    % update positions and binding
    currPos = nextPos;
    currBind = nextBind;

    % Check for binding/unbinding state changes.
    if currBind ~= 0 % enter this loop if the particle is currently bound
        binding_record(i+1)=1;
        probOff = koff*deltaT;
        if randomNumber < probOff
            nextBind = 0;  % move to unbound state
        else
            % First, pick a random nearby tether to attempt hopping to.
            nearbyIndices = findNearbyTethers(currPos,k,Ef,L,tether_locations,10);
            tetherIndex = datasample(nearbyIndices,1);
            
            % Calculate Delta G to hop between current tether and new
            % tether.
            DeltaG = (0.5*k*(wrapdistance(currPos,tether_locations(tetherIndex),L)^2 ...
                -wrapdistance(currPos,tether_locations(currBind),L)^2));
            
            % Calculate the probability of a hop.
            kHopCurrent = kHop*length(nearbyIndices)*exp(-DeltaG/2);
            probHop = kHopCurrent*deltaT;
            
            % Hop if needed, otherwise stay bound to original tether.
            if randomNumber < probOff +probHop
                nextBind = tetherIndex;
                hopCount = hopCount+1;
                if (probOff+probHop) > 1
                    hopOverageCount = hopOverageCount+1;
                end
            else
                nextBind = currBind;
            end
        end
        
    elseif currBind==0 % particle is unbound.  This loop attempts binding.
        % First, pick a random nearby tether to attempt binding to.
        nearbyIndices = findNearbyTethers(currPos,k,Ef,L,tether_locations,10);
        tetherIndex = datasample(nearbyIndices,1);
        
        % Calculate Delta G for this position and tether.
        DeltaG = -Ef+0.5*k*wrapdistance(currPos,tether_locations(tetherIndex),L)^2;
        % Calculate probability of binding to nearest tether:
        konCurrent = koff*length(nearbyIndices)*exp(-DeltaG);
        onProb = konCurrent*deltaT;
        %if onProb > 1
            %disp('onProb greater than one');
            %disp(['konCurrent = ' num2str(konCurrent)]);
            %disp(['deltaT = ' num2str(deltaT)]);
        %end
    
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
        dispFromCenter = wrapdisplacement(currPos,tether_locations(nextBind),L);
        %dispFromCenter = 0; % remove after debugging!
        % incorporate a term based on spring force.
        nextPos = currPos-D*k*dispFromCenter*deltaT + step;
    end

    % recording
    if mod(i, recsteps) == 0 
      x(jrec,1) = currPos;
      x(jrec,2) = currBind;
      jrec = jrec+1;
    end
end

% update last time
i = i + 1;
currPos = nextPos;
currBind = nextBind;
% recording
if mod(i, recsteps) == 0 
  x(jrec,1) = currPos;
  x(jrec,2) = currBind;
end

if plot_flag
    close all
    subplot(2,2,1)
    plot_time = timesteps;
    histogram(nonzeros(x(:,2)))
    title('tether locations')
    subplot(2,2,2)
    plot(x(:,1))
    title('position vs time')
    subplot(2,2,4)
    plot(nonzeros(x(:,2)))
    title('tether locations vs time')
    subplot(2,2,3)
    histfit(x(:,1), 100)
    title('histogram of locations')
end
catch err
  fprintf('%s',err.getReport('extended') );
  keyboard
end



