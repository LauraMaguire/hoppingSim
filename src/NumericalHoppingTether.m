function [ x, tether_locations, Ecurrent,distances] = NumericalHoppingTether( params, plot_flag )
try
 
% This code runs the simulation with a continuous model, taking in only
% non-dimensional parameters.

if nargin
    L = params.L;
    D = params.D;
    deltaT = params.deltaT;
    timesteps = params.timesteps;
    k = params.k;
    c = params.c;
    hop_probability = params.hop_probability; 
    r0 = params.r0;  
    Ef = params.Ef; 
    right_probability = params.right_probability; 
else
    N = 10^3; % Number of lattice sites
    timesteps = 10^7; % Total number of timesteps.
    k = 0.1; % spring constant in terms of lattice spacing and KbT
    % so that potential energy is 0.5*k x^2 (k = k/kB T)
    c = 0.05; % fraction lattice points with tether attachment sites.
    hop_probability = 0.1; % percent of time you attempt a hop
    r0 = 0.1; % this is basically kon for local density = 1 probability of binding.
    Ef = 100; % in KT
    right_probability = 0.5; % in case I want to include drift
    plot_flag = 1;
end

% Mak a tether vector with randomly-spaced tethers from a continuous
% uniform distribution.

% Set number of tethers:
M = round(L*c);
% Make a sorted list of random tethers:
tether_locations = L*sort(rand(M,1));
% 

% Need three types of steps, to unbind, to hop, or to just diffuse. It
% seems like diffusion is a different type of step.  So there should be
% some probability

% for now just do one particle.  x(i,1) = position, x(i,2) is well number (0 if
% unbound.
x = zeros(timesteps,2);
x(1,:) = [L/2 0]; % start at the center.

% intializing arrays for some diagnostic variables
Ecurrent = zeros(timesteps,1); % can remove when lifetime problem is solved
distances = zeros(timesteps,1);

for i=1:timesteps
    staterand = rand;    
    % Check for binding/unbinding state changes
    if x(i,2) ~= 0 % enter this loop is the particle is currently bound
        % LM: wrapdistance calculates how far particle is from center of
        % well, taking periodic boundary conditions into account.
        % energy_current = energy at current time
        dist = (wrapdistance(x(i,1),tether_locations(x(i,2)),L));
        distances(i) = dist;
        energy_current = 0.5*k*dist^2; % calculate the energy in this tether
        %disp(['Step ' num2str(i) ', tethered now, staterrand = ' num2str(staterrand)]);
        if staterand < hop_probability; % enter this loop if attempting a hop - need to come back and fix
            %disp(['R = ' num2str(staterrand) ', hp = ' num2str(hop_probability) ' , attempting hop']);
            % find the closest adjacent tether Could also randomly pick
            % right or left?  Would it be faster to call wrap only once and
            % have another variable?
            % test_tethers needs to test only "reasonably close" tethers -
            % not sure yet how to set that up
            test_tethers = [wrapdistance(tether_locations(wrap(x(i,2)+1,M)),x(i,1), N), wrap(x(i,2)+1,M); ...
                            wrapdistance(tether_locations(wrap(x(i,2)-1,M)),x(i,1), N), wrap(x(i,2)-1,M)];
            testDist = min(test_tethers(:,1));
            minDist = min(testDist);
            testIndex = find(testDist==minDist);
            index = datasample(testIndex,1);
%             [mindistance, test_index] = min(test_tethers(:,1));
%             index = test_tethers(test_index,2);
            %energy if you were attached to the closest adjacent teather
            % LM: Calculate the particle's energy if it were attached to
            % its nearest-neighbor tether instead of its current tether
            energy_nearest = 0.5*k*(mindistance)^2; 
            delta_energy = (energy_nearest - energy_current);
            if delta_energy < 0 % always accept move if it would decrease energy
                x(i+1,2) = index; % unbind from current tether, attach to new tether
            elseif  rand < exp(-delta_energy) % accept moves to higher energy with probability e^-energy;
                x(i+1,2) = index;  % unbind from current tether, attach to new tether
            else % rejected the move, so stay where you are
                x(i+1,2) = x(i,2);
            end
            
        elseif staterand < hop_probability + r0 % enter this loop if attempting to unbind
            %disp(['R = ' num2str(staterrand) ', hp+br = ' num2str(hop_probability+binding_rate) ' , attempting to unbind']);
            delta_energy = Ef-energy_current;
            Ecurrent(i) = energy_current; % can remove when lifetime problem is solved
            %disp(['Delta energy is ' num2str(delta_energy)]);
            unbindingRand = rand;
            if delta_energy < 0 % always accept move if it would decrease energy
                x(i+1,2) = 0; % move to unbound state
                %disp('DE < 0, unbinding!')
            elseif  unbindingRand < exp(-delta_energy) % accept moves to higher energy with probability e^-energy;
                x(i+1,2) = 0;  % move to unbound state
                %disp(['R = ' num2str(unbindingRand) '. exp(-Delta E) = ' num2str(exp(-delta_energy)) ' Unbinding!'])
            else % rejected the move, so stay where you are
                x(i+1,2) = x(i,2);
                %disp(['R = ' num2str(unbindingRand) '. exp(-Delta E) = ' num2str(exp(-delta_energy)) '. Staying bound!'])
            end
        else % if you don't try to hop or unbind, remain bound to the same tether
            x(i+1,2) = x(i,2);
            %disp('staying bound to same tether');
        end
        
        
        
    elseif x(i,2)==0 % particle is unbound.  This loop attempts binding.
        if staterand < r0
            %[mindistance, index] = min(wrapdistance(tether_locations,x(i,1),N)); % finds closest tether
            distList = wrapdistance(tether_locations,x(i,1),L);
            mindistance = min(distList); % finds closest tether 
            minIndexList = find(distList==mindistance);
            index = datasample(minIndexList,1); % randomly choose one index

            %accept with probability that depends on the binding energy
             %check the potential energy of the current location
            energy_nearest = 0.5*k*(mindistance)^2;
            delta_energy = energy_nearest - Ef;
            if delta_energy < 0 % always go down in energy
                x(i+1,2) = index;
                %disp(['Binding. tether location is ' num2str(x(i+1,2))]);
            elseif  rand < exp(-delta_energy) % accept moves to higher energy with probability e^-energy;
                x(i+1,2) = index;
                %disp(['Binding. tether location is ' num2str(x(i+1,2))]);
            else % rejected the move, so stay where you are
                x(i+1,2) = x(i,2);
            end
        else % reject move, so stay where you are. 
            x(i+1,2) = x(i,2);
        end
        
    end
    
    
    % Done with state changes - now deal with diffusion
    % set arbitrary gamma for now and assume F = - kx
    g = 0.1;
    sigma = 4*D*deltaT;
    step = 1/(2*sigma)*normrnd(0,sigma); %pick step size from a gaussian distribution
    
    if x(i+1,2) == 0 % unbound, accept move
        if rand < right_probability % attempt move to the right.
            x(i+1,1) = x(i,1)+ step;
        else % or to the left.
            x(i+1,1) = x(i,1)-step;
        end
    else % bound, test energy to see whether to accept move.
        % find displacement from center of well
        dispFromCenter = wrapdisplacement(x(i,1),tether_locations(x(i+1,2)),L);
        if rand < right_probability % attempt move to the right.
            x(i+1,1) = x(i,1)-g*k*dispFromCenter*deltaT + step;
        else % or to the left.
            x(i+1,1) = x(i,1)-g*k*dispFromCenter*deltaT-step;
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



