function [ x, tether_locations, Ecurrent,distances] = NumericalHoppingTether( params, plot_flag )
try
% To run: [x, tether_locations] = NumericalHoppingTether;

% To do:
% Make units make sense.
% periodic boundary conditions aren't right, need to make the tethers be
% periodic not the actual positions. 
% initialize random number generator
% Check if I've done the right thing with the binding rate constant -  Not sure how to deal with the rate constant here. 


    % Simulation and Physical Parameters and setup
if nargin
    N = params.N;
    timesteps = params.timesteps;
    k = params.k;
    c = params.c;
    hop_probability = params.hop_probability; 
    binding_rate = params.binding_rate;  
    binding_energy = params.binding_energy; 
    right_probability = params.right_probability; 
else
    N = 10^3; % Number of lattice sites
    timesteps = 10^7; % Total number of timesteps.
    k = 0.1; % spring constant in terms of lattice spacing and KbT
    % so that potential energy is 0.5*k x^2 (k = k/kB T)
    c = 0.05; % fraction lattice points with tether attachment sites.
    hop_probability = 0.1; % percent of time you attempt a hop
    binding_rate = 0.1; % this is basically kon for local density = 1 probability of binding.
    binding_energy = 100; % in KT
    right_probability = 0.5; % in case I want to include drift
    plot_flag = 1;
end

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


% Some checks
% if hop_probability+binding_rate > 1
%     fprintf('hop_probability+binding_rate > 1')
% end
if M==0
    pringf('must have at least one tether')
    M = 1;
end



% Need three types of steps, to unbind, to hop, or to just diffuse. It
% seems like diffusion is a different type of step.  So there should be
% some probability

% for now just do one particle.  x(i,1) = position, x(i,2) is well number (0 if
% unbound.
x = zeros(timesteps,2);
x(1,:) = [N/2 0]; % start at the center.
Ecurrent = zeros(timesteps,1); % can remove when lifetime problem is solved
distances = zeros(timesteps,1);

for i=1:timesteps
    staterand = rand;    
    % Check for binding/unbinding state changes
    if x(i,2) ~= 0 % enter this loop is the particle is currently bound
        % LM: wrapdistance calculates how far particle is from center of
        % well, taking periodic boundary conditions into account.
        % energy_current = energy at current time
        dist = (wrapdistance(x(i,1),tether_locations(x(i,2)),N));
        distances(i) = dist;
        energy_current = 0.5*k*dist^2; % calculate the energy in this tether
        %energy_current = 0.5*k*abs(x(i,1)-tether_locations(x(i,2)))^2; % energy, nothing fancy with wrapping
        %disp(['Step ' num2str(i) ', tethered now, staterrand = ' num2str(staterrand)]);
        if staterand < hop_probability; % enter this loop if attempting a hop
            %disp(['R = ' num2str(staterrand) ', hp = ' num2str(hop_probability) ' , attempting hop']);
            % find the closest adjacent tether Could also randomly pick
            % right or left?  Would it be faster to call wrap only once and
            % have another variable?
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
            
        elseif staterand < hop_probability + binding_rate % enter this loop if attempting to unbind
            %disp(['R = ' num2str(staterrand) ', hp+br = ' num2str(hop_probability+binding_rate) ' , attempting to unbind']);
            delta_energy = binding_energy-energy_current;
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
        if staterand < binding_rate
            %[mindistance, index] = min(wrapdistance(tether_locations,x(i,1),N)); % finds closest tether
            distList = wrapdistance(tether_locations,x(i,1),N);
            mindistance = min(distList); % finds closest tether 
            minIndexList = find(distList==mindistance);
            index = datasample(minIndexList,1); % randomly choose one index

            %accept with probability that depends on the binding energy
             %check the potential energy of the current location
            energy_nearest = 0.5*k*(mindistance)^2;
            delta_energy = energy_nearest - binding_energy;
            if delta_energy < 0 % always go down in energy
                x(i+1,2) = index;
                disp(['Binding. tether location is ' num2str(x(i+1,2))]);
            elseif  rand < exp(-delta_energy) % accept moves to higher energy with probability e^-energy;
                x(i+1,2) = index;
                disp(['Binding. tether location is ' num2str(x(i+1,2))]);
            else % rejected the move, so stay where you are
                x(i+1,2) = x(i,2);
            end
        else % reject move, so stay where you are. 
            x(i+1,2) = x(i,2);
        end
        
    end
    
    
    % Done with state changes - now deal with diffusion
    if rand < right_probability % attempt move to the right.
        test_position = x(i,1)+1;
    else % or to the left.
        test_position = x(i,1)-1;
    end
    
    if x(i+1,2) == 0 % unbound, accept move
        x(i+1,1) = test_position;
    else % bound, test energy to see whether to accept move.
        %check the potential energy of the current location
        % note, the state has already changed, so the tether_location is
        % x(i+1,2), the position is changed here, so x(i,1)
        energy_current = 0.5*k*wrapdistance(x(i,1),tether_locations(x(i+1,2)), N)^2;
        %and the energy of the test position
        energy_test_position = 0.5*k*wrapdistance(test_position, tether_locations(x(i+1,2)), N)^2;
        delta_energy = energy_test_position - energy_current;
        
        if delta_energy < 0 % always go down in energy
            x(i+1,1) = test_position;
        elseif  rand < exp(-delta_energy) % accept moves to higher energy with probability e^-energy;
            x(i+1,1) = test_position;
        else % rejected the move, so stay where you are
            x(i+1,1) = x(i,1);
        end
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



