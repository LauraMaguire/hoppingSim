function [boundOrNot, bindingRecord, lifetime] = LifetimeCalculator(x)
% Calculate average bound lifetime after a single run.  Input the x-vector
% that comes out of NumericalHoppingTether.  Outputs: (1) a vector
% running over all timesteps and containing a 1 if the particle is bound
% and a 0 if not, (2) a vector listing the length of each binding
% event (length of record vector will vary depending on the number of
% binding events) and (3) lifetime, the average length of a binding event.

% remove position information from x vectors, leaving only tether info
tetherInfo = squeeze(x(:,2).');
% initialize vector to contain a 1 if bound and 0 if not
boundOrNot = zeros(1, size(x,1));
% Loop over all timesteps to fill in boundOrNot
for step=1:size(x,1)
    tetherLocation = tetherInfo(step);
    if tetherLocation==0
        boundOrNot(step) = 0;
    else
        boundOrNot(step) = 1;
    end
end

% initialize counter to track number of binding events
eventCounter = 0;
% initialize counter to keep a running count of the length of each event
runningLengthCounter = 0;

% Loop over all but the first timestep to find the number and length of
% binding events.  Note that particle is always unbound on the first
% timestep.
for step=2:size(x,1)
    % add to the running count if the particle is bound
    if boundOrNot(step)==1
        runningLengthCounter = runningLengthCounter+1;
    end
    % enter this loop if ending a binding event (go from bound to unbound)
    if (boundOrNot(step)==0) && (boundOrNot(step-1)==1)
        % add an event to the counter
        eventCounter = eventCounter+1;
        % transfer the running count into the permanent record
        bindingRecord(eventCounter) = runningLengthCounter;
        % reset the running count
        runningLengthCounter = 0;
    end
end
% After looping through timesteps, set bindingRecord to zero if no
% unbinding events ever occurred.
if eventCounter==0
    bindingRecord=0;
end

% calculate the mean lifetime
lifetime = mean(bindingRecord);

end
