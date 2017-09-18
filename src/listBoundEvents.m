function [boundOrNot, bindingRecord] = listBoundEvents(x)
% Calculates macroscopic kon in units of per molarity per time.

% remove position information from x vectors, leaving only tether info
tetherInfo = squeeze(x(:,2).');
% initialize vector
boundOrNot = zeros(1, size(x,1));
% Loop over all timesteps to fill in unboundOrNot
for step=1:size(x,1)
    tetherLocation = tetherInfo(step);
    if tetherLocation==0 % fill in a zero if unbound
        boundOrNot(step) = 0;
    else % fill in a one if bound
        boundOrNot(step) = 1;
    end
end

% initialize counter to track number of binding events
eventCounter = 0;
% initialize counter to keep a running count of the length of each event
runningLengthCounter = 0;

% Loop over all but the first timestep to find the number and length of
% binding events.  Note that particle is always unbound on the first
% timestep. Ignore first timestep; negligible overall.
for step=2:size(x,1)
    % add to the running count if the particle is unbound
    if boundOrNot(step)==1
        runningLengthCounter = runningLengthCounter+1;
    end
    % enter this loop if beginning a binding event (go from unbound to
    % bound)
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
% binding events ever occurred.
if eventCounter==0
    bindingRecord=0;
end

end