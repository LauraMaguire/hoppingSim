function [boundList, unboundList] = listBindingEvents(bindingRecord)

t = length(bindingRecord);
boundList = zeros(t/10,1);
unboundList = zeros(t/10,1);

% initialize counter to track number of binding events
eventCounter = 0;
% initialize counter to keep a running count of the length of each event
runningLengthCounter = 0;

% Loop over all but the first timestep to find the number and length of
% binding events.  Note that particle is always unbound on the first
% timestep. Ignore first timestep; negligible overall.
for step=2:t
    % add to the running count if the particle is unbound
    if bindingRecord(step)==1
        runningLengthCounter = runningLengthCounter+1;
    end
    % enter this loop if beginning a binding event (go from unbound to
    % bound)
    if (bindingRecord(step)==0) && (bindingRecord(step-1)==1)
        % add an event to the counter
        eventCounter = eventCounter+1;
        % transfer the running count into the permanent record
        boundList(eventCounter) = runningLengthCounter;
        % reset the running count
        runningLengthCounter = 0;
    end
end

boundList = nonzeros(boundList);

% Do it again for unbound events.
eventCounter = 0;
runningLengthCounter = 0;
for step=2:t
    if bindingRecord(step)==0
        runningLengthCounter = runningLengthCounter+1;
    end
    if (bindingRecord(step)==1) && (bindingRecord(step-1)==0)
        % add an event to the counter
        eventCounter = eventCounter+1;
        % transfer the running count into the permanent record
        unboundList(eventCounter) = runningLengthCounter;
        % reset the running count
        runningLengthCounter = 0;
    end
end

unboundList = nonzeros(unboundList);

end
