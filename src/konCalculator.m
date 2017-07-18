function [kon,unboundOrNot, bindingRecord] = konCalculator(x)
% Calculates macroscopic kon in units of per molarity per time.

% remove position information from x vectors, leaving only tether info
tetherInfo = squeeze(x(:,2).');
% initialize vector
unboundOrNot = zeros(1, size(x,1));
% Loop over all timesteps to fill in unboundOrNot
for step=1:size(x,1)
    tetherLocation = tetherInfo(step);
    if tetherLocation==0 % fill in a one if unbound
        unboundOrNot(step) = 1;
    else % fill in a zero if bound
        unboundOrNot(step) = 0;
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
    if unboundOrNot(step)==1
        runningLengthCounter = runningLengthCounter+1;
    end
    % enter this loop if beginning a binding event (go from unbound to
    % bound)
    if (unboundOrNot(step)==0) && (unboundOrNot(step-1)==1)
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
if length(bindingRecord)>20 % fit to an exponential and extract time constant
    h = histogram(bindingRecord,round(length(bindingRecord/20)));
    x = h.BinWidth*([1:length(h.Values)]-0.5);
    y = h.Values;
    [fit, ~] = ExpFit(x,y);
    kon = -fit.b;
else
    kon = 0;
end

end