function [unboundOrNot, unboundRecord] = listUnboundEvents(x)
% Calculates macroscopic kon in units of per molarity per time.
numRec = size(x,1);
% remove position information from x vectors, leaving only tether info
tetherInfo = squeeze(x(:,2).');
% initialize vector
unboundOrNot = zeros(1, numRec);
unboundRecord = zeros(1, numRec);
% Set all unbound points to zeros
unboundOrNot( tetherInfo == 0 ) = 1;
% Loop over all but the first timestep to find the number and length of
% binding events.  Note that particle is always unbound on the first
% so we can't tell how long it's been bound. Ignore it.
for step=2:numRec
  if unboundOrNot(step)==1 && unboundOrNot(step-1)==0
    for stepLater = step+1:numRec
      % add to the running count if the particle is unbound
      % enter this loop if beginning a binding event (go from unbound to
      % bound)
      if (unboundOrNot(stepLater)==0) || (stepLater == numRec)
        % transfer the running count into the permanent record
        unboundRecord(step) = stepLater - step;
        break
      end
    end % loop over later times
  end
end % loop over binding events
end

