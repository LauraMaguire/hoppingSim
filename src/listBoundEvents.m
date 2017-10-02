function [boundOrNot, bindingRecord] = listBoundEvents(x)
% Calculates macroscopic kon in units of per molarity per time.
numRec = size(x,1);
% remove position information from x vectors, leaving only tether info
tetherInfo = squeeze(x(:,2).');
% initialize vector
boundOrNot = ones(1, numRec);
bindingRecord = zeros(1, numRec);
% Set all unbound points to zeros
boundOrNot( tetherInfo == 0 ) = 0;
% Loop over all but the first timestep to find the number and length of
% binding events.  Note that particle is always unbound on the first
% timestep. Ignore first timestep; negligible overall.
for step=2:numRec
  if boundOrNot(step)==1 && boundOrNot(step-1)==0
    for stepLater = 3:numRec
      % add to the running count if the particle is unbound
      % enter this loop if beginning a binding event (go from unbound to
      % bound)
      if (boundOrNot(stepLater)==0) || (stepLater == numRec)
        % transfer the running count into the permanent record
        bindingRecord(step) = stepLater - step;
        break
      end
    end % loop over later times
  end
end % loop over binding events
end
