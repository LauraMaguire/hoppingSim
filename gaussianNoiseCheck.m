% Check computeMSD with some gaussian noise
plotSummary = 1;
plotStd = 1;
% inputs
aveVal = 0;
sig = 1;
dim = 1;
numRuns = 5;  % number of runs
tSteps = 1e3; % num of time 'steps'
tWindow = 50;
maxptsMsd = 100;
useStart = 2;
quadFlag = 0;
% add paths
addpath('src')
% theory
xAveTheory = aveVal;
x2AveTheory = sig ^ 2 + aveVal .^ 2;
x4AveTheory = 3 .* sig ^4 + 6 * sig ^ 2 * aveVal ^2 + aveVal .^ 4;
msdTheory =  2 * sig ^ 2;
% calculate steps
tVec = 1:tSteps;
x = normrnd( aveVal / sqrt(dim), sig, [numRuns dim tSteps] );
x2 = x .^ 2;
% run msd on all data
[msdInfo, dtime] = computeMSD( x, maxptsMsd, quadFlag, useStart );
msdAve = msdInfo(:,1)';
% num msd for each run seperately
msdIndiv = zeros( numRuns, tSteps-1 );
msdSig = zeros( numRuns, tSteps-1 );
for ii = 1:numRuns
  msdTemp = computeMSD( x(ii,:,:), maxptsMsd, quadFlag, useStart );
  msdIndiv(ii,:) = msdTemp(:,1)';
  msdSig(ii,:) = msdTemp(:,2)';
end

% calculate simulated data average
xAveCalc = mean( x(:) );
x2AveCalc = mean( x2(:) );
x4AveCalc = mean( x2(:) .^ 2  );
msdCalc =  mean( msdAve(:) );

% Print how things went
fprintf('Theory <x> = %f Calc <x> = %f\n',...
  xAveTheory, xAveCalc)
fprintf('Theory <x>^2 = %f Calc <x>^2 = %f\n',...
  x2AveTheory, x2AveCalc)
fprintf('Theory <x>^4 = %f Calc <x>^4 = %f\n',...
  x4AveTheory, x4AveCalc)
fprintf('Theory msd = %f Calc msd = %f\n',...
  msdTheory, msdCalc)
% Average over a window
numWind = floor( tSteps / tWindow );
xWind  = zeros( numRuns, numWind );
x2Wind = zeros( numRuns, numWind);
msdWind = zeros( numRuns, numWind-1);
msdAveWind = zeros( 1, numWind-1);
msdStdWind = zeros( 1, numWind-1);
tWind = (tWindow+1) / 2 + tVec(1:tWindow:end);
dtWind = tWind(1:end-1);
for ii = 1:numWind
  ind = (1:tWindow) + (ii-1) * tWindow;
  xWind(:,ii) = mean( x(:,ind), 2 );
  x2Wind(:,ii) = mean( x2(:,ind), 2 );
  if ii < numWind
    msdWind(:,ii) = mean( msdIndiv(:,ind), 2 );
    msdAveWind(ii) = mean( msdAve(ind), 2 );
    msdStdWind(ii) = std( msdAve(ind), 0, 2 );
  end
end
% calculate variance and average amongst runs
msdWindAveRuns = mean( msdWind, 1 );
msdWindStdRuns = std( msdWind, 0, 1 );
msdWindStdRunsAveAll = mean( msdWindStdRuns );
msdWindStdRunsAveStart = mean( msdWindStdRuns( 1:round(numWind/2) ) );
msdWindStdRunsAveEnd = mean( msdWindStdRuns( round(numWind/2):end ) );
fprintf('msdStdRunsStart: %f, msdStdRunsAve: %f msdStdRunsEnd: %f\n',...
  msdWindStdRunsAveStart, msdWindStdRunsAveAll, msdWindStdRunsAveEnd);
% plot it

if dim == 1
  if plotSummary
    % if 1 dim, reshape
    x1d = squeeze(x);
    xSqr1d = squeeze(x2);
    % histogram
    numBins = 100;
    fig = figure();
    fig.WindowStyle = 'normal';
    subplot(3,3,1);
    hist( x1d(:), numBins );
    xlabel('$$ x $$ ');
    ylabel('counts');
    title('$$ x $$ data')
    subplot(3,3,2);
    hist( xSqr1d(:), numBins );
    xlabel('$$ x^2 $$ ');
    ylabel('counts');
    title('$$ x ^ 2 $$ data')
    subplot(3,3,3);
    hist( msdAve(:), numBins );
    xlabel('$$ x^2 $$ ');
    ylabel('counts');
    title('msdCompute')
    
    % no average
    % raw data
    subplot(3,3,4)
    plot( tVec, x1d )
    xlabel('$$ t $$');
    ylabel('$$ x $$');
    title('$$ x $$ data')
    % x^2 data
    subplot(3,3,5)
    plot( tVec, xSqr1d )
    xlabel('$$ t $$');
    ylabel('$$ x^2 $$');
    title('$$ x^2 $$ data')
    % msd
    subplot(3,3,6)
    plot( dtime, msdIndiv )
    hold
    plot( dtime, msdAve(:,1) )
    xlabel('$$ t $$');
    ylabel('$$ \langle x ^ 2 \rangle $$');
    title('msd Calc')
    
    % window ave
    % raw data
    subplot(3,3,7)
    plot( tWind, xWind  )
    xlabel('$$ t $$');
    ylabel('$$ x $$');
    title('window ave $$ x $$ data')
    % x^2 data
    subplot(3,3,8)
    plot( tWind, x2Wind  )
    xlabel('$$ t $$');
    ylabel('$$ x^2 $$');
    title('window ave $$ x^2 $$ data')
    % msd
    subplot(3,3,9)
    plot( dtWind, msdWind )
    hold
    errorbar( dtWind, msdAveWind, msdStdWind,'--' )
    errorbar( dtWind, msdWindAveRuns, msdWindStdRuns,'--' )
    xlabel('$$ t $$');
    ylabel('$$ \langle x ^ 2 \rangle $$');
    title('msd Calc')
  end
  if plotStd
    figure()
    plot( dtWind, msdStdWind, dtWind, msdWindStdRuns  )
    legend('all', 'runs')
  end
end