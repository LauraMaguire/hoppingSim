% [output]= findHorztlAsymp( t, msd, msdErr )
%   Description: Find the horizontal asymptote, maximum slope, and
%   intercept (if it can)
%   inputs: 
%   t: time 
%   msd: mean square displacement
%   msdErr: uncertainty in msd
function [output]= findHorztlAsymp(t,msd,errMsd)
% debug
debugMe = 0;
% make sure everything is right size and a column vector
if sum( size(t) ~= size(msd) ) || sum( size(t) ~= size(errMsd) )...
    || sum( size(msd) ~= size(errMsd) ) 
  error('inputs are different sizes')
end
[r, ~] = size(t);
if r == 1
  t = t'; msd = msd'; errMsd = errMsd';
end
% tuneable parameters (fit constraints)
binFac = 1; % tune number of bins
numBinBulk = 4;
bulkBinStart = 2; % what bulk bin to start looking at averages
dataSpreadThreshold = 0.1; % used for finding minimum slope in bulk bin
finalSlopeCutoff = -0.01; % only count final slopese if they are below this (b.c of noise)
% if the end slope is small, say it's steady
% make sure there isn't strange early slope behavior and late unsteady behavior
thresholdEarly= 0.005;
decidingFactor = 0.8;
thresholdSmallEndSlope = 0.02;
thresholdHugeEndSlope = 0.3;
% threshold for middle slope
thresholdMid = 0.03;
% Don't count noisy end points or start for steady state and max slope
startInd = 1;
endInd = 3;
% add paths
addpath('./src')
% Vector is noisey. Bin it to smooth it out
binNum = ceil( binFac * ( log( t(end) ) - log( t(1) ) ) );
spaceLog = round( logspace( 0, log10( length(t) ), binNum + 1 ) );
spaceLog = unique( spaceLog );
binNum = length( spaceLog ) - 1;
% Store slope and center of each bin
slopeBin = zeros( binNum, 1);
yinter = zeros( binNum, 1);
centerValyLog = zeros( binNum, 1);
centerValx= zeros( binNum, 1);
centerValxLog = zeros( binNum, 1);
aveBin = zeros( binNum, 1);
stdBin = zeros( binNum, 1);
numPnts = zeros( binNum, 1 );
binLength = zeros( binNum, 1 );
% Do a linear fit of data between points in a bin
for ii = 1:binNum
  indStart =  spaceLog(ii) ;
  indEnd = spaceLog(ii+1) ;
  yTempLog = log( msd( indStart:indEnd ) ./ t( indStart:indEnd ) );
  xTemp =  t( indStart:indEnd );
  xTempLog =  log( xTemp );
  errTemp = errMsd ( indStart:indEnd ) ./ ( msd(indStart:indEnd) );
  w = 1 ./ ( errTemp .^ 2 );
  ptsTemp = length(yTempLog);
  pfit = fit( xTempLog, yTempLog, 'poly1', 'weights',  w);
  slopeBin(ii) = pfit.p1;
  yinter(ii) =  pfit.p2;
  centerValyLog(ii) = pfit.p1 .* xTempLog( round( ptsTemp / 2 ) ) + pfit.p2;
  centerValxLog(ii) = mean(xTempLog);
  centerValx(ii) =  mean( xTemp );
  numPnts(ii) = ptsTemp;
  binLength(ii) = xTempLog(end) - xTempLog(1);
  aveBin(ii) = mean( yTempLog );
  stdBin(ii) = std( yTempLog );
end
% Do a bulk binnning to make sure you don't start analysis in crap
% skip first bin
binBulkSize = floor( binNum / numBinBulk  ); % number of bulk bins (used to average over noise)
slopeBinBulk = zeros( numBinBulk, 1 );
for ii = 1 : numBinBulk
  ind = binBulkSize * (ii - 1) + 1;
  indStart =  spaceLog(ind) ;
  indEnd = spaceLog(ind+binBulkSize) ;
  yTempLog = log( msd( indStart:indEnd ) ./ t( indStart:indEnd ) );
  xTempLog =  log( t( indStart:indEnd ) );
  errTemp = errMsd ( indStart:indEnd ) ./ ( msd(indStart:indEnd) );
  w = 1 ./ ( errTemp .^ 2 );
  pfit = fit( xTempLog, yTempLog, 'poly1', 'weights',  w);
  slopeBinBulk(ii) = pfit.p1;
end
% find mins
% Get spead of data. If it's small, ignore noisey end
minData = min ( msd./t );
maxData = max ( msd./t );
meanData = mean( msd ./ t);
spread = ( maxData - minData ) ./ meanData;
% find minimun slope of binner bins
if spread < dataSpreadThreshold
  [~, minBulkInd] = min( abs( slopeBinBulk(bulkBinStart:end-1) ) );
else
  [~, minBulkInd] = min( abs( slopeBinBulk(bulkBinStart:end) ) );
end
% Translate that to finer bin sizes
[~, minIndInBulkBin] = min( abs( slopeBin( ...
  ( minBulkInd ) * binBulkSize  + 1 : (minBulkInd+1)  * binBulkSize ) ) );
indMostZeroSlope =  minBulkInd * binBulkSize  + minIndInBulkBin;
%Store the slope of the first/last bin
slopeStart = slopeBin(1);
slopeEnd = slopeBin(end);
% Start at one and only go to indMinSlope unless it is diverging
% finalSlopes =  mean( slopeBin( end - endInd  : end ) );
finalSlopes =  mean( slopeBin( indMostZeroSlope : end ) );
if mean( finalSlopes ) < finalSlopeCutoff
  [slopeMostNeg, indSlopeMostNeg] = min( slopeBin(startInd : end) );
else
  [slopeMostNeg, indSlopeMostNeg] = min( slopeBin(startInd : indMostZeroSlope) );
end
indSlopeMostNeg =  indSlopeMostNeg + startInd - 1;
middleIndBulk = ceil( length( slopeBinBulk ) / 2 );
middleIndBins = round( length( slopeBin ) / 2 );
earlyBulkSlope = mean( slopeBinBulk( 1:middleIndBulk) );
lateBulkSlope = mean( slopeBinBulk( middleIndBulk:end ) );
dX =  sum( binLength( middleIndBulk:end) );
% if late slope if below a critical value, assume it's zero
critSlopeSmall = log(  ( 1 - thresholdSmallEndSlope / 2  ) ./ ...
  ( 1 + thresholdSmallEndSlope / 2 ) ) ./  dX;
% if late slope is too large a critical value, No steady'
critSlopeLate = log(  ( 1 - thresholdHugeEndSlope / 2  ) ./ ...
  ( 1 + thresholdHugeEndSlope / 2 ) ) ./  dX;
if abs( finalSlopes )  < abs( critSlopeSmall )
  lateBulkSlope = 0;
end
% set too big slope flag to zero
lateSlopeTooBig = 0;
if abs( finalSlopes )  > abs( critSlopeLate )
  lateSlopeTooBig = 1;
end
% make sure late slope isn't too big
if ( abs( decidingFactor * lateBulkSlope ) > abs( earlyBulkSlope ) ) && ...
    ( indMostZeroSlope < middleIndBins ) && abs( earlyBulkSlope ) > thresholdEarly
  lateSlopeTooBig = 1;
end
% Check that it got close to steady state to count: an a relative uncertainty of less than 10 %
deltaY = binLength(indMostZeroSlope) .* slopeBin(indMostZeroSlope);
relUncertaintyMin = 2 * abs( (1 - exp( deltaY ) ) ./ ( 1 + exp( -deltaY) ) );
% Check if the last few slopes are worse than the middle
% make sure uncertain in asymp from the smallest slope bin isn't too high
dX = sum( binLength(end-endInd : end) ) ;
deltaY = finalSlopes * dX;
relUncertaintyEnd = 2 * abs( (1 - exp( deltaY ) ) ./ ( 1 + exp( -deltaY) ) );
thresholdEnd = 0.1;
% determine steady state
if (relUncertaintyMin > thresholdMid) || (relUncertaintyEnd > thresholdEnd) ...
    ||  lateSlopeTooBig
  steadyState = 0;
else
  steadyState = 1;
end
% See if data equilbrated before calculating slope was useful
earlyAnom = 0;
if steadyState
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  countedBins = indMostZeroSlope:binNum ;
  indStart =  spaceLog(indMostZeroSlope) ;
  %indEnd = spaceLog(indMostZeroSlope+1);
  indEnd = spaceLog(end);
  vec = log( msd( indStart:indEnd ) ./ t( indStart:indEnd ) );
  err = errMsd( indStart:indEnd ) ./ msd( indStart:indEnd ) ;
  w = 1 ./ (err.^2);
  nPnts = length( vec );
  % calculate sums for averages
  sumVec = sum(vec .* w );
  sumW = sum( w );
  hAsymp = sumVec ./ sumW;
  sumstd2num = sum( ( vec - hAsymp ) .^ 2 .* w );
  sigh = sqrt( sumstd2num ./ ( ( nPnts - 1 ) ./ (nPnts) .*  sumW ) );
  % Calculate D
  D = exp(hAsymp);
  Dsig = sigh .* D;
  % Calculate early max slope to find the anomalous time. If it's too flat,
  % set anomalous time last center val in range
  ind2check = min( intersect( find( centerValyLog < hAsymp + sigh  ),  find( centerValyLog > hAsymp - sigh  ) ) );
  % Calculate t asym from intecept
  asympInter = ( hAsymp -  yinter(indSlopeMostNeg) ) ./ slopeMostNeg ;
  asympInterSig =  -sigh ./ slopeMostNeg ;
  tAnom = exp( asympInter );%
  tAnomSig = asympInterSig .* exp( asympInter );
  % make sure your tanom isn't larger than min. slope bin or a bin counted
  % in the average
  if ind2check <=  indSlopeMostNeg || centerValx(ind2check) < tAnom
    tAnom = centerValx( ind2check) ;
    tAnomSig = 0 ;
    earlyAnom = 1;
  end
  %%%%%%%%%%%%%%%%%%%%%%
else % not reached steady state
  D = NaN;
  Dsig = NaN;
  hAsymp = -Inf;
  sigh = NaN;
  tAnom = Inf;
  tAnomSig = NaN;
  countedBins = 0;
end
% Compile output
output.steadyState = steadyState;
output.earlyAnom = earlyAnom;
output.D = D;
output.Dsig = Dsig;
output.hAsymp = hAsymp;
output.hSig = sigh;
output.slopeEnd = slopeEnd;
output.slopeStart = slopeStart;
output.slopeMostNeg = slopeMostNeg;
output.yinterMostNeg = yinter(indSlopeMostNeg);
output.tAnom = tAnom;
output.tAnomSig = tAnomSig;
output.slopeBin = slopeBin;
output.centerBin = centerValyLog;
output.centerTimeBin = centerValx;
output.binLength = binLength;
output.upperBound = min(centerValyLog);
output.yinter = yinter;
output.countedBins= countedBins;
if debugMe
  keyboard
end
