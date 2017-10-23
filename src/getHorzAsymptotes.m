% [asymInfo, binInfo] = getHorzAsymptotes( x, y, errY, threshold, numBins,  ...
%    binSpacingType, plotFlag )
%
% outputs:
%   asymInfo: struct cantaining info on all the found horz. asymptotes
%   binInfo: To calculate asymptotes, data is binned. This contains bin info
%
% inputs: 
%   x: independent data
%   y: dependent data
%   errY: uncertainty in the dependent data
%   threshold: threshold value for calling a bin as a horz. asymptote
%   numBins: number of bins to divide data into
%   binSpacingType: type of bin spacing. either equidistant in indice number or x value
%   plotFlag: plotting flag
%
% plotinfo:
%   blue solid: data
%   black solid: max slope
%   red solid: horizontal asymptotes
%   green dash: anomalous time (intercept of max slope and asymptote)
%   faint black dot: bin divider
%   faint red dot: horizontal bin border
%%

function [asymInfo, binInfo] = getHorzAsymptotes( x, y, errY, threshold, numBins,  ...
  binSpacingType, plotFlag )

% break x up uniformly number of points
nPtns = length(x);
% set some defaults
if nargin == 2
  errY = [];
  threshold = 0.1;
  numBins = 10;
  binSpacingType = [];
  plotFlag = 0;
elseif nargin == 3
  threshold = 0.1;
  numBins = 10;
  binSpacingType = [];
  plotFlag = 0;
elseif nargin == 4
  numBins = 10;
  binSpacingType = [];
  plotFlag = 0;
elseif nargin == 5
  binSpacingType = [];
  plotFlag = 0;
elseif nargin == 6
  plotFlag = 0;
end

%if no weights, set them to one
if isempty(errY)
  w = ones(1, nPtns );
else
  w = 1 ./ errY .^ 2;
end
% make sure x,y,w are column vectors
if isrow(x); x = x.'; end
if isrow(y); y = y.'; end
if isrow(w); w = w.'; end

% linear
% dataBinInds = round( linspace( 1, nPtns, numBins+1) );
if isempty( binSpacingType ) || strcmp( binSpacingType, 'eqDeltaDist' )
  binSpacingType = 'eqDeltaDist';
  equalSpacesXVal = linspace( x(1), x(end), numBins+1 );
   dataBinInds = zeros( 1, numBins + 1 );
for ii = 1:numBins + 1
  [~,ind] = min( abs( x - equalSpacesXVal(ii) ) );
  dataBinInds(ii) = ind;
end
elseif strcmp( binSpacingType, 'eqDeltaInd' )
dataBinInds = round( linspace( 1, nPtns, numBins+1) );
else
  error('Cannot understand desired bin spacing')
end

% allocate
slopeBin = zeros(1, numBins);
yinterBin = zeros(1, numBins); 
xCenterBin = zeros(1, numBins);
numPtnsBin = zeros(1, numBins);
centerYValBin = zeros(1, numBins);
deltaYValBin = zeros(1, numBins);
thresTest = zeros(1, numBins);
% loop over bins
for ii = 1:numBins
  indStart = dataBinInds(ii);
  indEnd = dataBinInds(ii+1);
  xTemp = x( indStart:indEnd );
  yTemp = y( indStart:indEnd );
  wTemp = w( indStart:indEnd );
  pfit = fit( xTemp, yTemp, 'poly1', 'weights',  wTemp);
  slopeBin(ii) = pfit.p1;
  yinterBin(ii) =  pfit.p2;
  xCenterBin(ii) = ( xTemp(end) + xTemp(1) ) / 2;
  numPtnsBin(ii) = length(xTemp);
  centerYValBin(ii) = pfit.p1 .* xCenterBin(ii) + pfit.p2;
  deltaYValBin(ii) = pfit.p1 .* ( xTemp(end) - xTemp(1) );
  thresTest(ii) = abs( deltaYValBin(ii) ./ centerYValBin(ii) );
end

% get max slope and intercept
[maxSlopeSlopeVal, maxSlopeBin] = min( slopeBin );
maxSlopeInterceptVal = yinterBin( maxSlopeBin );

%find flat inds. If they are adjacent, plot them all together
acceptedFlatBins = thresTest < threshold;
failedInds = find( acceptedFlatBins == 0 );
maxFlats = length( failedInds ) + 1;
% add start and end
failedInds = [ 0 failedInds numBins+1 ];

flatInds = cell( 1, maxFlats  );
for ii = 1:length( failedInds )-1
  indRange = failedInds(ii)+1:failedInds(ii+1)-1;
  flatInds{ii} = indRange;
end
% get rid of empty cells
flatInds = flatInds( ~cellfun('isempty',flatInds) );
  
% now fit over range of flat data
numAsymptotes = length( flatInds );
slopeAsym = zeros(1, numAsymptotes);
yInterAsym = zeros(1, numAsymptotes); 
xCenterAsym = zeros(1, numAsymptotes);
numPtnsAsym = zeros(1, numAsymptotes);
centerYValAsym = zeros(1, numAsymptotes);
deltaYValAsym = zeros(1, numAsymptotes);
maxSlopeAsymIntercept =  zeros(1, numAsymptotes);
asymBinEnd = zeros(1, numAsymptotes); 
asymBinStart = zeros(1, numAsymptotes); 
aveW = zeros(1, numAsymptotes); 
stdW = zeros(1, numAsymptotes); 
aveUw = zeros(1, numAsymptotes); 
stdUw = zeros(1, numAsymptotes); 
for ii = 1:numAsymptotes
  % get data
  bins2count = flatInds{ii};
  indStart = dataBinInds( bins2count(1) );
  indEnd = dataBinInds( bins2count(end) + 1 );
  asymBinEnd(ii) = bins2count(end) + 1;
  asymBinStart(ii) = bins2count(1);
  xTemp = x( indStart:indEnd );
  yTemp = y( indStart:indEnd );
  wTemp = w( indStart:indEnd );
  nPntsTemp = length(xTemp);
  % fit it
  pfit = fit( xTemp, yTemp, 'poly1', 'weights',  wTemp);
  slopeAsym(ii) = pfit.p1;
  yInterAsym(ii) =  pfit.p2;
  % center and intercept values
  xCenterAsym(ii) = ( xTemp(end) + xTemp(1) ) / 2;
  numPtnsAsym(ii) = nPntsTemp ;
  centerYValAsym(ii) = pfit.p1 .* xCenterAsym(ii) + pfit.p2;
  deltaYValAsym(ii) = abs( pfit.p1 .* ( xTemp(end) - xTemp(1) ) );
  maxSlopeAsymIntercept(ii) = ( centerYValAsym(ii) - maxSlopeInterceptVal )...
    / maxSlopeSlopeVal;
  % averages
  sumVec = sum( yTemp .* wTemp );
  sumW = sum( wTemp );
  aveWtemp =  sumVec / sumW;
  sumstd2num = sum( ( yTemp - aveWtemp ) .^ 2 .* wTemp );
  stdWtemp  = sqrt( sumstd2num ./ ( ( nPntsTemp - 1 ) ./ (nPntsTemp) .*  sumW ) );
  aveW(ii) = aveWtemp;
  stdW(ii) = stdWtemp;
  aveUw(ii) = mean( yTemp );
  stdUw(ii) = std( yTemp );
end
% Store data
% asymptote info
asymInfo.numAsymptotes = numAsymptotes;
asymInfo.flatInds = flatInds;
asymInfo.maxNegSlope = maxSlopeSlopeVal;
asymInfo.aveW = aveW;
asymInfo.stdW = stdW;
asymInfo.aveUw = aveUw;
asymInfo.stdUw = stdUw;
asymInfo.slopeAsym = slopeAsym;
asymInfo.xCenterAsym  = xCenterAsym;
asymInfo.numPtnsAsym = numPtnsAsym;
asymInfo.centerYValAsym = centerYValAsym;
asymInfo.deltaYValAsym = deltaYValAsym;
asymInfo.maxSlopeAsymXIntercept  = maxSlopeAsymIntercept;
% bin info
binInfo.numBins = numBins;
binInfo.binSpacingType = binSpacingType;
binInfo.slopeBin = slopeBin;
binInfo.xCenterBin = xCenterBin;
binInfo.centerYValBin = centerYValBin;
binInfo.deltaYValBin = deltaYValBin;
binInfo.acceptedFlatBins = acceptedFlatBins;
%% Plotting %%
if plotFlag
% plot
figure()
plot( x, y );
ax = gca;
axisYLim =  ax.YLim;
ax.XLim = [ x(1) x(end) ];
hold
% max slope line
plot( x, maxSlopeSlopeVal * x + maxSlopeInterceptVal,'k-');
ax.YLim = axisYLim;
% plot bins
slopeLarge = 1000 .* max(y);
transparVal = 0.25;
for ii = 1:numBins-1
  indEnd = dataBinInds(ii+1);
  if isempty( intersect( ii+1, [asymBinEnd asymBinStart] ) )
    p = plot( x, slopeLarge .* ( x - x(indEnd) ), 'k:' );
  else
    p = plot( x, slopeLarge .* ( x - x(indEnd) ), 'r:' );
  end
  p.Color(4) = transparVal;
end
% plot asymptote regions
onesVec = ones( 1, nPtns );
for ii = 1:numAsymptotes
  % asymptote
  p = plot( x, centerYValAsym(ii) .* onesVec,'-' );
  p.Color = [0.75 0 0];
  % Their intercept
  p = plot( x, slopeLarge .* ( x - maxSlopeAsymIntercept(ii) ), '--'  );
  p.Color = [0 0.5 0];
end
title('Asymptotes, bins, max slope, and intercepts!');
ylabel('$$ y $$');
xlabel('$$ x $$');
hold off
end