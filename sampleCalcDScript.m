% sampleCalcDscript( numBins, threshold )
% Description: sample script that takes a sample data set, 
% 'horztlAsympTestData' and calculate the diffusion coefficient. 
%
% sample call:
% [diffInfo, asymInfo, binInfo] = sampleCalcDScript( 10, 0.1 )
%
% 
% outputs:
%   diffInfo: struct containing diffusion coefficient, standard deviation, 
%     the minimum anomalous scaling coefficient, and the anomalous time
%   asymInfo: struct cantaining info on all the found horz. asymptotes
%   binInfo: To calculate asymptotes, data is binned. This contains bin info
%
% inputs: 
%   numBins: number of bins to divide data into
%   threshold: threshold value for calling a bin as a horz. asymptote
%     If the relative change of the independent variable over a bin is 
%     less than threshold, say that that bin is flat. 
%     Smaller threshold == stricter asymptote criteria.
% 

function [diffInfo, asymInfo, binInfo] = sampleCalcDScript( numBins, threshold )

if nargin == 0
  numBins = 10;
  threshold = 0.1;
elseif nargin == 1
  threshold = 0.1;
end
% load data
load('horztlAsympTestData');

% test find horizontal asymptote finds the horizontal asymptote for the input data.
% Therefore, we need to take the logs outside of the main function. I've build a 
% wrapper for this, but first, let's see what calling the main function does.
x = log( dtime );
y = log( msd ./ dtime );
% Calculate error from error propagation (I believe this is right). However,
% I believe 'err' is the std of the data and not the std of the mean. They differ
% by a factor of 1 / sqrt(N). This wouldn't matter if N was constant, but the 
% late time data has a smaller N due tot eh time windows. 
% I am pretty sure we should be using the std of the mean and not the std. 
% We should discuss if you disagree.
errY = err ./ msd; 
% (I think it should be this, but I dont have the data) errY = err ./ msd ./ sqrt(N)

% plotting flag. I'd turn this on to convince yourself that the code is working 
% of if you think it's not.
plotFlag = 1;
% binning type. Seperate bins by equidistance indices or values. 
% This should be value for us. But feel free to look at the other type, 'eqDeltaInd'
binSpacingType = 'eqDeltaDist';
% run find getHorzAsymptotes. It get all horizontal asymptotes, even early time ones!
[asymInfo, binInfo] = getHorzAsymptotes( x, y, errY, threshold, ...
  numBins, binSpacingType, plotFlag );
% display it
disp(binInfo);
disp(asymInfo);

% Now, run the wrapper. This is most likely what you will be calling. It will actually calculate
% the diffusion coefficients after running getHorzAsymptotes. It basically does the above code
numPntsMsd = ones( size(msd) ); % This should be the number of point recorded in average.
plotFlag = 1; % plot option again
verbose = 1; % verbose flag. print if it found a diffusion coeff or not
errMean = err ./ sqrt( numPntsMsd );
[diffInfo] = getDfromMsdData( dtime, msd, errMean, ...
  threshold, numBins, plotFlag, verbose );
% display it
disp(diffInfo);
