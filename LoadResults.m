function r = LoadResults()
% LoadResults loads all hopping simulation output files in the current
% folder.  The folder should only contain hopping output files.

% Input: none.

% Output: r, a structure containing the fields described below. Fields may
% change as the simulation evolves.

%% Find all Matlab files in current folder; initialize r structure

d=dir('*.mat');
l = length(d);
r = struct();

%% Initialize arrays

% Cell arrays (to hold variable-length vectors)
r.filename       = cell(l,1); % filenames
r.msd            = cell(l,1); % msd vectors (nm^2)
r.errMean        = cell(l,1); % standard error in msd (nm^2)
r.dtime          = cell(l,1); % time axis vectors (us)

% Arrays of doubles
r.deltat         = zeros(1,l); % timestep (us)
r.lc             = zeros(1,l); % tether contour length (nm)
r.rhop           = zeros(1,l); % hopping rate parameter (dimensionless)
r.hopFreq        = zeros(l,1); % frequency of hopping while bound
r.hopOverageFreq = zeros(l,1); % frequency of hopping overage
r.khop           = zeros(1,l); % hopping rate (hops/us while bound)

%% Loop over all output files

for k=1:l
    
    disp(['Loading file ' num2str(k) ' of ' num2str(l) '.']);
    
    % Load hopping simulation output and record filename.
    fname               = d(k).name;
    data                = load(fname);
    r.filename{k}       = fname;

    % Record time axis, MSD, and error in MSD.
    r.dtime{k}          = data.results.dtime;
    r.msd{k}            = data.results.meanMSD;
    r.errMean{k}        = data.results.meanErr;

    % Record experimental input parameters.
    r.deltat(k)         = data.paramOut.deltaT;
    
    if isfield(data, 'paramOut.rHop')
    r.rhop(k)           = data.paramOut.rHop;
    else
    r.rhop(k)           = data.paramOut.kHop;
    end
    
    r.lc(k)             = data.paramOut.lc;
    
    % Record experimental outputs.
    r.hopFreq(k)        = data.results.hopFreq;
    r.hopOverageFreq(k) = data.results.hopOverageFreq;
    
    % Calculate hopping rate (hops/us while bound)
    r.khop(k)           = r.hopFreq(k)./r.deltat(k);

    % This section runs if binding AND unbinding were enabled.
    % Additional data can be calculated in that case.
    if (data.paramOut.bindFlag && data.paramOut.unbindFlag)
        
    % Initialize additional arrays
    r.kon               = zeros(1,l); % actual kon (us^-1 uM^-1)
    r.koff              = zeros(1,l); % actual koff (us^-1)
    r.kd                = zeros(1,l); % actual KD (M)
    r.onOverageCount    = cell(l,1);  % on-rate overage count
    
    % Fill arrays where possible
    r.kon(k)            = data.results.konCalc;
    r.koff(k)           = data.results.koffCalc;
    r.kd(k)             = r.koff(k)/r.kon(k);
    r.onOverageCount{k} = data.results.onOverageCount;
    end

end

disp('Done loading files.');

end