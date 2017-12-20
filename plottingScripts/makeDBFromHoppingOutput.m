function [DB,DBerr,kHopList,koffList,DBarray] ...
    = makeDBFromHoppingOutput(f,koffList)
% This scripts generates DB data using the outputs of the hopping
% simulation. Navigate to the folder in which the output files are located
% before running.  Folder should contain output files in ascending order of
% kHop (rHop). Assumes that lc is same for all files, and dtime and length
% of msd vector are same for all files.

% This script also generates sample plots of MSD and rho_MSD(t) for several
% values of khop and koff.  The plotting section assumes there are four
% files with rhop = 0,0.001,0.01,0.1, and that koffList =
% logspace(-3,-1,30).

% Inputs: f, the fraction of data to use in calculations.  f=0.9 uses first
%           90% of data and is usually a good value.
%         koffList, the list of koff values to use. Typical is
%           logspace(-3,-1,30).
% Output: DB, an array containing bound diffusion coefficients.  First
%           index is khop, second is koff.
%         DBerr, an array containing standard error in DB.  Same
%           organization as DB.
%         kHopList, list of kHop values
%         koffList, list of koff values (same as input)

%% User inputs
SetFigureDefaults(18,2); % first argument is default font size;
% second argument is default line width.

%% Load output files
r = LoadResults();

%% Determine contour length lc and rhop list from output files
lc = r.lc(1);
kHopList = r.khop;

%% Set up arrays for msd and error; define a time axis

disp('Loading MSD vectors.');

n = length(r.msd{1});
msdList = zeros(length(r.filename),f*n);
errList = zeros(length(r.filename),f*n);
for i=1:length(r.filename)
    %s = smooth(r.msd{i},1e2);
    s= r.msd{i};
    serr = r.errMean{i};
    msdList(i,:) = s(1:f*end);
    errList(i,:) = serr(1:f*end);
end
dtime = r.dtime{i}(1:f*end);

%% Set up array for lifetime distribution

disp('Generating lifetime distributions.');

distList = zeros(length(koffList),length(dtime));
for koffIndex=1:length(koffList)
    for tt=1:length(dtime)
        % rho = koff*e^(-koff*t)
        distList(koffIndex,tt) = koffList(koffIndex)*...
            exp(-koffList(koffIndex)*dtime(tt));
    end
end

%% Define the integrands
% The integrand is msd(rhop,t)*rho(koff,t).
% The lifetime integral goes in the denominator of the final expression and
% its integrand is rho(koff,t)*t.
% The error integrand is for calculating error in DB.

disp('Generating rho integrand.');

% Initialize arrays
integrand = zeros(length(kHopList),length(koffList),length(dtime));
ltintegrand = zeros(length(koffList),length(dtime));
errintegrand = zeros(length(kHopList),length(koffList),length(dtime));
% Loop over all koff, rhop, and t
for koffIndex = 1:length(koffList)
    for tt = 1:length(dtime)
        ltintegrand(koffIndex,tt) = dtime(tt)*distList(koffIndex,tt);
        for hopIndex =1:length(kHopList)
            integrand(hopIndex,koffIndex,tt) = ...
                msdList(hopIndex,tt)*distList(koffIndex,tt);
            errintegrand(hopIndex,koffIndex,tt) = ...
                errList(hopIndex,tt)*distList(koffIndex,tt);
        end
    end
end

%% Numerically integrate all integrands over time

disp('Numerically integrating.');

integral = sum(integrand,3);
lifetime = sum(ltintegrand,2);
err = sum(errintegrand,3);

%% Calculate bound diffusion coefficient and error
% In 1D, DB = integral / (2*lifetime).

disp('Calculating DB.');

% Initialize arrays
DB = zeros(length(kHopList),length(koffList));
DBerr = zeros(length(kHopList),length(koffList));
% Loop over all values of koff
for koffIndex = 1:length(koffList)
    DB(:,koffIndex) = integral(:,koffIndex)./(2*lifetime(koffIndex));
    DBerr(:,koffIndex) = err(:,koffIndex)./(2*lifetime(koffIndex));
end
%% Plot results with errorbars.

% Calculate expected results for rhop = 0:
y = (1.*koffList.*1.*lc)./(3.*1+koffList.*1.*lc);

figure
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% Convert x-axis from koff to KD (kon = 1e-3, diffusion-limited)
semilogx(koffList/1e-3,y,'k-');
hold all
for i=1:length(kHopList)
    errorbar(koffList/1e-3,DB(i,:)',DBerr(i,:)','o');
    %loglog(koffList,d(i,:)','o');
end

h = legend(horzcat('Tether Model',string(kHopList)));
ht = get(h,'Title');
set(ht,'String','$k_\mathrm{hop}$ ($\mu$s$^{-1}$)')
xlabel('Dissociation constant $K_D$ ($\mu$ M)');
ylabel('Bound diffusion ratio $D_B/D_F$');

%% Plot integrands
figure
plot(dtime,squeeze(integrand(1,1,:)));
hold all
plot(dtime,squeeze(integrand(1,15,:)));
plot(dtime,squeeze(integrand(1,30,:)));
axis([0,600,0,0.8]);
l = legend(string([koffList(1) koffList(15) koffList(30)]));
%l = legend({'$10^{-3}$','$10^{-2}$','$10^{-1}$'});
xlabel('Time ($\mu$s)');
ylabel('$\rho_{MSD}(t)$ (nm$^2$)');
v = get(l,'title');
set(v,'String','$k_{off}$ ($\mu$s$^{-1}$ $\mu$M$^{-1}$)')

%% Plot MSDs
figure
for i=1:length(kHopList)
loglog(dtime,squeeze(msdList(i,:)));
hold all
end
hold off
%l = legend({'$k_\mathrm{hop} = 0$ $\mu$s$^{-1}$','0.004 $\mu$s$^{-1}$','0.04 $\mu$s$^{-1}$' '0.4 $\mu$s$^{-1}$'});
l = legend(string(kHopList));
xlabel('Time ($\mu$s)');
ylabel('$\rho_{MSD}(t)$ (nm$^2$)');
v = get(l,'title');
set(v,'String','$k_\mathrm{hop}$ ($\mu$s$^{-1}$)')

%% Put results into a table that Mike will use to calculate selectivity

% Make one table per dataset (i.e. one table per Lc value). ALL VALUES ARE
% IN SI WITH NO PREFIXES.
% Five columns:
%   1) DF (m^2/s)
%   2) DB (m^2/s)
%   3) KD (M)
%   4) khop (s^-1)
%   5) error flag: 0 if DB value is true value, 1 if DB value is DB+DBerr,
%   -1 if DB value is DB-DBerr (this is a system to numerically propagate
%   error in DB into error in selectivity)
%
% Rows are grouped into blocks so that each large block represents one khop
% value, and smaller sub-blocks within give DB, DB+err, and DB-err.
%
% Mike Stefferson has a script that calculates selectivity using the
% resulting array.  Then he sends that result back to me and I use the
% script makeHoppingData.m to reformat it into something easy to use.
% Finally I send the makeHoppingData.m output back to Mike and he plots it
% in our standard format.  It's not that efficient but it works.

disp('Creating table for selectivity calculation.');

% Initialize an array of the correct size.
n = length(kHopList)*length(koffList);
m = length(koffList);
DBarray = zeros(n*3,5);

% Set DF = 1 nm^2/us = 1e-12 m^2/s for all values.
DBarray(:,1) = 1e-12;

% Fill in DB values as well as DB +/- err values (convert to m^2/s)
trueVal = [];
plusErr = [];
minuErr = [];
for i=1:length(kHopList)
    trueVal = horzcat(trueVal,DB(i,:));
    plusErr = horzcat(plusErr,DB(i,:)+DBerr(i,:));
    minuErr = horzcat(minuErr,DB(i,:)-DBerr(i,:));
end
DBarray(1:n,2) = 1e-12*trueVal;
DBarray(n+1:2*n,2) = 1e-12*plusErr;
DBarray(2*n+1:3*n,2) = 1e-12*minuErr;

% Fill in KD values (convert from koff using kon = 1e-3, then convert to M)
kdList = (koffList/1e-3)*1e-6;
DBarray(:,3) = repmat(kdList,1,3*length(kHopList))';

% Fill in error flags
DBarray(:,5) = vertcat(zeros(n,1),ones(n,1),-1*ones(n,1))';

% Fill in kHop values
k = [];
for i=1:length(kHopList)
    k = vertcat(k,kHopList(i)*ones(m,1));
end
DBarray(:,4) = repmat(k,3,1);

disp('Finished.');

end