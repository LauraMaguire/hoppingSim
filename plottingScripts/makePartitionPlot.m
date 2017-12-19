function [] = makePartitionPlot()
% This function plots the concentration profiles across a Nup-filled medium
% for two different boundary conditions: equilibrium (in which the TF
% concentration is equal in both reservoirs) and gradient (in which T(0) =
% 0 and T(L) = TL).  It shows the partitioning of TF into the Nup-filled
% medium.

% There are no inputs to the function, but the system parameters need to be
% set in the first section.

% There are no outputs, but a plot is generated.
%% Set by hand: input parameters for linear flux solver.

% This is for the gradient BC case.
params = struct();
params.DF = 1; % Free diffusion coefficient (nm^2/us)
params.DB = 1; % Bound diffusion coefficient (nm^2/us)
params.AB = 1; % Concentration of free TF at boundary (uM)
params.Nt = 1e3; % Total Nup concentation (uM)
params.ll = 100; % Tether contour length (nm) (not used since DB is fixed)
params.L = 100; % Length of pore (nm)
params.kon = 1e-3; % on-rate (us^-1 uM^-1) (diffusion-limited)
params.koff = 1e-2; % off-rate (us^-1)
params.x = 0.01*params.L*(1:100); % positions within gel (nm)

% Calculate free (A) and bound (C) concentrations of transport factor
[A,C] = ACSubNum(params,1,0);

%plot(A+C);
%% Make A+C traces

% Calculate A+C in equilibrium case
ACequil = ((params.AB*params.kon*params.Nt/params.koff)+params.AB)*ones(1,100);
% Create A+C trace in gradient case
ACgrad = fliplr(A+C);

%% Set up left and right reservoir regions

ACleft = ones(1,100);
ACright = zeros(1,100);
ACequil = horzcat(ACleft,ACequil, ACright);
ACgrad = horzcat(ACleft,ACgrad, ACright);

%% Create an x-axis that puts the edge of the pore at x=0
x = (-100:199);

%% Plot both A+C traces
figure
plot(x,ACequil)
hold all
plot(x,ACgrad, '--')
hold off
xlabel('Position x (nm)');
ylabel('Concentration T(x) + C(x) ($\mu$M)');
legend('Equilibrium', 'Transport');

end
