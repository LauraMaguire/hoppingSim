function [] = makePartitionPlot()
% This function plots the concentration profiles across a Nup-filled medium
% for two different boundary conditions: equilibrium (in which the TF
% concentration is equal in both reservoirs) and gradient (in which T(0) =
% 0 and T(L) = TL).  It shows the partitioning of TF into the Nup-filled
% medium.

% There are no inputs to the function, but the system parameters need to be
% set in the first section.
%%
params = struct();
params.DF = 1;
params.DB = 1;
params.AB = 1;
params.Nt = 1e3;
params.ll = 100;
params.L = 100;
params.kon = 1e-3;
params.koff = 1e-2;
params.x = 0.01*params.L*(1:100);

[A,C] = ACSubNum(params,1,0);

plot(A+C);
%%
ACequil = (params.AB*params.kon*params.Nt/params.koff)*ones(1,100);
ACgrad = fliplr(A+C);

%%
ACleft = ones(1,100);
ACright = zeros(1,100);
ACequil = horzcat(ACleft,ACequil, ACright);
ACgrad = horzcat(ACleft,ACgrad, ACright);

%%
x = (-100:199);

%%
plot(x,ACequil)
hold all
plot(x,ACgrad, '--')
hold off
xlabel('Position x (nm)');
ylabel('Concentration T(x) + C(x) ($\mu$M)');
legend('Equilibrium', 'Transport');

end
