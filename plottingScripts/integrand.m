close all
%%

plot(dtime,squeeze(integrand(1,1,:)));
hold all
plot(dtime,squeeze(integrand(1,15,:)));
plot(dtime,squeeze(integrand(1,30,:)));
hold off
axis([0,600,0,0.8]);
l = legend({'$10^{-3}$','$10^{-2}$','$10^{-1}$'});
xlabel('Time ($\mu$s)');
ylabel('$\rho_{MSD}(t)$ (nm$^2$)');
v = get(l,'title');
set(v,'String','$k_{off}$ ($\mu$s$^{-1}$ $\mu$M$^{-1}$)')
%%
plot(dtime,squeeze(distList(1,:)));
hold all
plot(dtime,squeeze(distList(15,:)));
plot(dtime,squeeze(distList(30,:)));
plot(dtime,smooth(squeeze(msdList(1,:)),100));
plot(dtime,smooth(squeeze(msdList(2,:)),100));
plot(dtime,smooth(squeeze(msdList(3,:)),100));
plot(dtime,smooth(squeeze(msdList(4,:)),100));
hold off
%axis([0,600,0,0.1]);
l = legend({'$10^{-3}$','$10^{-2}$','$10^{-1}$'});
xlabel('Time ($\mu$s)');
ylabel('$\rho_{MSD}(t)$ (nm$^2$)');
v = get(l,'title');
set(v,'String','$k_{off}$ ($\mu$s$^{-1}$ $\mu$M$^{-1}$)')