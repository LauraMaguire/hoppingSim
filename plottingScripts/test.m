DF = 1;
lp = 1;
lc = 100;
koffList = logspace(-6,2,100);
alph = zeros(length(koffList),length(results.dtime));
for koffIndex=1:length(koffList)
   alph(koffIndex,:) = (1 - exp(-2*DF*(3/(2*lc*lp))*results.dtime))/(3/(2*lc*lp));
    %alph(koffIndex,:) = (1 - exp(-2*DF*paramOut.k*dtime))/paramOut.k;
end
%%
DF = 1;
lp = 1;
lc = 100;
alph = (1 - exp(-2*DF*(3/(2*lc*lp))*results.dtime))/(3/(2*lc*lp));

%%
figure
loglog(results.dtime,results.meanMSD,'o');
%loglog(msdList(1,:));
hold all
loglog(results.dtime,alph)
hold off

%%
y = (1.*koffList.*1.*lc)./(3.*1+koffList.*1.*lc);
figure
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
loglog(koffList,y,'k-');
hold all
for i=1:length(khopList)
    errorbar(koffList,d(i,:)',derr(i,:)','o')
    errorbar(koffList,d_normal(i,:)',derr_normal(i,:)','o')
end
legend({'Tether Model','0','0.004', '0.04', '0.4'})
clear i
hold off

%%
DF = 1;
lp = 1;
lc = 100;
alph = (1 - exp(-2*DF*(3/(2*lc*lp))*r.dtime{1}))/(3/(2*lc*lp));

%%
figure
for i=1:length(r.filename)
loglog(r.dtime{i},r.msd{i},'o');
%loglog(msdList(1,:));
hold all
end
loglog(r.dtime{1},alph)
hold off