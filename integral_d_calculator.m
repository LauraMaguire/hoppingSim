
%%
f = 0.9;
lc = 100;
khopList = [0,0.001,0.01,0.1];
msdList = zeros(length(r.filename),f*1e5);
for i=1:length(r.filename)
    %s = smooth(r.msd{i},1e2);
    s= r.msd{i};
    msdList(i,:) = s(1:f*end);
end
dtime = r.dtime{i}(1:f*end);
clear i
%%
%koffList = [1e-4 1e-3 1e-2 1e-1];
koffList = logspace(-6,2,100);
distList = zeros(length(koffList),length(dtime));
for koffIndex=1:length(koffList)
    for tt=1:length(dtime)
        distList(koffIndex,tt) = koffList(koffIndex)*exp(-koffList(koffIndex)*dtime(tt));
    end
end
clear koff tt
%%
integrand = zeros(length(khopList),length(koffList),length(dtime));
ltintegrand = zeros(length(koffList),length(dtime));
for koffIndex = 1:length(koffList)
    for tt = 1:length(dtime)
        ltintegrand(koffIndex,tt) = dtime(tt)*distList(koffIndex,tt);
        for khopIndex =1:length(khopList)
            integrand(khopIndex,koffIndex,tt) = ...
            msdList(khopIndex,tt)*distList(koffIndex,tt);
        end
    end
end

clear khopIndex koffIndex tt
%%
integral = sum(integrand,3);
lifetime = sum(ltintegrand,2);
%%
d = zeros(length(khopList),length(koffList));
for koffIndex = 1:length(koffList)
    d(:,koffIndex) = koffList(koffIndex)*integral(:,koffIndex)/2;
    d(:,koffIndex) = integral(:,koffIndex)./(2*lifetime(koffIndex));
end
%%
y = (1.*koffList.*1.*lc)./(3.*1+koffList.*1.*lc);
figure
semilogx(koffList,y,'k-');
hold all
semilogx(koffList,d(:,:)')
legend({'Tether Model','0','0.004', '0.04', '0.4'})
hold off