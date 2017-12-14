
%% Plots log-log of instantaneous D and adds newD entry to r
for i=1:length(r.filename)
    dt = r.dtime{i}; 
    logD = r.msd{i}./(2*dt);
    errLogD = r.errMean{i}./(2*dt);
    %plot(dt,r.errMean{i}./(2*dt));
    %loglog(dt,logD);
    hold all;
    tt = find(dt == 1e3);
    r.newD(i) = logD(tt);
    r.newDerr(i) = errLogD(tt);

end

hold off;

%%
r.ratio = r.khop./r.koff;

%%
for i=1:length(r.filename)
    r.boundLifetime(i) = 1/r.koff(i);
    r.hopLifetime(i) = 1/r.khop(i);
    r.wellLifetime(i) = (2*100*1)/(3*1);
end

%%
for i=1:length(r.filename)
    r.diff{i} = r.msd{i}./(2*r.dtime{i});
    p(i) = loglog(r.dtime{i},r.diff{i});
    hold all;
    lifetimeIndex(i) = datasample(find(abs(r.dtime{i}-r.boundLifetime(i)) <0.5),1);
    l(i) = loglog(r.boundLifetime(i), r.diff{i}(lifetimeIndex(i)),'ko','MarkerSize',15);
    if not(isinf((r.hopLifetime(i))))
        hoppingIndex(i) = datasample(find(abs(r.dtime{i}-r.hopLifetime(i)) <10),1);
        hp(i) = loglog(r.hopLifetime(i), r.diff{i}(hoppingIndex(i)),'kx','MarkerSize',15);
    end
    wellIndex(i) = datasample(find(abs(r.dtime{i}-r.wellLifetime(i)) <0.5),1);
    w(i) = loglog(r.wellLifetime(i), r.diff{i}(wellIndex(i)),'ko','MarkerSize',5,'MarkerFaceColor','k');
end
line(1e3*ones(1,1000),0.001*(1:1000),'color','k');
hold off
h = legend([p,w(1),hp(1),l(1)],'0.4','0.04','0.004','0','well exploration lifetime','hopping lifetime','binding lifetime');
%h = legend({'0.4','0.04','0.004','0'});
xlabel('ln(t/\mus)');
ylabel('ln((<x^2>/2 t)/(nm^2/\mus))');

set(get(h,'title'),'String','k_{hop} (\mus^-1)')
clear h
%% Plot DB as a function of KD for analytic model, several tether lengths
lp =1;
df =1;
lc = [10, 100, 500, 1000, 1e6];
koff = logspace(-5,0);
kd = 1e3*koff;
db = zeros(length(koff),length(lc));
for lcIndex=1:length(lc)
    for koffIndex = 1:length(koff)
        db(koffIndex,lcIndex) = (koff(koffIndex)*lc(lcIndex)*lp*df)/...
            (3*df+koff(koffIndex)*lc(lcIndex)*lp);
    end
        
end
clear koffIndex lcIndex
%%
SetFigureDefaults(20,3);
semilogx(kd,db);
xlabel('$$K_D$$ ($\mu$M)');
ylabel('$$D_B/D_F$$');
h = legend({'10','100','500','1000','$$\infty$$'});
xlim([1e-2 1e3]);
ylim([0 1]);
tx = '$$l_c$$ (nm)';
%set(get(tx,'Interpreter'),'latex');
set(get(h,'title'),'String',tx)
clear h
pbaspect([1 1 1])