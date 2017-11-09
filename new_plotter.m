test = results.meanMSD - 2*results.pfCalc*results.dtime;

%%
y = tr3(1:100);
x = 1:100;

%%
 p1 =      0.1939;  %(0.1936, 0.1941)
 p2 =     0.08254;  %(0.06954, 0.09553)

 
line = p1*1:10e5 + p2;
plot(line);

%%
hold off;
loglog(tr3);hold all; loglog(0.2*(1:1e5));
loglog(0.2*(1:1e5).^0.95);

%%
figure;
loglog(tr3./(2*results.dtime));
hold on;
loglog(tr5./(2*results.dtime));
loglog((0.2*(1:1e6).^0.85)./(2*results.dtime));

%%
%disp('hello');
n =2;
w = nan(length(r.filename),n);
for i=1:length(r.filename)
    %disp('loop');
    dt = r.dtime{i}; 
    logD = r.msd{i}./(2*dt);
    errLogD = r.errMean{i}./(2*dt);
    %plot(dt,r.errMean{i}./(2*dt));
    %loglog(dt,logD);
    hold all;
    tt = find(dt == 1e3);
%     for j=1:n
%         %disp(num2str(j));
%         errtest = errLogD(tt+n/2-j);
%         w(i,n+1-j) = 1/errtest^2;
%     end
%     r.newD(i) = sum(w(i,:).*logD(tt-floor(n/2):tt+floor(n/2-1)))/sum(w(i,:));
%     r.newDerr(i) = 1/sqrt(sum(w(i,:)));
    r.newD(i) = logD(tt);
    r.newDerr(i) = errLogD(tt);
%     r.newD(i) = mean(logD(tt-5:tt+5));
%     r.newDerr(i) = std(logD(tt-5:tt+5))/sqrt(10-1);
end
% h = legend({'0.4','0.04','0.004','0'});
% xlabel('ln(t/\mus)');
% ylabel('ln((<x^2>/2 t)/(nm^2/\mus))');

% set(get(h,'title'),'String','k_{hop} (\mus^-1)')
% clear h

hold off;
% disp('goodbye');

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
