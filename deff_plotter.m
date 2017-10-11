

%%
figure
%h = errorbar(r.kd,r.dPost,r.dErr, 'ko');
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
hold all

hopValues = unique(r.khop);
hop = cell(1,length(hopValues));
for i=1:length(hopValues)
    hop{i} = find(r.khop == hopValues(i));
    leg{i} = num2str(hopValues(i));
    errorbar(r.kd(hop{i}),r.dPost(hop{i}),r.dErr(hop{i}),'o','MarkerSize',10);
end
leg{i+1} = 'Tether Model';

xlabel('KD (uM)');
ylabel('Deff / Dfree');

lc = 100;
xx = logspace(-3,1);
y = (1.*(xx*1e-3).*1.*lc)./(3.*1+(xx*1e-3).*1.*lc);
plot(xx,y,'b-');
axis([9e-4 1e1 0 1]);
l = legend(leg,'Location','southwest');
v = get(l,'title');
set(v,'String','kHop')
%%
figure
%h = errorbar(r.kd,r.dPost,r.dErr, 'ko');
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
hold all

hopValues = unique(r.khop);
hop = cell(1,length(hopValues));
for i=1:length(hopValues)
    hop{i} = find(r.khop == hopValues(i));
    leg{i} = num2str(hopValues(i));
    errorbar(r.koff(hop{i}),r.dBound(hop{i}),r.dErr(hop{i}),'o');
end
leg{i+1} = 'Tether Model';
xlabel('koff (us^-1)');
ylabel('Dbound / Dfree');

xx = logspace(-5,-2);
lc = 100;
y = (1.*xx.*1.*lc)./(3.*1+xx.*1.*lc);
semilogx(xx,y,'b-');

l = legend(leg,'Location','southeast');
v = get(l,'title');
set(v,'String','kHop')


%%
figure
%h = errorbar(r.kd,r.dPost,r.dErr, 'ko');
set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
hold all

hopValues = unique(r.khop);
hop = cell(1,length(hopValues));
for i=1:length(hopValues)
    hop{i} = find(r.khop == hopValues(i));
    leg{i} = num2str(hopValues(i));
    sel = r.bindFlux./r.nonbindFlux;
    err = r.bFluxErr./r.nonbindFlux;
    errorbar(r.kd(hop{i}),sel(hop{i}),err(hop{i}),'o','MarkerSize',10);
end

axis([1e-3 1e2 0 80]);
semilogx(xx/1e-3,y./nb,'b-');
l = legend(leg,'Location','southwest');
v = get(l,'title');
set(v,'String','kHop')
xlabel('KD (uM)');
ylabel('Selectivity');


%%
clear y
xx=logspace(-6,-1,30);
p = struct();
p.AB = 20;
p.L = 100;
p.Nt = 1e3;
p.ll=500;
p.DF = 1;
p.kon = 1e-3;

disp('Entering loop');
for j=1:length(xx)
    p.koff = xx(j);
    disp(num2str(j));
    [y(j), nb(j)] = subNum(p,1);
end
