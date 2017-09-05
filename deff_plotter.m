figure
hold all
plot(t1.results.dtime, t1.results.meanMSD);
plot(t2.results.dtime, t2.results.meanMSD);
plot(t3.results.dtime, t3.results.meanMSD);
plot(t4.results.dtime, t4.results.meanMSD);
%plot(t5.results.dtime, t5.results.meanMSD);
% plot(tr3.results.dtime(1:end), tr3.results.meanMSD(1:end));
% plot(tr4.results.dtime(1:end), tr4.results.meanMSD(1:end));
%legend('0.004','0.0044','0.02','0.04','0.01');
%plot(tr5.results.dtime(1:end/4), tr5.results.meanMSD(1:end/4));
%plot(tr6.results.dtime(1:end/8), tr6.results.meanMSD(1:end/8));

%%
figure
hold all
plot(t1.results.dtime(1:length(t1.results.Deff)), t1.results.Deff);
plot(t2.results.dtime(1:length(t2.results.Deff)), t2.results.Deff);
plot(t3.results.dtime(1:length(t3.results.Deff)), t3.results.Deff);
plot(t4.results.dtime(1:length(t4.results.Deff)), t4.results.Deff);
%plot(t5.results.dtime(1:length(t5.results.Deff)), t5.results.Deff);
% plot(tr3.results.dtime(1:end/2), tr3.results.Deff(1:end)/4);
% plot(tr4.results.dtime(1:end/2), tr4.results.Deff(1:end)/8);
legend('0.005','0.01','0.02','0.04');
%plot(tr5.results.dtime(1:end/8), tr5.results.Deff(1:end/4));
%plot(tr6.results.dtime(1:end/32), tr6.results.Deff(1:end/16));
%%
figure
hold all
plot(t6.results.dtime, t6.results.meanMSD);
plot(t7.results.dtime, t7.results.meanMSD);
plot(t8.results.dtime, t8.results.meanMSD);

% plot(tr3.results.dtime(1:end), tr3.results.meanMSD(1:end));
% plot(tr4.results.dtime(1:end), tr4.results.meanMSD(1:end));
legend('0.01','0.005','0.02');
%plot(tr5.results.dtime(1:end/4), tr5.results.meanMSD(1:end/4));
%plot(tr6.results.dtime(1:end/8), tr6.results.meanMSD(1:end/8));

%%

d=dir('*');
l = length(d);
deff = cell(1,l-3);
msd = cell(1,l-3);
dt = zeros(1,l-3);
kon = zeros(1,l-3);
koff = zeros(1,l-3);
deffCalc = zeros(1,l-3);
kd = zeros(1,l-3);
for k=4:l
    fname=d(k-3).name;
    data=load(fname);
    deff{k-3} = data.results.Deff;
    msd{k-3} = data.results.meanMSD;
    dt(k-3) = data.paramOut.deltaT;
    kon(k-3) = data.results.konCalc;
    koff(k-3) = data.results.koffCalc;
    deffCalc(k-3) = data.DeffCalc;
    kd(k-3) = koff(k-3)/kon(k-3);
end

%%

for k=1:8
    l = length(deff{k});
    tax = 1:l;
    plot(dt(k)*tax,deff{k})
    hold all
end
legend('show');

%%
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,deff{8});
plot(y);

%%
s = diff(deff{8});
y = s/dt(8);
plot(y);
hold all
y2 = filter(b,a,y);
plot(y2);

%%
xx = linspace(0,0.4);
plot(deffCalc,dRes,'bo');
hold all
plot(deffCalc(1:11),dRes(1:11),'ko');
plot(xx,xx,'r-');

%%
xx = linspace(0,0.4);
plot(deffCalc,dRes,'bo');
hold all
plot(deffCalc([3,7,11,15]),dRes([3,7,11,15]),'ko');
plot(xx,xx,'r-');

%%
xx = logspace(-3,0);
h = errorbar(deffCalc,dRes,dErr, 'ko');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%loglog(deffCalc,dRes,'ko');
hold all
plot(xx,xx,'r-');
%%
xx = logspace(-3,0);
h = errorbar(kd,dRes,dErr, 'ko');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%loglog(deffCalc,dRes,'ko');
hold all
%plot(x,x,'r-');

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
    errorbar(r.kd(hop{i}),r.dPost(hop{i}),r.dErr(hop{i}),'o');
end
l = legend(leg,'Location','southwest');
v = get(l,'title');
set(v,'String','kHop')
xlabel('KD (uM)');
ylabel('Deff / Dfree');

xx = logspace(-2,1);
y = (1.*xx.*1.*100)./(3.*1+xx.*1.*100);
semilogx(xx,y);

%%
figure
%h = errorbar(r.kd,r.dPost,r.dErr, 'ko');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold all

hopValues = unique(r.khop);
hop = cell(1,length(hopValues));
for i=1:length(hopValues)
    hop{i} = find(r.khop == hopValues(i));
    leg{i} = num2str(hopValues(i));
    errorbar(r.koff(hop{i}),r.dBound(hop{i}),r.dErr(hop{i}),'o');
end
leg{i+1} = 'Analytic Model';
xlabel('koff (us^-1)');
ylabel('Dbound / Dfree');

xx = logspace(-5,-2);
y = (1.*xx.*1.*100)./(3.*1+xx.*1.*100);
semilogx(xx,y);

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
    sel = results.bindFlux./results.nonbindFlux;
    errorbar(r.koff(hop{i}),sel(hop{i}),r.dErr(hop{i}),'o');
end

semilogx(xx,y./nb);
l = legend(leg,'Location','southwest');
v = get(l,'title');
set(v,'String','kHop')
xlabel('koff (us^-1)');
ylabel('Selectivity');


%%

xx=logspace(-5,-1,30);
p = struct();
p.AB = 20;
p.L = 100;
p.Nt = 1e3;
p.ll=100;
p.DF = 1;
p.kon = 1e-3;

disp('Entering loop');
for j=1:length(xx)
    p.koff = xx(j);
    disp(num2str(j));
    [y(j), nb(j)] = subNum(p,1);
end
