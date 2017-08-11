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
x = linspace(0,0.4);
plot(deffCalc,dRes,'bo');
hold all
plot(deffCalc(1:11),dRes(1:11),'ko');
plot(x,x,'r-');

%%
x = linspace(0,0.4);
plot(deffCalc,dRes,'bo');
hold all
plot(deffCalc([3,7,11,15]),dRes([3,7,11,15]),'ko');
plot(x,x,'r-');