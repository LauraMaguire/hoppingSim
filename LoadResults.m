function [msd,deff,dRes, dErr, deffCalc,dt,kon,koff,kd] = LoadResults()

d=dir('*');
l = length(d);
deff = cell(1,l-3);
msd = cell(1,l-3);
dt = zeros(1,l-3);
kon = zeros(1,l-3);
koff = zeros(1,l-3);
deffCalc = zeros(1,l-3);
kd = zeros(1,l-3);
dRes = zeros(1,l-3);
dErr = zeros(1,l-3);
for k=4:l
    fname=d(k).name;
    data=load(fname);
    deff{k-3} = data.results.Deff;
    msd{k-3} = data.results.meanMSD;
    dt(k-3) = data.paramOut.deltaT;
    kon(k-3) = data.results.konCalc;
    koff(k-3) = data.results.koffCalc;
    deffCalc(k-3) = data.DeffCalc;
    kd(k-3) = koff(k-3)/kon(k-3);
    [dRes(k-3), dErr(k-3)] = estimateDeff(dt(k-3)*1:length(deff{k-3}),deff{k-3});
end

end