function [msd,deff,dRes, dErr, deffCalc,dt,kon,koff,kd] = LoadResults()
z=2;

d=dir('*');
l = length(d);
deff = cell(1,l-z);
msd = cell(1,l-z);
dt = zeros(1,l-z);
kon = zeros(1,l-z);
koff = zeros(1,l-z);
deffCalc = zeros(1,l-z);
kd = zeros(1,l-z);
dRes = zeros(1,l-z);
dErr = zeros(1,l-z);
for k=(z+1):l
    fname=d(k).name;
    data=load(fname);
    deff{k-z} = data.results.Deff;
    msd{k-z} = data.results.meanMSD;
    dt(k-z) = data.paramOut.deltaT;
    kon(k-z) = data.results.konCalc;
    koff(k-z) = data.results.koffCalc;
    deffCalc(k-z) = data.DeffCalc;
    kd(k-z) = koff(k-z)/kon(k-z);
    [dRes(k-z), dErr(k-z)] = estimateDeff(dt(k-z)*1:length(deff{k-z}),deff{k-z});
end

end