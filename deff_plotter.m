d= dir('*');
%%
for k=4:length(d)
fname=d(k).name;
 
datatemp{k}=load(fname);

end

%%
for k=4:15
    kon(k)=datatemp{k}.results.konCalc;
    koff(k)=datatemp{k}.results.koffCalc;
    kd(k)=koff(k)/kon(k);
    DeffCalc(k) = datatemp{k}.DeffCalc;
    figure
    plot(datatemp{k}.results.Deff);
end


%%
x = logspace(-3,1);
loglog(DeffCalc(5:end),DeffObs(5:end),'bo');
hold all
loglog(x,x,'r-');
xlabel('Calculated Deff');
ylabel('Observed Deff');
