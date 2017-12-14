%% Plot DB as a function of KD for analytic model, several tether lengths
function [kd,db] = makeTetherDBs()
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

semilogx(kd,db);
end