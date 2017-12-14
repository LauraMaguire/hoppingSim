function out = makeHoppingData(selectivity)

out = struct();
%% 
out.kdVec = logspace(-6,-4,30);
out.kHopVec = [0,0.001,0.01,0.1];
%%
out.nuData = zeros(30,4);
out.nuErrData = zeros(30,4);
out.selData = zeros(30,4);
out.selErrPlusData = zeros(30,4);
out.selErrMinusData = zeros(30,4);
%% Fill in actual values for nu and selectivity
for khopIndex=1:4
    % set start and finish indices based on kHop
    fin = 30*khopIndex;
    start = fin - 29;
    % fill in nu data
    out.nuData(:,khopIndex) = selectivity.paramLoad(start:fin,2)./...
    selectivity.paramLoad(start:fin,1);
    % fill in selectivity data
    out.selData(:,khopIndex) = selectivity.val(start:fin);
end

%% Fill in values for error in nu
for khopIndex=1:4
    % set start and finish indices based on kHop
    fin = 120 + 30*khopIndex;
    start = fin - 29;
    % back-calculate error in DB
    DBerr = selectivity.paramLoad(start:fin,2) - ...
        selectivity.paramLoad((start-120):(fin-120),2);
    % fill in selectivity data
    out.nuErrData(:,khopIndex) = DBerr./selectivity.paramLoad(start:fin,1);
end

%% Fill in values for positive error in selectivity
for khopIndex=1:4
    % set start and finish indices based on kHop
    fin = 120 + 30*khopIndex;
    start = fin - 29;
    % fill in selectivity data
    out.selErrPlusData(:,khopIndex) = selectivity.val(start:fin);
end

%% Fill in values for negative error in selectivity
for khopIndex=1:4
    % set start and finish indices based on kHop
    fin = 240 + 30*khopIndex;
    start = fin - 29;
    % fill in selectivity data
    out.selErrMinusData(:,khopIndex) = selectivity.val(start:fin);
end

end
