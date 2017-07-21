function BFSum = sumNearbyBFs(x, k, L, Ef,tetherIndex, allTethers)

BFSum = 0;
for i=1:length(tetherIndex)
    dist = wrapdistance(x, allTethers(tetherIndex(i)),L);
    BFSum = BFSum + exp(Ef-k*dist^2/2);
end
end