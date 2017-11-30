function startLocations = findWellCenters(x, tl, Ef,k,threshold)
startLocations = [];
w = 2*sqrt(2*Ef/k);
for i=1:length(tl)
    testPositions = find(x(:,2) == i);
    sp = find(abs(x(testPositions,1)-tl(i)) < threshold*w);
    startLocations = vertcat(startLocations,testPositions(sp));
end
end