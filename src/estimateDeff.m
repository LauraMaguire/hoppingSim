function [deffEst,deffErr] = estimateDeff(dtime, deff)
    start = round(0.5*length(deff));
    %t = dtime(1:end/2);
    figure
    plot(dtime,deff);
    hold all
    plot(dtime(start:end),deff(start:end));
    deffEst = nanmean(deff(start:end));
    deffErr = nanstd(deff(start:end));
end