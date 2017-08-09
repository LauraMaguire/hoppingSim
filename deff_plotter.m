figure
hold all
plot(t1.results.dtime, t1.results.meanMSD);
plot(t2.results.dtime, t2.results.meanMSD);
plot(t3.results.dtime, t3.results.meanMSD);
plot(t4.results.dtime, t4.results.meanMSD);
plot(t5.results.dtime, t5.results.meanMSD);
% plot(tr3.results.dtime(1:end), tr3.results.meanMSD(1:end));
% plot(tr4.results.dtime(1:end), tr4.results.meanMSD(1:end));
legend('0.004','0.0044','0.02','0.04','0.01');
%plot(tr5.results.dtime(1:end/4), tr5.results.meanMSD(1:end/4));
%plot(tr6.results.dtime(1:end/8), tr6.results.meanMSD(1:end/8));

%%
figure
hold all
plot(t1.results.dtime(1:length(t1.results.Deff)), t1.results.Deff);
plot(t2.results.dtime(1:length(t2.results.Deff)), t2.results.Deff);
plot(t3.results.dtime(1:length(t3.results.Deff)), t3.results.Deff);
plot(t3.results.dtime(1:length(t4.results.Deff)), t4.results.Deff);
plot(t5.results.dtime(1:length(t5.results.Deff)), t5.results.Deff);
% plot(tr3.results.dtime(1:end/2), tr3.results.Deff(1:end)/4);
% plot(tr4.results.dtime(1:end/2), tr4.results.Deff(1:end)/8);
legend('0.005','0.01','0.02','0.001');
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