%%
figure
%h = errorbar(r.kd,r.dPost,r.dErr, 'ko');
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
hold all

tetherValues = unique(x(:,2));
th = cell(1,length(tetherValues));
for i=1:length(tetherValues)
    th{i} = find(x(:,2) == tetherValues(i));
    %leg{i} = num2str(round(mean(r.khop(hop{i})),3));
    plot(th{i},x(th{i},1),'o');
end
% leg{i+1} = 'Tether Model';
% 
% xlabel('K_D (\muM)');
% ylabel('D_B / D_F');
% 
% lc = 100;
% xx = logspace(-3,2);
% y = (1.*(xx*1e-3).*1.*lc)./(3.*1+(xx*1e-3).*1.*lc);
% plot(xx,y,'b-');
% %axis([9e-4 1e1 0 1]);
% l = legend(leg,'Location','southwest');
% v = get(l,'title');
% set(v,'String','k_{hop} (\mus^{-1})')