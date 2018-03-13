figure
for unit = 1:10
scatter(weeksTracked{unit},unit*ones(1,numel(weeksTracked{unit})),'filled')
hold on
end
xlabel('Weeks post surgery','FontWeight','bold')
ylabel('Unit Number','FontWeight','bold')
axis([0 16 0 11])
title('Weeks that a certain unit is detected','FontWeight','bold')