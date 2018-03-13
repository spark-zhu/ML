% checking along clusters..give me an overall amplitude measure. sum of top
% three channel peak to peak. 
clear AloneSummary AllSummary
for i =1 :numel(AloneUnitList)
         pList = P2PForFeature{AloneUnitList(i)};
    sortedList = sort(pList,'descend');
           AloneSummary(i)= sum(sortedList(1:3))
end
    

for i =1 :numel(P2PForFeature)
         pList = P2PForFeature{i};
    sortedList = sort(pList,'descend');
           AllSummary(i)= sum(sortedList(1:3));
end



PercentageSummary=histc(AloneSummary/3,20:20:200)./histc(AllSummary/3,20:20:200)

