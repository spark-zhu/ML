
% binranges = 0:150; % 150 is randomly selected to be as the common peak. 
% [bincounts] = histc(dis,binranges);
maxCutOff = 150;
Num_cluster_gt2units=zeros(1,maxCutOff);
Optimal_avg_num_unit_in_one_clu_ge3=zeros(1,maxCutOff);
Max_Num_unit_in_one_cluster=zeros(1,maxCutOff);

for cutOff = 1:maxCutOff 
    cutOff/maxCutOff 
BreakThrough = zeros(numel(DataForFeature),2);
BreakThrough(:,1)=1:numel(DataForFeature);
for i=1:numel(DataForFeature)
if BreakThrough(i,2)==0
BreakThrough(i,2)=max(BreakThrough(:,2))+1;


CurrentListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));
for j=1:numel(DataForFeature)
if min(MatrixUsed(find(BreakThrough(:,2)==BreakThrough(i,2)),j))<= cutOff
BreakThrough(j,2)=BreakThrough(i,2);
end
end
NewListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));



while numel(NewListForThisCluster)~=numel(CurrentListForThisCluster)
    CurrentListForThisCluster=NewListForThisCluster;
for j=1:numel(DataForFeature)
if min(MatrixUsed(find(BreakThrough(:,2)==BreakThrough(i,2)),j))<= cutOff
BreakThrough(j,2)=BreakThrough(i,2);
end
end
NewListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));
end





else
    continue;
end

end


sortList = sort(BreakThrough(:,2));
sortListCount =0;
for i=1:max(sortList)
sortListCount(i)=histc(sortList,i);
end 
Num_cluster_gt2units(cutOff)= numel(sortListCount(sortListCount>1)); % inverse u shape expected, should have optimal. 
% Med_Num_unit_in_one_cluster_g2(cutOff)= median(sortListCount(sortListCount>1))
Optimal_avg_num_unit_in_one_clu_ge3(cutOff)=  mean(sortListCount(sortListCount>2)); % shouldn't be too high or too low
Max_Num_unit_in_one_cluster(cutOff)= max(sortListCount); % shouldn't be unresonably high, monotonic increase. 
end
figure 
plot(Num_cluster_gt2units)
hold on 
plot(Optimal_avg_num_unit_in_one_clu_ge3*10,'r')
hold on
plot(Max_Num_unit_in_one_cluster,'k')
axis([0 50 0 100])
legend('Num cluster gt2units','Optimal avg num unit in one clu ge3*10','Max Num unit in one cluster')