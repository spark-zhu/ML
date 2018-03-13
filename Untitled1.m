
BreakThrough = zeros(numel(DataForFeature),2);
BreakThrough(:,1)=1:numel(DataForFeature);
cutOff = 40; 
for i=1:numel(DataForFeature)
    i/(numel(DataForFeature)-1)
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
