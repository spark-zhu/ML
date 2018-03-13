BreakThrough = zeros(numel(DataForFeature),2);
BreakThrough(:,1)=1:numel(DataForFeature);
for i=1:numel(DataForFeature)-1
    i/(numel(DataForFeature)-1)
if BreakThrough(i,2)==0
BreakThrough(i,2)=max(BreakThrough(i,2))+1;
for j=i+1:numel(DataForFeature)
if mean(PearsonMatrix(find(BreakThrough(:,2)==BreakThrough(i,2)),j))>0.99
BreakThrough(j,2)=BreakThrough(i,2);
end
end

else
    continue;
end

end