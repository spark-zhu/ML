% In this we want to investigat the usage of timestamps to affect merging
% of 2 cluster decision. Specifically, if spikes from 2 clusters always
% shows up at the same time(similar time), we would have a higher
% confidence that they are actually from the same neuron and hence the
% corresponding cluster should be merged. 

% we need the cluster and Time input from weighted_center.m

lag=20;% for each spike, we scan 20 adjacent time samples to look for another spike
% we do this for each cluster. 
cluster(cluster==0)=1;
lagMatrix=zeros(numel(unique(cluster)));
m=unique(cluster);
for i=1:numel(m) % cluster number
    list=find(cluster==m(i));
    i
    for j=1:numel(list) % index number counting the jth instance(spike) of cluster i
    
         timeList=(MyTimes>MyTimes(list(j)))&(MyTimes<=MyTimes(list(j))+lag);
         
         if sum(timeList)>0
            selected=cluster(timeList);
            for k=1:numel(selected)
            lagMatrix(i,selected(k))=lagMatrix(i,selected(k))+1;
            end
         end
        
    end
end

% give histogram plot for 32 clusters
for cl=1:numel(m)
figure 
bar((lagMatrix(cl,:)))
end

% probability estimation
h=histogram(cluster,'BinMethod','integers');
prier=h.Values;
prier=prier/sum(prier);
cl=5;
figure 
bar((lagMatrix(cl,:)),'y');
hold on;
plot(prier20*sum(lagMatrix(cl,:)),'r')
hold on;
plot(prier*sum(lagMatrix(cl,:)),'b')

prier20=sum(lagMatrix)/sum(sum(lagMatrix));
plot(prier);hold on;plot(prier20,'r')


p1Matrix=PlagMatrix;
for i=1:32
    for j=1:32
P2Matrix(i,j)=lagMatrix(i,j)/h.Values(i)/h.Values(j);
    end
end

for i=1:32
    
P3Matrix(:,i)=lagMatrix(:,i)/h.Values(i);
    
end


