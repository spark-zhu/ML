repeat=(diff(MyTimes)==0);% 4357 instance
repeat(end+1)=0; % which item in cluster is repeated. 
list1=find(repeat==1); % first item in the repeating pair
list2=list1+1; % second item in the repeating pair 
KlusterRepeat(:,1)=cluster(list1); % convert index to actual kluster number 
KlusterRepeat(:,2)=cluster(list2); % Note the array length drop from cluster to KlusterRepeat


% can the same cluster repeat iteself? 
selfRepeat=(KlusterRepeat(:,2)-KlusterRepeat(:,1))==0;
selfRepeat=find(selfRepeat==1); % give all row index in KlusaterRepeat.
% back convert to the time sample. 
y(:,1)=list1(selfRepeat);
y(:,2)=list2(selfRepeat);
% verification
Z(:,1)=cluster(list1(selfRepeat)); % verfied, all elements in the pair have same cluster number ,and same time,which is unclear why the elements are just counted as one spike)
Z(:,2)=cluster(list2(selfRepeat));
Z(:,3)=MyTimes(list1(selfRepeat))-MyTimes(list1(selfRepeat));
% hypothesis, the pair has different mask vector, that's why they are
% identified as 2 spikes from the same kluster, but detected on 2 different
% electrodes( not connected to each other , so not counted as 1 spike)
m1=mask(:,:,y(1,1));
m2=mask(:,:,y(1,2));
scatter(1:54,m1(2,:));hold on;scatter(1:54,m2(2,:),'r')
% verified. so two spikes could be at the same time, classified to be from
% the same kluster ( same neuron ) but seen on 2 different place of
% electrode grids. Perefectly resonable. 
% instancey y(1,1) is largely detected on ch1, y(1,2) is largely detected
% on ch2