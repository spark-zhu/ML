AllData = AllData';
lookUpTable = zeros(length(AllData{1}.list),4);
lookUpTable(:,1)=1;
lookUpTable(:,3)=AllData{1}.list;
lookUpTable(:,2)=1:length(AllData{1}.list);
lookUpTable(:,4)=WPS(1);

for i=2:size(AllData,1)
lookUpTable(end+1:end+length(AllData{i}.list),1)=i;
lookUpTable(end-length(AllData{i}.list)+1:end,2)=(1:length(AllData{i}.list))';
lookUpTable(end-length(AllData{i}.list)+1:end,4)=WPS(i);
lookUpTable(end-length(AllData{i}.list)+1:end,3)=AllData{i}.list;

end








% 
% 
% lookUpTable = zeros(length(AllData{1}.selectedClusters),4);
% lookUpTable(:,1)=1;
% lookUpTable(:,2)=1:length(AllData{1}.selectedClusters);
% lookUpTable(:,4)=WPS(1);
% 
% for i=2:size(AllData,1)
% lookUpTable(end+1:end+length(AllData{i}.selectedClusters),1)=i;
% lookUpTable(end-length(AllData{i}.selectedClusters)+1:end,2)=(1:length(AllData{i}.selectedClusters))';
% lookUpTable(end-length(AllData{i}.selectedClusters)+1:end,4)=WPS(i);
% end
% 
% lookUpTable(:,3)=1:size(lookUpTable,1);


% creating features. 
%appending
DataForFeature = AllData{1}.waveform;
for i=2:size(AllData,1)
DataForFeature=[DataForFeature;AllData{i}.waveform];
end


ValleyAcrossTime = AllData{1}.ValleyAcrossTime;
for i=2:size(AllData,1)
ValleyAcrossTime=[ValleyAcrossTime;AllData{i}.ValleyAcrossTime];
end


StdForFeature = AllData{1}.waveformStd;
for i=2:size(AllData,1)
StdForFeature=[StdForFeature;AllData{i}.waveformStd];
end


% MaskFeatureAll = AllData{1}.Intensity;
% for i=2:size(AllData,1)
% MaskFeatureAll=[MaskFeatureAll;AllData{i}.Intensity];
% end

if numel(AllData{1}.time)== numel(AllData{1}.waveform)
TimeFeatureAll = AllData{1}.time;
else
TimeFeatureAll = AllData{1}.time(~cellfun(@isempty,AllData{1}.time));    
end
for i=2:size(AllData,1)
    if numel(AllData{i}.time)== numel(AllData{i}.waveform)
TimeFeatureAll=[TimeFeatureAll; AllData{i}.time];
else
TimeFeatureAll=[TimeFeatureAll; AllData{i}.time(~cellfun(@isempty,AllData{i}.time))];    
end

end


try
WCForFeature = AllData{1}.WC;
for i=2:size(AllData,1)
WCForFeature=[WCForFeature;AllData{i}.WC];
end
catch
    LocForFeature = AllData{1}.Location;
for i=2:size(AllData,1)
    LocForFeature=[LocForFeature;AllData{i}.Location];
end
end
% Location Adjustment
LocAdj = cell2mat(LocForFeature);
LocAdj(:,3)=abs(LocAdj(:,3));
LocForFeature = mat2cell(LocAdj,ones(1,numel(LocForFeature)), [3]);
% P2P_Plot
P2PForFeature = AllData{1}.P2P;
P2P_summary_across_sessions=cell(1,numel(AllData));
P2P_summary_across_sessions{1} = max(cell2mat(P2PForFeature)');
for i=2:size(AllData,1)
P2PForFeature=[P2PForFeature;AllData{i}.P2P];
P2P_summary_across_sessions{i} = max(cell2mat(AllData{i}.P2P)');
end

max_P2P = cellfun(@max, P2P_summary_across_sessions,'UniformOutput',0);
mean_P2P = cellfun(@mean, P2P_summary_across_sessions,'UniformOutput',0);
median_P2P= cellfun(@median, P2P_summary_across_sessions,'UniformOutput',0);
std_P2P =  cellfun(@std, P2P_summary_across_sessions,'UniformOutput',0);
figure
errorbar(DPS,cell2mat(mean_P2P),cell2mat(std_P2P))
hold on;
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Amplitude (microV)')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
ylim([0 350])
% plot(DPS,median_P2P,'r')
hold on;
plot(DPS,cell2mat(max_P2P),'k')

%Unit Count Plot 
isSingle = zeros(1,numel(P2PForFeature));

for gp=1:numel(P2PForFeature)
    
   oneGroup = gp;
    
countsSave=zeros([size(0:1:50),numel(oneGroup)]); 
counts_MSave=zeros([size(0:20:1000),numel(oneGroup)]); 
counts_LSave=zeros([size(0:1000:20000),numel(oneGroup)]); 
Max_Ins_5s_FRSave = zeros(numel(oneGroup),1);
FR_avgSave = zeros(numel(oneGroup),1);

for unit=1:numel(oneGroup)
FiringTimeForThisUnit = TimeFeatureAll{oneGroup(unit)};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
FR_avg = numel(FiringTimeForThisUnit)/720;


countsSave(:,:,unit)=counts; 
counts_MSave(:,:,unit)=counts_M; 
counts_LSave(:,:,unit)=counts_L;  
Max_Ins_5s_FRSave(unit)=Max_Ins_5s_FR;
FR_avgSave (unit)= FR_avg;
end

GrptimeFeatures.counts=mean(countsSave,3);
GrptimeFeatures.counts_M=mean(counts_MSave,3);
GrptimeFeatures.counts_L=mean(counts_LSave,3);
GrptimeFeatures.Max_Ins_5s_FR=mean(Max_Ins_5s_FRSave);
GrptimeFeatures.FR_avg=mean(FR_avgSave);

counts = GrptimeFeatures.counts;
if counts(1)<3 % which means on average (across all sessions) , you see in 12 mins, less than 1 occasion, it fires within 1ms. 
    isSingle(gp)=1;
end
end
SingleCount = zeros(size(AllData));
MultiCount = zeros(size(AllData));

Single_fire_Count = zeros(size(AllData));
Multi_fire_Count = zeros(size(AllData));
FT_count  = cellfun(@numel,TimeFeatureAll);
for ses=1:numel(AllData)
SingleCount(ses) = sum(isSingle(find(lookUpTable(:,1)==ses)));
MultiCount(ses) = numel(AllData{ses}.list)-SingleCount(ses);
Single_fire_Count(ses) = FT_count(find(lookUpTable(:,1)==ses))'*isSingle(find(lookUpTable(:,1)==ses))';
Multi_fire_Count(ses) = FT_count(find(lookUpTable(:,1)==ses))'*~isSingle(find(lookUpTable(:,1)==ses))';
end
unit_count = [SingleCount MultiCount];
unit_FT_count = [Single_fire_Count  Multi_fire_Count];

figure
bar(DPS,unit_count,'stacked')
xlabel('Days')
ylabel('Unit Count')
% title('Single-Units')
legend('single','Multi')
set(gca,'FontSize',12,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')

unit_count=cellfun(@(x) numel(x.list),AllData);

figure
bar(DPS,unit_FT_count,'stacked')
xlabel('Days')
ylabel('Unit FT Count')
% title('Single-Units')
legend('single','Multi')
set(gca,'FontSize',12,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')

% create data structure
store.AllData = AllData;
% store.MaskFeatureAll = MaskFeatureAll;
store.DataForFeature=DataForFeature;
store.TimeFeatureAll=TimeFeatureAll;
store.lookUpTable=lookUpTable;
store.DPS=DPS;
store.files=files;
store.LocForFeature=LocForFeature;
store.StdForFeature=StdForFeature;
store.P2PForFeature=P2PForFeature;
store.ValleyAcrossTime=ValleyAcrossTime;
store1=store;
% save('store','store1','-v7.3');