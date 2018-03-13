% generate initial sortingFiles, create fake refine.xlsx
clear all
load('Y')
xlsContent = cell(numel(order),2);
xlsContent(:,1)  = arrayfun(@num2cell,order) ;
for i=1:numel(order)
xlsContent{i,2}  ='S';
end
xlswrite('Refine.xlsx',xlsContent);
%%
weighted_center_Refine
%%
mkdir('FR')
cd('FR')
%%
clear all
clc
files=dir('*.mat');
% importData into cells.
AllData = cell(size(files,1),1);
for i=1:size(AllData,1)
AllData{i}=importdata(files(i).name);

try
ReOrder(i) = AllData{i}.session;
catch
    AllData{i}.session = i;
    ReOrder(i) = AllData{i}.session;
end
    
end
[x,y]=sort(ReOrder);
AllData=AllData(y);
% Analysis Scheme 1, Easy, Lossy, USE the least amount of channels per day when comparing across different days. 
% Finding Sets Intersection 
% selectedChannels = AllData{1}.Dat_V_Map(:,2);
% for i=2:size(AllData,1)
% selectedChannels = intersect(selectedChannels,AllData{i}.Dat_V_Map(:,2));
% end
% 
% x(:,2)=selectedChannels;
% x(:,1)=1:size(x,1);
% Dat_V_Map=x;
% 
% for i=1:size(AllData,1)
% AllData{i}.s1waveform=AllData{i}.waveform;
% GiveS1_Ch = @(AgreedMapEle) find(AllData{i}.Dat_V_Map(:,2)==AgreedMapEle) ;
% selectedRows = arrayfun(GiveS1_Ch,selectedChannels);
% AllData{i,1}.s1waveform=cellfun(@(x) x(selectedRows,:), AllData{i,1}.s1waveform,'UniformOutput',0);
% end

% creating session - unit - pooled sample point - WPS (week post surgery) link
prompt='surgery Date in mmddyyyy: ';
surgeryDate = input(prompt,'s');
startNum = datenum(surgeryDate,'mmddyyyy');
DPS = zeros(1,numel(AllData));
for i=1:numel(AllData)
DPS(i) = datenum(AllData{i}.date,'mmddyyyy')-startNum;
WPS(i) = floor(DPS(i)/7);
end


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

StdForFeature = AllData{1}.waveformStd;
for i=2:size(AllData,1)
StdForFeature=[StdForFeature;AllData{i}.waveformStd];
end


MaskFeatureAll = AllData{1}.Intensity;
for i=2:size(AllData,1)
MaskFeatureAll=[MaskFeatureAll;AllData{i}.Intensity];
end

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
%P2P_Plot
P2PForFeature = AllData{1}.P2P;
P2P_summary_across_sessions=cell(1,numel(AllData));
P2P_summary_across_sessions{1} = max(cell2mat(P2PForFeature)');
for i=2:size(AllData,1)
P2PForFeature=[P2PForFeature;AllData{i}.P2P];
P2P_summary_across_sessions{i} = max(cell2mat(AllData{i}.P2P)');
end

max_P2P = cellfun(@max, P2P_summary_across_sessions);
mean_P2P = cellfun(@mean, P2P_summary_across_sessions);
median_P2P= cellfun(@median, P2P_summary_across_sessions);
std_P2P =  cellfun(@std, P2P_summary_across_sessions);

errorbar(DPS,mean_P2P,std_P2P)
hold on;
plot(DPS,median_P2P,'r')
hold on;
plot(DPS,max_P2P,'k')


ValleyAcrossTime=AllData{1}.ValleyAcrossTime;
% Unit Count Plot 
unit_count=cellfun(@(x) numel(x.list),AllData);
figure
plot(1:numel(AllData),unit_count);
% create data structure
store.AllData = AllData;
store.MaskFeatureAll = MaskFeatureAll;
store.DataForFeature=DataForFeature;
store.TimeFeatureAll=TimeFeatureAll;
store.lookUpTable=lookUpTable;
store.DPS=DPS;
store.files=files;
store.LocForFeature=LocForFeature;
store.StdForFeature=StdForFeature;
store.ValleyAcrossTime=ValleyAcrossTime;
store1=store;


save('store','store1');
%%
clear all
load('result-00.mat')
%%
validList = result.store.lookUpTable(result.AloneUnitList,3)
[~,~,raw] = xlsread('Refine.xlsx');
for i=1:size(raw,1)
if ~ismember(raw{i,1},validList)
    raw{i,2}='N';
end
end
xlswrite('Refine.xlsx',raw);
%%
weighted_center_Refine