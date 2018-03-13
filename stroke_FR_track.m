%% We break the stroke data in to three pieces(in time variable)


load('10182017_tracking.mat')
resultPre = result;
resultStroke = result;
resultPost = result;

resultPre.date='10152017';
resultStroke.date='10162017';
resultPost.date='10172017';

t1=[1200 2000]*20E3;
t2=[2300 3100]*20E3;
t3=[3200 4000]*20E3;

resultPre.time=cellfun(@(x) x(x>t1(1)&x<t1(2)),result.time,'UniformOutput',0);
resultStroke.time=cellfun(@(x) x(x>t2(1)&x<t2(2)),result.time,'UniformOutput',0);
resultPost.time=cellfun(@(x) x(x>t3(1)&x<t3(2)),result.time,'UniformOutput',0);


FR_Stroke=0;
FR_Stroke(1)= numel(cell2mat(resultPre.time'))/diff(t1/20E3)/size(result.Dat_V_Map,1);
FR_Stroke(2)= numel(cell2mat(resultStroke.time'))/diff(t2/20E3)/size(result.Dat_V_Map,1);
FR_Stroke(3)= numel(cell2mat(resultPost.time'))/diff(t3/20E3)/size(result.Dat_V_Map,1);

%%
% This code work on the assumption that selectd Cluster varibale is
% manually changed to select clusters. 


clear all
clc
files=dir('*.mat');
% importData into cells.
AllData = cell(size(files,1),1);
for i=1:size(AllData,1)
AllData{i} = importdata(files(i).name) ;

if numel(AllData{i}.selectedClusters)~=numel(AllData{i}.list)
    [C,ia] = setdiff(AllData{i}.list,AllData{i}.selectedClusters);
    AllData{i}.waveform(ia)=[];
    AllData{i}.waveformStd(ia)=[];
    AllData{i}.P2P(ia)=[];
    AllData{i}.Location(ia)=[];
    AllData{i}.ValleyAcrossTime(ia)=[];
    AllData{i}.time(ia)=[];
    AllData{i}.Avg_FR(ia)=[];
    AllData{i}.Max_Ins_5s_FR(ia)=[];
    AllData{i}.list(ia)=[];
end


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

max_P2P = cellfun(@max, P2P_summary_across_sessions);
mean_P2P = cellfun(@mean, P2P_summary_across_sessions);
median_P2P= cellfun(@median, P2P_summary_across_sessions);
std_P2P =  cellfun(@std, P2P_summary_across_sessions);
figure
errorbar(DPS,mean_P2P,std_P2P)
hold on;
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Amplitude (microV)')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
ylim([0 350])
% plot(DPS,median_P2P,'r')
hold on;
plot(DPS,max_P2P,'k')

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
save('store','store1');



%%
FR_final = [FR(1) FR_Stroke FR(2:end)'];
plot(FR_final,'LineWidth',2)
dayLabel ={'-1','Pr','S','Po','0','1','2','3','5','7','9','11'};
xticks(1:12)

xticklabels(dayLabel)
ylabel('Hz')
%% we redistribute and plot for each finger. 
 load('10182017_tracking-1.mat')
resultPre = result;
resultStroke = result;
resultPost = result;

resultPre.date='10152017';
resultStroke.date='10162017';
resultPost.date='10172017';

t1=[1200 2000]*20E3;
t2=[2300 3100]*20E3;
t3=[3200 4000]*20E3;



    load('D:\Box Sync\0919 stroke\Finger_map.mat');
Ch_Map= Ch_Map_new-15;
PeakCh   = cellfun(@(x)  find(x==max(x)),result.P2P) ;






Global_Time_Store = cell(4,8); % stored all the time stamps organized per channel  
Global_Unit_Store = cell(4,8); % store all the raw unit number(cluster label) per channel. 

for row=1:4
    for col=1:8
        unitOfThatCh = result.list(PeakCh==Ch_Map(row,col));
        validUnit = intersect(result.selectedClusters,unitOfThatCh);

    if ~isempty(validUnit)
    Global_Unit_Store{row,col} = validUnit;
     TimeStore = cell(1,numel(validUnit));    
    for unit=1:numel(validUnit)
    
    
 sample1 = result.time{result.list==validUnit(unit)};
% tt=0:0.1:time(end);
% allNum = Unit1Time;
% [counts] = histc(allNum,tt);
% % counts=counts/numel(group);
% counts=counts*10;

     TimeStore{unit}=sample1;       
            
            
    end
      Global_Time_Store{row,col}=TimeStore;
    end
    
    end
end
numOfElePerRow = zeros(4,1);
TimeStampPerRow = cell(4,1);
FR_Stroke=zeros(4,3);
for r=1:4
numOfElePerRow (r)= sum(cell2mat(cellfun(@(x) ~isempty(x),Global_Unit_Store(r,:),'UniformOutput',0)));
TimeStampPerRow{r}= cellfun(@(x) cell2mat(x),Global_Time_Store(r,:),'UniformOutput',0);
FR_Stroke(r,1)=numel(cell2mat(cellfun(@(x) x(x>t1(1)&x<t1(2)),TimeStampPerRow{r},'UniformOutput',0)))/numOfElePerRow (r)/diff(t1/20E3);
FR_Stroke(r,2)=numel(cell2mat(cellfun(@(x) x(x>t2(1)&x<t2(2)),TimeStampPerRow{r},'UniformOutput',0)))/numOfElePerRow (r)/diff(t2/20E3);
FR_Stroke(r,3)=numel(cell2mat(cellfun(@(x) x(x>t3(1)&x<t3(2)),TimeStampPerRow{r},'UniformOutput',0)))/numOfElePerRow (r)/diff(t3/20E3);
end

%FR = (Single_fire_Count+Multi_fire_Count)./cellfun(@(x) max(cell2mat(x.time')/20E3),AllData)./cellfun(@(x) size(x.Dat_V_Map,1),AllData);
%%
load('store.mat')
load('D:\Box Sync\0919 stroke\Finger_map.mat');
Ch_Map= Ch_Map_new-15;
for day=1:numel(store1.AllData)
    result=store1.AllData{day};
    PeakCh   = cellfun(@(x)  find(x==max(x)),result.P2P) ;






Global_Time_Store = cell(4,8); % stored all the time stamps organized per channel  
Global_Unit_Store = cell(4,8); % store all the raw unit number(cluster label) per channel. 

for row=1:4
    for col=1:8
        unitOfThatCh = result.list(PeakCh==Ch_Map(row,col));
        validUnit = intersect(result.selectedClusters,unitOfThatCh);

    if ~isempty(validUnit)
    Global_Unit_Store{row,col} = validUnit;
     TimeStore = cell(1,numel(validUnit));    
    for unit=1:numel(validUnit)
    
    
 sample1 = result.time{result.list==validUnit(unit)};
% tt=0:0.1:time(end);
% allNum = Unit1Time;
% [counts] = histc(allNum,tt);
% % counts=counts/numel(group);
% counts=counts*10;

     TimeStore{unit}=sample1;       
            
            
    end
      Global_Time_Store{row,col}=TimeStore;
    end
    
    end
end
numOfElePerRow = zeros(4,1);
TimeStampPerRow = cell(4,1);

for r=1:4
numOfElePerRow (r)= sum(cell2mat(cellfun(@(x) ~isempty(x),Global_Unit_Store(r,:),'UniformOutput',0)));
TimeStampPerRow{r}= cellfun(@(x) cell2mat(x),Global_Time_Store(r,:),'UniformOutput',0);
FR(r,day)=numel(cell2mat(cellfun(@(x) x ,TimeStampPerRow{r},'UniformOutput',0)))/numOfElePerRow (r)/max(cell2mat(result.time')/20E3);
if isnan(FR(r,day))
    FR(r,day)=0;
end
end


% FR = (Single_fire_Count+Multi_fire_Count)./cellfun(@(x) max(cell2mat(x.time')/20E3),AllData)./cellfun(@(x) size(x.Dat_V_Map,1),AllData);



end

%% FR_final = [FR(:,1) FR_Stroke_per_finger FR(:,2:end)];
h=figure
for i =1:4
subplot(4,1,i)
plot(FR(i,:),'LineWidth',2)
% dayLabel ={'-1','Pr','S','Po','0','1','2','3','5','7','9','11','13','22'};
% xticks(1:12)
% 
% xticklabels(dayLabel)
% ylabel('Hz')
% ylim([0 35])
end


