clear all
clc
files=dir('*.mat');
%% importData into cells.
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
NumEleFeature = size(AllData{1}.Dat_V_Map,1);
for i=2:size(AllData,1)
NumEleFeature=[NumEleFeature; size(AllData{i}.Dat_V_Map,1)];
end

DataForFeature = AllData{1}.waveform;
for i=2:size(AllData,1)
DataForFeature=[DataForFeature;AllData{i}.waveform];
end

FRForFeature = AllData{1}.Avg_FR;
for i=2:size(AllData,1)
FRForFeature=[FRForFeature;AllData{i}.Avg_FR];
end

SNRForFeature = (AllData{1}.SNR);
for i=2:size(AllData,1)
SNRForFeature=[SNRForFeature;(AllData{i}.SNR)];
end

Max5FRForFeature = AllData{1}.Max_Ins_5s_FR;
for i=2:size(AllData,1)
Max5FRForFeature=[Max5FRForFeature;AllData{i}.Max_Ins_5s_FR];
end

CVForFeature = AllData{1}.peakCV;
for i=2:size(AllData,1)
CVForFeature=[CVForFeature AllData{i}.peakCV];
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


PeakChList = zeros(numel(P2PForFeature),1);
for i=1:numel(P2PForFeature)
pl = P2PForFeature{i};
plmax=find(pl==max(pl));
PeakChList(i)=plmax(1);
end

FWHM_summary = zeros(numel(P2PForFeature),1);
valley2pkTime_summary = zeros(numel(P2PForFeature),1);

for i=1:numel(P2PForFeature)
FWHM_summary(i) = FWHM(DataForFeature{i}(PeakChList(i),:),1);
valley2pkTime_summary(i) = P2P_time_sample(DataForFeature{i}(PeakChList(i),:),1);
end

clear PeakWaveform PeakWaveform_scale
for i=1:numel(DataForFeature)
PeakWaveform(i,:) = DataForFeature{i}(PeakChList(i),:);
PeakWaveform_scale(i,:)=(PeakWaveform(i,:)-min(PeakWaveform(i,:)))/range(PeakWaveform(i,:));
end

for i=1:numel(DataForFeature)
    i
[~, I] = sort(P2PForFeature{i},'desc');
for top = 1:3
[minV, minI] = min(DataForFeature{i}(I(top),:));
minI=minI(1);
[minstd] = min(StdForFeature{i}(I(top),minI)); 
top3CV(i,top) = minstd/abs(minV);
end
end


clear p2pDist p2pFit 

for i=1:numel(DataForFeature)
    i
AmpList = sort(P2PForFeature{i},'desc');
p2pDist(i,:)= AmpList(1:20);
[fitresult, gof] = createFit_lin7to20(6:19, p2pDist(i,7:20));
p2pFit(i,:) = [fitresult.p1 fitresult.p2];
%  [fitresult, gof]=createFit_offset_exp(0:19,p2pDist(i,:))
%  p2pFit1(i,:) = [fitresult.a fitresult.b fitresult.m fitresult.c];
end
Lucky6 = p2pDist(:,10)./p2pDist(:,1);
DeviateFromLinear = p2pFit(:,2)./p2pDist(:,1);

boxplot(DeviateFromLinear)
errorbar(0:19,mean(p2pDist(p2pDist(:,1)<350,:)),std(p2pDist(p2pDist(:,1)<350,:)))
hold on
plot(0:19,median(p2pDist))
H1 = mean(p2pDist(:,1:6)');
boxplot(H1) % thres 

H2 = mean(p2pDist(:,7:12)');
H3=H2.*(H1<=231);

boxplot(H2) % thres H1 231
t=1:1:28;
v2pT=histc(valley2pkTime_summary,t);
% best_kmeans(score)
% eva=evalclusters(score,'kmean','gap','KList',[7:13])
[~,score,latent]  = pca(PeakWaveform_scale,'NumComponents',5);
% [idx,C] = kmeans(score,36); 
%   figure
%   clear sample_wav
% for clu=1:numel(unique(idx))
%   subplot(6,6,clu)
%     plot(mean(PeakWaveform_scale(idx==clu,:)))
%     sample_wav(clu,:) = mean(PeakWaveform_scale(idx==clu,:));
% title(num2str(sum(idx==clu)))
% end
% sample_wav([15 19],:)=[];
% 
% plot(mean(PeakWaveform(valley2pkTime_summary<=15,:)))
% sampleWaveform(1,:) = mean(PeakWaveform(valley2pkTime_summary<=15,:));
% sampleWaveform(2,:) = mean(PeakWaveform(valley2pkTime_summary>15,:));
% sampleWaveform(3,:) = mean(PeakWaveform);
% sampleWaveform(4,:) = median(PeakWaveform);
% sampleWaveform(5,:) = wav2(:,2)';
% RHO = corr(sample_wav',PeakWaveform');
% RHO_max = max(RHO);
% RHO_thre= 0;

% for wav = 1:size(sample_wav,1)
%     wav
%     for waveform = 1:size(PeakWaveform,1)
%         
% [cr,lgs]= xcorr(sample_wav(wav,:)',PeakWaveform_scale(waveform,:)',5,'coeff');
% XHO(wav,waveform)=max(cr);
%     end
% end
% XHO_max = max(XHO);
% histogram(RHO_max,100)
% load('H:/bi-weekly-result/p2pFit.mat')
% for i =1:4
%     p2pFit(:,i)=p2pFit(:,i).*(p2pDist(:,1)<350);
% end

% rearranging p2pFit
% for fit=1:size(p2pFit,1)
% if p2pFit(fit,3)>p2pFit(fit,1)
%    temp =  p2pFit(fit,1);
%    p2pFit(fit,1)=p2pFit(fit,3);
%    p2pFit(fit,3)=temp;
%    temp =  p2pFit(fit,2);
%    p2pFit(fit,2)=p2pFit(fit,4);
%    p2pFit(fit,4)=temp;
% end
% end
% mean(p2pFit(fit,2)) 
%%
thres(:,1) = (p2pDist(:,1)<=450); % first condition in amplitdue. 350 for 0109
thres(:,2) = FRForFeature>=0.01; % 
thres(:,3) = ~(top3CV(:,2)>=0.6 & top3CV(:,3)>=0.6 & top3CV(:,1)>=0.6);
thres(:,4) = valley2pkTime_summary>=3 & valley2pkTime_summary<=17;
thres(:,5) = FWHM_summary >=3 &  FWHM_summary <=15;  
% thres(:,6) = ((LocAdj(:,1)<=322.5).*(LocAdj(:,1)>=-87.5).*(LocAdj(:,2)>=-87.5).*(LocAdj(:,2)<=202.5).*(LocAdj(:,3)<=100));
thres(:,6) = ((LocAdj(:,1)<=465).*(LocAdj(:,1)>=-85).*(LocAdj(:,2)>=-85).*(LocAdj(:,2)<=265).*(LocAdj(:,3)<=100));

thres(:,7) = DeviateFromLinear<=0.6;

remain_index = prod(thres');


LocThre = ~((LocAdj(:,1)<=272.5).*(LocAdj(:,1)>=-37.5).*(LocAdj(:,2)>=--37.5).*(LocAdj(:,2)<=152.5).*(LocAdj(:,3)<=50));
TimeFast = -(1./p2pFit(:,2));
TimeSlow = -(1./p2pFit(:,4));
NoiseData=DataForFeature(remain_index==0);
NoiseData2=p2pDist(TimeSlow>-500 & TimeSlow<0,:);
RHO_noise1= thres(remain_index==0,:);%DeviateFromLinear(DeviateFromLinear>0.6);%top3CV(top3CV(:,2)>=0.5 & top3CV(:,3)>=0.5 & top3CV(:,1)>=0.5,:);%DeviateFromLinear(DeviateFromLinear>0.4 & DeviateFromLinear<0.5); %p2pDist(p2pDist(:,1)>300,1);%p2pDist(p2pDist(:,1)>186,1);%DeviateFromLinear(DeviateFromLinear>0.55 & p2pDist(:,1)<350);


% RHO_noise1 =p2pFit(TimeFast<0,1);
% RHO_noise2 =p2pFit(TimeFast<0,2);
% RHO_noise3 =p2pFit(TimeFast<0,3);
% RHO_noise4 =p2pFit(TimeFast<0,4);
% 
load('ch_map_pink.mat')
%load('Ch_Map_20161207_right_New.mat')

Ch_Map = Ch_Map_new-15;
for wav=1:numel(NoiseData)
    h=figure

for ch=1:32
[X_off Y_off]=find(Ch_Map==ch);
X_off=X_off*50;
Y_off=Y_off*100;
plot([1:41]+X_off,Y_off+NoiseData{wav}(ch,:),'b');
hold on
end
  title(num2str(RHO_noise1(wav,:)));
pause(1.5)
% title([num2str(RHO_noise1(wav)) ' ' num2str(RHO_noise2(wav)) ' ' num2str(RHO_noise3(wav)) ' ' num2str(RHO_noise4(wav)) ])


print(h,'-dpng',[num2str(wav) '.png'])

end

VIF=cell(1,numel(DataForFeature));
for unit=1:numel(DataForFeature)
waveform = DataForFeature{unit};
zeroRow = sum(waveform')==0;
waveform(zeroRow,:)=[];
    R0 = corrcoef(waveform'); % correlation matrix
VIF{unit}=diag(inv(R0))';
end
mean_VIF = cellfun(@mean, VIF);
std_VIF = cellfun(@std, VIF);
min_VIF = cellfun(@min, VIF);
max_VIF = cellfun(@max, VIF);

% save('sampleWaveforms','sampleWaveform')
plot(mean(PeakWaveform(FWHM_summary==5,:)))
%%   After this point, remain_index will play a key role in every aspect. I really hate this but there 
FWHM_across_ses = cell(numel(AllData),1);
valley2pkTime_across_ses = cell(numel(AllData),1);
for i=1:size(AllData,1)
FWHM_across_ses{i} = FWHM_summary((lookUpTable(:,1)==i)&(remain_index==1)');
valley2pkTime_across_ses{i} = valley2pkTime_summary(lookUpTable(:,1)==i & remain_index'==1);
end

mean_FWHM = cellfun(@mean, FWHM_across_ses);
std_FWHM = cellfun(@std, FWHM_across_ses);

mean_p2ptime = cellfun(@mean, valley2pkTime_across_ses);
std_p2ptime = cellfun(@std, valley2pkTime_across_ses);

figure
tx=1:numel(DPS);
xlimit = [155 175];
% tx=9:40;

errorbar(DPS(tx),mean_FWHM(tx)*0.05,std_FWHM(tx)*0.05,'lineWidth',2)
hold on;
set(gca,'FontSize',14,'fontWeight','bold')
xlabel('Days')
ylabel('Time (ms)')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
errorbar(DPS(tx),mean_p2ptime(tx)*0.05,std_p2ptime(tx)*0.05,'lineWidth',2)
legend('FWHM','v2pTime')
xlim(xlimit)

p2plist = max(cell2mat(P2PForFeature)');
for i=1:size(AllData,1)
P2P_summary_across_sessions{i} = p2plist(lookUpTable(:,1)== i & remain_index'==1);
end


max_P2P = cellfun(@max, P2P_summary_across_sessions);
mean_P2P = cellfun(@mean, P2P_summary_across_sessions);
median_P2P= cellfun(@median, P2P_summary_across_sessions);
std_P2P =  cellfun(@std, P2P_summary_across_sessions);
figure
errorbar(DPS(tx),mean_P2P(tx),std_P2P(tx),'LineWidth',2)
hold on;
set(gca,'FontSize',14,'fontWeight','bold')
xlabel('Days')
ylabel('Amplitude (microV)')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
ylim([0 450])
% plot(DPS,median_P2P,'r')
hold on;
plot(DPS(tx),max_P2P(tx),'k--')
xlim(xlimit)
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

try
countsSave(:,:,unit)=counts; 
counts_MSave(:,:,unit)=counts_M; 
counts_LSave(:,:,unit)=counts_L;  
catch
end
Max_Ins_5s_FRSave(unit)=Max_Ins_5s_FR;
FR_avgSave (unit)= FR_avg;
end

GrptimeFeatures.counts=mean(countsSave,3);
GrptimeFeatures.counts_M=mean(counts_MSave,3);
GrptimeFeatures.counts_L=mean(counts_LSave,3);
GrptimeFeatures.Max_Ins_5s_FR=mean(Max_Ins_5s_FRSave);
GrptimeFeatures.FR_avg=mean(FR_avgSave);

counts = GrptimeFeatures.counts;
if (counts(1)+counts(2))<=0.01*numel(FiringTimeForThisUnit) % which means on average (across all sessions) , you see in 12 mins, less than 1 occasion, it fires within 1ms. 
    isSingle(gp)=1;
end
end
SingleCount = zeros(size(AllData));
MultiCount = zeros(size(AllData));

Single_fire_Count = zeros(size(AllData));
Multi_fire_Count = zeros(size(AllData));
FT_count  = cellfun(@numel,TimeFeatureAll);

for ses=1:numel(AllData)
SingleCount(ses) = sum(isSingle' &  (lookUpTable(:,1)==ses) & remain_index'==1);
MultiCount(ses) = sum((~isSingle') & (lookUpTable(:,1)==ses) & remain_index'==1);
Single_fire_Count(ses) = sum((remain_index==1)'.*FT_count.*(lookUpTable(:,1)==ses).*isSingle');
Multi_fire_Count(ses) = sum((remain_index==1)'.*FT_count.*(lookUpTable(:,1)==ses).*~isSingle');
recording_time(ses) = max(cell2mat(AllData{ses}.time'))/20E3;
end
unit_count = [SingleCount MultiCount];
unit_FT_count = [Single_fire_Count./recording_time'  Multi_fire_Count./recording_time'];

figure
bar(DPS(tx),unit_count(tx,:),'stacked')
xlabel('Days')
ylabel('Unit Count')
% title('Single-Units')
legend('single','Multi')
set(gca,'FontSize',14,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
xlim(xlimit)
figure
bar(DPS(tx),unit_FT_count(tx,:),'stacked')
xlabel('Days')
ylabel('Total Units Firing Count/second')
% title('Single-Units')
legend('single','Multi')
set(gca,'FontSize',14,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
xlim(xlimit)
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
store.FWHM_summary=FWHM_summary;
store.valley2pkTime_summary = valley2pkTime_summary;
store.PeakWaveform=PeakWaveform;
store.PeakWaveform_scale=PeakWaveform_scale;
store.unit_FT_count=unit_FT_count;
store.remain_index = remain_index;
store.unit_count=unit_count;
store.is_single = isSingle;

store.max_P2P = max_P2P;
store.mean_P2P = mean_P2P;
store.std_P2P =  std_P2P;

store.mean_FWHM = mean_FWHM;
store.std_FWHM =std_FWHM;

store.mean_p2ptime =mean_p2ptime;
store.std_p2ptime = std_p2ptime ;


store1=store;
save('store','store1','-v7.3');
dateTime=cell(1,1);
%% importDate continuous data type convert to random date
%copy from excel recording data and time

subFolderList = dir ; 

subFolderName = {subFolderList.name};
desiredFolder = cellfun(@(y) y(1)=='r' , subFolderName );
subFolderList=subFolderList(desiredFolder);
% subFolderList([end])=[] 
%%
FolderList=subFolderList;

DPS = 0;
AllData = cell(1,1);
subList_pre = 0;
subList_now = 0;
for Folder =1:numel(FolderList)
    Folder
    cd(FolderList(Folder).name)
    Datestr = FolderList(Folder).name
    cd (Datestr([end-4:end-3 end-1:end end-9:end-6]))
    Datestr=Datestr([end-4:end-3 end-1:end end-9:end-6]);
  mainFolder = pwd;
  
   subFolderList = dir('datasets/'); 
     subFolderName = {subFolderList.name};
     desiredFolder = cellfun(@(y) y(1)=='d' , subFolderName );
     subFolderList=subFolderList(desiredFolder);
     subList_now=numel(subFolderList);
     smallOffset = datevec(between(datetime(dateTime{Folder},'Format','yyyy-MM-dd_HH-mm-ss'),datetime(dateTime{Folder+1},'Format','yyyy-MM-dd_HH-mm-ss')));
     if (smallOffset(4)+smallOffset(3))~=0
         
     smallOffset = smallOffset(4)+smallOffset(3)*24- subList_pre +1;
     subList_pre=subList_now
     else
         smallOffset=smallOffset(4);
         subList_pre=subList_now
     end
     DPS
     for subFolder =1:numel(subFolderList)-1
     DPS(end+1)=DPS(end)+1+smallOffset;
     smallOffset=0;
          dsName = ['ds' num2str(subFolder)];

          files=dir(['datasets\' dsName '\Refine_mountain\' Datestr '_tracking.mat']);
% importData into cells.


AllData{end+1}=importdata( [files.folder '\' files.name]);

     end 
      cd ..
      cd ..
end

DPS(1) = [];
AllData(1)=[];
WPS=DPS;


%% We can do a plot at this point
% load('/Users/Spark/Xie_software/tracking analysis/Ch_Map_20161207_right.mat')
% Dat_V_Map(:,2)=selectedChannels;
% Dat_V_Map(:,1)=1:size(selectedChannels,1);
% 
% for i=1:size(Ch_Map,1)
%     for j=1:size(Ch_Map,2)
%     if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
%     Ch_Map(i,j)=0;
%     end
%     end
% end
% 
% pdfCellArray=cell(1,length(DataForFeature));
% for clu=1:length(DataForFeature)
%     pdfCellArray{clu}=[num2str(clu) '.pdf'];
% MeanPlot=DataForFeature{clu};
% % 4-by-8 Channel Plot showing correlation by visual inspection
% tOffset=(size(DataForFeature{i},2)-1)/2;
% h=figure;
% t=[-tOffset/20e3:1/20e3:tOffset/20e3];
% for row=1:4
%     for col=1:8
%     subplot(4,8,(row-1)*8+col)
%     if Ch_Map(row,col)~=0
%     temp=find(Dat_V_Map(:,2)==Ch_Map(row,col));
%     plot(t,MeanPlot(temp,:),'LineWidth',2)
%     axis([t(1) t(end) min(min(MeanPlot))-10 max(max(MeanPlot))+10])
%     end
%     end
% end
% title(['sample-' num2str(clu)  ' day-' num2str(lookUpTable(clu,1))])
% print(h,'-dpdf',pdfCellArray{clu})
% end
%% serialization 
% DataForFeature=cellfun(@(x)reshape(x',numel(x),1),DataForFeature,'UniformOutput',0); %NOTE !!! the reshape is not trival, do double check you are getting what you expected.
% 
%% purely finding the largest waveform/ Peak Amplitude Histogram.
% Feature_LargestWaveform=cell(size(DataForFeature)); % naively ignoring positive peaks. 
% Feature_PA_Histogram=cell(size(DataForFeature));
% Feature_Min_each_Channel=cell(size(DataForFeature));
% Feature_Max_each_Channel=cell(size(DataForFeature));
% Feature_P2P_each_Channel=cell(size(DataForFeature));
% for i=1:numel(DataForFeature)
%     Feature_Min_each_Channel{i}=min(DataForFeature{i}');% what you get is the min from each channel.
%     Feature_Max_each_Channel{i}=max(DataForFeature{i}');
%     Feature_P2P_each_Channel{i}=Feature_Max_each_Channel{i}-Feature_Min_each_Channel{i};
%     Feature_P2P_each_Channel_sqrd{i}= Feature_P2P_each_Channel{i}.^2;
%     % find the global min, may not be useful themselves, but is useful to
% % extract the largest waveform. 
% [c1,i1]=min(DataForFeature{i});
% [c2,i2]=min(c1);
% gMin=c2;
% ChMin=i1(i2);
% tMin= i2;
% Feature_LargestWaveform{i} = DataForFeature{i}(ChMin,:);
% % find the global max. 
% [c1,i1]=max(DataForFeature{i});
% [c2,i2]=max(c1);
% gMax=c2;
% ChMax=i1(i2);
% tMax= i2;
% % why we calculate gEx, just to check for interneurons. 
% if abs(gMax)>=abs(gMin)
% gEx=gMax;
% else 
%     gEx=gMin;
% end
% 
% Feature_PA_Histogram{i}=[gMax,gMin,gEx];
% 
% ChMinList(i)=ChMin;
% end
% check for alignment
% figure
% for i=1:numel(DataForFeature)
% plot(Feature_LargestWaveform{i})
% hold on;
% end
% about to misaligning the global minimum by about 1 point. 

%% prepare for PCA analysis. 
% clear PCA_features
% for i=1:numel(DataForFeature)
% PCA_features(i,:)=Feature_P2P_each_Channel{i};
% end
% 
% [~,score,latent,tsquared,explained,mu] = pca(PCA_features,'Economy',0);
% % we could locate weighted center location. 
% [X,Y] = meshgrid([1 2 3 4 5 6 7 8],[1 2 3 4]); % Y gives row,  X gives col. 
% load('/Users/Spark/Xie_software/tracking analysis/Ch_Map_20161207_right.mat')
% 
% Ch_Map_2 = Ch_Map; 
% for i=1:size(Ch_Map,1)
%     for j=1:size(Ch_Map,2)
%     if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
%     Ch_Map_2(i,j)=0;
%     else 
%     Ch_Map_2(i,j)=find(Dat_V_Map(:,2)==Ch_Map(i,j));
%     end
%     end
% end
% for i=1:size(Ch_Map,1)
%     for j=1:size(Ch_Map,2)
%     if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
%     Ch_Map(i,j)=0;
%     end
%     end
% end
% Feature_for_WC=Feature_P2P_each_Channel;
% Feature_ReMap=cell(size(Feature_for_WC));
% for m=1:numel(Feature_for_WC)
% temp=zeros(4,8);
% 
% for i=1:4
%     for j=1:8
%         if Ch_Map_2(i,j)~=0
%         temp(i,j)=Feature_for_WC{m}(Ch_Map_2(i,j));
%         end
%     end
% end
% Feature_ReMap{m}=temp;
% end
% WC_row=0;WC_col=0;
% for i=1:numel(Feature_ReMap)
% WC_row(i)=sum(sum(Feature_ReMap{i}.*Y))/sum(sum(Feature_ReMap{i}));
% WC_col(i)=sum(sum(Feature_ReMap{i}.*X))/sum(sum(Feature_ReMap{i}));
% end
% % location information. Min Channel Location. ch that has the largest
% % negative waveform. 
% for i=1:numel(ChMinList)
% m=find(Ch_Map==selectedChannels(ChMinList(i)));
% row=mod(m,4);
% if row==0
%     row=4;
% end
% col=ceil(m/4);
% rowList(i)=row;
% colList(i)=col;
% end
% % zero mean everthing location series. 
% rowList=rowList-2.5;
% WC_row=WC_row-2.5;
% colList=colList-4.5;
% WC_col=WC_col-4.5;
% % note that std([1 2 3 4])=1.291,      std([1 2 3 4 5 6 7 8])=2.4495
% % note that std(PCA1)=43.18  std(PCA2)=19.25
% % adjustable gain. 
% gRow = (43+19)/2/1.291;
% gCol = (43+19)/2/2.4495;
% 
% rowList=rowList*gRow;
% colList=colList*gCol;
% WC_col=WC_col*gCol;
% WC_row=WC_row*gRow;
% 
%% PCA 3D plotting or PC_1_2 + one dimentional location info
%Third_Axis = score(:,3);Third_Name='PC3';
%Third_Axis = rowList; Third_Name='rowList';
% Third_Axis = WC_col;Third_Name='WC_col';
% Third_Axis = WC_row; Third_Name='WC_row';
% figure
% colorCode=['r','g','b','k'];
% for i = 1:numel(unique(lookUpTable(:,1)))
% list=find( lookUpTable(:,1)==i);
% scatter3(score(list,1),score(list,2),Third_Axis(list),20,colorCode(i),'filled')
% hold on;
% end
% xlabel('P1')
% ylabel('P2')
% zlabel(Third_Name)
% legend('1','2','3','4')

%% Clustering Feature Preparation. 
% Clu_Feat=zeros(length(Feature_LargestWaveform),4);
% Clu_Feat(:,1:2) = score(:,1:2);
% Clu_Feat(:,3:4) =  [WC_row' WC_col'];
% Clu_Feat=score(:,1:5);
% Clu_Feat=cell2mat(Feature_P2P_each_Channel');
% eva=evalclusters(Clu_Feat,'kmean','gap','KList',[2:25],'Distance','cityblock')
% % best_kmeans gives 10 as optimal, gap statistics gives 13 as optimal. 
% % we would probably do 10, 11, 12 ,13 as four group of plots. This is info
% % before Jan 18. On Jan19 citiblock gives 15 ?P2P correction for ch
% % weighted center. 
% KList=[13 15 18];
% clear idx
% for k=1:numel(KList)
% idx(:,k) = kmeans(Clu_Feat,KList(k),'replicate',5,'Distance','cityblock');
% end

%% plotting
% colorCode=['r','g','b','k'];
% for gra=KList
% mkdir(num2str(gra))
% cd(num2str(gra))
% cluID=idx(:,find(KList==gra));
% tOffset=20;
% t=[-tOffset/20e3:1/20e3:tOffset/20e3];
% pdfCellArray=cell(1,length(numel(unique(cluID))));
% 
% for clu=1:numel(unique(cluID))
%  pdfCellArray{clu}=[num2str(clu) '.pdf'];
%     h=figure;
% list=find(cluID==clu);
% 
% for row=1:4
%     for col=1:8
%     subplot(4,8,(row-1)*8+col)
%     if Ch_Map(row,col)~=0
%     temp=find(Dat_V_Map(:,2)==Ch_Map(row,col));
%     
%      for unit=1:numel(list)
%      plot(t,DataForFeature{list(unit)}(temp,:),colorCode(lookUpTable(find(lookUpTable(:,3)==list(unit)),1)))
%      hold on;
%      end
% 
%     
%     axis([t(1) t(end) min(min(DataForFeature{list(unit)}))-10 max(max(DataForFeature{list(unit)}))+10])
%     
%     
%     
%     end
%     end
% end
% print(h,'-dpdf',pdfCellArray{clu})
% save('pdfNames','pdfCellArray')
% end
% append_pdfs(['K=' num2str(gra) '.pdf'],pdfCellArray{:})
% save('IDX','cluID')
% cd ..
% close all
% pause(2);
% gra
% end
%% Do Pearson Correlation Coefficient Analysis.
% PearsonMatrix = zeros(numel(DataForFeature));
% MSEMatrix = zeros(numel(DataForFeature));
% SpearsMatrix = zeros(numel(DataForFeature));
% count=0;
% total=nchoosek(numel(DataForFeature),2);
% 
%  fprintf(1,'progress: ')
% for i=1:numel(DataForFeature)-1
%     for j=i+1:numel(DataForFeature)
%     count=count+1;
%     fprintf(1,'%4.2f',count/total); 
%     MaskThre=0.16;
%     select = (MaskFeatureAll{i}>MaskThre)|(MaskFeatureAll{j}>MaskThre);
% %     select = select > 0.1;
%     if sum(select>0)
%     selW1 = DataForFeature{i}(logical(select),:);
%     selW2 = DataForFeature{j}(logical(select),:);
%     corLis = zeros(1,size(selW1,1));
%     corLisSp = zeros(1,size(selW1,1));
%     MSELis = zeros(1,size(selW1,1));
%     for n=1:numel(corLis)
%     %corLis(n)=corr(selW1(n,:)',selW2(n,:)');
%     corLisSp(n)=corr(selW1(n,:)',selW2(n,:)','type','Spearman');
%     MSELis(n)=mean((selW1(n,:)-selW2(n,:)).^2);
%     end
%     
%     pTp1 =P2PForFeature{i}(logical(select));
%     pTp2 =P2PForFeature{j}(logical(select));
%     
%     Wa = sigmf(pTp1,[0.1 65]);
%     Wb = sigmf(pTp2,[0.1 65]);
%     W3 = -0.9*sigmf(min(pTp1,pTp2),[0.05 110])+1;
%     W1 = max(Wa,Wb); 
%     corLisSp(corLisSp<0.1)=0.1;% remove negative parts
%     corLisSp(isnan(corLisSp))=0.1;% remove again;
%     
%     W2 = (corLisSp).^(3/4);  % so long as it generally agrees it's fine.
%     % W2(W2>=5)=5;  % limit of shape difference penalization
%     Weight = W1;  % if any one of the pairs being compared have high P2P, this comparison gains weight. 
% 
% MSELis= MSELis.*W3;    
% MSELis= MSELis./W2;
%      
%      
%      
%     Wa = (MaskFeatureAll{i}(logical(select)));
%     Wb = (MaskFeatureAll{j}(logical(select)));
%     W2 = (1 - abs(Wa-Wb));
%     W2(W2<MaskThre)=MaskThre;  % If the mask value differs significantly, that's very bad as well, this penalization capped at 6.25 for thres=0.16;
%      W2(isnan(W2))=MaskThre;   % correct a second time. 
%     MSELis= MSELis./W2;
% 
%      
%      
%      
%     
%  
%     MSEMatrix(i,j) = sum(MSELis.*Weight)/sum(Weight);
%     else
%     PearsonMatrix(i,j) = NaN;
%     SpearsMatrix(i,j) = NaN;
%     MSEMatrix(i,j) = NaN;
%     end
%     fprintf(1,'\b\b\b\b')
%     end
% end
% fprintf('\n')
% 
% MatrixUsed=MSEMatrix;
% % RHO = corr((PCA_features').^3);
% dis = reshape(MatrixUsed,1,numel(DataForFeature)*numel(DataForFeature));
% dis = dis(dis>0);
% figure
% hist(dis,10000)
% t=[-10/20e3:1/20e3:30/20e3];
% MatrixUsed=MatrixUsed+MatrixUsed';
% uisave
% %%
% 
% 
% 
% BreakThrough = zeros(numel(DataForFeature),2);
% BreakThrough(:,1)=1:numel(DataForFeature);
% cutOff = 22; 
% for i=1:numel(DataForFeature)
%     i/(numel(DataForFeature)-1)
% if BreakThrough(i,2)==0
% BreakThrough(i,2)=max(BreakThrough(:,2))+1;
% 
% 
% CurrentListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));
% for j=1:numel(DataForFeature)
% if min(MatrixUsed(find(BreakThrough(:,2)==BreakThrough(i,2)),j))<= cutOff
% BreakThrough(j,2)=BreakThrough(i,2);
% end
% end
% NewListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));
% 
% 
% 
% while numel(NewListForThisCluster)~=numel(CurrentListForThisCluster)
%     CurrentListForThisCluster=NewListForThisCluster;
% for j=1:numel(DataForFeature)
% if min(MatrixUsed(find(BreakThrough(:,2)==BreakThrough(i,2)),j))<= cutOff
% BreakThrough(j,2)=BreakThrough(i,2);
% end
% end
% NewListForThisCluster = find(BreakThrough(:,2)==BreakThrough(i,2));
% end
% 
% 
% 
% 
% 
% else
%     continue;
% end
% 
% end
% figure
% hist(BreakThrough(:,end), max(BreakThrough(:,2)))
% % NumOfClu(thres/5-1) = max(BreakThrough(:,2));
% 
% % end
% %% Plotting All Clusters...If they contains more than 2 Units
% cluID=BreakThrough(:,2);
% 
% pdfCellArray=cell(1,length(numel(unique(cluID))));
% load('C:\Users\xie.lab.ws1\Documents\MATLAB\Hanlin_tracking\recording 01-13-17_right\rgb')
% for clu=1:numel(unique(cluID))
%  pdfCellArray{clu}=[num2str(clu) '.png'];
% 
% list=find(cluID==clu);
% if numel(list)>1
%         h=figure;
% for row=1:4
%     for col=1:8
%     subplot(4,8,(row-1)*8+col)
%     for unit=1:numel(list)
%         Ch_Map=AllData{lookUpTable(find(lookUpTable(:,3)==list(unit)),1)}.Ch_Map;
%         temp=Ch_Map(row,col)-15;
%     if Ch_Map(row,col)>0 && MaskFeatureAll{list(unit)}(temp)>0.16
%     
%     plot(t,DataForFeature{list(unit)}(temp,:),'Color',rgb(lookUpTable(find(lookUpTable(:,3)==list(unit)),1),2:4)/255)
%     hold on;
%     end
% 
%     axis([t(1) t(end) -100 25])
%     end
%     end
% end
% set(gcf, 'Position', get(0,'Screensize'));
% print(h,'-dpng',pdfCellArray{clu})
% end
% end
% save('pdfNames','pdfCellArray')
