% import optical analysis data? go into the reference folder of the dataset
% we need the time data ready
% ROI results in csv format
% intan output waveform data in xls 97-04 format
% we need variable p from align_elec_optic.m
% align_elec_optic
read_Intan_RHD2000_file
y=[amplifier_channels.native_order];
selection = (y<48 & y >15);
amplifier_channels=amplifier_channels(selection);

path=pwd;
k1 = findstr('TSeries',path);
k2 = findstr('\Ref',path);
datasetName=path(k1:k2-1);
%%
y=importdata('Results010.csv');
ROI_raw=y.data;
ROI_raw=ROI_raw(:,2:end);
numOfROI=size(ROI_raw,2)/4; % find out how many ROIs are included in this file.
ROI=zeros(size(ROI_raw,1),numOfROI);

for i=1:numOfROI
ROI(:,i)=ROI_raw(:,4*i-2); % this might need to be changed. 
F0trace = ROI(1:600,i);
F0sort=sort(F0trace);
F0=mean(F0sort(1:60));
ROI(:,i) = (ROI(:,i)-F0)/F0;
end
%% plottig all the Imaging trace
hold on
 plot((1:numel(w1))/oF,w1/max(w1)*5+20,'b');
hold on
plot((1:numel(w2))/oF,w2/max(w2)*5+25,'r');
hold on

for roi = 7%1:size(ROI,2)
Opti = ROI(:,roi);

plot(realTime,Opti+30,'k') %+(roi-1)*10
hold on
end
plot([240 240],[30,40],'Color','k','lineWidth',3)
text(240,31,'0');
text(240,40,'10');
hold on
plot([245 245],[20,25],'Color','b','lineWidth',3)
text(248,20,'0');
text(248,25,num2str(max(w1)));

hold on
plot([248 248],[25,30],'Color','r','lineWidth',3)
text(251,26,'0');
text(251,30,num2str(max(w2)));

plot([0 50],[45,45],'Color','r','lineWidth',3)
text(0,46,'0');
text(50,46,'50');


set(gca,'YTickLabel',[]) 
ylabel('df/F')
xlabel('Time(s)')
set(gca,'FontSize',12,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')

%%
try
load('Time.mat') % load timing data.
catch
    load('realTime.mat') % load timing data.

end
ElecList=dir('*.xls'); % get all electrical recording spike timing data file spec.
ElecRec=cell(size(ElecList,1),1);
RecContent=cell(size(ElecList,1),1);

for i=1:size(ElecList,1) 
book=xlsread(ElecList(i).name); % if this step fails, try open the xls file and save as 97-04 xls file.
Unit_time=book(2:end,2:3); % unit number in the second col, time in the third col
Unit_time(:,2)=Unit_time(:,2);%*p(1)+p(2);
waveform=book(2:end,7:end); % col 4:6 is PCA components. 
ElecRec{i}=Unit_time;
RecContent{i}=waveform;
% determine how many units are there.
unique(Unit_time(:,1))
end
% clear book y ROI_raw

% mean-std normalization.
% for i=1:numOfROI
%  ROI(:,i)=(ROI(:,i)-mean(ROI(:,i)))/std(ROI(:,i));
% end
%%
%  channels = unique(book(:,1)); %[11 12 16];
figure
for channel = 1:numel(channels)
    ch = channels(channel)
selected = book(:,1)== ch;
unitNum = book(selected,2);
if min(unitNum)>1
    unitNum=unitNum-1;
end 
waveform = book(selected,7:46);
time=book(selected,3);


FirUnitIndex = find(unitNum==1);
SecUnitIndex = find(unitNum==2);
try 
    ThirdUnitIndex = find(unitNum==3);
catch
end

% 
% subplot(1,2,2)
FirUnitWaveform = waveform(FirUnitIndex,:);
SecUnitWaveform = waveform(SecUnitIndex,:);

try 
    ThirdUnitWaveform = waveform(ThirdUnitIndex,:);
catch
end

t = (0:-1+size(FirUnitWaveform,2))/20;


 
Unit1Time = time(FirUnitIndex);
% group = cell(1,numel(RealtimeStim));
% for i = 1:numel(group)
% group{i} = Unit1Time((Unit1Time>=RealtimeStim(i)- before)&(Unit1Time<RealtimeStim(i)+after));
% end


Unit2Time = time(SecUnitIndex);
% group2 =  cell(1,numel(RealtimeStim));
% for i = 1:numel(group2)
% group2{i} = Unit2Time((Unit2Time>=RealtimeStim(i)- before)&(Unit2Time<=RealtimeStim(i)+after));
% end


try
    Unit3Time = time(ThirdUnitIndex);
% group3 =  cell(1,numel(RealtimeStim));
% for i = 1:numel(group3)
% group3{i} = Unit3Time((Unit3Time>=RealtimeStim(i)- before)&(Unit3Time<=RealtimeStim(i)+after));
% end
catch
    
end
for pt = 1:numel(Unit1Time)
plot([Unit1Time(pt) Unit1Time(pt)],[(channel-1)*3 (channel-1)*3+0.8],'b','linewidth',0.1)
hold on
end
for pt = 1:numel(Unit2Time)
plot([Unit2Time(pt) Unit2Time(pt)],[(channel-1)*3+1 (channel-1)*3+1.8],'r','linewidth',0.1)
hold on
end
end



yticks([0:3:24])
yticklabels({'1','2','3','4','5','6','7','8','9'})
% set(gca,'YTickLabel',[]) 
ylabel('Ch')
xlabel('Time(s)')
set(gca,'FontSize',12,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
%% plotting 
% clearvars -except realTime ROI datasetName ElecRec numOfROI
% colorList=['r','g','k','m','y']; % Assume maximum 5 units can be detected on those electrodes. 
% MarkerList=['*','o','x','p','h']; % Assume maximum 5 electrodes 
% for i=1:numOfROI  % create one figure for one ROI
% figure
% 
% h=plot(realTime,ROI(:,i),'b'); %ploting optical recording
% set(h,'LineWidth',1)
% legend(h,'Optical')
% title([datasetName ' ROI ' num2str(i)])
% xlabel('time(s)')
% ylabel('N std above mean(fluorescence Intensity)')
% hold on; 
% %plotting electrical recording spike timing data, one type of marker per
% %electrode, one color per unit. 
% 
% 
% 
% for electrode=1:numel(ElecRec)
% Unit_time=ElecRec{electrode};
% for rec=1:size(Unit_time,1)
%     % offset each electrical recording trail by 0.2*std of Optical recording. 
% h=scatter(Unit_time(rec,2),1.8+0.2*electrode,10,colorList(Unit_time(rec,1)),MarkerList(electrode));
% hold on; 
% end
%  %ElecList(electrode).name
% end
% 
% end
%% PLotting New
% channels = unique(book(:,1)); %[11 12 16];
for channel = 6%1:numel(channels)
    ch = channels(channel);
selected = book(:,1)== ch;
unitNum = book(selected,2);
if min(unitNum)>1
    unitNum=unitNum-1;
end 
waveform = book(selected,7:46);
time=book(selected,3);


FirUnitIndex = find(unitNum==1);
SecUnitIndex = find(unitNum==2);
try 
    ThirdUnitIndex = find(unitNum==3);
catch
end

% 
% subplot(1,2,2)
FirUnitWaveform = waveform(FirUnitIndex,:);
SecUnitWaveform = waveform(SecUnitIndex,:);

try 
    ThirdUnitWaveform = waveform(ThirdUnitIndex,:);
catch
end

t = (0:-1+size(FirUnitWaveform,2))/20;

% errorbar (t,mean(FirUnitWaveform),std(FirUnitWaveform),'Color',[.7 .7 .7])
% hold on

%% temp change color
% temp =SecUnitWaveform;
% SecUnitWaveform=FirUnitWaveform;
% FirUnitWaveform=temp;
%%





% plot(t,mean(FirUnitWaveform),'b','lineWidth',3)


% hold on 
% errorbar (t,mean(SecUnitWaveform),std(SecUnitWaveform),'Color',[.7 .7 .7])
% hold on








% plot(t,mean(SecUnitWaveform),'g','lineWidth',3)

% try
%     hold on 
%     plot(t,mean(ThirdUnitWaveform ),'r','lineWidth',3)
% 
% catch
% end


% xlabel('t(ms)')
% ylabel('Amplitude(uV)')

% stimWave(:,1)=stimWave(:,1)-stimWave(1,1);
% stimWave(:,1)=stimWave(:,1)*1000;
% plot(stimWave(:,1),stimWave(:,2))
% plot(stimWave(:,1),stimWave(:,2))
% plot(stimWave(:,1),stimWave(:,2),'LineWidth',2)
% 

 % find out t-critical.

 
Unit1Time = time(FirUnitIndex);
% group = cell(1,numel(RealtimeStim));
% for i = 1:numel(group)
% group{i} = Unit1Time((Unit1Time>=RealtimeStim(i)- before)&(Unit1Time<RealtimeStim(i)+after));
% end


Unit2Time = time(SecUnitIndex);
% group2 =  cell(1,numel(RealtimeStim));
% for i = 1:numel(group2)
% group2{i} = Unit2Time((Unit2Time>=RealtimeStim(i)- before)&(Unit2Time<=RealtimeStim(i)+after));
% end


try
    Unit3Time = time(ThirdUnitIndex);
% group3 =  cell(1,numel(RealtimeStim));
% for i = 1:numel(group3)
% group3{i} = Unit3Time((Unit3Time>=RealtimeStim(i)- before)&(Unit3Time<=RealtimeStim(i)+after));
% end
catch
    
end


% selectedTrials = 1:numel(group);  %[ 2 3 4 6 8 9 ];
%%
%% temp change color
% clear temp
% temp =group2;
% group2=group;
% group=temp;
%%
% subplot(1,3,1)
% for i = 1:numel(selectedTrials)
%    scatter(group{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),i*ones(1,numel(group{selectedTrials(i)})),1,'d','filled','MarkeredgeColor','b','MarkerfaceColor','r') 
%    hold on; 
%    scatter(group2{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),0.1+i*ones(1,numel(group2{selectedTrials(i)})),1,'o','filled','MarkeredgeColor','g','MarkerfaceColor','g') 
%    hold on; 
%    try
%        scatter(group3{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),0.2+i*ones(1,numel(group3{selectedTrials(i)})),1,'+','filled','MarkeredgeColor','r','MarkerfaceColor','b') 
%    hold on; 
%    catch
%        
%    end
% end
% plot([0 0],[0,i],'g','LineWidth',2)
% axis([0 20 0 5])
% xlabel('t(s)')
% ylabel('Trial(#)')
title(['Channel ' amplifier_channels(ch).custom_channel_name])

% print(h,'-depsc',[num2str(ch) '.eps'])



% for i = 1:numel(selectedTrials)
% group{selectedTrials(i)}=group{selectedTrials(i)}-RealtimeStim(selectedTrials(i));
% group2{selectedTrials(i)}=group2{selectedTrials(i)}-RealtimeStim(selectedTrials(i));
% 
% 
% 
%    try
%  group3{selectedTrials(i)}=group3{selectedTrials(i)}-RealtimeStim(selectedTrials(i));
% 
%    catch
%        
%    end
% end


windowLen = 0.5*20E3;
v= 1/windowLen*20E3*ones(1,windowLen);
subplot(1,2,1)
oF = 100;

tt=0:0.1:time(end);
allNum = Unit1Time;
[counts] = histc(allNum,tt);
% counts=counts/numel(group);
counts=counts*10;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;

counts1=counts/max(counts);

%running average
sample1 = floor(Unit1Time*20E3);
binVersion=zeros(size(1:sample1(end)));
binVersion(sample1)=1;
w = conv(binVersion,v,'same'); 
w = downsample(w,20E3/oF);

optiTime = round(realTime*oF);
sROI = spline(optiTime,ROI',1:optiTime (end))';
if numel(w)>=size(sROI,1) % append to have same length 
sROI(end+1:numel(w),:)=repmat(median(sROI),numel(w)-size(sROI,1),1);
else
    w(end+1:size(sROI,1))=0;
end
% plot(tt,counts,'b')
% text(tt(Index),1.5,num2str(MaxRate))
% hold on;
xcorrList = zeros(size(ROI,2),2);
for roi=1:size(ROI,2)
[r,lags]= xcorr(sROI(:,roi),w,15*oF,'coeff'); % should loop through all columns of r. 
[Y,I] = max(r);
xcorrList(roi,1)=Y;
xcorrList(roi,2)=lags(I)/oF;
end

xcorrList1=xcorrList;
w1=w;
sROI1=sROI;






allNum = Unit2Time;
[counts] = histc(allNum,tt);
% counts=counts/numel(group2);
counts=counts*10;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;
counts2=counts/max(counts);

% plot(tt,counts2+2,'g')
% text(tt(Index),3.5,num2str(MaxRate))
% 
% hold on;
try 
sample1 = floor(Unit2Time*20E3);
binVersion=zeros(size(1:sample1(end)));
binVersion(sample1)=1;
w = conv(binVersion,v,'same'); 
w = downsample(w,20E3/oF);

optiTime = round(realTime*oF);
sROI = spline(optiTime,ROI',1:optiTime (end))';
if numel(w)>=size(sROI,1) % append to have same length 
sROI(end+1:numel(w),:)=repmat(median(sROI),numel(w)-size(sROI,1),1);
else
    w(end+1:size(sROI,1))=0;optiTime
end
% plot(tt,counts,'b')
% text(tt(Index),1.5,num2str(MaxRate))
% hold on;
xcorrList = zeros(size(ROI,2),2);
for roi=1:size(ROI,2)
[r,lags]= xcorr(sROI(:,roi),w,15*oF,'coeff'); % should loop through all columns of r. 
[Y,I] = max(r);
xcorrList(roi,1)=Y;
xcorrList(roi,2)=lags(I)/oF;
end

xcorrList2=xcorrList;
w2=w;
sROI2=sROI;

catch
end



try 
allNum = Unit3Time;
[counts] = histc(allNum,tt);
% counts=counts/numel(group3);
counts=counts*10;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;

counts3=counts/max(counts);

% plot(tt,counts3+4,'r')
% text(tt(Index),5.5,num2str(MaxRate))
sample1 = floor(Unit3Time*20E3);
binVersion=zeros(size(1:sample1(end)));
binVersion(sample1)=1;
w = conv(binVersion,v,'same'); 
w = downsample(w,20E3/oF);

optiTime = round(realTime*oF);
sROI = spline(optiTime,ROI',1:optiTime (end))';
if numel(w)>=size(sROI,1) % append to have same length 
sROI(end+1:numel(w),:)=repmat(median(sROI),numel(w)-size(sROI,1),1);
else
    w(end+1:size(sROI,1))=0;
end
% plot(tt,counts,'b')
% text(tt(Index),1.5,num2str(MaxRate))
% hold on;
xcorrList = zeros(size(ROI,2),2);
for roi=1:size(ROI,2)
[r,lags]= xcorr(sROI(:,roi),w,15*oF,'coeff'); % should loop through all columns of r. 
[Y,I] = max(r);
xcorrList(roi,1)=Y;
xcorrList(roi,2)=lags(I)/oF;
end

xcorrList3=xcorrList;
w3=w;
sROI3=sROI;
catch
    
end
% hold on
% plot([0 0],[0,6],'k','LineWidth',1)
% 
% axis([-30 30 0 6])
% h=figure;
% set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
% for roi =1:size(ROI,2)
% Opti = ROI(:,roi);
% Opti = Opti-min(Opti);
% Opti = Opti/max(Opti);
% 
% plot(realTime,Opti+roi,'k')
% hold on;
% scatter(tt,counts1+roi,3,'filled','MarkerFaceColor','b');
% hold on;
%  plot(tt,counts1+roi,'b--');
% 
% hold on;
%  plot(tt,counts2+roi,'g--');
%  hold on;
% scatter(tt,counts2+roi,3,'filled','MarkerFaceColor','g');
% 
% hold on;
% try
%     
%     hold on;
%  plot(tt,counts3+roi,'r--');
%  hold on;
% scatter(tt,counts3+roi,3,'filled','MarkerFaceColor','r');
% catch
% end
% end

h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
for roi =1:size(ROI,2)
Opti = ROI(:,roi);
Opti = Opti-min(Opti);
% Opti =sqrt(Opti);
Opti = Opti/max(Opti);

plot(realTime,Opti+roi,'k')
% hold on;
% scatter(tt,counts1+roi,3,'filled','MarkerFaceColor','b');
hold on;
 plot((1:numel(w1))/oF,w1/max(w1)+roi,'b');
 hold on;
text(numel(w1)/oF,roi+0.5,[num2str(xcorrList1(roi,1)) '   ' num2str(xcorrList1(roi,2))],'Color','b')
hold on;
try
 plot((1:numel(w2))/oF,w2/max(w2)+roi,'g');
 hold on;
 text(numel(w2)/oF+25,roi+0.5,[num2str(xcorrList2(roi,1)) '   ' num2str(xcorrList2(roi,2))],'Color','g')
catch
end
%  hold on;
% scatter(tt,counts2+roi,3,'filled','MarkerFaceColor','g');

hold on;
try
    
    hold on;
 plot((1:numel(w3))/oF,w3/max(w3)+roi,'r');
 hold on;
  text(numel(w3)/oF+50,roi+0.5,[num2str(xcorrList3(roi,1)) '  ' num2str(xcorrList3(roi,2))],'Color','r')

%  hold on;
% scatter(tt,counts3+roi,3,'filled','MarkerFaceColor','r');
catch
end
end

set(gcf, 'Position', get(0,'Screensize'));

name='0703-0829';
title(['Channel ' amplifier_channels(ch).custom_channel_name])

print(h,'-dpng',[name '-' amplifier_channels(ch).custom_channel_name '.png'])
savefig(h,[name '-' amplifier_channels(ch).custom_channel_name])
set(h,'PaperPosition',[0 0 53 27])
print(h,'-dpdf',[amplifier_channels(ch).custom_channel_name '.pdf'])

















end
close all