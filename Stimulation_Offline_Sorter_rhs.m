%%
name='stimulation20uA22180222145757';
Book1=eval(name);
mkdir(name)
cd(name)
max(data2')

%%

stim_ch=20;  
stim = data2(stim_ch,:);
try 
    fs = frequency_parameters.amplifier_sample_rate;
catch
    fs = 30E3;
end

plot(1/fs:1/fs:numel(stim)/fs,stim)
%%
minThres = 19;
stimOn= find(stim>minThres);
timeTotal = 1/fs:1/fs:numel(stim)/fs;
timeStimOn=timeTotal(stimOn);
plot(diff(timeStimOn))
%%
minTime =0.9;
tselecte = diff(timeStimOn)>minTime;
tselecte(end+1)=true;

RealtimeStim = timeStimOn(tselecte);
stem(RealtimeStim)

%%                     
 after = 0.5; % in seconds, should be bigger than one but smaller than stimulation period. 
 before = 0.5;
channels = unique(Book1(:,1)); %[11 12 16];
for channel = 1:numel(channels)
    ch = channels(channel);
selected = Book1(:,1)== ch;
unitNum = Book1(selected,2);
waveform = Book1(selected,7:46);
time=Book1(selected,3);


FirUnitIndex = find(unitNum==1);
SecUnitIndex = find(unitNum==2);
try 
    ThirdUnitIndex = find(unitNum==3);
catch
end

h=figure;
subplot(1,3,2)
FirUnitWaveform = waveform(FirUnitIndex,:);
SecUnitWaveform = waveform(SecUnitIndex,:);

try 
    ThirdUnitWaveform = waveform(ThirdUnitIndex,:);
catch
end

t = (0:-1+size(FirUnitWaveform,2))/30;

% errorbar (t,mean(FirUnitWaveform),std(FirUnitWaveform),'Color',[.7 .7 .7])
% hold on

%% temp change color
% temp =SecUnitWaveform;
% SecUnitWaveform=FirUnitWaveform;
% FirUnitWaveform=temp;
%%





plot(t,mean(FirUnitWaveform),'b','lineWidth',3)


hold on 
% errorbar (t,mean(SecUnitWaveform),std(SecUnitWaveform),'Color',[.7 .7 .7])
% hold on








plot(t,mean(SecUnitWaveform),'g','lineWidth',3)

try
    hold on 
    plot(t,mean(ThirdUnitWaveform ),'r','lineWidth',3)

catch
end


xlabel('t(ms)')
ylabel('Amplitude(uV)')

% stimWave(:,1)=stimWave(:,1)-stimWave(1,1);
% stimWave(:,1)=stimWave(:,1)*1000;
% plot(stimWave(:,1),stimWave(:,2))
% plot(stimWave(:,1),stimWave(:,2))
% plot(stimWave(:,1),stimWave(:,2),'LineWidth',2)
% 

 % find out t-critical.

 
Unit1Time = time(FirUnitIndex);
group = cell(1,numel(RealtimeStim));
for i = 1:numel(group)
group{i} = Unit1Time((Unit1Time>=RealtimeStim(i)- before)&(Unit1Time<RealtimeStim(i)+after));
end


Unit2Time = time(SecUnitIndex);
group2 =  cell(1,numel(RealtimeStim));
for i = 1:numel(group2)
group2{i} = Unit2Time((Unit2Time>=RealtimeStim(i)- before)&(Unit2Time<=RealtimeStim(i)+after));
end


try
    Unit3Time = time(ThirdUnitIndex);
group3 =  cell(1,numel(RealtimeStim));
for i = 1:numel(group3)
group3{i} = Unit3Time((Unit3Time>=RealtimeStim(i)- before)&(Unit3Time<=RealtimeStim(i)+after));
end
catch
    
end


selectedTrials = 1:numel(group);  %[ 2 3 4 6 8 9 ];
%%
%% temp change color
% clear temp
% temp =group2;
% group2=group;
% group=temp;
%%
subplot(1,3,1)
for i = 1:numel(selectedTrials)
   scatter(group{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),i*ones(1,numel(group{selectedTrials(i)})),1,'d','filled','MarkeredgeColor','b','MarkerfaceColor','r') 
   hold on; 
   scatter(group2{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),0.1+i*ones(1,numel(group2{selectedTrials(i)})),1,'o','filled','MarkeredgeColor','g','MarkerfaceColor','g') 
   hold on; 
   try
       scatter(group3{selectedTrials(i)}-RealtimeStim(selectedTrials(i)),0.2+i*ones(1,numel(group3{selectedTrials(i)})),1,'+','filled','MarkeredgeColor','r','MarkerfaceColor','b') 
   hold on; 
   catch
       
   end
end
plot([0 0],[0,i],'g','LineWidth',2)
% axis([0 20 0 5])
xlabel('t(s)')
ylabel('Trial(#)')
title(['Channel ' amplifier_channels(ch).custom_channel_name])

% print(h,'-depsc',[num2str(ch) '.eps'])



for i = 1:numel(selectedTrials)
group{selectedTrials(i)}=group{selectedTrials(i)}-RealtimeStim(selectedTrials(i));
group2{selectedTrials(i)}=group2{selectedTrials(i)}-RealtimeStim(selectedTrials(i));



   try
 group3{selectedTrials(i)}=group3{selectedTrials(i)}-RealtimeStim(selectedTrials(i));

   catch
       
   end
end



subplot(1,3,3)

tt=-0.5:0.05:0.5;
allNum = cell2mat(group');
[counts] = histc(allNum,tt);
counts=counts/numel(group);
counts=counts*2;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;

counts=counts/max(counts);

plot(tt,counts,'b')
text(tt(Index),1.5,num2str(MaxRate))
hold on;

allNum = cell2mat(group2');
[counts] = histc(allNum,tt);
counts=counts/numel(group2);
counts=counts*2;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;
counts2=counts/max(counts);

plot(tt,counts2+2,'g')
text(tt(Index),3.5,num2str(MaxRate))

hold on;

try 
allNum = cell2mat(group3');
[counts] = histc(allNum,tt);
counts=counts/numel(group3);
counts=counts*2;
[MaxRate, Index] = max(counts);
MaxRate=floor(MaxRate*100)/100;

counts3=counts/max(counts);

plot(tt,counts3+4,'r')
text(tt(Index),5.5,num2str(MaxRate))

catch
    
end
hold on
plot([0 0],[0,6],'k','LineWidth',1)

axis([-0.5 0.5 0 1.5])

title(name(end-9:end-6))
print(h,'-dpng',[name(end-9:end-6) '-' num2str(ch) '.png'])
savefig(h,[name(end-9:end-6) '-' num2str(ch)])
print(h,'-dpdf',[num2str(ch) '.pdf'])

















end
% namelist = dir('*.pdf');
% append_pdfs([name '.pdf'],namelist.name)
%% curve

%%
ISI_in_MS_bins = diff(double(Unit2Time))*1000; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
figure
ax=axes('Units','centimeters','position', [0.6 0.6 5 4.5]);
        set(ax,'FontSize',6)
        if ~isempty(counts)
bar(0.5:1:50.5,counts)
xlim([0 51])
        end   
 
 ax=axes('Units','centimeters','position', [6.5 0.6 5 4.5]);
        set(ax,'FontSize',6)
if ~isempty(counts_M)
bar(10:20:1010,counts_M)
xlim([0 1020])
end
        

    ax=axes('Units','centimeters','position', [12 0.6 5 4.5]);
    set(ax,'FontSize',6)
if ~isempty(counts_L)
bar(0.5:1:20.5,counts_L)
xlim([0 21])
end

