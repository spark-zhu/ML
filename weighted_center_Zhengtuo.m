% ch_num= length(amplifier_channels);
% figure(1);clf; co=get(gca,'colororder');  
% offset= max(max(amplifier_data))*0.5;
% 
% d = fdesign.highpass('Fst,Fp,Ast,Ap',0.015,0.025,60,1);
%    Hd = design(d,'equiripple');
%    
% %    for i=1:ch_num
% i=19;
%        datafilt(i,:)  = filter(Hd,amplifier_data(i,:));
%        figure;
%        plot(t_amplifier(1:1447200), datafilt(i,1:1447200)+(i-1)*300,'-','color',co(mod(i,7)+1,:)); hold on; 
%        text(max(t_amplifier)+2,(i-1)*offset,sprintf(amplifier_channels(i).native_channel_name),'color',co(mod(i,7)+1,:));
% %    end
%    
clear 
clc
DatFile='wrapping_right.dat';

%% 
nChansInRawFile=26;
nChansInDatFile=26;
tBefore=20;
tAfter=20;

% MyTimes=h5read('awake_8-21-16.kwik','/channel_groups/0/spikes/time_samples');
% MyTimes=h5read('hybrid_10sec.kwik','/channel_groups/0/spikes/time_samples');

FileInfo = dir(DatFile);
Source = memmapfile(DatFile, 'Format', {'int16', [nChansInRawFile, (FileInfo.bytes/nChansInDatFile/2)], 'x'});
% FullSpikes = zeros(nChansInDatFile, nTimeSamplesPerWaveform, nSpikes);
% for i=1:nSpikes
%      FullSpikes(:,:,i) = (Source.Data.x(1:nChansInDatFile,MyTimes(i)-tBefore:MyTimes(i)+tAfter),1);
% end

cluster=h5read('awake_8-21-16.kwik','/channel_groups/0/spikes/clusters/main');
% cluster=h5read('hybrid_10sec.kwik','/channel_groups/0/spikes/clusters/main');
load('Ch_Map_and_DatVMap.mat')
MyTimes(end)=[];
cluster(end)=[];
list=unique(cluster);
for clu=1:numel(list)
    
p=find(cluster==list(clu));
%%%%%%%%%%%%%%%%%%%%%%%%%%  extract waveforms
clear unit MeanPlot 
for k=1:length(p);
    unit(:,:,k)=Source.Data.x(1:nChansInDatFile,MyTimes(p(k))-tBefore:MyTimes(p(k))+tAfter);
%    plot(unit(5,:));% THIS SEEMS LIKE STRANGE, WHY GET DATA FROM ALL
%    CHANNELS. 
%    hold on
end
MeanPlot=mean(unit,3);
%% 4-by-8 Channel Plot showing correlation by visual inspection
for i=1:26
MeanPlot(i,:)=MeanPlot(i,:)-MeanPlot(i,1);
end
h=figure;
t=[-tBefore/20e3:1/20e3:tAfter/20e3];
for row=1:4
    for col=1:8
    subplot(4,8,(row-1)*8+col)
    if Ch_Map(row,col)~=0
    temp=find(Dat_V_Map(:,2)==Ch_Map(row,col));
    plot(t,MeanPlot(temp,:),'LineWidth',2)
    axis([t(1) t(end) min(min(MeanPlot))-10 max(max(MeanPlot))+10])
    end
    end
end
print(h,'-depsc',[num2str(list(clu)) '.eps'])
end
%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  plot single cluster for everyone channel
for m=1:length(unit(:,1,1))
for l=1:length(unit(1,:,1))
mean_spikes(m,l)=mean(unit(m,l,:));     %%%%%%%%%%%% channel #
end
plot(mean_spikes(m,:));      %SO AGAIN THE POINT IS WHY IT'S REVELANT TO
% PLOT ALL CHANNELS, NOT ALL ARE CAPTURING THE CURRENT SPIKE RIGHT.
% unit 18*81*7742   (7742)they are actually from the same unit. 
% all samples are averaged in time, not in channel. 
hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files=dir('unit*');
% clust=[4,5,9,12,19,26,28,29]; 8-24-16    %% WHAT ARE THOSE
% clust=[14,17,26,28,30,40]; 8-21-16
% clust=[2,4,5,6,7,8,13,14,16,17,20,21,26,27,28,29,30,31,32]
for n=1:length(files);                       %%%%%%%%%%%  n=spike #
    unit=load(files(n).name);
    unit=unit.unit;
    for m=1:length(unit(:,1,1))     %%%%%%%%%%%% m = channel #
     for l=1:length(unit(1,:,1));
      mean_spikes(m,l,n)= mean(unit(m,l,:));     %%%%%%%%%%%% channel #
     end 
% plot(mean_spikes(m,:));
% hold on
end
end
%%%%%%%%%%%%%%%%%%%%%%%% aligh the waveforms  SO THAT THEY ALL START AT
%%%%%%%%%%%%%%%%%%%%%%%% ZERO
for n=1:length(files);   
    for m=1:18;
      mean_spikes_alighed(m,:,n)=mean_spikes(m,:,n)-mean_spikes(m,1,n);
end    
end



%% 
figure 
for n=1:length(files);  % why channel 16. 
plot(mean_spikes_alighed(16,:,n));
hold on
end

figure     % plotting from one channel, the 7 units selected.
plot(mean_spikes_alighed(17,:,1),'color',[0 1 0.74117648601532]); %dark green
hold on

plot(mean_spikes_alighed(17,:,2),'color',[0 0.447058826684952 0.74117648601532]);%light blue
hold on

plot(mean_spikes_alighed(17,:,7),'color',[0.929411768913269 0.694117665290833 0.125490203499794]); %yello
hold on

plot(mean_spikes_alighed(17,:,3),'color',[0.46 0.674 0.188]);%light green
hold on

plot(mean_spikes_alighed(17,:,4),'color',[0 0 0.68]);%purple
hold on

plot(mean_spikes_alighed(17,:,6),'color',[0.85 0.35 0.1]);%orange
hold on

plot(mean_spikes_alighed(17,:,5),'color',[0 1 1]);%cyan
%%

for n=1:length(files);
for m=11:17;    
p2p(n,m)=max(mean_spikes_alighed(m,:,n))-min(mean_spikes_alighed(m,:,n));
if p2p(n,m)<35  % select a threshold for peak to peak. 
    p2p(n,m)=0;
end
end
end

figure 
for m=1:18; % plot at the 2nd unit, all channels. 
plot(mean_spikes_alighed(m,:,2));
hold on
end

%figure
%plot(mean_spikes_alighed(15,:));

s1=mean_spikes(10,:,1);
s2=mean_spikes(10,:,2);
% s1 and s2 never used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%locate the neurons 

%%%%%8-24
c10=[0,90,0];c13=[30,90,0];c17=[60,90,0];c5=[90,90,0];c9=[0,60,0];c12=[30,60,0];c16=[60,60,0];c6=[90,60,0];c2=[120,60,0];c8=[0,30,0];c11=[30,30,0];
c15=[60,30,0];c7=[90,30,0];c3=[120,30,0];c14=[30,0,0];c18=[60,0,0];c4=[90,0,0];


%%%%%%%%%%%%%%%for 6 electrodes,11,12,13,15,16,17
for ii=1:length(files);      % ii represents unit,center of gravity of 6 electrodes. 
xy(ii,:)=(c11*p2p(ii,11)+c12*p2p(ii,12)+c13*p2p(ii,13)+c15*p2p(ii,15)+c16*p2p(ii,16)+c17*p2p(ii,17))/(p2p(ii,17)+p2p(ii,16)+p2p(ii,15)+p2p(ii,13)+p2p(ii,12)+p2p(ii,11));

end

x=xy(:,1);
y=xy(:,2);

     electrode=[1,3;1,2;1,1;2,3;2,2;2,1]
     electrode=electrode.*30;
xx=electrode(:,1);
yy=electrode(:,2);

figure
plot(xx,yy,'MarkerSize',60,'Marker','square','LineStyle','none')
xlim([-15,114]);
ylim([-30,114]);
hold on

axis off 
axis equal
%ploting electrodes finsihed. Start plotting neuron positions relative to
%electrodes.

           plot(x(1),y(1),'o','markersize',10,'MarkerFaceColor',[0 1 0.74117648601532],'MarkerEdgeColor',[0 1 0.74117648601532]);
            hold on
            
            plot(x(2),y(2),'o','markersize',10,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532]);
            hold on
             plot(x(7),y(7),'o','markersize',10,'MarkerFaceColor',[0.929411768913269 0.694117665290833 0.125490203499794],'MarkerEdgeColor',[0.929411768913269 0.694117665290833 0.125490203499794]);
            hold on
            plot(x(3),y(3),'o','markersize',10,'MarkerFaceColor',[0.46 0.674 0.188],'MarkerEdgeColor',[0.46 0.674 0.188]);
            hold on
          plot(x(4),y(4),'o','markersize',10,'MarkerFaceColor',[0 0 0.68],'MarkerEdgeColor',[0 0 0.74]);
hold on

            plot(x(5),y(5),'o','markersize',10,'MarkerFaceColor',[0.85 0.35 0.1],'MarkerEdgeColor',[0.85 0.35 0.1]);
           hold on
           
              plot(x(6),y(6),'o','markersize',10,'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1]);


% x=[3.6,0.4,2.58,3.27,3.37,2.65,1.3,2.12,0.9,0.68,0.32,2.42];
% y=[1.24,1.5,2,2.06,1.95,2.12,1.3,1.79,2.6,2.75,1.22,1.84];
% figure
% for q=1:12
%     plot(x(q),y(q),'o','markersize',10)
%     hold on
% end

%%%%%8-25
a1=(c17*4+c16*5+c6*4+c15*4)/(4+5+4+4);
a2=(c10*5+c13*5.5+c17*4.5+c9*3.5+c12*4)/(5+5.5+4.5+3.5+4);
a3=(c17*6+c5*5.5+c16*6+c6*5+c15*3+c7*3)/(6+5.5+6+5+3+3);
a4=(c5*4+c6*4+c2*5+c3*3)/(4+4+5+3);
a5=(c17*6+c5*9+c16*5+c6*7+c2*5.5+c15*3+c7*3+c3*2)/(6+9+5+7+5.5+3+3+2);
a6=(c10*5+c13*3.5+c9*4+c12*4+c8*4+c11*3)/(5+3.5+4+4+4+3);
a7=(c12*5+c16*5+c11*4+c15*5)/(5+5+4+5);
a8=(c13*6.5+c17*12+c5*10+c12*6+c16*10.5+c6*9.5+c15*6+c7*5.5+c14*2.5+c18*3+c4*2.5)/(6.5+12+10+6+10.5+9.5+6+5.5+2.5+3+2.5);
a9=(c10*3.5+c13*4+c17*3.5+c9*2.5+c12*3+c16*3)/(3.5+4+3.5+2.5+3+3);
a10=(c17*3+c5*4.5+c16*3.5+c6*5+c2*5.5+c7*4+c3*5)/(3+4.5+3.5+5+5.5+4+5);

x=[a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1),a9(1),a10(1)];
y=[a1(2),a2(2),a3(2),a4(2),a5(2),a6(2),a7(2),a8(2),a9(2),a10(2)];
figure
for q=1:10
    plot(x(q),y(q),'o','markersize',10)
    hold on
end
%%%%%8-27
a1=(c17*13+c5*16+c16*13+c6*19+c2*21+c15*10+c7*11+c3*17)/(13+16+13+19+21+10+11+17+9+11);
a2=(c13*12+c17*13+c9*12+c12*20+c16*19+c8*11+c11*18+c15*17)/(12+13+12+20+19+11+18+17);
a3=(c13*9+c17*11+c5*10+c16*12+c6*11.5+c15*8)/(9+11+10+12+11.5+8);
a4=(c17*6+c5*10+c16*7+c6*10+c2*12+c3*7)/(6+10+7+10+12+7);
a5=(c10*12+c13*13+c17*12+c5*11+c9*9+c12*11+c16*11+c6*11)/(12+13+12+11+9+11+11+11);
a6=(c17*17+c5*28+c16*16+c6*23+c2*17+c15*8+c7*9+c3*11)/(17+28+16+23+17+8+9+11);
a7=(c17*4+c5*3+c12*10+c16*16+c6*12+c11*7+c15*13)/(4+3+10+16+12+7+13);
a8=(c10*14+c13*10+c9*12+c12*10+c8*10)/(14+10+12+10+10);
a9=(c13*6.5+c17*12+c5*10+c12*6+c16*10.5+c6*9.5+c15*6+c7*5.5+c14*2.5+c18*3+c4*2.5)/(6.5+12+10+6+10.5+9.5+6+5.5+2.5+3+2.5);;
a10=(c10*10+c13*12+c17*8+c9*8+c12*9+c16*8)/(12+10+8+8+9+8);
a11=(c9*3.5+c8*7+c11*5)/(3.5+7+5);

x=[a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1),a9(1),a10(1)];
y=[a1(2),a2(2),a3(2),a4(2),a5(2),a6(2),a7(2),a8(2),a9(2),a10(2)];
figure
for q=1:10
    plot(x(q),y(q),'o','markersize',10)
    hold on
end
%%%%%%%%%%%%%%%%%%%%plot neuron locations




