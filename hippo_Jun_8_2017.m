%% Hippocampus Analysis Script

% we want to see several stuff, what stuff are there ? 
% 1. theta oscillation 5-15Hz roughly, not sure. Buzaki paper is
% 6-10Hz,awake
% 2. ripples, @ 110-200 Hz, awake rest or non-rem sleep.
% 3. gamma frequency 30-120Hz awake and alert. 

% sharp wave and ripple are of time scalee 100ms. 
% read_Intan_RHD2000_file
clear all
close all
prompt='Folder name: ';
Folderstr=input(prompt,'s');
prompt='Channel_Map=0109?: 1 for y, 0 for 0524, 2 for 0807, 3 for 0809';
ChMapNum=input(prompt);


path1=pwd;
prompt=([path1([end-4:end-3 end-1:end end-9:end-6]) '?']);
choice=input(prompt,'s');
if strcmp(choice,'1')
str = path1([end-4:end-3 end-1:end end-9:end-6]);
else
str = choice  ;  
end


prompt='CMR?: ';
CMR=input(prompt);
DIR=dir('*.rhd');
%%
for i= 10:12%1:numel(DIR)
read_Intan_RHD2000_file_special(DIR(i).name);

if i==10
data=amplifier_data;
else
    data=[data amplifier_data];
end

% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
end
% Do Channel rejection after ch 47, before ch16
y=[amplifier_channels.native_order];imp=[amplifier_channels.electrode_impedance_magnitude];

selection = (y<48 & y >15)&(imp<3E6);

amplifier_channels=amplifier_channels(selection);
data=data(selection,:);


x(:,2)=[amplifier_channels.native_order]';
x(:,1)=1:size(x,1);
Dat_V_Map=x;
save('Dat_V_Map','Dat_V_Map')

switch ChMapNum
   case 1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 0
   load('ch_map_pink.mat');
   case 2
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_large.mat');
   case 3
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_Small.mat');
    case 4
    load('D:\0919 stroke\Finger_map.mat');
        
end
Ch_Map=Ch_Map_new;

for i=1:size(Ch_Map,1)
for j=1:size(Ch_Map,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j), 1))
        Ch_Map(i,j)=0;
    end
end
end

Ch_Map_2 = Ch_Map; 
for i=1:size(Ch_Map,1)
for j=1:size(Ch_Map,2)
if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
Ch_Map_2(i,j)=0;
else 
Ch_Map_2(i,j)=find(Dat_V_Map(:,2)==Ch_Map(i,j));
end
end
end

amplifier_data=data; 
[b,a]=butter(4,[100 200]/10000,'bandpass');
fprintf(1,'Filtering progress: ')
for i=1:size(amplifier_data,1)
    fprintf(1,'%4.2f',i/size(amplifier_data,1));   
amplifier_data(i,:) = filtfilt(b,a,amplifier_data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')
if CMR ==1
amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);
end


amplifier_data_10=data; 
[b,a]=butter(8,[1 50]/100,'bandpass');
fprintf(1,'Filtering progress: ')
amplifier_data_10=downsample(amplifier_data_10',100);
for i=1:size(amplifier_data_10,1)
    fprintf(1,'%4.2f',i/size(amplifier_data_10,1));   
amplifier_data_10(i,:) = filter(b,a,amplifier_data_10(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')
amplifier_data_10=amplifier_data_10';
clear amplifier_data_1
for ch=1:24
    vq = repmat(amplifier_data_10(ch,:),100,1);
    vq = reshape(vq,1,numel(vq(:)));
amplifier_data_1(ch,:)=vq;
end
amplifier_data_10=amplifier_data_1;
%% Nov202017 Modification
% Target_Sig = sum(amplifier_data.^2);

for ch=1:size(data,1)
    ch
    Target_Sig=amplifier_data(ch,:).^2;
    v= 1/1000*ones(1,1000);
Target_Sig = conv(Target_Sig,v,'same'); 

%     nor_tar_sig  = (Target_Sig-mean(Target_Sig))/std(Target_Sig);%(log(Target_Sig) - mean(log(Target_Sig)))/std(log(Target_Sig));
%     [x,y]=ecdf(nor_tar_sig);
% 
%  plot(y,x)
% hold on
%  xlim([5 15])
%   ylim([0.99 1])
stdThre=5;
thres = mean(Target_Sig )+stdThre*std(Target_Sig);
[pks ,locs] = findpeaks(Target_Sig,'MinPeakHeight',thres ,'MinPeakDistance',50E-3*20E3);
tBefore=0.25*20E3;
tAfter=0.25*20E3;
t=(-5000:5000)/20;
valid = locs>(tBefore) &  locs<(size(amplifier_data,2)-tAfter);
locs=locs(valid);

subsetNum = 36;
if numel(locs)<36
subsetNum = numel(locs);

end
subset = randperm(numel(locs));
locs=locs(subset(1:subsetNum));
figure 
for row=1:6
    for col=1:6
        sample = (row-1)*6+col;
        if sample>numel(locs)
        break
        end
    subplot(6,6,(row-1)*6+col)
    sig_segment = amplifier_data_10(ch,locs(sample)-tBefore:locs(sample)+tAfter); 
    sig_segment_fil = amplifier_data(ch,locs(sample)-tBefore:locs(sample)+tAfter)/std(amplifier_data(ch,:)); 
yyaxis left
    plot(t,3*sig_segment-sig_segment(5000)*3+200);
    ylabel('3X  uV')
    axis([-200 200 -1000 1000])
% axis([125.1 125.75 -1000 1000])
% yyaxis right
hold on 
% [b,a]=butter(4,[100 200]./10E3,'bandpass');
yyaxis right
plot(t,sig_segment_fil-5);
axis([-200 200 -15 15])
ylabel('N std')
hold on



    end
end
title(num2str(numel(pks)))
set(gcf, 'Position', get(0,'Screensize'));
Ele(ch)=numel(pks);
saveas(gcf,num2str(ch),'png')
saveas(gcf,num2str(ch),'fig')

% hist(amplifier_data(ch,:),100);
% hold on; 
% title(num2str(3*std(amplifier_data(ch,:))))
end





%%
std(Target_Sig)
stdThre = 3.75;
thres = exp(mean(log(Target_Sig))+stdThre*std(log(Target_Sig)));
[pks ,locs] = findpeaks(Target_Sig,'MinPeakHeight',thres ,'MinPeakDistance',50E-3*20E3);    
plot(Target_Sig)
hold on;
plot([1 numel(Target_Sig)],[thres  thres],'g')
% for i=1:31
%     figure
%     plot(amplifier_data(i,:))
% end

tBefore=0.25*20E3;
tAfter=0.25*20E3;
valid = locs>(tBefore) &  locs<(size(amplifier_data,2)-tAfter);
locs=locs(valid);
pks=pks(valid);

tBefore_Xcor=0.025*20E3;
tAfter_Xcor=0.025*20E3;


valid = locs>(tBefore) &  locs<(size(data,2)-tAfter);
locs=locs(valid);
%%
mkdir(Folderstr)
cd(Folderstr)

save('stdThre','stdThre')
pdfCellArray=cell(1,size(locs,1));

clear powerNormalised

PowerStore = cell(1,numel(locs));
RawSigStore = cell(1,numel(locs));
Optimal_Lag_Store = cell(numel(locs),1);

Lag_Matrix = cell2mat(Optimal_Lag_Store);
Lag_Ref_Ch = 9;
for sample=1:numel(locs)
    sample
pdfCellArray{sample}=[num2str(sample) '.pdf'];

    Peak_Std = (log(Target_Sig(locs(sample)))- mean(log(Target_Sig)))/std(log(Target_Sig));


    sig_segment = data(:,locs(sample)-tBefore:locs(sample)+tAfter); 
    sig_segment_fil = amplifier_data(:,locs(sample)-tBefore:locs(sample)+tAfter); 

    Sig_for_Xcor = amplifier_data(:,locs(sample)-tBefore_Xcor:locs(sample)+tAfter_Xcor);
    clear Optimal_lag
    for ch = 1:size(Sig_for_Xcor,1)
    
    [r,lags] = xcorr(Sig_for_Xcor(Lag_Ref_Ch,:),Sig_for_Xcor(ch,:));    
    Optimal_lag(ch)=lags(r==max(r));
    end
        Optimal_Lag_Store{sample} = Optimal_lag;
        
        
    PlotFigure = 1;
    if PlotFigure==1
        
        

    fs   = 20E3;
window = 100E-3*fs; % 100ms seconds of data 
overlap = floor(window/10*9);
for ch = 1:size(sig_segment,1)
 signal = sig_segment(ch,:);%-mean(amplifier_data);
[s,f,t] = spectrogram(signal,window,overlap,5:200,20E3);
power=s.*conj(s);
clear  powerNormalised2 powerNormalised1 powerNormalised3
for fre = 1:numel(f)
    
    powerNormalised1(fre,:,ch) =  power(fre,:)./max(power(fre,:));

end
for fre = 1:numel(f)
    
    powerNormalised3(fre,:,ch) =  (power(fre,:)-mean(power(fre,:)))./std(power(fre,:));

end
for tp = 1:numel(t)
    
    powerNormalised2(:,tp,ch) =  power(:,tp)./max(power(:,tp));

end
powerNormalised(:,:,ch)=powerNormalised3(:,:,ch);%+powerNormalised2(:,:,ch);
end
    PowerStore{sample} = powerNormalised;
    RawSigStore{sample}= sig_segment;
    
    
% newPower = PowerStore{1};
% for i=2:numel(locs)
% newPower=newPower+PowerStore{i};
% end
% powerNormalised=newPower;
% powerNormalised=powerNormalised/numel(locs);    
%     
% 
% newSig = RawSigStore{1};
% for i=2:numel(locs)
% newSig=newSig+RawSigStore{i};
% end
% sig_segment=newSig;
% sig_segment=sig_segment/numel(locs);   



h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
set(h, 'Visible', 'off');
Ox=0.1;
Oy=0.2;
gap=0.1;
Gw=1.7;
Gh=4.9;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [Y(row,col) X(row,col)  Gh Gw]);
        set(ax,'FontSize',6)
%         set(ax,'FontWeight','bold')

        if Ch_Map_2(row,col)~=0
            temp=Ch_Map_2(row,col);

% clear t;
%          t=125:1/20E3:125.75;
%          yyaxis left
% plot(t,0.75*amplifier_data(temp,125*20E3:125.75*20E3)+500);
% axis([125.1 125.75 -1000 1000])
% yyaxis right
% [b,a]=butter(4,[100 200]./10E3,'bandpass');
% plot(t,0.75*filter(b,a,amplifier_data(temp,125*20E3:125.75*20E3))-100);
% axis([125.1 125.75 -200 100])

imagesc(powerNormalised(:,:,temp))
% imagesc(mean(powerNormalised,3))
ax = gca;
ax.YDir = 'normal';
% title(['ch' num2str(ch)])
colormap(jet)
% axis([5000 5030 5 200])
% caxis([-3 3])
 caxis([-0.5 max(max(max(powerNormalised(:,:,:))))])
% title( ['Ch' num2str(Ch_Map_new(row,col) )])
        else
            text(0,0,num2str(Peak_Std))
            axis([-1 1 -1 1])
        end
if row~=1 || col~=1
    axis off 
      set(gca,'xticklabel',[])
end

    end
end


% h=figure;
% set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
Ox=13.8;
Oy=0.2;
gap=0.1;
Gw=1.7;
Gh=4.9;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [Y(row,col) X(row,col)  Gh Gw]);
        set(ax,'FontSize',6)
%         set(ax,'FontWeight','bold')

        if Ch_Map_2(row,col)~=0
            temp=Ch_Map_2(row,col);

            
            
            
% clear t;
%          t=125:1/20E3:125.75;
%          yyaxis left
plot(2*sig_segment(temp,:)+500);
% axis([125.1 125.75 -1000 1000])
% yyaxis right
hold on 
% [b,a]=butter(4,[100 200]./10E3,'bandpass');
plot(2*sig_segment_fil(temp,:)-100);
hold on
axis([0 size(sig_segment,2) -1000 1000])
plot([size(sig_segment,2) size(sig_segment,2)],[0 2000],'lineWidth',2)





% 
% histogram(Lag_Matrix(:,temp),'BinMethod','integers')
% axis([-10 10 0 189])
% 
% 




% axis([125.1 125.75 -200 100])

% imagesc(powerNormalised(:,:,temp))
% ax = gca;
% ax.YDir = 'normal';
% title(['ch' num2str(ch)])
% colormap(jet)
% axis([5000 5030 5 200])

% title( ['Ch' num2str(Ch_Map_new(row,col) )])
        end
% if row~=1 || col~=1
    axis off 
      set(gca,'xticklabel',[])
% end

    end
end


set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')




set(gcf, 'Position', get(0,'Screensize'));
set(h,'PaperPosition',[0 0 53 27])
% print(h,'-dpdf',[num2str(sample) '.pdf'])

    end
end
close all
% namelist = dir('*.pdf');
% append_pdfs([str '.pdf'],pdfCellArray{:})
% 
%% Dr.Xie Cross overlay plots of identified neurons.   
%  
color=jet(120);
color=color(5:104,:);
firstRow=[0 0 0];
color=[firstRow;color];

%listUnit = [2 3 4 5 7 12]; % for 0606
listUnit = [19 37 34 9 14 4]; % for 0616

Lia = ismember(list,listUnit);
Cselected5=WCList(Lia,:);
Ox=3;
Oy=9;
gap=0.3;
Gw=1.5;
Gh=1.5;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);
% Weighted Center Calculation
% WC_Y=sum(sum(Weighted_P2P_Remap.*(Y+1)))/sum(sum(Weighted_P2P_Remap));
% WC_X=sum(sum(Weighted_P2P_Remap.*(X+1)))/sum(sum(Weighted_P2P_Remap));
% MUSIC Location Calculation
figure
for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [Y(row,col) X(row,col) Gw Gh]);
        set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
        
        set(ax,'Box','on')

        set(ax,'LineWidth',2)
%         plot(0,0)
%         T = text(-0.2,0,num2str(floor(abs(P2P_Remap(row,col)))));
%         set(T,'FontSize',6)
%         
%         set(T,'Color',[0 0 1])
        
        set(ax,'xtick',[])
        set(ax,'ytick',[])
        set(ax,'XColor',[0.7 0.7 0.7])
        set(ax,'YColor',[0.7 0.7 0.7])
        
%         if SelectedLocation(row,col)==1
%         set(ax,'XColor',[1 0 0])
%         set(ax,'YColor',[1 0 0])
%         end
        
        end
end



for center= 1:size(Cselected5,1)
ax=axes('Units','centimeters','position', [Cselected5(center,2)/25*1.5-0.5+Oy Cselected5(center,1)/25*1.5+Ox-0.5 1 1]);
ZSize = 60;
RDistance = abs(Cselected5(center,3));
RDistance(RDistance>60)=60;
NorDistance = RDistance/60*100; % full color specified for 2 spacing. 30um
scatter(0,0,30+(ZSize)^2,'+','LineWidth',2,'MarkerEdgeColor',color(101-floor(NorDistance),:))
% fx=floor((WC_X-Ox)/1.5*25*10)/10;
% fy=floor((WC_Y-Oy)/1.5*25*10)/10;
% txt1 = ['       ' num2str(floor(Cselected5(center,3))) ];
% T=text(0,0,txt1);
% set(T,'FontWeight','bold')
% set(T,'FontSize',5)
set(ax,'Visible','off')
end


%% New Code For Lag Plots.


pdfCellArray=cell(1,size(locs,1));

clear powerNormalised

locs1=locs;  locs =randi([1 numel(Target_Sig)],1, numel(locs1));

 
 PowerStore = cell(1,numel(locs));
Optimal_Lag_Store = cell(numel(locs),1);
ReferenceCh = 8;
for sample=1:numel(locs)
    sample/numel(locs)

pdfCellArray{sample}=[num2str(sample) '.pdf'];

    Peak_Std = (log(Target_Sig(locs(sample)))- mean(log(Target_Sig)))/std(log(Target_Sig));


    sig_segment = data(:,locs(sample)-tBefore:locs(sample)+tAfter); 
    sig_segment_fil = amplifier_data(:,locs(sample)-tBefore:locs(sample)+tAfter); 

    Sig_for_Xcor = amplifier_data(:,locs(sample)-tBefore_Xcor:locs(sample)+tAfter_Xcor);
    clear Optimal_lag
    for ch = 1:size(Sig_for_Xcor,1)
    
    [r,lags] = xcorr(Sig_for_Xcor(ReferenceCh,:),Sig_for_Xcor(ch,:));    
    Optimal_lag(ch)=lags(r==max(r));
    end
        Optimal_Lag_Store{sample} = Optimal_lag;
        
        
    PlotFigure = 0;

end
%%
close all
Lag_Matrix = cell2mat(Optimal_Lag_Store);


h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')

Ox=10.8;
Oy=0.2;
gap=0.1;
Gw=1.7;
Gh=4.9;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [Y(row,col) X(row,col)  Gh Gw]);
        set(ax,'FontSize',6)
%         set(ax,'FontWeight','bold')

        if Ch_Map_2(row,col)~=0
            temp=Ch_Map_2(row,col);

            
    
 histogram(Lag_Matrix(:,temp),'BinMethod','integers')
 axis([-10 10 0 189])



% title( ['Ch' num2str(Ch_Map_new(row,col) )])
        end
% if row~=1 || col~=1

      set(gca,'yticklabel',[])
% end

    end
end