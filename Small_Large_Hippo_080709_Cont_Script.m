load('D:\0517-hippo\ch_map_pink.mat')
switch ChMapNum
   case 1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 0
   load('D:\0517-hippo\ch_map_pink.mat');
   case 2
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_large.mat');
   case 3
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_Small.mat');
        
end
%%
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
    
    


h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
set(h, 'Visible', 'off');
Ox=0.1;
Oy=0.2;
gap=0.1;
Gw=1.5;
Gh=4.9;
numGx = 16;
numGy = 2;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:2
    for col=1:16

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
Ox=0.1;
Oy=10.5;
gap=0.1;
Gw=1.5;
Gh=4.9;
numGx = 16;
numGy = 2;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:2
    for col=1:16

        ax=axes('Units','centimeters','position', [Y(row,col) X(row,col)  Gh Gw]);
        set(ax,'FontSize',6)
%         set(ax,'FontWeight','bold')

        if Ch_Map_2(row,col)~=0
            temp=Ch_Map_2(row,col);

            
            
            
% clear t;
%          t=125:1/20E3:125.75;
%          yyaxis left
plot(10*sig_segment(temp,:)+500);
% axis([125.1 125.75 -1000 1000])
% yyaxis right
hold on 
% [b,a]=butter(4,[100 200]./10E3,'bandpass');
plot(20*sig_segment_fil(temp,:)-100);
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
print(h,'-dpdf',[num2str(sample) '.pdf'])

    end
end
% close all
namelist = dir('*.pdf');
append_pdfs([str '.pdf'],pdfCellArray{:})


newPower = PowerStore{1};
for i=2:numel(locs)
newPower=newPower+PowerStore{i};
end
powerNormalised=newPower;
powerNormalised=powerNormalised/numel(locs);    
    

newSig = RawSigStore{1};
for i=2:numel(locs)
newSig=newSig+RawSigStore{i};
end
sig_segment=newSig;
sig_segment=sig_segment/numel(locs);   
save('imp_data','powerNormalised','sig_segment');
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
ReferenceCh = 13;
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

Ox=0;
Oy=0.2;
gap=0.1;
Gw=1.5;
Gh=4.5;
numGx = 16;
numGy = 2;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

for row=1:2
    for col=1:16

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