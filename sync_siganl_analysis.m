%% sync_siganl_analysis 
pulse =  RedLEDArea;%  selected pulse variable. % ADCdata;%GreenLEDArea;%
close all 
figure
plot(pulse)
DifPulse  = diff(pulse); 
figure
plot(diff(pulse))
thres = 450; % input desired threshold.
findpeaks(DifPulse,'minPeakHeight',thres,'MinPeakDistance',100*2);
clear locs pks
[pks locs] = findpeaks(DifPulse,'minPeakHeight',thres,'MinPeakDistance',100*2);
locs = locs+1;%(compensate for diff function loss)
% locs(87)=[]
figure
plot(diff(locs))
%locs(find(diff(locs)<1537)+1)=[];
%%
Intan.locs=locs;
Green.locs = locs;
Red.locs = locs;
%%
Green.locs(:,2)=Intan.locs;
Red.locs(:,2)=Intan.locs;
Green.All(:,1)=Green.locs(1,1):Green.locs(end,1);
vq = interp1(Green.locs(:,1),Intan.locs,Green.All);
Green.All(:,2)=vq;

Red.All(:,1)=Red.locs(1,1):Red.locs(end,1);
vq = interp1(Red.locs(:,1),Intan.locs,Red.All);
Red.All(:,2)=vq;
%%
windowLen = 0.05*20E3; % window length in sample for Firing rate calculation
v= 1/windowLen*20E3*ones(1,windowLen);
 oF = 100; % downsample fold. 
ignoredUnitNumber = []; % noise unit eliminated. 
%        load('D:\0517-hippo\ch_map_pink.mat');

% Ch_Map= Ch_Map_new-15;
% 
figure
count = 0 ;
FRstore = cell(1,1);
binVersion=zeros(size(1:max(cell2mat(result.time'))));
for unit=1:numel(result.selectedClusters)
    unit
sample1 = result.time{result.list==result.selectedClusters(unit)};
binVersion(:)=0;
binVersion(sample1)=1;
w = conv(binVersion,v,'same'); 
w = downsample(w,oF);
FRstore{unit}=w;
end
ttFRsample = downsample(0:max(cell2mat(result.time')),oF) ;

%%
% GreenCenter(147184:end,:)=[];
% RedCenter(147442:end,:)=[];
while(sum(GreenCenter(:,2)==0)>0||sum(GreenCenter(:,1)==0)>0)
for i=2:numel(GreenCenter(:,2))-1
    i
    if GreenCenter(i,2)==0
    GreenCenter(i,2)=(GreenCenter(i-1,2)+GreenCenter(i+1,2))/2;
    end
    
    if GreenCenter(i,1)==0
    GreenCenter(i,1)=(GreenCenter(i-1,1)+GreenCenter(i+1,1))/2;
    end
    
end
end

while(sum(RedCenter(:,2)==0)>0||sum(RedCenter(:,1)==0)>0)
for i=2:numel(RedCenter(:,2))-1
    if RedCenter(i,2)==0
    RedCenter(i,2)=(RedCenter(i-1,2)+RedCenter(i+1,2))/2;
    end
    
    if RedCenter(i,1)==0
    RedCenter(i,1)=(RedCenter(i-1,1)+RedCenter(i+1,1))/2;
    end
    
end
end
% uisave
% while(sum(isnan(GreenCenter(:,2)))>0||sum(isnan(GreenCenter(:,1)))>0)
% for i=2:numel(GreenCenter(:,2))-1
%     if isnan(GreenCenter(i,2))
%     GreenCenter(i,2)=(GreenCenter(i-1,2)+GreenCenter(i+1,2))/2;
%     end
%     
%     if isnan(GreenCenter(i,1))
%     GreenCenter(i,1)=(GreenCenter(i-1,1)+GreenCenter(i+1,1))/2;
%     end
%     
% end
% 
%  if isnan(GreenCenter(1,2))
%     GreenCenter(1,2)=GreenCenter(2,2);
%  end
%     
%  if isnan(GreenCenter(end,2))
%     GreenCenter(end,2)=GreenCenter(end-1,2);
%  end
%     
% if isnan(GreenCenter(1,1))
%     GreenCenter(1,1)=GreenCenter(2,1);
%  end
%     
%  if isnan(GreenCenter(end,1))
%     GreenCenter(end,1)=GreenCenter(end-1,1);
%     end
%  
%     
% 
% 
% end


%%
figure

% Kinetics = diff(GreenCenter(Green.All(:,1),2));
K1 = GreenCenter(Green.All(:,1),2); % Y
K2 = GreenCenter(Green.All(:,1),1); % X
Kinetics =K2;%sqrt(diff([K1(1);K1]).^2+diff([K2(1);K2]).^2);% 2*abs(diff([K2(1);K2]));%

for unit=1:numel(FRstore)
    offset = 200;
    
plot(Green.All(:,2)/20E3, Kinetics+offset*(unit-1),'b')
hold on;
plot(ttFRsample/20E3,FRstore{unit}+offset*(unit-1),'r')
hold on;
Tracker = interp1(Green.All(:,2)/20E3,Kinetics,ttFRsample/20E3,'linear' ,0);

TrackerSpeed = abs(interp1(Green.All(:,2)/20E3,diff([Kinetics(1); Kinetics]),ttFRsample/20E3,'linear' ,0));

TrackerAcc = abs(interp1(Green.All(:,2)/20E3,diff(diff([Kinetics(1); Kinetics(1); Kinetics])),ttFRsample/20E3,'linear' ,0));

CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(Tracker))',Tracker(~isnan(Tracker))');

% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerSpeed))',TrackerSpeed(~isnan(TrackerSpeed))');
% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerAcc))',TrackerAcc(~isnan(TrackerAcc))');
% 

text(5600,offset*(unit-1)+100,[num2str(CorrelationUnit(unit)*100) '%'])

end

%%
figure

% Kinetics = diff(GreenCenter(Green.All(:,1),2));
K1 = RedCenter(Red.All(:,1),2);
K2 = RedCenter(Red.All(:,1),1);
Kinetics =K2;%2*abs(diff([K2(1);K2]));%sqrt(diff([K1(1);K1]).^2+diff([K2(1);K2]).^2);

for unit=1:numel(FRstore)
    offset = 200;
    
h=plot(Red.All(:,2)/20E3, Kinetics+offset*(unit-1),'b');
% alpha(h,0.5)
hold on;
plot(ttFRsample/20E3,FRstore{unit}+offset*(unit-1),'r')
hold on;
Tracker = interp1(Red.All(:,2)/20E3,Kinetics,ttFRsample/20E3,'linear' ,0);

TrackerSpeed = abs(interp1(Red.All(:,2)/20E3,diff([Kinetics(1); Kinetics]),ttFRsample/20E3,'linear' ,0));

TrackerAcc = abs(interp1(Red.All(:,2)/20E3,diff(diff([Kinetics(1); Kinetics(1); Kinetics])),ttFRsample/20E3,'linear' ,0));

CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(Tracker))',Tracker(~isnan(Tracker))');

% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerSpeed))',TrackerSpeed(~isnan(TrackerSpeed))');
% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerAcc))',TrackerAcc(~isnan(TrackerAcc))');
% 

text(4200,offset*(unit-1)+100,[num2str(CorrelationUnit(unit)*100) '%'])

end
%%
figure


K1 = RedCenter(Red.All(:,1),2);
K2 = RedCenter(Red.All(:,1),1);
Kinetics = 2*abs(diff([K2(1);K2]));%K2;%sqrt(diff([K1(1);K1]).^2+diff([K2(1);K2]).^2);
% Kinetics = abs(diff([K2(1);K2]));



% windowLen = 5 % window length in sample for Firing rate calculation
% v= 1/windowLen*20E3*ones(1,windowLen);
% w = conv(binVersion,v,'same'); 
% 

clear CorrelationUnit

for unit=1:numel(FRstore)
    offset = 200;
    
plot(Red.All(:,2)/20E3, Kinetics+offset*(unit-1),'b')
hold on;
plot(ttFRsample/20E3,FRstore{unit}+offset*(unit-1),'r')
hold on;
Tracker = interp1(Red.All(:,2)/20E3,Kinetics,ttFRsample/20E3,'linear' ,0);

TrackerSpeed = abs(interp1(Red.All(:,2)/20E3,diff([Kinetics(1); Kinetics]),ttFRsample/20E3,'linear' ,0));

TrackerAcc = abs(interp1(Red.All(:,2)/20E3,diff(diff([Kinetics(1); Kinetics(1); Kinetics])),ttFRsample/20E3,'linear' ,0));

CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(Tracker))',Tracker(~isnan(Tracker))');

% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerSpeed))',TrackerSpeed(~isnan(TrackerSpeed))');
% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerAcc))',TrackerAcc(~isnan(TrackerAcc))');
% 

text(5600,offset*(unit-1)+100,[num2str(CorrelationUnit(unit)*100) '%'])

end
%% X/Y Location Tunning Curve Plot 
figure

% Kinetics = diff(GreenCenter(Green.All(:,1),2));
K1 = GreenCenter(Green.All(:,1),2);
K2 = GreenCenter(Green.All(:,1),1);
Kinetics =K1;%sqrt(diff([K1(1);K1]).^2+diff([K2(1);K2]).^2);

for unit=1:numel(FRstore)
     offset = 50;
%     
% plot(Green.All(:,2)/20E3, Kinetics+offset*(unit-1),'b')
% hold on;
% plot(ttFRsample/20E3,FRstore{unit}+offset*(unit-1),'r')
% hold on;
 Tracker = interp1(Green.All(:,2)/20E3,Kinetics,ttFRsample/20E3,'linear' ,0);
 
 Tracker = round(Tracker);
 clear FRperTrack
 for bin = min(Tracker):1:max(Tracker)
 FRperTrack(bin-min(Tracker)+1)=mean(FRstore{unit}(Tracker==bin));
 end
 
 
 plot(min(Tracker):1:max(Tracker),FRperTrack+offset*(unit-1))
 hold on
% TrackerSpeed = abs(interp1(Green.All(:,2)/20E3,diff([Kinetics(1); Kinetics]),ttFRsample/20E3,'linear' ,0));
% 
% TrackerAcc = abs(interp1(Green.All(:,2)/20E3,diff(diff([Kinetics(1); Kinetics(1); Kinetics])),ttFRsample/20E3,'linear' ,0));
% 
% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(Tracker))',Tracker(~isnan(Tracker))');
% 
% % CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerSpeed))',TrackerSpeed(~isnan(TrackerSpeed))');
% % CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerAcc))',TrackerAcc(~isnan(TrackerAcc))');
% % 
% 
% text(1500,offset*(unit-1)+100,[num2str(CorrelationUnit(unit)*100) '%'])
% 
end


%%


K1 = RedCenter(Red.All(:,1),2);
K2 = RedCenter(Red.All(:,1),1);
Kinetics =sqrt(diff([K1(1);K1]).^2+diff([K2(1);K2]).^2);
Kinetics = K2;%abs(diff([K2(1);K2]));

% windowLen = 5 % window length in sample for Firing rate calculation
% v= 1/windowLen*20E3*ones(1,windowLen);
% w = conv(binVersion,v,'same'); 
%
clear CorrelationUnit

for unit=1:numel(FRstore)
 
    

Tracker = interp1(Red.All(:,2)/20E3,Kinetics,ttFRsample/20E3,'linear' ,0);

TrackerSpeed = abs(interp1(Red.All(:,2)/20E3,diff([Kinetics(1); Kinetics]),ttFRsample/20E3,'linear' ,0));

TrackerAcc = abs(interp1(Red.All(:,2)/20E3,diff(diff([Kinetics(1); Kinetics(1); Kinetics])),ttFRsample/20E3,'linear' ,0));

CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(Tracker))',Tracker(~isnan(Tracker))');

% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerSpeed))',TrackerSpeed(~isnan(TrackerSpeed))');
% CorrelationUnit(unit)=corr(FRstore{unit}(~isnan(TrackerAcc))',TrackerAcc(~isnan(TrackerAcc))');
% 


end

h=figure;
   Ox=15;    
Oy=15;
gap=20;
        Gw=30;
        Gh=30;        
    
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        
      for r=1:4
        for c=1:8
%         scatter(X(r,c),Y(r,c),22000,[0 0 0],'s','LineWidth',2)
        pAx = X(r,c)-12.5;
        pAy = Y(r,c)-12.5;
        pCx = X(r,c)+12.5;
        pCy = Y(r,c)+12.5;
        plot([pAx pCx],[pAy pAy],'k','LineWidth',1)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',1)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',1)
        end
      end
    
      hold on
        
        scaleCoeff = 50/50;

clear  Cselected5

[aAndB aInds bInds] = d_uIntersect(result.list,result.selectedClusters);
Cselected5=cell2mat(result.Location);
Cselected5=Cselected5(aInds,:);%noise rejection
% Cselected5=Cselected5(1:unit,:);
rgb=redblue(40);
for center= 1:size(Cselected5,1)
% ax=axes('Units','centimeters','position', [Cselected5(center,1)*scaleCoeff-0.5+Ox Cselected5(center,2)*scaleCoeff+Oy-0.5 1 1]);
ZSize = 60;
% RDistance = sqrt(sum((Cselected5(center,:)-[highX highY 0]).^2));
% RDistance(RDistance>60)=60;
% NorDistance = RDistance/60*100; % full color specified for 2 spacing. 30um
% if center~=size(Cselected5,1)
% scatter(Cselected5(center,1),Cselected5(center,2),2000,'.','LineWidth',2,'MarkerEdgeColor',rgb(floor(CorrelationUnit(center)*100+20),:))
scatter(Cselected5(center,1),Cselected5(center,2),500,'filled','LineWidth',1,'MarkerFaceColor',rgb(min(floor(CorrelationUnit(center)*100+20),size(rgb,1)),:),'MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',1)

% else
% scatter(0,0,300,'*','LineWidth',2,'MarkerEdgeColor',rgb(unit,:))    
% end
% fx=floor((WC_X-Ox)*scaleCoeff*10)/10;
% fy=floor((WC_Y-Oy)*scaleCoeff*10)/10;
% 

% txt1 = ['       ' num2str(floor(Cselected5(center,3))) ];
% T=text(0,0,txt1);
% set(T,'FontWeight','bold')
% set(T,'FontSize',5)

% set(ax,'Visible','off')
hold on;
end
axis([-40 440 -40 200])


m=gca;
m.DataAspectRatio=[120 120 1];
colorbar
colormap(rgb);
caxis([-0.2 0.2])