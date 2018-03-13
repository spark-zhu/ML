
% Xue_List = zeros(5,2);
% load('D:\Hanlin_tracking\recording 01-13-17_right\rgb')
% Import    UnitToBeDrawn -by- Session.Unit.Week , manaul input from Book2
XueCollection=cell(1,4);%4unit to be drawn;Fill in values before proceed.
Weeks = 1:11; % weeks to be monitored. Or Days to be monitored, or different units to be monitored 
% Note if a session number is zero, it means that at that week, a unit is
% not found.
XueCollection{1,5}=[];
%%
clear PrincessTable
for Unit=5%1:4   % use any one 1:4 or only the fifth one or create a sixth one 
for week = 1:numel(Weeks) % you can also treat this as different units. 
ses = XueCollection{Unit}(1,week);
if ses~=0
unitName = XueCollection{Unit}(2,week);
unitIndex = find(AllData{ses}.selectedClusters == unitName);
PrincessTable(Unit,week)= lookUpTable( lookUpTable(:,1)==ses & lookUpTable(:,2)==unitIndex, 3); 
end
end
end
%%
colorlist = ['r','g','b','y'];
clear colorlist
colorlist(1,:)=rgb(1,:);
colorlist(2,:)=rgb(2,:);
colorlist(3,:)=rgb(7,:);
colorlist(4,:)=rgb(9,:);
for UNIT = 1:4
list = PrincessTable(UNIT,:);
 h=figure
        Ox=1.6;
        Oy=21;
        gap=0.1;
        Gw=1.8;
        Gh=4.5;
        numGx = numel(Weeks);
        numGy = 1;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        
       % select a channel first 
       
       store = zeros(32,41,numel(list));
       for i=1:numel(list)
           if list(i)~=0
           store(:,:,i)=DataForFeature{list(i)};
           end
       end
       MeanStore=mean(store,3);
       AvgP2P=max(MeanStore')-min(MeanStore');
       selectedCh = find(AvgP2P==max(AvgP2P));
       range = [min(min(store(selectedCh,:,:))) max(max(store(selectedCh,:,:)))];
       
        for row=1:1
        for col=1:numel(Weeks)
        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')
        
        if list(col)~=0
        plot(store(selectedCh,:,col),'LineWidth',3,'Color',colorlist(UNIT,:))
        
        end
        axis([1 41 -180 50])
        box off
        title(['Day ' num2str(Weeks(col)*7)],'FontWeight','bold','FontSize',12)
        if col~=1
            axis off
        end
        if col==1
            ylabel('Amplitude (uV)','FontWeight','bold','FontSize',12)
             set(ax,'XTick',[1 41])
            set(ax,'xticklabels',{'0','2'})
            xlabel('t (ms)','FontWeight','bold','FontSize',12)
        end
        end
        end
        set(gcf, 'Position', get(0,'Screensize'));
 set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')
set(h,'PaperPosition',[0 0 40 27])
print(h,'-dpng',num2str(UNIT))
end
       
%% Plotting amplitude distribution 

for UNIT = 5:5
h=figure;
Ox=12.5    ;    
Oy=12.5;
        gap=5;
        Gw=25;
        Gh=25;
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
        list = PrincessTable(UNIT,:);
        list(list==0)=[];
        P2PStore = cell2mat(P2PForFeature(list));
        clear MeanStore
        for mI = 1:numel(list)
        MeanStore(:,:,mI)=DataForFeature{list(mI)};
        end
        MeanWaveform = mean(MeanStore,3);
        % create the variable SelectedLocation=zeros(4,8) and input some
        % ones to indicate locations used for localization.
        LocationList = reshape(((Ch_Map_new-15).*SelectedLocation),32,1);
        XList=reshape(X.*SelectedLocation,32,1);
        XList = XList(XList>0);
        YList = YList(YList>0);
        ZList = zeros(size(XList));
        YList=reshape(Y.*SelectedLocation,32,1);
        ampselection = LocationList(LocationList>0)
        ZAMP = MeanWaveform(ampselection,:);
        [U,S,V] = svd(ZAMP);
        En=U(:,2:end);
        
        f = @(x)parameterfun(x,XList,YList,ZList,En);
        [x,fval] = fminunc(f,[-100 -100 0])
        
        for ch=1:32
            if sum(P2PStore(:,ch))==0
                P2Pmean(ch)=0;
            else
        P2Pmean(ch) = mean(P2PStore(P2PStore(:,ch)>0,ch));
        
            end
        end
        
        for row =1:4
            for col =1:8
            Z(row,col) =   P2Pmean(Ch_Map_new(row,col)-15);
            end
        end
        createFit_2D(X,Y,Z)
       
        m=colorbar;
ylabel(m, 'Peak To Peak Amplitude','FontSize',10)
        daspect([105 104 1])


        axis off 
print(h,'-dpng',[ 'Mask' num2str(UNIT)])
end




%% PLotting across different time
for Month = [1 5 9]



 list = PrincessTable(:,Month);

        h=figure;
        set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
        Ox=1.6;
        Oy=1.6;
        gap=0.6;
        Gw=3;
        Gh=4.5;
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        t=[-10/20e3:1/20e3:30/20e3];

        for row=1:4
            for col=1:8
        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')
                for unit=1:numel(list)
                    Ch_Map=AllData{lookUpTable(find(lookUpTable(:,3)==list(unit)),1)}.Ch_Map;
                    temp=Ch_Map(row,col)-15;
                    if Ch_Map(row,col)>0 && MaskFeatureAll{list(unit)}(temp)>=0
                        
                        plot(t,DataForFeature{list(unit)}(temp,:),'Color',colorlist(unit,:),'LineWidth',3)
                        hold on;
                    end
                    axis([t(1) t(end) -150 50])
                    
                    
                end
                if row~=4 | col~=1
                set(ax,'XTick',[])
                set(ax,'YTick',[])
                end
                if row==4 & col==1
             ylabel('Amplitude (uV)','FontWeight','bold','FontSize',10)
             set(ax,'XTick',[t(1) t(end)])
            set(ax,'xticklabels',{'0','2'})
            xlabel('t (ms)','FontWeight','bold','FontSize',12)
                end
            end
        end  

 set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')
set(h,'PaperPosition',[0 0 40 27])
print(h,'-dpng',[ 'Week' num2str(Month)],'-loose')
end
%%

%% Plotting ISI

% Sessions = 4:7; Actually Let's just do this for everyday.  
figure
for UNIT = 1:4
    TS_thatUnit_thatWeek = cell(1,numel(Weeks));
HIST_thatUnit_thatDay = zeros(numel(Weeks),25);

for week = 1:numel(Weeks)
list = PrincessTable(UNIT,week);
if list>0
TS_thatUnit_thatDay{week}=unique(cell2mat(TimeFeatureAll(list)));
ISI_in_MS_bins = diff(double(TS_thatUnit_thatDay{week}))/20;
[counts,centers] = hist(ISI_in_MS_bins,10:20:510);
counts=counts(1:25);
centers=centers(1:25);
HIST_thatUnit_thatDay(week,:)=counts/sum(counts);
end
end
subplot(2,2,UNIT)

bar3(HIST_thatUnit_thatDay',0.5)
view(-77,25)
daspect([9 25 1])
pbaspect auto
set(gca,'GridLineStyle','none')

if UNIT==4
set(gca,'YTick',[5 15 25])
set(gca,'yticklabels',{'100','300','500'})
set(gca,'XTick',[1 4 7 9])
set(gca,'xticklabels',{'6','9','12','14'})

ylabel('ms','FontWeight','bold','FontSize',8)
xlabel('week','FontWeight','bold','FontSize',8)
zlabel('Frequency','FontWeight','bold','FontSize',8)
end
if UNIT~=4
set(gca,'YTick',[])
set(gca,'XTick',[])
end
end




%    
%     
% 
%      
%         clear Days
%         
%         for unit=1:numel(list)
%         Days(unit)=lookUpTable(find(lookUpTable(:,3)==list(unit)),1);
%         end
%         Days = unique(Days);
%         
%         Ox=0.6;
%         Oy=22;
%         gap=0.6;
%         Gw=18;
%         Gh=3;
%         ax=axes('Units','centimeters','position', [Ox Oy Gw Gh]);
%         set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
%         for day=1:numel(Days)
%             scatter(Days(day),1,250, rgb(Days(day),2:4)/255,'filled')
%             hold on;
%         end
%         axis([1 numel(AllData) 0 2])
%         set(gcf, 'Position', get(0,'Screensize'));
%                 set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')
% set(h,'PaperPosition',[0 0 40 27])
% print(h,'-dpdf',pdfCellArray{clu},'-loose')
%% Location Distribution Plot


clear color       
color=jet(49);
% color=color(5:104,:);
firstRow=[0 0 0];
color=[firstRow;color];
color=flip(color,1);

% PrincessTable(1:4,:)=[];
totalList = PrincessTable;


rgb = linspecer(numel(totalList));
% plotting 
h=figure;
        set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
        Ox=0;
        Oy=1.6;
        Gw=50;
        Gh=24;
   
        ax=axes('Units','centimeters','position', [Ox Oy Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')

        
        
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
    for unit = 1:numel(totalList)
        unit
        P2P = P2PForFeature{totalList(unit)};
        try
       Wx = NeuronLocationList(unit,1);
        Wy = NeuronLocationList(unit,2);
        Wz = NeuronLocationList(unit,3);

        catch
            
        [sortedIntensity IntensityOrder]=sort(P2P,'descend');
Ch_Map_2 = Ch_Map_new-15;
for Location =1:numel(IntensityOrder)
Row(Location) = mod(find(Ch_Map_2==IntensityOrder(Location)),4);
if Row(Location) == 0
Row(Location)=4;
end
Col(Location) = ceil(find(Ch_Map_2==IntensityOrder(Location))/4);
end
% ROUND 1
Graph = zeros(4,8);
InitialSelection = IntensityOrder(1:6);
for Location =1:6
Graph(Row(Location),Col(Location))=1;
end
CC = bwconncomp(Graph,4);
clear CountList
if numel(CC.PixelIdxList)>1
for listE =1:numel(CC.PixelIdxList)
CountList(listE) = numel(CC.PixelIdxList{listE});
end
[sortValue sortI]=sort(CountList);

    ToBeCleared=Ch_Map_2(CC.PixelIdxList{sortI(1)});
for del=1:numel(ToBeCleared)
    Row(IntensityOrder==ToBeCleared(del))=[];
    Col(IntensityOrder==ToBeCleared(del))=[];
IntensityOrder(IntensityOrder==ToBeCleared(del))=[];
end
end
% ROUND 2
clear CountList
Graph = zeros(4,8);
InitialSelection = IntensityOrder(1:6);
for Location =1:6
Graph(Row(Location),Col(Location))=1;
end
CC = bwconncomp(Graph,4);


for listE =1:numel(CC.PixelIdxList)
CountList(listE) = numel(CC.PixelIdxList{listE});
end
[sortValue sortI]=sort(CountList);
if numel(CC.PixelIdxList)>1
    ToBeCleared=Ch_Map_2(CC.PixelIdxList{sortI(1)});
for del=1:numel(ToBeCleared)
    Row(IntensityOrder==ToBeCleared(del))=[];
    Col(IntensityOrder==ToBeCleared(del))=[];
IntensityOrder(IntensityOrder==ToBeCleared(del))=[];
end

end
%
SelectedElements = Ch_Map_2(CC.PixelIdxList{sortI(end)});
if numel(P2P(SelectedElements)>25)<3
if numel(SelectedElements)>4
    clear sortValue sortI
[sortValue sortI]=sort(P2P(SelectedElements));
SelectedElements=SelectedElements(sortI(end-3:end));
end
end

clear Row Col Graph
for Location =1:numel(SelectedElements)
Row(Location) = mod(find(Ch_Map_2==SelectedElements(Location)),4);
if Row(Location) == 0
Row(Location)=4;
end
Col(Location) = ceil(find(Ch_Map_2==SelectedElements(Location))/4);
end
% 
SelectedLocation = zeros(4,8);
for Location =1:numel(SelectedElements)
SelectedLocation(Row(Location),Col(Location))=1;
end

%Prepare for estimation
%
Ox=12.5;    
Oy=12.5;
        gap=5;
        Gw=25;
        Gh=25;
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
%         load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')

        %MeanWaveform = mean(MeanStore,3);
        % create the variable SelectedLocation=zeros(4,8) and input some
        % ones to indicate locations used for localization.
        LocationList = reshape(((Ch_Map_2).*SelectedLocation),32,1);
        XList=reshape(X.*SelectedLocation,32,1);
        YList=reshape(Y.*SelectedLocation,32,1);
        XList = XList(XList>0);
        YList = YList(YList>0);
        ZList = zeros(size(XList));
        Loc2List =  reshape(((Ch_Map_2).*SelectedLocation),32,1);
        ampselection = Loc2List(Loc2List>0);
%         for amp=1:numel(ampselection)
%         ampselection(amp)=find(Dat_V_Map(:,2)==ampselection(amp));
%         end
        MeanPlot = DataForFeature{totalList(unit)};
        ZAMP = MeanPlot(ampselection,:);
        [U,S,V] = svd(ZAMP);
        En=U(:,2:end);
        clear f
        f = @(x)parameterfun(x,XList,YList,ZList,En);
        
% WC_Y=sum(sum(Weighted_P2P_Remap.*(Y+1)))/sum(sum(Weighted_P2P_Remap));
% WC_X=sum(sum(Weighted_P2P_Remap.*(X+1)))/sum(sum(Weighted_P2P_Remap));
% MUSIC Location Calculation

options.Display='off';
tic
[gX,gY,gZ]=meshgrid([-108 -36 -12 -4 4 12 36 108],[-108 -36 -12 -4 4 12 36 108],[-108 -36 -12 -4 4 12 36 108]);
ch=find(P2P==max(P2P));

highX = X(find(Ch_Map_2==ch));
highY = Y(find(Ch_Map_2==ch));
gX=gX+highX;
gY=gY+highY;
Full_Coordinates = cell(8,8,8);
FVAL = zeros(8,8,8);
for cx=1:8
    for cy=1:8
       
        for cz=1:8
        [MUSIC,fval] = fminunc(f,double([gX(cx,cy,cz) gY(cx,cy,cz) gZ(cx,cy,cz)]),options);
        Full_Coordinates{cx,cy,cz}=MUSIC;
        FVAL(cx,cy,cz)=fval;
        end
    end
end
toc
FVAL_list = reshape(FVAL,512,1);
Full_Coordinates_List = cell2mat(reshape(Full_Coordinates,512,1));

% Find Global Minimum.
OverallMinCorD = Full_Coordinates_List(FVAL_list == min(FVAL_list),:);



        Wx = OverallMinCorD(1,1);
        Wy = OverallMinCorD(1,2);
        Wz = abs(OverallMinCorD(1,3));
        if Wz> size(color,1)
            Wz = size(color,1)-1;
        end
        NeuronLocationList(unit,:)=[Wx Wy Wz];
    end
        maxP2P= max(P2PForFeature{totalList(unit)});
if Wz> size(color,1)
            Wz = size(color,1)-1;
end

scatter(Wx,Wy,400,'+','LineWidth',2,'MarkerEdgeColor',color(floor(abs(Wz)+1),:));
text(Wx,Wy-10,num2str(unit+9));
% if unit==1
% colorbar('location','SouthOutside')   
% end        
%     scatter(double(Wx),double(Wy),maxP2P^2/100,rgb(unit,:),'.','LineWidth',4)
%     text(double(Wx),double(Wy)+0.5,[num2str(floor(P2P)) 'uV'],'FontSize',8)
    hold on
    [C ,I] = max(P2PForFeature{totalList(unit)});
    wf= DataForFeature{totalList(unit)}(I,:);
    t=-10:30;
    hold on;
    timeConstant=1.7;
    ampConstant = 6;
    timeShift = -5;
    ampShift = 10;
    plot(t/timeConstant+Wx+timeShift,wf/ampConstant+Wy+ampShift,'Color',rgb(unit,:),'LineWidth',4)
    hold on; 
        plot([t(1)/timeConstant+Wx+timeShift t(1)/timeConstant+Wx+timeShift],[min(wf)/ampConstant+Wy+ampShift (min(wf)+50)/ampConstant+Wy+ampShift],'Color',rgb(unit,:),'LineWidth',3.5)
    hold on; 

    end
     
        
        Ox=12.5;
        Oy=12.5;
        gap=5;
        Gw=25;
        Gh=25;
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
        plot([pAx pCx],[pAy pAy],'k','LineWidth',2)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',2)
        end
    end
    
    
    
    
      axis([-15 250 -15 130])  
ylabel('Micrometers','FontWeight','bold','FontSize',10)
xlabel('Micrometers','FontWeight','bold','FontSize',10)
title('Relative Location of Units Tracked on the electrode array','FontWeight','bold','FontSize',12)

        
        set(gcf, 'Position', get(0,'Screensize'));
        set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')
set(h,'PaperPosition',[0 0 50 27])
axis off
print(h,'-dpng',[ 'Week'],'-loose')
%% Plotting Weekly Summary 
result.store=store1;
DPS = result.store.DPS;
lookUpTable = result.store.lookUpTable;
P2PForFeature=result.store.P2PForFeature;
weeksTracked = result.Groups;
weeksTracked(1)=[];

% need two arraies, one showing the p2p of each group , the othe one marks
% single/ multi unit of each group. 
%%
TimeFeatureAll= result.store.TimeFeatureAll;
isSingle = zeros(1,numel(weeksTracked));
P2P = max(cell2mat(P2PForFeature)');
for gp=1:numel(weeksTracked)
    
    oneGroup = weeksTracked{gp};
    P2Pavg(gp) = mean(P2P(oneGroup));
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


isSingleUnit = zeros(1,size(lookUpTable,1));

for unit=1:numel(isSingleUnit)
    
FiringTimeForThisUnit = TimeFeatureAll{unit};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);
Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;
FR_avg = numel(FiringTimeForThisUnit)/720;
if counts(1)<3 % which means on average (across all sessions) , you see in 12 mins, less than 1 occasion, it fires within 1ms. 
    isSingleUnit(unit)=1;
end


end


%%
% convert session to WPS
week10AllUnits = find(lookUpTable(:,4)==17);
weeksTracked_P2P = weeksTracked;
for gp=1:numel(weeksTracked)
unitsFromThatGp = weeksTracked{gp};
unitsFromThatGp_P2P = P2P(weeksTracked{gp});
weeksTracked{gp} = arrayfun(@(x) lookUpTable(x,4), unitsFromThatGp);
[C,ia,ic]= unique(weeksTracked{gp});
weeksTracked{gp} = C;
weeksTracked_P2P{gp} = unitsFromThatGp_P2P(ia);
wk10= intersect(unitsFromThatGp,week10AllUnits);
if ~isempty(wk10)
week10units(gp) =wk10(1);
end
end

%%
select = cellfun(@(x) max(x)-min(x) , weeksTracked);
select = cellfun(@numel , weeksTracked);
% [a,b]=sort(select,'descend');
% weeksTracked=weeksTracked(b);
% select = cellfun(@numel , weeksTracked);
% 
% bk = weeksTracked;
weeksTracked = weeksTracked(select>4);
weeksTracked_P2P=weeksTracked_P2P(select>4);
isSingle=  isSingle(select>4);
P2Pavg = P2Pavg(select>4);
week10units=week10units(select>4);


% minColor = min(P2Pavg);
% maxColor = max(P2Pavg);
minColor = log(min(cell2mat(weeksTracked_P2P)));
maxColor = log(max(cell2mat(weeksTracked_P2P)));

colorInt=maxColor-minColor;

color=jet(120);
color=color(5:104,:);
firstRow=color(1,:);
color=[firstRow;color];
% 'Color',color(floor((value-minColor)/colorInt*100)+1,:)
Single_weeks_tracked = weeksTracked(logical(isSingle));
Single_P2Pavg =  P2Pavg(logical(isSingle));
Single_weeks_tracked_P2P = weeksTracked_P2P(logical(isSingle));
% sort according to color  P2P 
 [a,b]=sort(Single_P2Pavg);
Single_weeks_tracked=Single_weeks_tracked(b);
Single_P2Pavg=Single_P2Pavg(b);
Single_weeks_tracked_P2P=Single_weeks_tracked_P2P(b);


figure 
for i=1:numel(Single_weeks_tracked)
    days=Single_weeks_tracked{i}*7;
    P2P_list = Single_weeks_tracked_P2P{i};
    for j=1:numel(days)
%         if P2P_list(j)>200
%         P2P_list(j)=200;
%         end
    rectangle('Position',[days(j)-3,i-0.5,6,0.8],'FaceColor',color(floor((log(P2P_list(j))-minColor)/colorInt*100)+1,:)) 
    hold on
    end
% scatter(weeksTracked{i}*7,i*ones(size(weeksTracked{i})),200*ones(size(weeksTracked{i})),'filled','MarkerFaceColor',rgb(i,:))
% hold on
end
colormap(jet)
c=colorbar;
colormap jet 
c.Ticks=[0 0.2 0.4 0.6 0.8 1];
c.TickLabels=num2cell(floor(exp(c.Ticks*colorInt+minColor)));
set(gca,'XTick',42:7:91)
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Unit #')
title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
Single_weeks_tracked_save=Single_weeks_tracked;
Single_P2Pavg_save=Single_P2Pavg;
Single_P2P_save=Single_weeks_tracked_P2P;

%%
Single_weeks_tracked = weeksTracked(~logical(isSingle));
Single_P2Pavg = P2Pavg(~logical(isSingle));
Single_weeks_tracked_P2P = weeksTracked_P2P(logical(~isSingle));

% sort according to color  P2P 
 [a,b]=sort(Single_P2Pavg);
Single_weeks_tracked=Single_weeks_tracked(b);
Single_P2Pavg=Single_P2Pavg(b);
Single_weeks_tracked_P2P=Single_weeks_tracked_P2P(b);



figure 
for i=1:numel(Single_weeks_tracked)
    days=Single_weeks_tracked{i}*7;
        P2P_list = Single_weeks_tracked_P2P{i};

    for j=1:numel(days)
    rectangle('Position',[days(j)-3,i-0.5,6,0.8],'FaceColor',color(floor((log(P2P_list(j))-minColor)/colorInt*100)+1,:)) 
    hold on
    end
% scatter(weeksTracked{i}*7,i*ones(size(weeksTracked{i})),200*ones(size(weeksTracked{i})),'filled','MarkerFaceColor',rgb(i,:))
% hold on
end
colormap(jet)
c=colorbar;
colormap jet 
c.Ticks=[0 0.2 0.4 0.6 0.8 1];
c.TickLabels=num2cell(floor(exp(c.Ticks*colorInt+minColor)));

set(gca,'XTick',42:7:91)
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Unit #')
title('Multi-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
%% Heying
% Single_weeks_tracked_save=Single_weeks_tracked;
% Single_P2Pavg_save=Single_P2Pavg;

Single_weeks_tracked=[Single_weeks_tracked Single_weeks_tracked_save];
Single_P2Pavg=[Single_P2Pavg Single_P2Pavg_save];
Single_P2P = [Single_weeks_tracked_P2P Single_P2P_save];

figure 
for i=1:numel(Single_weeks_tracked)
    days=Single_weeks_tracked{i}*7;
    P2P_list = Single_P2P{i};
    for j=1:numel(days)
    rectangle('Position',[days(j)-3,i-0.5,6,0.8],'FaceColor',color(floor((log(P2P_list(j))-minColor)/colorInt*100)+1,:)) 
    hold on
    end
% scatter(weeksTracked{i}*7,i*ones(size(weeksTracked{i})),200*ones(size(weeksTracked{i})),'filled','MarkerFaceColor',rgb(i,:))
% hold on
end
colormap(jet)
c=colorbar;
colormap jet 
c.Ticks=[0 0.2 0.4 0.6 0.8 1];
c.TickLabels=num2cell(floor(exp(c.Ticks*colorInt+minColor)));

set(gca,'XTick',42:7:119)
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Unit #')
title('Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
%%
    figure
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
    week10units(week10units==0)=[];
    LocMat = cell2mat(result.store.LocForFeature(week10units));
    minColor = min(P2P(week10units));
     maxColor = max(P2P(week10units));
    colorInt=maxColor-minColor;

for Loc =1:numel(week10units)
       scatter(LocMat(Loc,1),LocMat(Loc,2),LocMat(Loc,3).^2+300,'*','LineWidth',4,'MarkerEdgeColor',color(floor((P2P(week10units(Loc))-minColor)/colorInt*100)+1,:))
       hold on;
end
axis([-40 440 -40 200])
figure
 for Loc = 1:numel(week10units)
 
 scatter3(LocMat(Loc,1),LocMat(Loc,2),LocMat(Loc,3),(P2P(week10units(Loc))/1).^2,'filled','MarkerFaceColor','green') 
hold on
 end
view(30, 60);

axis([-40 440 -40 200])
%% plotting all unit across weeks or all trackable units across sesssions. 
    minColor = log(min(P2P));
     maxColor = log(max(P2P));
    colorInt=maxColor-minColor;
color=jet(120);
color=color(5:104,:);
firstRow=color(1,:);
color=[firstRow;color];

for ses = 1:1% numel(unique(lookUpTable(:,1)))
 h=figure('Name',num2str(ses));
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
        plot([pAx pCx],[pAy pAy],'k','LineWidth',2)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',2)
        end
    end
    hold on 
    week10AllUnits = find(lookUpTable(:,1)==ses);
    week10units= week10AllUnits;
    P2P = max(cell2mat(P2PForFeature)');

    LocMat = cell2mat(result.store.LocForFeature(week10units));

for Loc =1:numel(week10units)
       scatter(LocMat(Loc,1),LocMat(Loc,2),mean(LocMat(:,3)).^2+200,'*','LineWidth',2.5,'MarkerEdgeColor',color(floor((log(P2P(week10units(Loc)))-minColor)/colorInt*100)+1,:))
       hold on;
end
axis([-40 440 -40 200])
set(gca,'DataAspectRatio',[120 120 1])
c=colorbar;
colormap jet 
c.Ticks=[0 0.2 0.4 0.6 0.8 1];
c.TickLabels=num2cell(floor(exp(c.Ticks*colorInt+minColor)));
title(['Day ' num2str(DPS(ses))])
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2,'Position',[0.08 0.04 0.9 0.9])
saveas(h,[h.Name '.png'])
end
%% PCA-RGB-PLot
    P2P = max(cell2mat(P2PForFeature)');

MaxWaveFormCollector = zeros(numel(P2P),41);
for wave=1:numel(P2P)
    maxP2Pch=find(P2P(wave)==P2PForFeature{wave});
    maxP2Pch=maxP2Pch(1);
MaxWaveFormCollector(wave,:)=result.store.DataForFeature{wave}(maxP2Pch,:);
MaxWaveFormCollector(wave,:)=(MaxWaveFormCollector(wave,:)-min(MaxWaveFormCollector(wave,:)))/P2P(wave);
end
[~,score,~]  = pca(MaxWaveFormCollector(remain_index==1,:),'NumComponents',3);

thre=0.95;
for pc=1:3
sortPC = sort(score(:,pc));
CutValue = sortPC(floor(size(score,1)*thre));
score(score(:,pc)>CutValue,pc)=CutValue;

CutValue = sortPC(floor(size(score,1)*(1-thre)));
score(score(:,pc)<CutValue,pc)=CutValue;

end

clear RGB
for pc=1:3
RGB(:,pc) = (score(:,pc)-min(score(:,pc)))/(max(score(:,pc))-min(score(:,pc)));
end
RGB_new = zeros(size(lookUpTable,1),3);
RGB_new(remain_index==1,:)=RGB;
%%
clear F
for ses = 1:numel(unique(lookUpTable(:,1)))
 h=figure('Name',num2str(ses));
 set(h,'Visible','off');
 ses
    
Ox=12.5;    
Oy=12.5;
        gap=5;
        Gw=25;
        Gh=25;
if ChMapNum~=1
Ox=15;    
Oy=15;
gap=20;
        Gw=30;
        Gh=30;
end
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
        plot([pAx pCx],[pAy pAy],'k','LineWidth',2)
        hold on
        plot([pCx pCx],[pAy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pCx],[pCy pCy],'k','LineWidth',2)
        hold on 
        plot([pAx pAx],[pAy pCy],'k','LineWidth',2)
        end
    end
    hold on 
    week10AllUnits = lookUpTable(:,1)==ses & remain_index'==1;
    week10units= find(week10AllUnits==1);
 
    LocMat = cell2mat(result.store.LocForFeature(week10units));
     
for Loc =1:sum(week10AllUnits)
       scatter(LocMat(Loc,1),LocMat(Loc,2), P2P(week10units(Loc)),'filled','LineWidth',2.5,'MarkerFaceColor',RGB_new(week10units(Loc),:),	'MarkerFaceAlpha',0.8)
       hold on;
end
centerX = 190;
centerY = 90;
LocForPCA_raw = [LocMat(:,1)-centerX LocMat(:,2)-centerY];

clear LocForPCA
LocForPCA=[];
for Loc = 1:sum(week10AllUnits)
    weight = floor(P2P(week10units(Loc))^2/10);
LocForPCA = [LocForPCA;repmat(LocForPCA_raw(Loc,:),weight,1)];
end
[coeff,score,latent,tsquared,explained,mu]= pca(LocForPCA);
Center = [centerX centerY];
clear PC1 PC2
PC1 = Center+coeff(:,1)'*explained(1)*5;
PC2 = Center+coeff(:,2)'*explained(2)*5;

hold on
plot([centerX PC1(1)],[centerY PC1(2)],'r','LineWidth',2)
hold on
plot([centerX PC2(1)],[centerY PC2(2)],'b','LineWidth',2)

% axis([-87.5 322.5 -87.5 202.5])
axis([-85 465 -85 265])

set(gca,'DataAspectRatio',[120 120 1])
c=colorbar;
% colormap jet 
% c.Ticks=[0 0.2 0.4 0.6 0.8 1];
% c.TickLabels=num2cell(floor(exp(c.Ticks*colorInt+minColor)));
title(['Day ' num2str(DPS(ses))])
set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2,'Position',[0.08 0.04 0.9 0.9])
% saveas(h,[h.Name '.png'])
F(ses) = getframe(h);    
close(h)
end
v = VideoWriter('PCA-unscale-clamp95.avi');
% movie(F)
v.FrameRate=1;
open(v)
writeVideo(v,F);
close(v)
X_color = (1-RGB(:,1)+RGB(:,2))/2;
Y_color = 1/sqrt(3)*((1-RGB(:,1)+RGB(:,2))/2+RGB(:,3));
%%
% figure
% for pt=1:size(RGB,1)
% scatter(X_color(pt),Y_color(pt),'filled','MarkerFaceColor',RGB(pt,:))
% hold on
% end
RGB_color_time = jet(size(RGB,1));
figure
for pt=1:size(RGB,1)
scatter3(RGB(pt,1),RGB(pt,2),RGB(pt,3),'filled','MarkerFaceColor',RGB_color_time(pt,:))
hold on
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')


%% Calculating boundary distance 
Locations = cell2mat(result.store.LocForFeature);
Boundary_Distance = zeros(1,size(Locations,1));
Boundary=[min(Gx) max(Gx) min(Gy) max(Gy)];
for i=1:numel(Boundary_Distance)
Boundary_Distance(i)=Electrode_Unit_distance(Locations(i,1:2), Boundary);
end
figure
[y,x] = ecdf(Boundary_Distance(isSingleUnit==1));
plot(x,y)

[y,x] = ecdf(Boundary_Distance(isSingleUnit==0));
hold on
plot(x,y,'r')
legend('single','multi')
xlabel('distance to eletrode boundary')

figure 
[counts] = histc(Boundary_Distance(isSingleUnit==1),-75:5:100);
plot(-72.5:5:102.5,counts/sum(counts));
[counts] = histc(Boundary_Distance(isSingleUnit==0),-75:5:100);
hold on
plot(-72.5:5:102.5,counts/sum(counts),'r');
legend('single','multi')
xlabel('distance to eletrode boundary')

figure
errorbar([1 2],[mean(Boundary_Distance(isSingleUnit==1)) mean(Boundary_Distance(isSingleUnit==0))],[std(Boundary_Distance(isSingleUnit==1)),std(Boundary_Distance(isSingleUnit==0))])
xlim([0 3])
%%
figure
for i=1:numel(weeksTracked)
    days=weeksTracked{i}*7;
    for j=1:numel(days)
    rectangle('Position',[days(j)-3,i,6,0.8],'FaceColor',color(floor((P2Pavg(i)-minColor)/colorInt*100)+1,:)) 
    hold on
    end
% scatter(weeksTracked{i}*7,i*ones(size(weeksTracked{i})),200*ones(size(weeksTracked{i})),'filled','MarkerFaceColor',rgb(i,:))
% hold on
end
colormap(jet)
axis([0 18*7 0 12])

% figure
% for i=1:numel(weeksTracked)
% scatter(weeksTracked{i}*7,i*ones(size(weeksTracked{i})),200*ones(size(weeksTracked{i})),'filled','MarkerFaceColor',rgb(i,:))
% hold on
% end
axis([0 18*7 0 12])
set(gca,'XTick',7:7:126) 
set(gca,'YTick',1:11)
ylabel('Unit Index','FontWeight','bold','FontSize',12)
xlabel('Days Post Surgery Tracked','FontWeight','bold','FontSize',11)
title('Summary of Tracked Units','FontWeight','bold','FontSize',11)
%% heying 

%% Plotting tracking of overall waveform across time, you could also change it to plot identified units in one cluster
% PrincessTable(1:4,:)=[];
totalList = week10AllUnits;
t=-10:30

rgb = linspecer(numel(totalList));
    AllData=result.store.AllData;
WkUnitList=week10AllUnits;     
    dis=0;
h=figure;
        set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
        Ox=0.3;
        Oy=0.3;
        gap=0.3;
        Gw=6;
        Gh=6;
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        
        for row=1:4
            for col=1:8
        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')
                for unit=1:numel(WkUnitList)
                    Ch_Map=AllData{lookUpTable(find(lookUpTable(:,3)==WkUnitList(unit)),1)}.Ch_Map;
                    temp=Ch_Map(row,col)-15;
                    if Ch_Map(row,col)>0 && MaskFeatureAll{WkUnitList(unit)}(temp)>0.1
                        
                        plot(t+(unit-1)*dis,DataForFeature{WkUnitList(unit)}(temp,:),'Color',rgb(unit,:),'LineWidth',2)
                        hold on;
                    end
                    axis([t(1) t(end)+(numel(WkUnitList)-1)*dis -150 50])
                    set(gca,'xticklabel',[])
                    set(gca,'yticklabel',[])
                    
                end
            end
        end

%% Signal to Noise estimate
clear all
close all
folderList=dir;
folderList(1:12)=[];
NoiseCollection=cell(1,1);
Legend=cell(1,1);
for folder=1:numel(folderList)
  Name=folderList(folder).name;
   if ~isempty(strfind(Name,'record'))
       loc=strfind(Name,'g');
       Legend{end+1}=folderList(folder).name(loc+1:end);
      cd(folderList(folder).name)
load('stdV.mat')
%AvgImpedance(AvgImpedance>2E6)=2E6;
NoiseCollection{end+1}=stdNoise;

    cd ..
   end
   
end
NoiseCollection=NoiseCollection(2:end);
Legend=Legend(2:end);

prompt='surgery Date in mmddyyyy: ';
surgeryDate = input(prompt,'s');
startNum = datenum(surgeryDate,'mmddyyyy');
DPS = zeros(1,numel(Legend));
for i=1:numel(Legend)
DPS(i) = datenum(Legend{i},'yyyy-mm-dd')-startNum;
WPS(i) = floor(DPS(i)/7);
end

mean_std_noise = cellfun(@mean, NoiseCollection);
std_std_noise = cellfun(@std, NoiseCollection);

figure
errorbar(DPS,mean_std_noise,std_std_noise)
hold on;
set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('All channels Std ')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
ylim([0 20])



SNR = (max_P2P.^2)./(3*mean_std_noise).^2;% from tracking oversampling new metrics 
plot(DPS,SNR);
opt.XLabel = 'Days Post Surgery'; % xlabel
opt.YLabel = 'SNR'; %ylabel
opt.XLim = [0, DPS(end)+5]; % set x axis limit
opt.YLim = [0, 80]; % set y axis limit

% Save? comment the following line if you do not want to save
opt.FileName = 'plotAxisLimit.jpg'; 

% create the plot
setPlotProp(opt);



opt.XLabel = 'Days Tracked'; % xlabel
opt.YLabel = 'number of neurons'; %ylabel
opt.XLim = [0, x(end)]; % set x axis limit
opt.YLim = [0, 12]; % set y axis limit

% Save? comment the following line if you do not want to save
opt.FileName = 'plotAxisLimit.jpg'; 

% create the plot
setPlotProp(opt);
    