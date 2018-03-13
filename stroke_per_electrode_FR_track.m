%%?
load('10182017_tracking-1.mat');
resultPre = result;
resultStroke = result;
resultPost = result;

resultPre.date='10152017';
resultStroke.date='10162017';
resultPost.date='10172017';

t1=[1200 2000]*20E3;
t2=[2300 3100]*20E3;
t3=[3200 4000]*20E3;


load('Finger_map.mat');
Ch_Map= Ch_Map_new-15;
PeakCh   = cellfun(@(x)  find(x==max(x)),result.P2P) ;





%%
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
load('store.mat');
load('D:\Box Sync\0919 stroke\Finger_map.mat');
Ch_Map= Ch_Map_new-15;
channel_FR = NaN(numel(store1.AllData),4,8);
daymark = [-1 0 1 2 3 5 7 9 11 13 16 22 29 36 43 50 57];
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
    
    channel_FR(day,row,col)=numel(cell2mat(Global_Time_Store{row,col}));
    else
       if result.Ch_Map_2(row, col)~=0
           channel_FR(day,row,col)=0;
       end
    end
    end
end
channel_FR(day,:,:)=channel_FR(day,:,:)/max(cell2mat(result.time')/20E3);
numOfElePerRow = zeros(4,1);
TimeStampPerRow = cell(4,1);

for r=1:4
 numOfElePerRow (r)= sum(~isnan(channel_FR(day,r,:)));
 TimeStampPerRow{r}= cellfun(@(x) cell2mat(x),Global_Time_Store(r,:),'UniformOutput',0);
 FR(r,day)=numel(cell2mat(cellfun(@(x) x ,TimeStampPerRow{r},'UniformOutput',0)))/numOfElePerRow (r)/max(cell2mat(result.time')/20E3);
% if isnan(FR(r,day))
%     FR(r,day)=0;
% end
figure(r);clf; title(sprintf('%d',r));
set(gcf,'position',[100+r*300,400,300 800]);
for col=1:8
    subplot(8,1,col);
    plot(daymark, channel_FR(:,r,col),'o-','markersize',5, 'linewidth',2);
    xlim([-2 50]);
    if col~=8
    set(gca,'fontsize',14,'xticklabel',[]);
    end
end

% FR = (Single_fire_Count+Multi_fire_Count)./cellfun(@(x) max(cell2mat(x.time')/20E3),AllData)./cellfun(@(x) size(x.Dat_V_Map,1),AllData);
end
end



h=figure
for i =1:4
subplot(4,1,i)
plot(daymark,FR(i,:),'LineWidth',2)
% dayLabel ={'-1','Pr','S','Po','0','1','2','3','5','7','9','11','13','22'};
% xticks(1:12)
% 
% xticklabels(dayLabel)
% ylabel('Hz')
% ylim([0 35])
end
