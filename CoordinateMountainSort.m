%% for open-ephys
clear all
% addpath(genpath('/home/hanlin/MatlabFiles/ML_code'))
DIR=dir('*.continuous');
prompt='Folder name: ';
Folderstr=input(prompt,'s');

prompt='CMR?: ';
CMR=input(prompt);

prompt='Hippo?: ';
Hip=input(prompt);

path1=pwd;
prompt=([path1([end-4:end-3 end-1:end end-9:end-6]) '?']);
choice=input(prompt,'s');
if strcmp(choice,'1')
str = path1([end-4:end-3 end-1:end end-9:end-6]);
else
str = choice  ;  
end

badCh = [20 27 33];
scale=0.195;
prompt='Channel_Map=: 1 for 0109, 0 for 0524/0620, 2 for 0807, 3 for 0809, 4 for 0919stroke';
ChMapNum=input(prompt);

%
%pause(18000)

[b,a]=butter(4,300/10000,'high');


%
timePerhour = 1; % in hour.
clear x
%
for i=1:numel(DIR)
    clear data
    try
data=load_open_ephys_data_faster(DIR(i).name, 'unscaledInt16');
    catch data=load_open_ephys_data(DIR(i).name);
    end
    if i==1
numOfHours = ceil(size(data,1)/timePerhour/3600/20E3);
for fol = 1:numOfHours
mkdir(num2str(fol))
end

samplePerBlock =timePerhour*3600*20E3;

    end
 
% NoisePerMin = 0;
% for Min = 1: floor(numel(data)/60/20E3)
%     Min
%   seg = data( (Min-1)*60*20E3+1:(Min)*60*20E3);
%   seg = filter(b,a,seg);
%   seg(abs(seg)>6*std(seg))=[];
%   
% NoisePerMin(Min)=std(seg);
% end



for fol = 1:numOfHours
    [i fol]
DataSeg = data((fol-1)*samplePerBlock+1:min(fol*samplePerBlock,numel(data)));



save([num2str(fol) '/' DIR(i).name(7:8)],'DataSeg','-v6')
end

% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
delete(DIR(i).name)
end
% x(:,1)=1:size(x,1);
% % manually reject several channels. 
% data([19 32 26]-15,:)=[];



%
clear data DataSeg

% numOfHours=0;
% 
% prompt='numOfHours: ';
% numOfHours=input(prompt);



% importData into cells.

for fol = 1:numOfHours
    fol
    cd(num2str(fol))
    clear amplifier_channels amplifier_data x 
    files=dir('*.mat');
    
    [b,a]=butter(4,300/10000,'high');
fprintf(1,'Filtering progress: ')

for i=1:numel(files)
    i
    if i==1
    dataS=double(importdata(files(i).name)');
    amplifier_data=zeros(32,size(dataS,2));
    amplifier_data(1,:)=filtfilt(b,a,dataS*scale);
    clear dataS
    end
    dataS = filtfilt(b,a,double(importdata(files(i).name)'*scale));
    amplifier_data(i,1:min(numel(dataS),size(amplifier_data,2))) = dataS(1:min(numel(dataS),size(amplifier_data,2)));
    clear dataS
end

amplifier_data(badCh-16,:)=[];

for i=1:numel(files)
    
    fileName = files(i).name;
    x(i,2)=str2double(fileName(1:2))-1;
end
 x(badCh-16,:)=[];
 x(:,1)=1:size(x,1);
Dat_V_Map=x;


 clear DataSeg;
% window_length = 10000;
% dend=floor(size(amplifier_data,2)/window_length)*window_length;
% amplifier_data=amplifier_data(:,1:dend);
% data_copy=amplifier_data;
% balcklist = [45 46 47];
% amplifier_data(balcklist-15,:)=[];
% amplifier_channels(balcklist-15)=[];


% data=amplifier_data(1,:);
% [b,a]=butter(4,300/10000,'high');

% med = median(amplifier_data);
% 
% data=med;
% data=amplifier_data(25,:);
%  fil=filtfilt(b,a,data);

if CMR ==1
amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);
elseif CMR==2
    list1=1:15;
    list2=16:31;
amplifier_data(list1,:) = amplifier_data(list1,:) - repmat(median(amplifier_data(list1,:)),size(amplifier_data(list1,:),1),1);
amplifier_data(list2,:) = amplifier_data(list2,:) - repmat(median(amplifier_data(list2,:)),size(amplifier_data(list2,:),1),1);
         
end


size(amplifier_data,2)/20e3/60

% find 12 mins around the mid point 
mkdir(Folderstr)
cd(Folderstr)


save('Dat_V_Map','Dat_V_Map')


amplifier_data=floor(amplifier_data*100); % 100 times bigger to avoid huge round off error;
% fid=fopen('test.dat','w');
% fwrite(fid,data,'int16');
% fclose(fid)
amplifier_data=int16(amplifier_data);
writemda(amplifier_data,[str '.mda'],'int16');
% m = matfile([str '.mat'],'Writable',true); %
% m.data=data;
clear amplifier_data
savePath=pwd;
createGeometryCSV;

save('param','CMR','ChMapNum')
cd ..
cd ..
end
%
for folder = 1:30
    cd(num2str(folder))
    delete('*.mat')
    cd ..
end
%% coordinate mountainSort

% excecute daemon command. 
unix('mp-daemon-start hanlin','-echo')
% get which mouse is this
% get which day it is 
% get what data it is and enter 
path1=pwd;
Datestr = path1([end-4:end-3 end-1:end end-9:end-6]);
     
% get all the hourly data and prepare to enter
     
     subFolderList = dir ; 
     subFolderName = {subFolderList.name};
     desiredFolder = cellfun(@(y) ~isempty(str2num(y)) , subFolderName );
     subFolderList=subFolderList(desiredFolder);
% create a folder called   mouse-day in the per day subfolder. 
mountainFolder = Datestr;
dayFolder = pwd;
mkdir(mountainFolder)
% copy pipeline files into this folder 
unix(['cp /home/hanlin/mountainlab/examples/xxxx-xxxx/pipelines.txt ./' mountainFolder])
unix(['cp /home/hanlin/mountainlab/examples/xxxx-xxxx/instructions.txt ./' mountainFolder])
unix(['cp /home/hanlin/mountainlab/examples/xxxx-xxxx/curation.script ./' mountainFolder])
fileID = fopen([ './' mountainFolder '/datasets.txt'],'w');
fclose(fileID);
mkdir([ './' mountainFolder '/datasets'])

 for subFolder =1:numel(subFolderList)
     dsName = ['ds' num2str(subFolder)];
 mkdir([ './' mountainFolder '/datasets/' dsName])
 unix(['mv ' [ './' num2str(subFolder) '/CMR/' Datestr '.mda']        [ ' ./' mountainFolder '/datasets/' dsName '/']])
 unix(['mv ' [ './' num2str(subFolder) '/CMR/geom.csv']        [ ' ./' mountainFolder '/datasets/' dsName '/']])
 unix(['mv ' [ './' num2str(subFolder) '/CMR/Dat_V_Map.mat']        [ ' ./' mountainFolder '/datasets/' dsName '/']])

 unix(['cp /home/hanlin/mountainlab/examples/xxxx-xxxx/datasets/ds1/params.json' [ ' ./' mountainFolder '/datasets/' dsName '/']])

fileID = fopen([ './' mountainFolder '/datasets.txt'],'a');
fprintf(fileID,[dsName ' datasets/' dsName '\n'])
fclose(fileID);
cd([ '' mountainFolder '/datasets/' dsName ''])
 unix(['prv-create' [ ' ' Datestr '.mda'  ' raw.mda.prv']])
 
cd ..
cd ..
cd ..
 end
     
cd(mountainFolder)
     
 for subFolder =1:numel(subFolderList)
          dsName = ['ds' num2str(subFolder)];

 unix(['kron-run msnt45 ' dsName],'-echo')
 
 end
   


%% Second Generation
 subFolderList = dir ; 

subFolderName = {subFolderList.name};
desiredFolder = cellfun(@(y) y(1)=='r' , subFolderName );
subFolderList=subFolderList(desiredFolder);
subFolderList([1])=[] 
FolderList=subFolderList;

for Folder =1:numel(FolderList)
    cd(FolderList(Folder).name)
    Datestr = FolderList(Folder).name
    cd (Datestr([end-4:end-3 end-1:end end-9:end-6]))
    Datestr=Datestr([end-4:end-3 end-1:end end-9:end-6]);
  mainFolder = pwd;
  
   subFolderList = dir ; 
     subFolderName = {subFolderList.name};
     desiredFolder = cellfun(@(y) ~isempty(str2num(y)) , subFolderName );
     subFolderList=subFolderList(desiredFolder);
     
 for subFolder =1:numel(subFolderList)
          dsName = ['ds' num2str(subFolder)];

 cd(['datasets/' dsName])



 warning off
% [~,~,raw] = xlsread('Refine.xlsx');
% load('Locations.mat')
% Sav_LocationList = list; 
% Ori_WCList = WCList;
% clear list WCList
ChMapNum=0;


% get str,dataFile, data from mda
str = Datestr;

   data = readmda([str '.mda']);
   SzData2 = size(data,2);
   SzData1 = size(data,1);

    load('Dat_V_Map.mat')


% Get CH_MAP
if ChMapNum ==1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
elseif ChMapNum==0
     load('D:\0517-hippo\ch_map_pink.mat');

elseif ChMapNum==4
    load('D:\Box Sync\0919 stroke\Finger_map.mat');
    elseif ChMapNum==3
    load('D:\0925 stroke\recording 2017-10-11\PCB_Finger.mat');
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



nChansInRawFile=size(Dat_V_Map,1);
nChansInDatFile=size(Dat_V_Map,1);
tBefore=10;
tAfter=30;


      ElecList=dir('*.csvn'); % get all electrical recording spike timing data file spec.
ImpedanceRec=zeros(32,size(ElecList,1));



for i=1:size(ElecList,1) 
book=importdata(ElecList(i).name); 

impedance=book.data(:,2);

if numel(impedance)>32
ImpedanceRec(:,i)=impedance(17:48);
else
    ImpedanceRec(:,i)=impedance;
end

end



if size(ImpedanceRec,2)~=1
AvgImpedance=median(ImpedanceRec')';
else
AvgImpedance=ImpedanceRec;    
end

%AvgImpedance(AvgImpedance>2E6)=2E6;

Ch_Map_Imp=Ch_Map_new-15;
for row=1:4
    for col=1:8
    newImpMatrix(row, col) =  AvgImpedance(Ch_Map_Imp(row,col));
    if newImpMatrix(row, col)>3E6
      newImpMatrix(row, col)= 3E6;  
    end
    end
end
newImpMatrix=newImpMatrix/1E6;


cd(mainFolder)
mdaName =['output/msnt45--' dsName '/firings.mda'];


A=readmda(mdaName);

cluster=A(3,:);
MyTimes=A(2,:);
cd(['datasets/' dsName])

color=jet(120);
color=color(5:104,:);
firstRow=[0 0 0];
color=[firstRow;color];

valid = MyTimes>(tBefore) &  MyTimes<(SzData2-tAfter);
MyTimes=MyTimes(valid);
cluster=cluster(valid);
% Mask=Mask(valid,:);
valid= (floor(cluster)==cluster)&(cluster<=numel(unique(cluster)));
MyTimes=MyTimes(valid);
cluster=cluster(valid);

list=unique(cluster);


pdfCellArray=cell(1,size(list,1));

WCList = zeros(numel(list),3);

% Prepare saving structure
mkdir('Refine_mountain')
cd ('Refine_mountain')
% oldList = cell2mat(raw(:,1));
% newList = zeros(size(oldList));

            result.date=str;
%             result.session = folderNum;
            result.Dat_V_Map=Dat_V_Map;
            result.Ch_Map=Ch_Map;
            result.Ch_Map_2=Ch_Map_2;
                        result.selectedClusters=list;

            result.waveform=cell(size(result.selectedClusters));
            result.waveformStd=cell(size(result.selectedClusters));
%             result.Intensity=cell(size(result.selectedClusters));
            result.P2P = cell(size(result.selectedClusters));
            result.Location = cell(size(result.selectedClusters));
            result.ValleyAcrossTime = cell(size(result.selectedClusters));
            result.time =cell(size(result.selectedClusters));
            result.Max_Ins_5s_FR=zeros(size(result.selectedClusters));
            result.Avg_FR=zeros(size(result.selectedClusters));
            result.newImpMatrix=newImpMatrix;
            result.AvgImpedance=AvgImpedance;
              
% 
SecValley = cell(1,numel(list));
clear newList
for clu=1:numel(list)%[3 57]
    
    proceed = 1;
    newList(clu)=1;

    
    if proceed == 1
    clu
pdfCellArray{clu}=[num2str(list(clu)) '.pdf'];
p=find(cluster==list(clu));
FiringTimeForThisUnit = MyTimes(p);
        result.time{clu}=MyTimes(p);

ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
[counts_M] = histc(ISI_in_MS_bins,0:20:1000);
[counts_L] = histc(ISI_in_MS_bins,0:1000:20000);


Sec_bins = double(FiringTimeForThisUnit)/20E3;
[counts_5s_1] = histc(Sec_bins,0:5:720);
[counts_5s_2] = histc(Sec_bins,2.5:5:722.5);
Max_Ins_5s_FR = max(max(counts_5s_1),max(counts_5s_2))/5;



clear unit %MeanPlot
for k=1:length(p);
    unit(:,:,k)=1/100*double(data(1:nChansInDatFile,MyTimes(p(k))-tBefore:MyTimes(p(k))+tAfter));
    %    plot(unit(5,:));% THIS SEEMS LIKE STRANGE, WHY GET DATA FROM ALL
    %    CHANNELS.
    %    hold on
    for ch=1:nChansInDatFile
         unit(ch,:,k)=unit(ch,:,k)-unit(ch,1,k);% first point aligned. 
    end
end

% ValleyAcrossTime = min(reshape(unit,[size(unit,1)*size(unit,2),size(unit,3)]));

MeanPlot=mean(unit,3);
StdPlot=std(unit,0,3);
% Mask_clu=Mask(cluster==list(clu),:);
% if size(Mask_clu,1)~=1
% Intensity=max(mean(Mask_clu),median(Mask_clu));
% else
% Intensity=Mask_clu;
% end    
P2P = max(MeanPlot')-min(MeanPlot');
ch=find(P2P==max(P2P));
[sortP2P,sP2P]=sort(P2P,'descend');

ValleyAcrossTime      = min(reshape(unit(ch,:,:),size(unit,2),size(unit,3)));
SecValleyAcrossTime = min(reshape(unit(sP2P(2),:,:),size(unit,2),size(unit,3)));
    SecValley{clu}=SecValleyAcrossTime;

Valley = min(MeanPlot');
 Weighted_P2P = P2P;
for i=1:4
for j=1:8
    if Ch_Map_2(i,j)~=0
    Weighted_P2P_Remap(i,j)=Weighted_P2P(Ch_Map_2(i,j));
    Valley_Remap(i,j)=floor(Valley(Ch_Map_2(i,j)));
    P2P_Remap(i,j)=P2P(Ch_Map_2(i,j));
    end
end
end
% Location Estimation Using the MUSIC method
% 

[sortedIntensity IntensityOrder]=sort(P2P,'descend');

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

% Prepare for estimation
%
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
%         load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
    
        %MeanWaveform = mean(MeanStore,3);
        % create the variable SelectedLocation=zeros(4,8) and input some
        % ones to indicate locations used for localization.
        LocationList = reshape(((Ch_Map_new-15).*SelectedLocation),32,1);
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
        
        ZAMP = MeanPlot(ampselection,:);
        [U,S,V] = svd(ZAMP);
        En=U(:,2:end);
        clear f
        f = @(x)parameterfun(x,XList,YList,ZList,En);
        
WC_Y=sum(sum(Weighted_P2P_Remap.*(Y+1)))/sum(sum(Weighted_P2P_Remap));
WC_X=sum(sum(Weighted_P2P_Remap.*(X+1)))/sum(sum(Weighted_P2P_Remap));
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



    if 2>1 % case 1 =, new unit  strSearchResult==1




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

WCList(clu,:)= OverallMinCorD(1,1:3);
    else

    end
    
    
    


t=[-tBefore/20e3:1/20e3:tAfter/20e3]*1000;
% 4-by-8 Channel Plot showing correlation by visual inspection
for i=1:size(MeanPlot,1)
    Valley =Valley -MeanPlot(i,1);
    MeanPlot(i,:)=MeanPlot(i,:)-MeanPlot(i,1);
end

h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
set(h, 'Visible', 'off');
Ox=3;
Oy=8.5;
gap=0.4;
Gw=1.5;
Gh=2.25;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);


for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
           title(num2str(newImpMatrix(row,col)))
        if Ch_Map(row,col)~=0
            temp=find(Dat_V_Map(:,2)==Ch_Map(row,col));
            errorbar(t,MeanPlot(temp,:),StdPlot(temp,:),'Color',[0.9 0.9 0.9])

            hold on;
            if min(min(MeanPlot))<-100 || max(max(MeanPlot))>25
                plot(t,MeanPlot(temp,:),'Color','r','LineWidth',2)
%                 axis([t(1) t(end) -150 50])
                axis([t(1) t(end) min(min(MeanPlot)) max(max(MeanPlot))])

            else
                plot(t,MeanPlot(temp,:),'Color','r','LineWidth',1)
%                 axis([t(1) t(end) -100 25])
                axis([t(1) t(end) min(min(MeanPlot)) max(max(MeanPlot))])

            end
            set(gca,'xticklabel',[])
                       title(num2str(newImpMatrix(row,col)))

        end


    end
end

  Ox=3;
 Oy=20.5;
gap=0.3;
Gw=1.5;
Gh=1.5;
if ChMapNum~=1

gap=0.8;
        Gw=1.2;
        Gh=1.2;
end
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);
% Weighted Center Calculation
WC_Y=sum(sum(Weighted_P2P_Remap.*(Y+1)))/sum(sum(Weighted_P2P_Remap));
WC_X=sum(sum(Weighted_P2P_Remap.*(X+1)))/sum(sum(Weighted_P2P_Remap));
% MUSIC Location Calculation

for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
        
        set(ax,'Box','on')

        set(ax,'LineWidth',2)
        plot(0,0)
        T = text(-0.2,0,num2str(floor(abs(P2P_Remap(row,col)))));
        set(T,'FontSize',6)
        
        set(T,'Color',[0 0 1])
        
        set(ax,'xtick',[])
        set(ax,'ytick',[])
        set(ax,'XColor',[0.7 0.7 0.7])
        set(ax,'YColor',[0.7 0.7 0.7])
        
        if SelectedLocation(row,col)==1
        set(ax,'XColor',[1 0 0])
        set(ax,'YColor',[1 0 0])
        end
        
        end
end

scaleCoeff = 1.5/25;
if ChMapNum~=1

scaleCoeff = 1.2/30;

end
clear  Cselected5
Cselected5=OverallMinCorD;
for center= 1:size(Cselected5,1)
ax=axes('Units','centimeters','position', [Cselected5(center,1)*scaleCoeff-0.5+Ox Cselected5(center,2)*scaleCoeff+Oy-0.5 1 1]);
ZSize = 60;
RDistance = sqrt(sum((Cselected5(center,:)-[highX highY 0]).^2));
RDistance(RDistance>60)=60;
NorDistance = RDistance/60*100; % full color specified for 2 spacing. 30um
scatter(0,0,30+(ZSize)^2,'+','LineWidth',2,'MarkerEdgeColor',color(101-floor(NorDistance),:))
fx=floor((WC_X-Ox)*scaleCoeff*10)/10;
fy=floor((WC_Y-Oy)*scaleCoeff*10)/10;


% txt1 = ['       ' num2str(floor(Cselected5(center,3))) ];
% T=text(0,0,txt1);
% set(T,'FontWeight','bold')
set(T,'FontSize',5)

set(ax,'Visible','off')
end
 ax=axes('Units','centimeters','position', [0.6 0.6 5 4.5]);
        set(ax,'FontSize',6)
        if ~isempty(counts)
bar(0.5:1:50.5,counts)
xlim([0 51])
        end

  ax=axes('Units','centimeters','position', [3 5.7 16 2.5]);
        set(ax,'FontSize',4)
 scatter(double(FiringTimeForThisUnit)/20E3/60,ValleyAcrossTime,5*ones(size(ValleyAcrossTime)),'+')
   
 
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

ax=axes('Units','centimeters','position', [17.5 0.6 5 4.5]);
        set(ax,'FontSize',6)
T=text(.1,.1,[str '-----' num2str(list(clu))  ]); %
T=text(.1,.2,[num2str(sum(cluster==list(clu))) '-----' num2str(numel(cluster)) ]); %
T=text(.1,.3,[num2str(sum(cluster==list(clu))/(SzData2/20E3)) '  spk/sec']); %
T=text(.1,.4,['Max5s = ' num2str(Max_Ins_5s_FR) '  spk/sec']); %



set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')




set(gcf, 'Position', get(0,'Screensize'));
set(h,'PaperPosition',[0 0 53 27])
print(h,'-dpdf',[num2str(list(clu)) '.pdf'])
% print(h,'-depsc',[num2str(list(clu)) '.eps'])




        StandardMeanPlot=zeros(32,size(MeanPlot,2));
        StandardStdPlot=zeros(32,size(MeanPlot,2));
%         StandardIntensity=zeros(1,32);
        StandardP2P=zeros(1,32);
        
        for ch=1:32
        if ~isempty(find((Dat_V_Map(:,2)-15)==ch))
        StandardMeanPlot(ch,:)=MeanPlot(find(Dat_V_Map(:,2)==ch+15),:);
        StandardStdPlot(ch,:)=StdPlot(find(Dat_V_Map(:,2)==ch+15),:);
%         StandardIntensity(ch)=Intensity(find(Dat_V_Map(:,2)==ch+15));
        StandardP2P(ch)=P2P(find(Dat_V_Map(:,2)==ch+15));
        end
        end
        
        result.waveform{clu}=StandardMeanPlot;
%         result.Intensity{clu}=StandardIntensity;
        result.P2P{clu}=StandardP2P;
        result.waveformStd{clu}=StandardStdPlot;
        result.Location{clu}=WCList(clu,:);
        result.Max_Ins_5s_FR(clu)= Max_Ins_5s_FR;
        result.Avg_FR(clu)= sum(cluster==list(clu))/(SzData2/20E3);
        result.ValleyAcrossTime{clu}=ValleyAcrossTime;

    end
end
clear data
        result.ValleyAcrossTime = result.ValleyAcrossTime(newList>0);
        result.waveform=result.waveform(newList>0);
%         result.Intensity=result.Intensity(newList>0);
        result.P2P=result.P2P(newList>0);
        result.waveformStd=result.waveformStd(newList>0);
        result.Location= result.Location(newList>0);
        result.Max_Ins_5s_FR= result.Max_Ins_5s_FR(newList>0);
        result.Avg_FR= result.Avg_FR(newList>0);
        result.time= result.time(newList>0);

list = list (newList>0);
WCList = WCList(newList>0,:);
pdfCellArray = pdfCellArray(newList>0);
newList = newList(newList>0);
            result.list=list;
            result.newList=newList;
            
            
            
            
% save([result.date '_tracking'],'result')
% save([mainFolderPath '\' result.date '_tracking'],'result')
save('Locations','list','newList','WCList')

        Gw=25;
        Gh=25;
        numGx = 8;
        numGy = 4;
        gap=5;
        
if ChMapNum~=1
Gw = 30;
Gh= 30;
        numGx = 8;
        numGy = 4;
        gap=20;
        


end
clear WCList1
WCList1(:,1)=WCList(:,1)-0.5*(numGx*Gw+(numGx-1)*gap);
WCList1(:,2)=WCList(:,2)-0.5*(numGy*Gh+(numGy-1)*gap);

[theta,rho] = cart2pol(WCList1(:,1),WCList1(:,2));
[Y,I] = sort(theta,'descend');
% alternatively we could use this TSPNN
userConfig = struct('xy',WCList1(:,1:2)); 
resultStruct = tsp_nn(userConfig);
I = resultStruct.optRoute;
pdfCellArray=pdfCellArray(I);
% newStr = [str(5:8) '-' str(1:2) '-' str(3:4)]
newStr =str;
%delete(pdfCellArray{:})
%
result.waveform=result.waveform(I)';
%         result.Intensity=result.Intensity(I);
        result.P2P=result.P2P(I)';
        result.waveformStd=result.waveformStd(I)';
        result.Location= result.Location(I)';
        result.Max_Ins_5s_FR= result.Max_Ins_5s_FR(I)';
        result.Avg_FR= result.Avg_FR(I)';
        result.time= result.time(I)';
        result.ValleyAcrossTime=result.ValleyAcrossTime(I)';


            result.list=list(I)';
            result.newList=newList(I)';
            %
save([result.date '_tracking'],'result')
% save([mainFolderPath '\' result.date '_tracking'],'result')

append_pdfs([newStr '.pdf'],pdfCellArray{:})
order = list(I);
save('Y','order');

close all

cd(mainFolder)
 end
 cd ..
  cd ..
end