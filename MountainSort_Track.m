% fileList = dir('*.mat');
% 
% day1sum = load(fileList(1).name);
% day2sum = load(fileList(2).name);
% day1sum=day1sum.result;
% day2sum=day2sum.result;
% 
% simMatrix = zeros(numel(day1sum.waveform),numel(day2sum.waveform));
% for i=1:size(simMatrix,1)
%     for j=1:size(simMatrix,2)
%      [i,j]   
%     simMatrix(i,j)= sum(sum((day1sum.waveform{i} -day2sum.waveform{j}).^2));
%     
%     end
% end
% Match_Matrix_row = zeros(numel(day1sum.waveform),numel(day2sum.waveform));
% for i=1:size(Match_Matrix_row,1)
%     Match_Matrix_row(i,find(simMatrix(i,:)==min(simMatrix(i,:))))=1;
% end
% 
% Match_Matrix_col = zeros(numel(day1sum.waveform),numel(day2sum.waveform));
% for j=1:size(Match_Matrix_col,2)
%     Match_Matrix_col(find(simMatrix(:,j)==min(simMatrix(:,j))),j)=1;
% end
% 
% final_match = Match_Matrix_col.*Match_Matrix_row;
% [X Y]=find(final_match ==1);
% load('ch_map_pink.mat')
% Ch_Map = Ch_Map_new-15;
% 
% for wav=1:numel(X)
%     h=figure
% 
% for ch=1:32
% [X_off Y_off]=find(Ch_Map==ch);
% X_off=X_off*50;
% Y_off=Y_off*100;
% plot([1:41]+X_off,Y_off+day1sum.waveform{X(wav)}(ch,:),'b');
% hold on
% plot([1:41]+X_off,Y_off+day2sum.waveform{Y(wav)}(ch,:),'r');
% end
% 
% 
% print(h,'-dpdf',[num2str(wav) '.pdf'])
% 
% end
% 
%%
DataForFeature = store1.DataForFeature;
P2PForFeature = store1.P2PForFeature;
%% Press Last 
for i=1:numel(DataForFeature)
    i
[AmpList,P_index] = sort(P2PForFeature{i},'desc');
weight = 1:numel(P_index);
weight = 1 - sigmf(weight',[0.6 10]);
finalweight = zeros(size(weight));
for w = 1:numel(weight)
finalweight(P_index(w))=weight(w);
end
DataForFeature{i}=repmat(finalweight,1,size(DataForFeature{i},2)).*DataForFeature{i};
DataForFeature{i}(DataForFeature{i}>0)=0;
DataForFeature{i}=DataForFeature{i}(:,2:2:end);
end
LocForFeature = cell2mat(store1.LocForFeature);
MSEMatrix = zeros(numel(DataForFeature));
LocMatrix = zeros(numel(DataForFeature));
lookUpTable = store1.lookUpTable;


% SpearsMatrix = zeros(numel(DataForFeature));
count=0;
total=nchoosek(numel(DataForFeature),2);

 fprintf(1,'progress: ')
for i=1:numel(DataForFeature)-1
    for j=i+1:numel(DataForFeature)
    count=count+1
%     s1= sum(sum((DataForFeature{i}(:,1:37) -DataForFeature{j}(:,3:39)).^2));
%     s2= sum(sum((DataForFeature{i}(:,2:38) -DataForFeature{j}(:,3:39)).^2));
%     s3= sum(sum((DataForFeature{i}(:,3:39) -DataForFeature{j}(:,3:39)).^2));
%     s4= sum(sum((DataForFeature{i}(:,4:40) -DataForFeature{j}(:,3:39)).^2));
%     s5= sum(sum((DataForFeature{i}(:,5:41) -DataForFeature{j}(:,3:39)).^2));
    slide = 5;
    square = 0;
    slide_start = slide+1;
    slide_end = size(DataForFeature{i},2)-slide ; 
    slide_plus = (slide_end - slide_start);
    S_array = zeros(1,slide*2+1);
    for s=1:numel(S_array)
        if square == 1
     S_array(s)= sum(sum((DataForFeature{i}(:,s:s+slide_plus) -DataForFeature{j}(:, slide_start:slide_end)).^2));
        
        else
     S_array(s)= sum(sum(abs((DataForFeature{i}(:,s:s+slide_plus) -DataForFeature{j}(:, slide_start:slide_end)))));
        end
        
     end

    MSEMatrix(i,j) = min(S_array);
    MSE_debug{i,j} = S_array;
    LocMatrix(i,j)= pdist([LocForFeature(i,1:2);LocForFeature(j,1:2)],'euclidean');
    end
end

%% Match For One Round
totalSession = numel(unique(lookUpTable(:,1)));
AllsesCount = []
AllMatchRecord = [] 
for thisSes = 1:totalSession-1
    nextSes = thisSes+1 
    thisSesCandidateList  = (lookUpTable(:,1)==thisSes  & remain_index' ) ;
    
    nextSesCandidateList  = (lookUpTable(:,1)==nextSes  & remain_index' ) ;
    AllsesCount = [AllsesCount numel(find(thisSesCandidateList>0))];
    if thisSes==totalSession-1
         AllsesCount = [AllsesCount numel(find(nextSesCandidateList>0))];
    end
simMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simMatrix(thisSesCandidateList,nextSesCandidateList)=MSEMatrix(thisSesCandidateList,nextSesCandidateList);
simLocMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simLocMatrix(thisSesCandidateList,nextSesCandidateList)=LocMatrix(thisSesCandidateList,nextSesCandidateList);

Match_Matrix_row = zeros(numel(DataForFeature),numel(DataForFeature));
consider = 3; 
for i=1:size(Match_Matrix_row,1)
    if min(simMatrix(i,:))~=Inf
    [listMSE LMSE_index]=sort(simMatrix(i,:));
    
    LocList = simLocMatrix(i,LMSE_index(1:consider));
    [LocMin LocMin_index ] = sort(LocList);
    
    selected = LMSE_index(LocMin_index(1));
        Match_Matrix_row(i,selected)=1;
    
    end
end


Match_Matrix_col = zeros(numel(DataForFeature),numel(DataForFeature));
for j=1:size(Match_Matrix_col,2)
    if min(simMatrix(:,j))~=Inf
    [listMSE LMSE_index]=sort(simMatrix(:,j));
    LocList = simLocMatrix(LMSE_index(1:consider),j);
    [LocMin LocMin_index ] = sort(LocList);
    selected = LMSE_index(LocMin_index(1));
    Match_Matrix_col(selected,j)=1;
    end
end

final_match = Match_Matrix_col.*Match_Matrix_row;
[X,Y]=find(final_match ==1);
AllMatchRecord=[AllMatchRecord; [X Y]];
end

totalPossiblePairs =  0 
totalPossibleUnit = 0  
for i=1:numel(AllsesCount)-1
totalPossiblePairs = totalPossiblePairs+ min(AllsesCount(i),AllsesCount(i+1))
if i == numel(AllsesCount)-1
totalPossibleUnit = totalPossiblePairs+ min(AllsesCount(i),AllsesCount(i+1))
end 
end
%% Match For Two Round
totalSession = numel(unique(lookUpTable(:,1)));

AllMatchRecord = [] 
for thisSes = 1:totalSession-1
    nextSes = thisSes+1 
    
    thisSesCandidateList  = (lookUpTable(:,1)==thisSes  & remain_index' ) ;
    
    nextSesCandidateList  = (lookUpTable(:,1)==nextSes  & remain_index' ) ;

simMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simMatrix(thisSesCandidateList,nextSesCandidateList)=MSEMatrix(thisSesCandidateList,nextSesCandidateList);
simLocMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simLocMatrix(thisSesCandidateList,nextSesCandidateList)=LocMatrix(thisSesCandidateList,nextSesCandidateList);

for round = 1:7
Match_Matrix_row = zeros(numel(DataForFeature),numel(DataForFeature));
consider = 3; 
for i=1:size(Match_Matrix_row,1)
    if min(simMatrix(i,:))~=Inf
    [listMSE LMSE_index]=sort(simMatrix(i,:));
    
    LocList = simLocMatrix(i,LMSE_index(1:consider));
    [LocMin LocMin_index ] = sort(LocList);
    
    selected = LMSE_index(LocMin_index(1));
        Match_Matrix_row(i,selected)=1;
    
    end
end


Match_Matrix_col = zeros(numel(DataForFeature),numel(DataForFeature));
for j=1:size(Match_Matrix_col,2)
    if min(simMatrix(:,j))~=Inf
    [listMSE LMSE_index]=sort(simMatrix(:,j));
    LocList = simLocMatrix(LMSE_index(1:consider),j);
    [LocMin LocMin_index ] = sort(LocList);
    selected = LMSE_index(LocMin_index(1));
    Match_Matrix_col(selected,j)=1;
    end
end

final_match = Match_Matrix_col.*Match_Matrix_row;
[X,Y]=find(final_match ==1);
AllMatchRecord=[AllMatchRecord; [X Y]];


simMatrix(X,:)=Inf;
simMatrix(:,Y)=Inf;

simLocMatrix(X,:)=Inf;
simLocMatrix(:,Y)=Inf;


end
end


%% Binary Linear Programming
totalSession = numel(unique(lookUpTable(:,1)));

AllMatchRecord = [] 
for thisSes = 1:totalSession-1
    nextSes = thisSes+1 
    
    thisSesCandidateList  = (lookUpTable(:,1)==thisSes  & remain_index' ) ;
    
    nextSesCandidateList  = (lookUpTable(:,1)==nextSes  & remain_index' ) ;

simMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simMatrix(thisSesCandidateList,nextSesCandidateList)=MSEMatrix(thisSesCandidateList,nextSesCandidateList);
simLocMatrix = Inf(numel(DataForFeature),numel(DataForFeature));
simLocMatrix(thisSesCandidateList,nextSesCandidateList)=LocMatrix(thisSesCandidateList,nextSesCandidateList);

n1=numel(find(thisSesCandidateList>0));
n2=numel(find(nextSesCandidateList>0));

f=zeros(1,n1*n2);
intcon = n1*n2;
UnitRecord = zeros(2,n1*n1);
A=zeros(n1,n1*n2);
b=ones(n1,1);
lb=zeros(n1*n2,1);
ub=ones(n1*n2,1);

numZeroRow  = 1;
for i=1:size(Match_Matrix_row,1)
    if min(simMatrix(i,:))~=Inf
    
        f(((numZeroRow-1)*n2+1):(numZeroRow*n2))= simMatrix(i,simMatrix(i,:)~=Inf);
        
        UnitRecord(1,((numZeroRow-1)*n2+1):(numZeroRow*n2)) = i; 
        UnitRecord(2,((numZeroRow-1)*n2+1):(numZeroRow*n2)) = find(simMatrix(i,:)~=Inf); 
        
        A(numZeroRow,((numZeroRow-1)*n2+1):(numZeroRow*n2))= 1;
             
   numZeroRow = numZeroRow +1;
    end
end

selected_index = intlinprog(-f,intcon,A,b,[],[],lb,ub);


final_match = Match_Matrix_col.*Match_Matrix_row;
[X,Y]=find(final_match ==1);
AllMatchRecord=[AllMatchRecord; [X Y]];
end

%%

Dist_Of_Matched_record = arrayfun(@(x) LocMatrix(AllMatchRecord(x,1),AllMatchRecord(x,2)), 1:size(AllMatchRecord,1));
MSE_Of_Matched_record = arrayfun(@(x) MSEMatrix(AllMatchRecord(x,1),AllMatchRecord(x,2)), 1:size(AllMatchRecord,1));
%%
selected_wav  = find(MSE_Of_Matched_record>200 & MSE_Of_Matched_record<=240);
selected_wav  = find(Dist_Of_Matched_record>30 & Dist_Of_Matched_record<=40);
selected_wav  = [50 102 135 233 574 788];% 1:820; %AllMatchRecord;
clear F OKMatrix
for wav =1:numel(selected_wav)
    wav
h=figure('OuterPosition',[500   300   800   600]);

MeanPlot =  DataForFeature{AllMatchRecord(selected_wav(wav),1)};
MeanPlot2 = DataForFeature{AllMatchRecord(selected_wav(wav),2)};

set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
 set(h, 'Visible', 'on');
Ox=1;
Oy=1;
gap=0.4;
Gw=1.5;
Gh=2.25;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

load('ch_map_pink.mat');
Ch_Map=Ch_Map_new-15;
for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
         
        if Ch_Map(row,col)~=0
            temp=Ch_Map(row,col);
    
           
                plot(MeanPlot(temp,:),'Color','r','LineWidth',2)
                hold on ;
                plot(MeanPlot2(temp,:),'Color','b','LineWidth',2)
                
               

           
            
            set(gca,'xticklabel',[])
             

        end
           ylim([min(min(min(MeanPlot)),min(min(MeanPlot2))) max(max(max(MeanPlot)),max(max(MeanPlot2)))])

    end
end
% AllMatchRecord(selected_wav(wav),1);
% AllMatchRecord(selected_wav(wav),2);
title( mat2str(floor(MSE_debug{AllMatchRecord(selected_wav(wav),1),AllMatchRecord(selected_wav(wav),2)}/10)))
 ax=axes('Units','centimeters','position', [5 12 1 1]);


 u1 = AllMatchRecord(selected_wav(wav),1);
 u2 = AllMatchRecord(selected_wav(wav),2);
 
 text(0,0,['Loc: ' num2str(floor(LocMatrix(u1,u2)))  '  MSE: '  num2str(floor(MSEMatrix(u1,u2))) '  Pair ' mat2str([selected_wav(wav) u1 u2])])
 axis off




 OK = input('Unit Ok?');
OKMatrix(selected_wav(wav))=OK;
F(wav) = getframe(h);   
close(h)
end
  

v = VideoWriter('strange_track.avi');
% movie(F)
v.FrameRate=1;
open(v)
writeVideo(v,F);
close(v)
%% Plot Training Set
scatter(Dist_Of_Matched_record(OKMatrix==1),MSE_Of_Matched_record(OKMatrix==1),'r.')
hold on
scatter(Dist_Of_Matched_record(OKMatrix==2),MSE_Of_Matched_record(OKMatrix==2),'b.')
hold on
scatter(Dist_Of_Matched_record(OKMatrix==3),MSE_Of_Matched_record(OKMatrix==3),'k.')
xlabel('MAE')
ylabel('Distance')
legend('Valid','Not Sure','Invalid')
%% Plot
P2PForFeature=store1.P2PForFeature;
P2P = max(cell2mat(P2PForFeature)');



% AllMatchRecord(Dist_Of_Matched_record>60,:)=[]; 
recordMatrix = NaN(size(AllMatchRecord,1),totalSession); % assume each unit can be matched across all sessions. 
for rec = 1:size(AllMatchRecord,1)
    
   [x,y]=find(recordMatrix== AllMatchRecord(rec,1));
   if isempty(x) % we don't see this record before, write it to the next non-empty line
   [xN,yN]=find(isnan(recordMatrix)~=1);
      if ~isempty(xN)
       recordMatrix(max(xN)+1, lookUpTable(AllMatchRecord(rec,1),1))=AllMatchRecord(rec,1);
       recordMatrix(max(xN)+1, lookUpTable(AllMatchRecord(rec,2),1))=AllMatchRecord(rec,2);
      else   
       recordMatrix(0+1, lookUpTable(AllMatchRecord(rec,1),1))=AllMatchRecord(rec,1);
       recordMatrix(0+1, lookUpTable(AllMatchRecord(rec,2),1))=AllMatchRecord(rec,2);
       end
   else   % we have seen this record before, write it after last occurance
       recordMatrix(x, lookUpTable(AllMatchRecord(rec,1),1))=AllMatchRecord(rec,1);
       recordMatrix(x, lookUpTable(AllMatchRecord(rec,2),1))=AllMatchRecord(rec,2);
   end
end
spreadMatrix  = recordMatrix;
spreadMatrix(isnan(min(spreadMatrix')),:)=[];
spreadMatrix(isnan(spreadMatrix))=-1000;
clear spreadMatrix_sort a b emp
for rec =1:size(spreadMatrix,1)
    emp(rec) = sum(spreadMatrix(rec,:)==-1000);
[a , b] =sort(emp) ;   

end
spreadMatrix_sort =spreadMatrix(b,:);
%%
spreadMat = spreadMatrix;
for i=1:size(spreadMat,1)
    for j=1:size(spreadMat,2)
    if spreadMatrix(i,j)==-1000
    spreadMat(i,j)=-100;
    else
        spreadMat(i,j) = P2P(spreadMatrix(i,j));
        
    end
    end
end
%% spreadMat = spreadMatrix;
for i=1:size(spreadMat,1)
    for j=1:size(spreadMat,2)
    if spreadMatrix_sort(i,j)==-1000
    spreadMat(i,j)=-100;
    else
        spreadMat(i,j) = P2P(spreadMatrix_sort(i,j));
        
    end
    end
end
%%
figure
imagesc(spreadMatrix)
set(gca,'Ydir','Normal')
colormap jet
figure
imagesc(spreadMatrix_sort)
set(gca,'Ydir','Normal')
colormap jet

figure
imagesc(spreadMat)
set(gca,'Ydir','Normal')
colormap jet
%% GiveVideos of All Tracked Record 
for record = 1:50%numel(spreadMatrix_sort)
clear F
selectedList = spreadMatrix_sort(record,spreadMatrix_sort(record,:)~=-1000);
for ses = 1:numel(selectedList)

h=figure('OuterPosition',[-5    35   714   540]);
MeanPlot = DataForFeature{selectedList(ses)};
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
 set(h, 'Visible', 'off');
Ox=1;
Oy=1;
gap=0.4;
Gw=1.5;
Gh=2.25;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

load('ch_map_pink.mat');
Ch_Map=Ch_Map_new-15;
for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')
         
        if Ch_Map(row,col)~=0
            temp=Ch_Map(row,col);
    
           
                plot(MeanPlot(temp,:),'Color','r','LineWidth',2)
                
               

           
            
            set(gca,'xticklabel',[])
             

        end
           axis( [ 1 41 min(min(MeanPlot)) max(max(MeanPlot))])

    end
end
 title([num2str(lookUpTable(selectedList(ses),1)) '---'  num2str(selectedList(ses)) ]  )
% saveas(h,[h.Name '.png'])
F(ses) = getframe(h);    
close(h)
end
v = VideoWriter([num2str(record) 'track.avi'] );
% movie(F)
v.FrameRate=1;
open(v)
writeVideo(v,F);
close(v)
end