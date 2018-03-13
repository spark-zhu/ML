% cross_correlation_time_stamp_based
% data assumed to be in allrecords-by-2 matrix
clear all 
close all
MyTimes=h5read('2017-01-1310mins.mat.kwik','/channel_groups/0/spikes/time_samples');
cluster=h5read('2017-01-1310mins.mat.kwik','/channel_groups/0/spikes/clusters/main');
data(:,1)=double(cluster);
data(:,2)=double(MyTimes)/20E3;
% clust=[]; 8-21-16 
clust=[37 14 26 40 28 38 17 30]; % This line needs to be adjusted accordingly for new dataset.
newdata=zeros(1,2);
for i=1:numel(data(:,1))
    if ismember(data(i,1),clust)==1
    newdata(end+1,:)=data(i,:);
    end
end
newdata=newdata(2:end,:);

%channel_list=unique(data(:,1)); version 1
 channel_list=clust; % version 2
 data=newdata;       % version 2
 %%
binSize=0.001; % unit in seconds, so this is 2ms. 
%
[aAndB aInds bInds] = d_uIntersect(result.list,result.selectedClusters);
track_data = result.time;
track_data=track_data(aInds,:);%noise rejection

track_data = cellfun(@(x) double(x)/20E3, track_data,'UniformOutput' ,0);
%% adjust size of the data ch to match the bin size;
bilateral_lag=0.03;
offsetSample=bilateral_lag/binSize;

minTime=min(cell2mat(track_data'))-bilateral_lag-binSize;
maxTime=max(cell2mat(track_data'))+bilateral_lag+binSize;

howmanyBins=floor((maxTime-minTime)/binSize);
binCountSummary=zeros(numel(track_data),howmanyBins);

%% go over all the records. 
for i=1:numel(track_data)
    binCountSummary(i,ceil((track_data{i}-minTime)/binSize))=1;
end

CrossCorrMatrix = cell(numel(track_data),numel(track_data));
CrossPeakMatrix = zeros(numel(track_data),numel(track_data));
ChiMatrix = ones(numel(track_data),numel(track_data));
for i=1:numel(track_data)
for j=i+1:numel(track_data)
    [i j]
CrossCorrMatrix{i,j}=xcorr(binCountSummary(i,:),binCountSummary(j,:),bilateral_lag/binSize);


try 
    var=CrossCorrMatrix{i,j};
    var=var(var>0);
CrossPeakMatrix(i,j)=max(var)/median(var);
catch
end


try 
var = CrossCorrMatrix{i,j};
ele = bilateral_lag/binSize;
left = sum(var(1:ele));
right = sum(var(ele+2:end));
if left>5&&right>5
total=left + right;
Expect = total/2;
if Expect>5
Chi = sum(([left right]-Expect).^2)/Expect;
ChiMatrix(i,j)=Chi;
end
end
catch
end


end
end
%%
Ans=ChiMatrix>15;

for i=1:numel(track_data)
for j=i+1:numel(track_data)
if Ans(i,j)
figure
bar(CrossCorrMatrix{i,j})
title(num2str(ChiMatrix(i,j)))
end
end
end

AllCompare = reshape(ChiMatrix,size(CrossPeakMatrix,1)*size(CrossPeakMatrix,2),1);
AllCompare=AllCompare(AllCompare~=0);
AllCompare(AllCompare>50)=50;
hist(AllCompare,100)
[x ,y]=ecdf(AllCompare);
plot(y,x)

for i=1:numel(track_data)
for j=i+1:numel(track_data)
    var=CrossCorrMatrix{i,j};
    var(ele+1)=[];
    var=var(var>0);
CrossPeakMatrix(i,j)=max(var)/median(var);
end
end



% CrossPeakMatrix=CrossPeakMatrix+CrossPeakMatrix';
AllCompare = reshape(CrossPeakMatrix,size(CrossPeakMatrix,1)*size(CrossPeakMatrix,2),1);
AllCompare=AllCompare(AllCompare~=0);
AllCompare=AllCompare(AllCompare<=10);
hist(AllCompare,100)
[x ,y]=ecdf(AllCompare);
plot(y,x)
NF_cutoff = 2.2;
%%


for i=1:numel(track_data)
for j=i+1:numel(track_data)
if CrossPeakMatrix(i,j)>NF_cutoff && max(CrossCorrMatrix{i,j})>60 % ChiMatrix(i,j)>15 
h=figure
bar(CrossCorrMatrix{i,j})
title([num2str(i) '--' num2str(j)])
% saveas(h,[num2str(i) '--' num2str(j)],'png')
end
end
end

timeLag_In_Samples=(-(bilateral_lag/binSize):(bilateral_lag/binSize))*binSize*1000;
%% Position Plot
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
rgb=jet(20);
for i=1:numel(track_data)
for j=i+1:numel(track_data)
if CrossPeakMatrix(i,j)>NF_cutoff && ChiMatrix(i,j)>15 && max(CrossCorrMatrix{i,j})>60
    
var = CrossCorrMatrix{i,j};
ele = bilateral_lag/binSize;
left = sum(var(1:ele));
right = sum(var(ele+2:end));

imbalance = max(left,right)/(left+right);
imbalance =max(left,right)/min(left,right);
    if right>left
davinci( 'arrow', 'X', [Cselected5(i,1) Cselected5(j,1)], 'Y', [Cselected5(i,2) Cselected5(j,2)],'FaceColor',rgb(ceil((imbalance)*5),:),'LineWidth',1.5,'Color',rgb(ceil((imbalance)*5),:),'Shaft.Width',2,'Head.Width',5,'Head.Sweep',0.5,'Shaft.Width',1) 
    else
        davinci( 'arrow', 'X', [Cselected5(j,1) Cselected5(i,1)], 'Y', [Cselected5(j,2) Cselected5(i,2)],'FaceColor',rgb(ceil((imbalance)*5),:),'LineWidth', 1.5,'Color',rgb(ceil((imbalance)*5),:),'Shaft.Width',2,'Head.Width',5,'Head.Sweep',0.5,'Shaft.Width',1) 

    end
hold on 
scatter(Cselected5(i,1),Cselected5(i,2),100,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',1)
hold on 
scatter(Cselected5(j,1),Cselected5(j,2),100,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',1)

end
end
end
colorbar
colormap(rgb);
caxis([0 4])

%%  Direct Connection 
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
rgb=jet(10);
for i=1:numel(track_data)
for j=i+1:numel(track_data)
if CrossPeakMatrix(i,j)>NF_cutoff && ChiMatrix(i,j)>15 && max(CrossCorrMatrix{i,j})>60
    
var = CrossCorrMatrix{i,j};
ele = bilateral_lag/binSize;
left = max(var(1:ele));
right = max(var(ele+2:end));
middle = var(ele+1);


imbalance = CrossPeakMatrix(i,j);
    if right==max([right left])
davinci( 'arrow', 'X', [Cselected5(i,1) Cselected5(j,1)], 'Y', [Cselected5(i,2) Cselected5(j,2)],'FaceColor',rgb(ceil((imbalance)),:),'LineWidth',1.5,'Color',rgb(ceil((imbalance)),:),'Shaft.Width',2,'Head.Width',5,'Head.Sweep',3,'Shaft.Width',1,'FaceAlpha',0.5) 
    delay = find(right==var(ele+2:end));
%     text(0.5*(Cselected5(i,1)+Cselected5(j,1)),0.5*(Cselected5(i,2)+Cselected5(j,2)),num2str(delay))
    end
    if left==max([right left])
        davinci( 'arrow', 'X', [Cselected5(j,1) Cselected5(i,1)], 'Y', [Cselected5(j,2) Cselected5(i,2)],'FaceColor',rgb(ceil((imbalance)),:),'LineWidth', 1.5,'Color',rgb(ceil((imbalance)),:),'Shaft.Width',2,'Head.Width',5,'Head.Sweep',3,'Shaft.Width',1,'FaceAlpha',0.5) 
    delay = find(left==var(ele:-1:1));
%     text(0.5*(Cselected5(i,1)+Cselected5(j,1)),0.5*(Cselected5(i,2)+Cselected5(j,2)),num2str(delay))
    end
    
%     if middle==max([right left middle])
%         davinci( 'arrow', 'X', [Cselected5(j,1) Cselected5(i,1)], 'Y', [Cselected5(j,2) Cselected5(i,2)],'ArrowType','double','FaceColor',rgb(ceil((imbalance)),:),'LineWidth', 1.5,'Color',rgb(ceil((imbalance)),:),'Shaft.Width',2,'Head.Width',5,'Head.Sweep',3,'Head.Length',4,'Shaft.Width',1) 
% 
%     end
hold on 
scatter(Cselected5(i,1),Cselected5(i,2),100,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',1)
text(Cselected5(i,1),Cselected5(i,2),num2str(i),'FontSize',14)

hold on 
scatter(Cselected5(j,1),Cselected5(j,2),100,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',1)
text(Cselected5(j,1),Cselected5(j,2),num2str(j),'FontSize',14)
end
end
end
colorbar
colormap(rgb);
caxis([0 10])
%%


rgb=redblue(40);
for center= 1:size(Cselected5,1)
% ax=axes('Units','centimeters','position', [Cselected5(center,1)*scaleCoeff-0.5+Ox Cselected5(center,2)*scaleCoeff+Oy-0.5 1 1]);

% RDistance = sqrt(sum((Cselected5(center,:)-[highX highY 0]).^2));
% RDistance(RDistance>60)=60;
% NorDistance = RDistance/60*100; % full color specified for 2 spacing. 30um
% if center~=size(Cselected5,1)
% scatter(Cselected5(center,1),Cselected5(center,2),2000,'.','LineWidth',2,'MarkerEdgeColor',rgb(floor(CorrelationUnit(center)*100+20),:))

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
%% plot
figure
for i=1:numel(track_data)
for j=1:numel(track_data)
    [i j]
subplot(numel(track_data),numel(track_data),(i-1)*numel(track_data)+j)
bar(timeLag_In_Samples,CrossCorrMatrix{i,j})
ax=gca;
ax.XTick = [-10 0 10];
% ylim([0 0.25])
% xlim([min(timeLag_In_Samples) max(timeLag_In_Samples)]);
% title([num2str(channel_list(i)) '  -  '  num2str(channel_list(j))])
hold on;
end
end

XCorReuslt = zeros(numel(track_data),numel(track_data));
for i=1:numel(track_data)
for j=i+1:numel(track_data)
    XCorReuslt(i,j)=max([CrossCorrMatrix{i,j} CrossCorrMatrix{j,i}]);
end
end


XCorReuslt = zeros(numel(track_data),numel(track_data));
for i=1:numel(track_data)
for j=1:numel(track_data)
    
    if i~=j
    XCorReuslt(i,j)=max(CrossCorrMatrix{i,j});
    end
    
    end
end
