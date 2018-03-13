% ForceMatch_Weekly , would use LookUpTable column 4 week label. 
%MatrixUsed=MatrixUsed+eye(size(MatrixUsed,1))*max(MatrixUsed(:));
%%
weeklabel     =    lookUpTable(:,4);
weeklist      = unique(weeklabel);
NumOfWeek     =   numel(weeklist);
clear ConversionMatrix
ConversionMatrix = zeros(sum(weeklabel==min(weeklabel)), NumOfWeek);
ConversionMatrix(1:sum(weeklabel==min(weeklabel)),1)= 1:sum(weeklabel==min(weeklabel));
    
for wk =1:NumOfWeek-1
for unitIndex = 1:size(ConversionMatrix,1)
    
    if ConversionMatrix(unitIndex,wk)~=0

TargetList = (MatrixUsed(ConversionMatrix(unitIndex,wk),:).*(weeklabel==weeklist(wk+1))');
TargetList(TargetList==0)=max(TargetList);
ConversionMatrix(unitIndex,wk+1) = find(TargetList==min(TargetList));   

    end

end 
leftover = setdiff(find(weeklabel==weeklist(wk+1)),ConversionMatrix(:,wk+1));

if ~isempty(leftover)
ConversionMatrix(end+1:end+numel(leftover),wk+1)=leftover;
end

end 

%% Reverse the label

cTableWk(:,1) = weeklist;
cTableWk(:,2) = flip(weeklist);
for label=1:numel(weeklabel)
weeklabel(label)=cTableWk(find(cTableWk(:,1)==weeklabel(label)),2);
end

%%
weeklist      = unique(weeklabel);
NumOfWeek     =   numel(weeklist);
clear ConversionMatrix
ConversionMatrix = zeros(sum(weeklabel==min(weeklabel)), NumOfWeek);
ConversionMatrix(1:sum(weeklabel==min(weeklabel)),1)= numel(weeklabel):-1:numel(weeklabel)-sum(weeklabel==min(weeklabel))+1;
    
for wk =1:NumOfWeek-1
for unitIndex = 1:size(ConversionMatrix,1)
    
    if ConversionMatrix(unitIndex,wk)~=0

TargetList = (MatrixUsed(ConversionMatrix(unitIndex,wk),:).*(weeklabel==weeklist(wk+1))');
TargetList(TargetList==0)=max(TargetList);
ConversionMatrix(unitIndex,wk+1) = find(TargetList==min(TargetList));   

    end

end 
leftover = setdiff(find(weeklabel==weeklist(wk+1)),ConversionMatrix(:,wk+1));

if ~isempty(leftover)
ConversionMatrix(end+1:end+numel(leftover),wk+1)=leftover;
end

end 


%%
close all
% using the unit with highest firing rate in that week, plot weekly
% tracking result 
% cluID=BreakThrough2(:,2);
t=[-10/20e3:1/20e3:30/20e3];
% pdfCellArray=cell(1,length(numel(unique(cluID))));
numOfWks = numel(unique(lookUpTable(:,4)));
rgb = linspecer(max(lookUpTable(:,4)),'sequential');%jet(max(lookUpTable(:,4))); 
TrackedWks = [] 
for clu=1:size(ConversionMatrix,1)
    pdfCellArray{clu}=[num2str(clu) '.png'];
    
    list=ConversionMatrix(clu,ConversionMatrix(clu,:)>0);
    wkNums = unique(lookUpTable(:,4));
    WkUnitList=[];
    for wk = 1:numOfWks
       WkUnitShortList=list(lookUpTable(list,4) == wkNums(wk));
       if ~isempty(WkUnitShortList)
       clear FRcount
       for us=1:numel(WkUnitShortList)
       FRcount(us)= numel(TimeFeatureAll{WkUnitShortList(us)})
       end
           WkUnitList(end+1)=  WkUnitShortList( FRcount==max(FRcount));
       end 
    end
    
    if numel(WkUnitList)>1
        h=figure;
        set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')
        Ox=0.6;
        Oy=0.6;
        gap=0.6;
        Gw=3;
        Gh=4.5;
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
                    if Ch_Map(row,col)>0 && MaskFeatureAll{WkUnitList(unit)}(temp)>0.01
                        
                        plot(t,DataForFeature{WkUnitList(unit)}(temp,:),'Color',rgb(lookUpTable(find(lookUpTable(:,3)==WkUnitList(unit)),4),:),'LineWidth',2)
                        hold on;
                    end
                    axis([t(1) t(end) -125 25])
                    
                    
                end
            end
        end
        clear WkNums
        
        for unit=1:numel(WkUnitList)
        WkNums(unit)=lookUpTable(find(lookUpTable(:,3)==WkUnitList(unit)),4);
        end
%         Days = unique(Days);
        TrackedWks(end+1)=numel(WkUnitList);
       
        
        totalList = 1000*ones(1,numel(ConversionMatrix(clu,:)));
        for total = 1:numel(totalList)-1
           if ConversionMatrix(clu,total)~=0
           totalList(total) = MatrixUsed(ConversionMatrix(clu,total),ConversionMatrix(clu,total+1));
           end
        end
        
        totalList=flip(totalList);
        totalList(1)=[];
        
        
        Ox=0.6;
        Oy=21;
        gap=0.5;
        Gw=2.4;
        Gh=4.5;
        numGx = max(wkNums);
        numGy = 1;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        
       % select a channel first 
       store = zeros(32,41,numel(WkUnitList));
       for i=1:numel(WkUnitList)
           store(:,:,i)=DataForFeature{WkUnitList(i)};
       end
       MeanStore=mean(store,3);
       AvgP2P=max(MeanStore')-min(MeanStore');
       selectedCh = find(AvgP2P==max(AvgP2P));
       range = [min(min(store(selectedCh,:,:))) max(max(store(selectedCh,:,:)))];
       
        for row=1:1
            for col=1:max(wkNums)
        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')
        
        unitsOfThatWk = WkUnitList(find(WkNums==col));
        if ~ isempty(unitsOfThatWk)
        store = zeros(32,41,numel(unitsOfThatWk));
       for i=1:numel(unitsOfThatWk)
           store(:,:,i)=DataForFeature{unitsOfThatWk(i)};
       end
       MeanStore=mean(store,3);
        plot(MeanStore(selectedCh,:),'Color',rgb(col,:),'LineWidth',3)
        [C,I] = min(MeanStore(selectedCh,:));
       text(I, C, num2str(ConversionMatrix(clu,col)))
        title([num2str(totalList(col))])
        axis([1 41 -180 50])
        end
        end
        end
        
        
%         Ox=21;
%         Oy=21;
%         gap=0.1;
%         Gw=10;
%         Gh=4.5;
%         numGx = 1;
%         numGy = 1;
%         Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
%         Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
%         [X,Y] = meshgrid(Gx,Gy);
%         ax=axes('Units','centimeters','position', [X Y Gw Gh]);
%         set(ax,'FontSize',8)
%         set(ax,'FontWeight','bold')

        
        
        
        
        
        
        set(gcf, 'Position', get(0,'Screensize'));
                set(h,'PaperUnits','centimeters','PaperPositionMode','Manual')
set(h,'PaperPosition',[0 0 40 27])

print(h,'-dpng',pdfCellArray{clu},'-loose')
    end
end
% DIR=dir('*.pdf');
% append_pdfs('result.pdf',DIR.name)
% delete(DIR.name)
close all
% save('pdfNames
clear track
for i=2:9
track(i-1)=histc(TrackedWks,i);
end
bar(2:9,track)

