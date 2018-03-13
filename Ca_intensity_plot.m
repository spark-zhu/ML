%% Get Absolute time for a 
% clear all
% fileName='TSeries-08312017-1543-014';% input this manually.
%numOfCol_XML=178858;% this number is for the main folder. % find out this manually.
numOfCol_XML=148568;
TimeString=TSeries090620171348017; %importfileXML_ColC(fileName,1,numOfCol_XML);
TimeString=cellfun(@char,TimeString,'UniformOutput',0);
realTime=0;
for i=1:numOfCol_XML
    k = findstr('absoluteTime',TimeString{i});
    if ~isempty(k)
        temp= findstr('"',TimeString{i});
        numStart=temp(1);numEnd=temp(2);
        realTime(end+1)=str2num(TimeString{i}(numStart+1:numEnd-1));
    end
end
figure
plot(diff(realTime))
realTime=realTime(2:end);
uisave('realTime')
%%
% totalFileNum=3;
% startline=2;
% endline=3251;
% filenameList=cell(totalFileNum,1);
% nameSystem='TM_MC_ROC';
% for i=1:totalFileNum
%     filenameList{i}=[nameSystem num2str(i)];
% end
% 
% for i=1:totalFileNum
% AreaM_Intensity(:,i) = importfile(filenameList{i},startline,endline);
% end
% AreaM_Int_mean=mean(AreaM_Intensity);
% AreaM_Int_std=std(AreaM_Intensity);
% 
% for i=1:totalFileNum
%     figure
%     plot(realTime,AreaM_Intensity(:,i),'b');hold on;
%     plot([realTime(1) realTime(end)],[AreaM_Int_mean(i)+2*AreaM_Int_std(i),AreaM_Int_mean(i)+2*AreaM_Int_std(i)],'g')
%     plot([realTime(1) realTime(end)],[AreaM_Int_mean(i)+3*AreaM_Int_std(i),AreaM_Int_mean(i)+3*AreaM_Int_std(i)],'r')
%     xlabel('time(s)')
%     ylabel('meanAreaIntensity(a.u.)')
%     legend('signal','mean+2std','mean+3std','Location','northwest')
%     title(['TMSA ROI  ' num2str(i)])
% end