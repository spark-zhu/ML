clear all
close all
folderList=dir;
ImpedanceSummary=zeros(32,1);
Legend=cell(1,1);
for folder=1:numel(folderList)
  Name=folderList(folder).name;
   if ~isempty(strfind(Name,'record'))
       loc=strfind(Name,'g');
       Legend{end+1}=folderList(folder).name(loc+1:end);
     cd(folderList(folder).name)
       ElecList=dir('*.csv'); % get all electrical recording spike timing data file spec.
ImpedanceRec=zeros(32,size(ElecList,1));



for i=1:size(ElecList,1) 
book=importdata(ElecList(i).name); 

impedance=book.data(:,2);
Enable=book.data(:,1);
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
ImpedanceSummary(:,end+1)=AvgImpedance;

     cd ..
   end
   
end
ImpedanceSummary=ImpedanceSummary(:,2:end);
figure
Legend=Legend(2:end);
plot(ImpedanceSummary)
legend(Legend)
path=pwd;
a=strfind(path,'\');
title(path(a(end)+1:end));
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
Ch_Map=Ch_Map_new-15;
Transformed_impedance_summary=cell(1,numel(Legend));
%
for ses=1:numel(Legend)
    clear newImpMatrix
for row=1:4
    for col=1:8
    newImpMatrix(row, col) =  ImpedanceSummary(Ch_Map(row,col),ses);
    if newImpMatrix(row, col)>2E6
      newImpMatrix(row, col)= 2E6;  
    end
    end
end
Transformed_impedance_summary{ses} = newImpMatrix;
end


figure 
for i=1:numel(Legend)
subplot(8,8,i)
imagesc(Transformed_impedance_summary{i});
caxis([250E3 2E6])
end
colormap(jet)

%% exclude all possible bad channels 
% ImpedanceSummary([1 2 8 28],:)=[];
% ImpedanceSummary([11 17],:)=[];

for ch= 1:numel(Legend)
  impCh =   ImpedanceSummary(:,ch);
 
meanImp(ch) = mean(impCh(impCh<2E6));
stdImp(ch) = std(impCh(impCh<2E6));
Nelectrode(ch) = sum(impCh<2E6);
if Nelectrode(ch)==0
   Nelectrode(ch)=NaN; 
end
end

prompt='surgery Date in mmddyyyy: ';
surgeryDate = input(prompt,'s');
startNum = datenum(surgeryDate,'mmddyyyy');
DPS = zeros(1,numel(Legend));
for i=1:numel(Legend)
DPS(i) = datenum(Legend{i},'yyyy-mm-dd')-startNum;
WPS(i) = floor(DPS(i)/7);
end


figure
%yyaxis left
errorbar(DPS(12:end),meanImp(12:end),stdImp(12:end),'Linewidth',2)

set(gca,'FontSize',12,'fontWeight','bold')
xlabel('Days')
ylabel('Impedance (ohms)')
% title('Single-Units')
set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')
 ylim([0 2E6])
%  hold on;
% yyaxis right
% scatter(DPS,Nelectrode,'filled')
% hold on
% plot(DPS,Nelectrode)
% % plot(DPS,median_P2P,'r')
% hold on;
% plot(DPS,max_P2P,'k')


%% Correlate Noise Level and Impedance. 
 NoiseFilSummary=cell(1,numel(Legend));
 NoiseRawSummary=cell(1,numel(Legend));
for folder=1:numel(folderList)
    folder
  Name=folderList(folder).name;
   if ~isempty(strfind(Name,'record'))
       loc=strfind(Name,'g');
       cd(folderList(folder).name)
       
   
DIR=dir('*.rhd');
 if ~isempty(DIR)
i=floor(numel(DIR)/2);
read_Intan_RHD2000_file_special(DIR(i).name);

data=amplifier_data;

% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters folderList syncFolder folder localFolder storageFolder folderName distributionOfLabor

% Do Channel rejection after ch 47, before ch16
y=[amplifier_channels.native_order];
selection = (y<48 & y >15);
amplifier_channels=amplifier_channels(selection);
data=data(selection,:);

% Commmon Reference Montage.

amplifier_data=data; clear data;
window_length = 10000;
dend=floor(size(amplifier_data,2)/window_length)*window_length;
amplifier_data=amplifier_data(:,1:dend);
% balcklist = [45 46 47];
% amplifier_data(balcklist-15,:)=[];
% amplifier_channels(balcklist-15)=[];
RawstdNoise  = std(amplifier_data');
data=amplifier_data(1,:);
[b,a]=butter(4,300/10000,'high');
[b,a]=butter(2,300/10000,'high');
fprintf(1,'Filtering progress: ')
for i=1:size(amplifier_data,1)
    fprintf(1,'%4.2f',i/size(amplifier_data,1));   
amplifier_data(i,:) = filtfilt(b,a,amplifier_data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')


stdNoise  = std(amplifier_data');
clear x
x(:,2)=[amplifier_channels.custom_order]';
      x(:,1)=1:size(x,1);
      Dat_V_Map=x;
      load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
Ch_Map=Ch_Map_new;
for i=1:size(Ch_Map,1)
    for j=1:size(Ch_Map,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
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

NoiseMapFil = zeros(4,8);
for i=1:size(NoiseMapFil,1)
    for j=1:size(NoiseMapFil,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
    NoiseMapFil(i,j)=0;
    else 
    NoiseMapFil(i,j)=stdNoise(Ch_Map_2(i,j));
    end
    end
end

NoiseMap = zeros(4,8);
for i=1:size(NoiseMap,1)
    for j=1:size(NoiseMap,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
    NoiseMap(i,j)=0;
    else 
    NoiseMap(i,j)=RawstdNoise(Ch_Map_2(i,j));
    end
    end
end



% save('stdV','stdNoise')
% for i=1:size(amplifier_data,1)
%  
% amplifier_data(i,:) = amplifier_data(i,:)/stdNoise(i)*mean(stdNoise);
% 
% end
% amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);


%  fprintf(1,'feature calculating progress: ')
% 
%  for i=1:dend/window_length
%  fprintf(1,'%4.2f',i/dend*window_length);   
%  temp=amplifier_data(:,(i-1)*window_length+1:i*window_length);
%  feat1(i)=median(mean((temp').^2));
%  feat2(i)=max(max(abs(temp)));
%  fprintf(1,'\b\b\b\b')
%  end
% fprintf('\n')
%  
% 
%  
%  mod_feat1=repmat(feat1,window_length,1);
%  align_feat1=reshape(mod_feat1,dend,1);
%  
%  mod_feat2=repmat(feat2,window_length,1);
%  align_feat2=reshape(mod_feat2,dend,1);
%  titlePath=pwd;
% figure('Name',titlePath) 
% hist(floor(feat1),max(floor(feat1))-min(floor(feat1))+1);
% title('MS_CRM')
% %selected 2 threshold. 87  and 124 
% 
% % figure
% % [x1,y1]=ecdf(feat1);
% % plot(y1,x1)
% % title('MS')
% 
% figure('Name',titlePath) 
% hist(floor(feat2),max(floor(feat2))-min(floor(feat2))+1);
% title('Max_CRM')
% selected 2 threshold. 87  and 124 

% figure
% [x2,y2]=ecdf(feat2);
% plot(y2,x2)
% title('Max')
 NoiseFilSummary{folder} = NoiseMapFil;
 NoiseRawSummary{folder} = NoiseMap;
 else
        NoiseFilSummary{folder} = [];
 NoiseRawSummary{folder} = [];
 end

   cd ..
   end
end
%%
NoiseFilSummaryReShape{folder} = cell(size(Legend));
 NoiseRawSummaryReShape{folder} = cell(size(Legend));
 Transformed_impedance_summaryReshape{folder}     = cell(size(Legend));
 
 
for folder=1:numel(Legend)
if isempty(NoiseFilSummary{folder})||isempty(NoiseRawSummary{folder})
NoiseFilSummary{folder} = zeros(4,8);
 NoiseRawSummary{folder} = zeros(4,8);
end
NoiseFilSummaryReShape{folder} = reshape(NoiseFilSummary{folder},1,32);
 NoiseRawSummaryReShape{folder} = reshape( NoiseRawSummary{folder} ,1,32);
 Transformed_impedance_summaryReshape{folder}     = reshape(Transformed_impedance_summary{folder},1,32);
 
end
NoiseFilSummaryReShape=NoiseFilSummaryReShape(1:35);
 NoiseRawSummaryReShape=NoiseRawSummaryReShape(1:35);
 Transformed_impedance_summaryReshape=Transformed_impedance_summaryReshapee(1:35);

figure
scatter(cell2mat(Transformed_impedance_summaryReshape),cell2mat(NoiseRawSummaryReShape))

Imp = cell2mat(Transformed_impedance_summaryReshape);
Raw = cell2mat(NoiseRawSummaryReShape);
Fil = cell2mat(NoiseFilSummaryReShape);


selected1 = Fil>0 & Fil < 40;
selected2 = Raw>0 & Raw < 200;
selected3 = Imp>6E5 & Imp<3E6;
FinalSel = selected1 & selected2 & selected3;

Imp = Imp(FinalSel);
Raw = Raw(FinalSel);
Fil = Fil(FinalSel);