clear all
close all
clear all
close all

% prompt='Folder name: ';
% Folderstr=input(prompt,'s');
% 
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



prompt='Channel_Map=: 1 for 0109, 0 for 0524/0620, 2 for 0807, 3 for 0809, 4 for 0919stroke';
ChMapNum=input(prompt);



DIR=dir('*.rhd');
for i=1:numel(DIR)
read_Intan_RHD2000_file_special(DIR(i).name);

if i==1
data=amplifier_data;
else
    data=[data amplifier_data];
end

% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
end
% Do Channel rejection after ch 47, before ch16
y=[amplifier_channels.native_order];imp=[amplifier_channels.electrode_impedance_magnitude];
selection = (y<48 & y >15)&(imp<2E6); 
if strcmp(str,'10252017')||strcmp(str,'10272017')||strcmp(str,'10292017')
selection(32:end)=0;
end
amplifier_channels=amplifier_channels(selection);
data=data(selection,:);
%


amplifier_data=data; 
window_length = 10000;
dend=floor(size(amplifier_data,2)/window_length)*window_length;
amplifier_data=amplifier_data(:,1:dend);
% data_copy=amplifier_data;
% balcklist = [45 46 47];
% amplifier_data(balcklist-15,:)=[];
% amplifier_channels(balcklist-15)=[];


% data=amplifier_data(1,:);
% [b,a]=butter(4,300/10000,'high');
[b,a]=butter(4,300/10000,'high');
fprintf(1,'Filtering progress: ')
for i=1:size(amplifier_data,1)
    fprintf(1,'%4.2f',i/size(amplifier_data,1));   
amplifier_data(i,:) = filtfilt(b,a,amplifier_data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')
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


%  stdNoise = std(amplifier_data');
%  for ch=1:size(amplifier_data,1)
%  sig=amplifier_data(ch,:);
%  sig(sig>6*stdNoise(ch))=[];
%  stdNoise(ch)=std(sig);
%  end
%  save('stdV','stdNoise')


 fprintf(1,'feature calculating progress: ')

 for i=1:dend/window_length
 fprintf(1,'%4.2f',i/dend*window_length);   
 temp=amplifier_data(:,(i-1)*window_length+1:i*window_length);
 feat1(i)=median(mean((temp').^2));
 feat2(i)=max(max(abs(temp)));
 fprintf(1,'\b\b\b\b')
 end
 clear temp
fprintf('\n')
 
%  feat1=mean(newfil.^2);
%  feat1(feat1>400)=400;
%  
%  newfil=reshape(data,20000,numel(data)/20000);
%  feat2=max(abs(newfil));
%  feat2(feat2>600)=600;
 
 mod_feat1=repmat(feat1,window_length,1);
 align_feat1=reshape(mod_feat1,dend,1);
 
 mod_feat2=repmat(feat2,window_length,1);
 align_feat2=reshape(mod_feat2,dend,1);
 titlePath=pwd;
figure('Name',titlePath) 
hist(floor(feat1),max(floor(feat1))-min(floor(feat1))+1);
title('MS')
%selected 2 threshold. 87  and 124 

% figure
% [x1,y1]=ecdf(feat1);
% plot(y1,x1)
% title('MS')

figure('Name',titlePath) 
hist(floor(feat2),max(floor(feat2))-min(floor(feat2))+1);
title('Max')


param.feat1=feat1;
param.feat2=feat2;

param.thre_Median_MS = 400;
param.Median_MS_Include = sum(feat1<param.thre_Median_MS)/numel(feat1);

param.thre_Max_abs = 400;
param.Max_abs_Include =  sum(feat2<param.thre_Max_abs)/numel(feat2);

param.total_Include = sum(feat1<param.thre_Median_MS&feat2<param.thre_Max_abs)/numel(feat1);

[param.Median_MS_Include param.Max_abs_Include param.total_Include]



data=data(:,logical(align_feat1<param.thre_Median_MS)&logical(align_feat2<param.thre_Max_abs));
clear mod_feat2 mod_feat1 align_feat1 align_feat2 

amplifier_data=data; 
[b,a]=butter(4,3/10000,'high');
fprintf(1,'Filtering progress: ')
for i=1:size(amplifier_data,1)
    fprintf(1,'%4.2f',i/size(amplifier_data,1));   
amplifier_data(i,:) = filtfilt(b,a,amplifier_data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')


if CMR ==1
amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);
elseif CMR==2
    list1=1:15;
    list2=16:31;
amplifier_data(list1,:) = amplifier_data(list1,:) - repmat(median(amplifier_data(list1,:)),size(amplifier_data(list1,:),1),1);
amplifier_data(list2,:) = amplifier_data(list2,:) - repmat(median(amplifier_data(list2,:)),size(amplifier_data(list2,:),1),1);
         
end
close all
plot(amplifier_data(1,:));
title(str)
prompt='thres: ';
thres=input(prompt);
selected = abs(amplifier_data(1,:))<thres;
amplifier_data=amplifier_data(:,selected);
dur=10;
try
amplifier_data=amplifier_data(:,floor(size(amplifier_data,2)/2)-dur*60*20e3+1:size(amplifier_data,2)/2+dur*60*20e3);
catch
    
end


 freqrange=[4 8;8 15;15 30;30 60;60 110;4 30;4 110];
 Fs = 20E3;
 window = 60; %unit in seconds.
 % assume non-overlapping windows. 
 window_length =  floor(window*Fs);
dend=floor(size(amplifier_data,2)/window_length)*window_length;
amplifier_data=amplifier_data(:,1:dend);
result = cell(size(freqrange,1),1);
% loop through per frequency range per channel
for range=1:numel(result)
   
    perChannelEachWindow = zeros(size(amplifier_data,1),size(amplifier_data,2)/window_length);
   for ch=1: size(amplifier_data,1)
   [range ch]
       perChannelEachWindow(ch,:)=bandpower(reshape(amplifier_data(ch,:),window_length,size(amplifier_data,2)/window_length),Fs,freqrange(range,:));
   end
   result{range}= perChannelEachWindow;
end
x(:,2)=[amplifier_channels.native_order]';
x(:,1)=1:size(x,1);
Dat_V_Map=x;
save(['D:\Box Sync\0919 stroke\' str 'LFP'],'result','Dat_V_Map');




%%
% clear all
% clc
% files=dir('*.mat');
% % importData into cells.
% AllData = cell(size(files,1),1);
% load('D:\Box Sync\0919 stroke\Finger_map.mat')
% Pcut=0.05; 
% for i=1:size(AllData,1)
% AllData{i}=importdata(files(i).name);
% Ch_Map=Ch_Map_new;
% 
%     DataSort  = cellfun(@(x) sort(x,2),AllData{i}.result,'UniformOutput' ,0);
%     SampleCutOff = floor((1-Pcut)*size(AllData{i}.result{1},2));
%     DataCut=cellfun(@(x) x(:,1:SampleCutOff),DataSort,'UniformOutput' ,0);
%     Datastore= cellfun(@(x) [mean(x');std(x')], DataCut,'UniformOutput' ,0);
%     DataReStore = cell(4,8);
% 
% Dat_V_Map=AllData{i}.Dat_V_Map;
% m=0;
% for m=1:size(Ch_Map,1)
% for j=1:size(Ch_Map,2)
%     if isempty(find(Dat_V_Map(:,2)==Ch_Map(m,j), 1))
%         Ch_Map(m,j)=0;
%     end
% end
% end
% 
% Ch_Map_2 = Ch_Map; 
% for m=1:size(Ch_Map,1)
% for j=1:size(Ch_Map,2)
% if isempty(find(Dat_V_Map(:,2)==Ch_Map(m,j)))
% Ch_Map_2(m,j)=0;
% else 
% Ch_Map_2(m,j)=find(Dat_V_Map(:,2)==Ch_Map(m,j));
% 
% DataReStore{m,j}= cellfun(@(x) x(:,Ch_Map_2(m,j)),Datastore,'UniformOutput' ,0);
% end
% end
% end
% 
% AllData{i}.Ch_Map_2=Ch_Map_2;
% AllData{i}.DataReStore=DataReStore;
% 
% 
% try
% ReOrder(i) = AllData{i}.session;
% catch
%     AllData{i}.session = i;
%     ReOrder(i) = AllData{i}.session;
% end
%     
% 
% 
%    
% 
% 
% end
% [x,y]=sort(ReOrder);
% AllData=AllData(y);
% 
% 
% Global_result = cell(size(AllData{1}.result));
% for band = 1:numel(Global_result)
% 
%     bandresult=cell(4,8);
%     
% for m=1:size(Ch_Map,1)
% for j=1:size(Ch_Map,2)
% track_result=NaN(2,numel(AllData));
%     for day=1:numel(AllData)
%     
%         if ~isempty(AllData{day}.DataReStore{m,j})
%         track_result(:,day)= AllData{day}.DataReStore{m,j}{band};
%         end
%         
%     end
%     bandresult{m,j}=track_result;
% end
% end
% Global_result{band}=bandresult;
% end
% %%
% freqrange=[4 8;8 15;15 30;30 60;60 110;4 30;4 110];
% dayLabel ={'-4','-1','0','1','2','3','5','7','9','11','13','16','22','29','36','43','50'};
% for band = 1:numel(Global_result)
% h=figure('Name',[num2str(freqrange(band,1)) '-' num2str(freqrange(band,2))]);
% 
% for m=1:size(Ch_Map,2)
% for j=1:size(Ch_Map,1)
% 
%     subplot(8,4,(m-1)*4+j)
%     errorbar(1:17,Global_result{band,1}{j,m}(1,:),Global_result{band,1}{j,m}(2,:),'lineWidth',1)
%     set(gca,'FontSize',6)
% xticks(1:16)
% xticklabels(dayLabel)
% xlim([0 17])
% % ylim([0 100])
% end
% end
% set(gcf, 'Position', get(0,'Screensize'));
% saveas(h,h.Name)
% saveas(h,[h.Name '.png'])
% end
% 
% 
