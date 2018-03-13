clear all
close all
%
files = dir

     subFolderList = dir ; 
     subFolderName = {subFolderList.name};
     desiredFolder = cellfun(@(y) y(1)=='r' , subFolderName );
     subFolderList=subFolderList(desiredFolder);
%      subFolderList=subFolderList(41:48);
%%
for folder=1:numel(subFolderList)
Folderstr='CMR';
HD_flag=2;
Filter=1;
CMR=1;
Hip=0;
dur=0;
ChMapNum=5;
global_save_path='D:\Box Sync\xie.lab.ws1\0131 mouse\3-4-2017';
cd(subFolderList(folder).name)
path1=pwd;
str = path1([end-4:end-3 end-1:end end-9:end-6]);

DIR=dir('*.rhd');
if isempty(DIR)
      DIR=dir('*.rhs');
end
clear data
for i=1:numel(DIR)
    try
read_Intan_RHD2000_file_special(DIR(i).name);
    catch
  read_Intan_RHS2000_file(DIR(i).name);

    end
if i==1
data=amplifier_data;

try
ADCdata=board_adc_data;
catch
end

try
data2=stim_data;
catch
end


else
    data=[data amplifier_data];
try
ADCdata=[ADCdata board_adc_data];
catch
end

try
data2=[data2 stim_data];
catch
end

end
% size(ADCdata)
% clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
end
% Do Channel rejection after ch 47, before ch16

% amplifier_channels=amplifier_channels(1:58);
% data=data(1:58,:);
y=[amplifier_channels.native_order];imp=[amplifier_channels.electrode_impedance_magnitude];

if HD_flag==1
    selection = (y<129 & y >1)&(imp<3E6); 
% selection = (imp<3E6);
elseif HD_flag==2
  selection = (imp<3E6);
else
   selection = (y<48 & y >15)&(imp<3E6);
end
% selection = 1:31;
amplifier_channels=amplifier_channels(selection);
Fs= frequency_parameters.amplifier_sample_rate;
clear x
if HD_flag~=2
x(:,2)=[amplifier_channels.native_order]';
else
x(:,2)=[amplifier_channels.native_order]'+1;
end

x(:,1)=1:size(x,1);
Dat_V_Map=x;

data=data(selection,:);




window_length = Fs/2;
dend=floor(size(data,2)/window_length)*window_length;
data=data(:,1:dend);
% data_copy=amplifier_data;
% balcklist = [45 46 47];
% amplifier_data(balcklist-15,:)=[];
% amplifier_channels(balcklist-15)=[];


% data=amplifier_data(1,:);
% [b,a]=butter(4,300/10000,'high');
if Filter==1
[b,a]=butter(4,300/(Fs/2),'high');
fprintf(1,'Filtering progress: ')
for i=1:size(data,1)
    fprintf(1,'%4.2f',i/size(data,1));   
data(i,:) = filtfilt(b,a,data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')
end
% med = median(amplifier_data);
% 
% data=med;
% data=amplifier_data(25,:);
%  fil=filtfilt(b,a,data);

if CMR ==1
data = data - repmat(median(data),size(data,1),1);
elseif CMR==2
    list1=1:15;
    list2=16:31;
data(list1,:) = data(list1,:) - repmat(median(data(list1,:)),size(data(list1,:),1),1);
data(list2,:) = data(list2,:) - repmat(median(data(list2,:)),size(data(list2,:),1),1);
         
end

try
 stdNoise = std(data');
catch
    clear stdNoise
for sample=1:size(data,1)
 sample/size(data,1)
    stdNoise(sample)= std(data(sample,:));
end 
end


 for ch=1:size(data,1)
 sig=data(ch,:);
 sig(sig>6*stdNoise(ch))=[];
 stdNoise(ch)=std(sig);
 end
 save('stdV','stdNoise')


 fprintf(1,'feature calculating progress: ')
clear feat1 feat2
 for i=1:dend/window_length
 fprintf(1,'%4.2f',i/dend*window_length);   
 temp=data(:,(i-1)*window_length+1:i*window_length);
 feat1(i)=median(mean((temp').^2));
 feat2(i)=max(max(abs(temp)));
 fprintf(1,'\b\b\b\b')
 end
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
%selected 2 threshold. 87  and 124 

% figure
% [x2,y2]=ecdf(feat2);
% plot(y2,x2)
% title('Max')
%
  ElecList=dir('*.csv'); % get all electrical recording spike timing data file spec.
ImpedanceRec=zeros(32,size(ElecList,1));



for i=1:size(ElecList,1) 
book=importdata(ElecList(i).name); 

impedance=book.data(:,2);

if numel(impedance)>32 && HD_flag~=1
ImpedanceRec(:,i)=impedance(17:48);
else
    try
    ImpedanceRec(:,i)=impedance;
    
    catch
    clear ImpedanceRec 
    ImpedanceRec(:,i)=impedance;
    end
    end

end

%what percentage of data does a each threhold includes? 
param.feat1=feat1;
param.feat2=feat2;

param.thre_Median_MS = 400;
param.Median_MS_Include = sum(feat1<param.thre_Median_MS)/numel(feat1);

param.thre_Max_abs = 400;
param.Max_abs_Include =  sum(feat2<param.thre_Max_abs)/numel(feat2);

param.total_Include = sum(feat1<param.thre_Median_MS&feat2<param.thre_Max_abs)/numel(feat1);

[param.Median_MS_Include param.Max_abs_Include param.total_Include]

% 
% 
%(:,logical(align_feat1<param.thre_Median_MS)&logical(align_feat2<param.thre_Max_abs));
size(data,2)/20e3/60

% find 12 mins around the mid point 

%


if dur>0
dur=6;
try
data=data(:,floor(size(data,2)/2)-dur*60*20e3+1:size(data,2)/2+dur*60*20e3);
catch
end
end

for sample=1:size(data,1)
 sample/size(data,1)
data(sample,:)=floor(data(sample,:)'*100);
end 
 % 100 times bigger to avoid huge round off error;
% fid=fopen('test.dat','w');
% fwrite(fid,data,'int16');
% fclose(fid)
data=int16(data);

amount=size(data,2);
mkdir([global_save_path str])
writemda(data,[global_save_path str '\' str '.mda'],'int16');
savePath = [global_save_path str ];
save([global_save_path str '\Dat_V_Map.mat'],'Dat_V_Map','amount')
try
save([global_save_path str '\STI_DATA'],'data2','-v7.3')
catch
end
% writemda(data(:,4E7:7E7),[str '.mda'],'int16');

% m = matfile([str '.mat'],'Writable',true); %
% m.data=data;
try
save('sync_Signal','ADCdata')
catch
end
try
    savePath=global_save_path;
createGeometryCSV;
catch
'Could Not Generate Ch Map, Be Advised'    
end
% save(str,'data','-v7.3')
%save('param','param','CMR','ChMapNum','ADCdata')
cd ..
close all
end
