clear all
close all
DIR=dir('*.rhs');
for i=1:numel(DIR)
read_Intan_RHS2000_file(DIR(i).name);

if i==1
data=amplifier_data;
else
    data=[data amplifier_data];
end

clearvars -except DIR data amplifier_data amplifier_channels frequency_parameters
end
window_length = 20000;
amplifier_data=data;
dend=floor(size(amplifier_data,2)/window_length)*window_length;
amplifier_data=amplifier_data(:,1:dend);

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
 fprintf(1,'feature calculating progress: ')

 for i=1:dend/window_length
 fprintf(1,'%4.2f',i/dend*window_length);   
 temp=amplifier_data(:,(i-1)*window_length+1:i*window_length);
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
%% what percentage of data does a each threhold includes? 
param.feat1=feat1;
param.feat2=feat2;

param.thre_Median_MS = 100;
param.Median_MS_Include = sum(feat1<param.thre_Median_MS)/numel(feat1);

param.thre_Max_abs = 260;
param.Max_abs_Include =  sum(feat2<param.thre_Max_abs)/numel(feat2);

param.total_Include = sum(feat1<param.thre_Median_MS&feat2<param.thre_Max_abs)/numel(feat1);

[param.Median_MS_Include param.Max_abs_Include param.total_Include]



data=amplifier_data(:,logical(align_feat1<param.thre_Median_MS)&logical(align_feat2<param.thre_Max_abs));
size(data,2)/20e3/60

%% find 12 mins around the mid point 
prompt='File name: ';
str = input(prompt,'s');
x(:,2)=[amplifier_channels.custom_order]';
x(:,2)=x(:,2)+1;
load('C:\Users\xie.lab.ws1\Desktop\RHD_TO_RHS_MAP.mat')
for i=1:numel(x(:,2))
x(i,2)=RHD_TO_RHS_MAP(find(RHD_TO_RHS_MAP(:,2)==x(i,2)),1);
end
x(:,1)=1:size(x,1);
Dat_V_Map=x;
save('Dat_V_Map','Dat_V_Map')
dur=6;
data=data(:,floor(size(data,2)/2)-dur*60*20e3+1:size(data,2)/2+dur*60*20e3);
data=floor(data'*100); % 100 times bigger to avoid huge round off error;
% fid=fopen('test.dat','w');
% fwrite(fid,data,'int16');
% fclose(fid)
data=int16(data);
save(str,'data')
save('param','param')
% copied the python file. 
[status,cmdout] = dos('copy "D:\Xue_device_ZhengtuoSurgery_Feb6_tracking\recording 2017-02-21\NoRestriction\mat2dat.py" mat2dat.py','-echo')
% swtich to dat.
[status,cmdout] = dos(['python mat2dat.py ' str '.mat'],'-echo')
% copy dat file for postive thresholding
parentPath = pwd;
mkdir('positive')
cd('positive')

[status,cmdout] = dos(['copy ' '"' parentPath '\' str '.mat.dat' '" ' str '.mat.dat'],'-echo')
cd(parentPath)
% create combine.prm, prepare for klusta
RHS_Get_ChMap_ChMap2_DAT2CH
% start Klusta and perform clustering
[status,cmdout] = dos('activate klusta & klusta combine.prm','-echo')
clearvars -except str
dataFile=str; % inheritated from the GiveMat_Acro_Noise.m 
load([dataFile '.mat'])
RHS_weighted_center
curDir = pwd;
cd('..')
load('savePath.mat')
cd(curDir);
[status,cmdout] = dos(['copy ' [newStr '.pdf '] savePath],'-echo')
close all 


clearvars -except str
dataFile=str; % inheritated from the GiveMat_Acro_Noise.m 
load([dataFile '.mat'])
cd('positive')
[status,cmdout] = dos('activate klusta & klusta combine.prm','-echo')
RHS_weighted_center
% plot everything
close all 






