clc
currentpath=pwd;
s1=strfind(currentpath,'recording');
if isempty(s1)
    folder_name = uigetdir;
    cd(folder_name)
    currentpath=pwd;
    s1=strfind(currentpath,'recording');
end
s2=strfind(currentpath,'\');
if numel(s2)>find(s2==(s1-1))
mainFolderPath = currentpath(1:s2(find(s2==(s1-1))+1)-1);
else
mainFolderPath = currentpath;
end


   if exist([mainFolderPath '\Dat_V_Map.mat'], 'file')
   load([mainFolderPath '\Dat_V_Map.mat'])
   elseif exist('Dat_V_Map.mat', 'file')
   load('Dat_V_Map.mat')
    else
      read_Intan_RHD2000_file
      y=[amplifier_channels.native_order];imp=[amplifier_channels.electrode_impedance_magnitude];

selection = (y<48 & y >15)&(imp<3E6);
% selection = 1:31;
amplifier_channels=amplifier_channels(selection);
      x(:,2)=[amplifier_channels.native_order]';
      x(:,1)=1:size(x,1);
      Dat_V_Map=x;
      save([mainFolderPath '\Dat_V_Map.mat'],'Dat_V_Map')
   end
   
   
switch ChMapNum
   case 1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 0
   load('D:\0517-hippo\ch_map_pink.mat');
   case 2
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_large.mat');
   case 3
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_Small.mat');
    case 4
    load('D:\Box Sync\0919 stroke\Finger_map.mat');
        
end

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
% clearvars -except Ch_Map Ch_Map_2 Dat_V_Map str folderList syncFolder folder localFolder storageFolder folderName distributionOfLabor
graph=ConnectedGraph(Ch_Map_2);
graph=graph(8:end);
numOfChs=num2str(size(Dat_V_Map,1));

textProbe = fileread('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\Klusta files\probeFile.prb');
newStr = [textProbe(1:321) graph textProbe(322:end)];
newStr = [newStr(1:155) numOfChs newStr(156:end)];
fileID = fopen('probe.prb','w');
fprintf(fileID, newStr);
fclose(fileID);


numOfChs=num2str(size(Dat_V_Map,1));
textPrm = fileread('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\Klusta files\combinedfiles.prm');
newStr = [textPrm(1:174) numOfChs textPrm(175:end)];

% [file, path, ~] = ...
%     uigetfile('*.dat', 'Select an DAT', 'MultiSelect', 'off');

DatName=dir('*.dat');
DatFile =DatName(1).name; %[path,file];%'wrapping_12_28.dat';
file=DatFile(1:end-4);
newStr = [newStr(1:19) file newStr(20:end)];
fileID = fopen('combine.prm','w');
fprintf(fileID, newStr);
fclose(fileID);

% 
% 
parentPath = pwd;
try
cd('positive')
textPrm = fileread('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\Klusta files\combinedfiles_Positive.prm');
newStr = [textPrm(1:174) numOfChs textPrm(175:end)];
DatName=dir('*.dat');
DatFile =DatName(1).name; %[path,file];%'wrapping_12_28.dat';
file=DatFile(1:end-4);
newStr = [newStr(1:19) file newStr(20:end)];
fileID = fopen('combine.prm','w');
fprintf(fileID, newStr);
fclose(fileID);


textProbe = fileread('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\Klusta files\probeFile.prb');
newStr = [textProbe(1:321) graph textProbe(322:end)];
newStr = [newStr(1:155) numOfChs newStr(156:end)];
fileID = fopen('probe.prb','w');
fprintf(fileID, newStr);
fclose(fileID);
catch
end


cd(parentPath)
% 
