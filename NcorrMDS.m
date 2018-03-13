%%
clear all
try ChMapNum
catch 
  prompt='Channel_Map=0109?: 1 for y, 0 for 0524, 2 for 0807, 3 for 0809';
ChMapNum=input(prompt);
end

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


path1=mainFolderPath;
dataFile=path1([end-4:end-3 end-1:end end-9:end-6]);


try
    load([dataFile '.mat'])
str=dataFile;


catch
    
try 
    cd .. 
        load([dataFile '.mat'])
        str=dataFile;
       

    cd(currentpath)

catch
   prompt=([path1([end-4:end-3 end-1:end end-9:end-6]) '?']);
choice=input(prompt,'s');
if strcmp(choice,'1')
dataFile = path1([end-4:end-3 end-1:end end-9:end-6]);
else
dataFile = choice  ;  
end
load(dataFile)
str=dataFile;

end
    

end

%%
TwoD=0;
data=double(data); 
R = zeros(size(data,2));
for row=1:size(R,1)-1
for col=row+1:size(R,1)
    [row col]
R(row,col)=max(xcorr(data(:,row),data(:,col),20, 'coeff'));
end
end
R=R+R';
for i=1:size(R,1)
    R(i,i)=1;
end
NLnR = -log(R)*100; % major assumption
 
if TwoD==0
 Y = mdscale(NLnR,3);
 x=Y(:,1);
 y=Y(:,2);
 z=Y(:,3);
elseif TwoD==1
     Y = mdscale(NLnR,2);
 x=Y(:,1);
 y=Y(:,2);

end
%% 

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
      x(:,2)=[amplifier_channels.custom_order]';
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
%% M
M=Ch_Map_2;
finalAnswer = cell(1,1);
finalAnswer{1} = [0,0];
for i=1:size(M,1)
    for j=1:size(M,2)
        if M(i,j)~=0
        a=[];
        b=[];
        c=[];
        d=[];
        e=[];
        
        if j+1<=size(M,2)
        a=M(i,j+1);
        end
        
        if i+1<=size(M,1)
        b=M(i+1,j);
        end
        
        if j-1>=1
        e=M(i,j-1);
        end
        
        if ~isempty(a)&& ~isempty(b)
            
        c=M(i+1,j+1);
            
        end
        
        if ~isempty(e)&& ~isempty(b)
          
            d=M(i+1,j-1);
    
        end
        
        if ~isempty(a)&&a~=0
          finalAnswer{end+1}=[M(i,j) ,a];
        end
%         if ~isempty(c)&&c~=0
%           finalAnswer{end+1}=[M(i,j) ,c];
%         end
        if ~isempty(b)&&b~=0
          finalAnswer{end+1}=[M(i,j) ,b];
        end
%         if ~isempty(d)&&d~=0
%          finalAnswer{end+1}=[M(i,j) ,d];
%         end
        end 
    end
end
finalAnswer(1)=[];
%%
if TwoD==0
figure
for i=1:numel(finalAnswer)
    pt1 = finalAnswer{i}(1);
        pt2 = finalAnswer{i}(2);

    plot3(x([pt1 pt2]),y([pt1 pt2]),z([pt1 pt2]),'b')
    hold on
end
txt = arrayfun(@num2cell, Dat_V_Map(:,2));
text(Y(:,1),Y(:,2),Y(:,3)+10,txt')
hold on
scatter3(Y(:,1),Y(:,2),Y(:,3),100,'filled')
elseif TwoD==1
figure
for i=1:numel(finalAnswer)
    pt1 = finalAnswer{i}(1);
        pt2 = finalAnswer{i}(2);

    plot(x([pt1 pt2]),y([pt1 pt2]),'b')
    hold on
end
txt = arrayfun(@num2cell, 1:size(R,1));
text(Y(:,1),Y(:,2)+5,txt')
hold on
scatter(Y(:,1),Y(:,2),100,'filled')
end