% create geometry.csv
% load('04232017.mat')
try Dat_V_Map;
catch
load('Dat_V_Map.mat')
end
% load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
switch ChMapNum
   case 1
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat');
   case 0
   load('ch_map_pink.mat');
   case 2
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_large.mat');
   case 3
    load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\NewMap_Hippo_Small.mat');
    case 4
    load('D:\Box Sync\0919 stroke\Finger_map.mat');
    case 5
    load('ch_map_Feb222018.mat');
        
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
% Ox=12.5;    
% Oy=12.5;
%         gap=5;
%         Gw=25;
%         Gh=25;
 if Hip~=1       
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
        
      CSVList = zeros(size(Dat_V_Map));
      for ele = 1:size(Dat_V_Map,1)
      [row,col]=find(Ch_Map_2==ele);
      CSVList(ele,1)=X(row,col);
      CSVList(ele,2)=Y(row,col);
      end
      csvwrite(['geom.csv'],CSVList);
 end
      
      %%
      if Hip==1
      Ox=15;    
Oy=15;
gap=20;
Gw=30;
Gh=30;
        
numGx = 16;
numGy = 2;
Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);
        
      CSVList = zeros(size(Dat_V_Map));
      for ele = 1:size(Dat_V_Map,1)
      [row,col]=find(Ch_Map_2==ele);
      CSVList(ele,1)=X(row,col);
      CSVList(ele,2)=Y(row,col);
      end
      csvwrite('geom.csv',CSVList);
      
      end
      
      