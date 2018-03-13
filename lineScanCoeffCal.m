% xmlGetLineScan

 
clear all
close all
folderList=dir;
ImpedanceSummary=zeros(1,3);
Legend=cell(1,1);
for folder=1:numel(folderList)
  Name=folderList(folder).name;
   if ~isempty(strfind(Name,'LineScan'))
       loc=strfind(Name,'44');
       Legend{end+1}=folderList(folder).name(loc+3:end);
     cd(folderList(folder).name)
       ElecList=dir('*.xml'); 
       
 theStruct=parseXML(ElecList(1).name);
 scanTime = str2double(theStruct.Children(4).Children(8).Children(6).Children(8).Attributes(2).Value);
 micronPerPixel = str2double(theStruct.Children(2).Children(36).Children(2).Attributes(2).Value);
 coefficients = micronPerPixel/scanTime;

ImpedanceSummary(end+1,:)=[scanTime micronPerPixel coefficients];

     cd ..
   end
   
end