     subFolderList = dir ; 
     subFolderName = {subFolderList.name};
     desiredFolder = cellfun(@(y) ~isempty(str2num(y)) , subFolderName );
     subFolderList=subFolderList(desiredFolder);
     subFolderName=subFolderName(desiredFolder);
     
     subFolderName=subFolderName([subFolderList.isdir]);
     subFolderList=subFolderList([subFolderList.isdir]);
     for folder =1:numel(subFolderList)
         cd(subFolderName{folder})
     delete('*.mat')
     cd ..
     end