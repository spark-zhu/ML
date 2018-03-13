
load('savePath.mat')

% 


FolderContent = dir;
candidateFolder = zeros(size(FolderContent));
for folder=1:numel(FolderContent)
candidateFolder (folder ) = ~isempty(strfind(FolderContent(folder).name,'record')) && FolderContent(folder).isdir;
end
folderList=FolderContent(logical(candidateFolder));


for folder = 46:50
folderName = folderList(folder).name;
cd(folderName)
str = [folderName(end-4:end-3) folderName(end-1:end) folderName(end-9:end-6)];
try
[status,cmdout] = dos(['copy ' [str '.pdf '] savePath],'-echo')
catch
end
cd .. 
end