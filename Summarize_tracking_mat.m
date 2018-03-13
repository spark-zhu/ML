%% pullout tracking.mat
clear all
mouseStr='mouse05242017';
files = dir;
filesName = {files.name};
desiredFolder = cellfun(@(y) y(end)=='7', filesName);
files=files(desiredFolder);
 files=files(8:36)
%%
for file=1:numel(files)
system(['copy E:\' mouseStr '\' files(file).name '\Refine_mountain\' files(file).name '_tracking.mat'  ' E:\'])
end