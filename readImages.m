stack=zeros(77,220,100,'uint16');
for i=1:100
stack(:,:,i)=imread('Substack (861-960).tif',i);
end

Avg=imread('AVG_TSeries-09072016-1402-004.tif');
end
rawList=dir('*.tif');


close all
figure
subplot(2,1,1)
imagesc(img)
colorbar
subplot(2,1,2)
imagesc(stack(:,:,1))
colorbar
    newStack=stack;
for i=1:100
newStack(:,:,i)=stack(:,:,i)/2+Avg;
end
for i=1:100
imwrite(newStack(:,:,i),['newstack' num2str(i) '.tif'],'tif');
end

fileList(1)=rawList(1);
% clear all ch2 data
for i=2:numel(rawList)
if ~isempty(findstr('Ch1',rawList(i).name))
    fileList(end+1)=rawList(i);
end
end