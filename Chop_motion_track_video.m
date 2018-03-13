% Chop motion track video
vidObj=VideoReader('motion-track-1211mouse.avi');
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

k = 1;
while hasFrame(vidObj)
    k
    s(k).cdata = readFrame(vidObj);
    k = k+1;
end

v = VideoWriter('newfile.avi');
open(v)
for file=2:10
    v = VideoWriter([num2str(file) '.avi']);
open(v)

writeVideo(v,s(((file-1)*180+1):min(1797,(file*180))))
close(v)
end