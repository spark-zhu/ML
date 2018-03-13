% input para1= time window length after eletrical event to find calcium
% event. e.g. 5 seconds. 
% input parameter 2= a threshold to determine the existence of cal event. 
% e.g.  mean + nth std. 
% loop will be go through all electrical event, outer loop, go through each
% channel, for each event, if immediately followed by a cal event, plot the
% electrical event in a subplot. The whole plot contains all subplots, each
% generated from one electrodes. in each subplot, each sorted unit should have its
% own color.
% input parameter 3= designated ROI cal trace. 
% input paprameter 4= puative corresponding electrical trace.

thre=3;   % unit is std.
win=5;% unit s
cal=ROI(:,3); % using ROI3's calcium trace for this task.
numOfelec=size(ElecRec,1);
figure
colorList=['r','g','k','m','y']; % Assume maximum 5 units can be detected on those electrodes. 

for elec=1:numOfelec
subplot(1,numOfelec,elec)
Rec=ElecRec{elec};
waveform=RecContent{elec};
for event=1:numel(Rec(:,1))
  elecTime=Rec(event,2);
  list= (realTime <= (elecTime+win) & realTime >= (elecTime));
  if max(cal(list))>=thre
  plot(waveform(event,:),colorList(Rec(event,1)))
  hold on;
  end
end
title(ElecList(elec).name)

end

