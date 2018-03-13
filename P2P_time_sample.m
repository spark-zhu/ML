function output = P2P_time_sample(waveform, negative)
% assume peak happens after valley 
if negative==1

    
    minSample = find(waveform==min(waveform));
    minSample=minSample(1);
    maxSample = find(waveform(minSample:end)==max(waveform(minSample:end)))+minSample-1;
    thres=(max(waveform(minSample:end))-min(waveform))*0.9+min(waveform);
    output = sum(waveform(minSample:maxSample)<=thres);
end
end