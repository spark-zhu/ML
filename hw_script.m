t=(0:numel(hF)-1)/20E3;
plot(t,hF)
xlabel('Time(s)')
ylabel('Voltage (uV)')
title('Raw Data')

[b,a]=butter(4,100/10000);
lP = filtfilt(b,a,hF);
plot(t,lP)
xlabel('Time(s)')
ylabel('Voltage (uV)')
title('Low Passed')

[b,a]=butter(4,300/10000,'high');
lP = filtfilt(b,a,hF);
plot(t,lP)
xlabel('Time(s)')
ylabel('Voltage (uV)')
title('High Passed')

[pks1,locs1] = findpeaks(-double(lP),'MINPEAKHEIGHT',60);
unitTime1=locs1;
store1 = zeros(numel(unitTime1),40);
figure
for i =1 : numel(unitTime1)
store1(i,:)=lP(unitTime1(i)-10:unitTime1(i)+29);
plot(0:1/20:1.95,store1(i,:))

hold on
end
xlabel('Time(ms)')
ylabel('Voltage (uV)')
title('Threshold at -60uV')


[pks1,locs2] = findpeaks(-double(lP),'MINPEAKHEIGHT',130);
unitTime2 = locs2; 
unitTime1 = setdiff(locs1,locs2);

store1 = zeros(numel(unitTime1),40);
figure
for i =1 : numel(unitTime1)
store1(i,:)=lP(unitTime1(i)-10:unitTime1(i)+29);
plot(0:1/20:1.95,store1(i,:))

hold on
end

store2 = zeros(numel(unitTime2),40);

for i =1 : numel(unitTime2)
store2(i,:)=lP(unitTime2(i)-10:unitTime2(i)+29);
plot(0:1/20:1.95,store2(i,:),'r')
hold on
end

xlabel('Time(ms)')
ylabel('Voltage (uV)')
title('Threshold at -60uV and -130uV')

scatter(unitTime1/20E3,ones(numel(unitTime1),1),'r.')
hold on
scatter(unitTime2/20E3,2*ones(numel(unitTime2),1),'b.')
axis([0 100 0 3])
xlabel('Time(s)')
title('Firing raster plot')

    hist(unitTime1/20E3,0.5:1:99.5)
    xlabel('Time(s)')
    ylabel('Firing Rate (Hz)')
    title('Unit 1 Firing Rate')
[timeLag_In_Samples,CrossCorrMatrix,binCountSummary]=crossCorr(unitTime1',unitTime2');
plot(timeLag_In_Samples,CrossCorrMatrix{1,2})