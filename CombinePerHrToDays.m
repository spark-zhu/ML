clear all
prompt='wait: ';
wait=input(prompt);
pause(wait)
% dataStore = zeros(28,82800*20E3,'int16');
for folder =1:17
    try
    folder 
    clear data

data=readmda(['ds' num2str(folder) '\12242017.mda'],0);
data=int16(data);
writemda(data,['ds' num2str(folder) '\12242017.mda'],'int16');

%     dataStore(:,(folder-1)*72E6+1:(folder)*72E6)=data;

    catch
        ['error at' num2str(folder)]
    end
end
%%
clear all
% dataStore = zeros(28,82800*20E3,'int16');
for folder =1:25
    try
    folder 
    clear data

data=readmda([num2str(folder) '\CMR\12142017.mda'],0);
data=int16(data);
writemda(data,[num2str(folder) '\CMR\12142017.mda'],'int16');

%     dataStore(:,(folder-1)*72E6+1:(folder)*72E6)=data;

    catch
        ['error at' num2str(folder)]
    end
end


%%
clear all
prompt='wait: ';
wait=input(prompt);
pause(wait)
dateStr='12202017';
finalFolder = 24;
data=readmda(['ds' num2str(finalFolder) '\' dateStr '.mda'],1);
dataStore = zeros(28,size(data,2)+3600*(finalFolder-1)*20E3,'int16');


for folder =1:finalFolder
    try
    folder 
    clear data

data=readmda(['ds' num2str(folder) '\' dateStr '.mda'],1);
% data=int16(data);
% writemda(data,['ds' num2str(folder) '\12182017.mda'],'int16');

    dataStore(:,(folder-1)*72E6+1:min((folder)*72E6,size(dataStore,2)))=data;

    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        ['error at' num2str(folder)]
    end
end
cd('E:/')
writemda(dataStore,[dateStr '.mda'],'int16');
clear all
%%
clear all
prompt='wait: ';
wait=input(prompt);
pause(wait)
dateStr='11272017';
finalFolder = 18;
data=readmda([num2str(finalFolder+19) '\CMR\' dateStr '.mda'],1);
dataStore = zeros(29,size(data,2)+3600*(finalFolder-1)*20E3,'int16');


for folder =1:finalFolder
    try
    folder 
    clear data

data=readmda([num2str(folder+19) '\CMR\' dateStr '.mda'],1);
% data=int16(data);
% writemda(data,['ds' num2str(folder) '\12182017.mda'],'int16');

    dataStore(:,(folder-1)*72E6+1:min((folder)*72E6,size(dataStore,2)))=data;

    catch e
        fprintf(1,'There was an error! The message was:\n%s',e.message);
        ['error at' num2str(folder)]
    end
end
cd('F:/')
 dateStr='11282017';
writemda(dataStore,[dateStr '.mda'],'int16');
clear all