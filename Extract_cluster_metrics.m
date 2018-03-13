clear all
mouseStr='mouse0620';
files = dir('*.mat');
filesName = {files.name};
desiredFolder = cellfun(@(y) y(end)=='t', filesName);

files=files(desiredFolder);
%%
for file=1:numel(files)
file
     

Cluster_result = loadjson(['cluster_metrics' files(file).name(1:8) '.json']);
Sorting_result = load(files(file).name);
Cluster_result = Cluster_result.clusters;
clu_label = cellfun(@(x) x.label,Cluster_result)
Cluster_result1=Cluster_result;
Cluster_result=cell(1,max(clu_label));
for clu=1:numel(clu_label)
Cluster_result{clu_label(clu)}=Cluster_result1{clu};
end
Sorting_result =Sorting_result.result;
Cluster_result = Cluster_result(Sorting_result.list);
Sorting_result.isolation=cellfun(@(x) x.metrics.isolation, Cluster_result)';
Sorting_result.noise_overlap=cellfun(@(x) x.metrics.noise_overlap, Cluster_result)';
Sorting_result.SNR=cellfun(@(x) x.metrics.peak_snr, Cluster_result,'UniformOutput',0 )';
Sorting_result.bursting_parent=cellfun(@(x) x.metrics.bursting_parent, Cluster_result)';
Sorting_result.overlap_clu=cellfun(@(x) x.metrics.overlap_cluster, Cluster_result)';

for unit=1:numel(Sorting_result.list)
wav = Sorting_result.waveform{unit};
wav_std = Sorting_result.waveformStd{unit};

P2P = min(wav');
ch=find(P2P==min(P2P));
ch=ch(1);
tp=find(wav(ch,:)==min(P2P));
tp=tp(1);
Sorting_result.peakCV(unit)=abs(wav_std(ch,tp)/min(P2P));

FiringTimeForThisUnit = Sorting_result.time{unit};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
maxTime = max(cell2mat(Sorting_result.time'))/20E3;
try
Sorting_result.multi1(unit)= counts(1);
catch
    try
Sorting_result.multi2(unit)= counts(2);
    catch
    end
end
% we allow one spike per min.
end
result = Sorting_result;
save(files(file).name,'result')
end
%%
for file=1:numel(files)
unix(['cp ' files(file).name '/Refine_mountain/' files(file).name  '_tracking.mat' ' /media/D1/0109tracking'])
end
%%
Datestr='12182017'
% subFolderList = dir(['datasets/']); 
% subFolderName = {subFolderList.name};
% desiredFolder = cellfun(@(y) y(1)=='d' , subFolderName);
% 
% subFolderList=subFolderList(desiredFolder);
     
% for subFolder =1:numel(subFolderList)
%           dsName = ['ds' num2str(subFolder)];

Cluster_result = loadjson('cluster_metrics.json');
Sorting_result = load([Datestr '_tracking.mat']);
Cluster_result = Cluster_result.clusters;
Sorting_result =Sorting_result.result;
Cluster_result = Cluster_result(Sorting_result.list);
Sorting_result.isolation=cellfun(@(x) x.metrics.isolation, Cluster_result)';
Sorting_result.noise_overlap=cellfun(@(x) x.metrics.noise_overlap, Cluster_result)';
Sorting_result.SNR=cellfun(@(x) x.metrics.peak_snr, Cluster_result)';
Sorting_result.bursting_parent=cellfun(@(x) x.metrics.bursting_parent, Cluster_result)';

for unit=1:numel(Sorting_result.list)
wav = Sorting_result.waveform{unit};
wav_std = Sorting_result.waveformStd{unit};

P2P = min(wav');
ch=find(P2P==min(P2P));
tp=find(wav(ch,:)==min(P2P));
Sorting_result.peakCV(unit)=abs(wav_std(ch,tp)/min(P2P));

FiringTimeForThisUnit = Sorting_result.time{unit};
ISI_in_MS_bins = diff(double(FiringTimeForThisUnit))/20; % Intervals in miliseconds. 
[counts] = histc(ISI_in_MS_bins,0:1:50);
maxTime = max(cell2mat(Sorting_result.time'))/20E3;
Sorting_result.multi(unit)= counts(1)+counts(2); % we allow one spike per min.
end

result = Sorting_result;
% Sorting_result = load(['datasets/' dsName '/Refine_mountain/' Datestr '_tracking.mat']);
save([Datestr '_tracking.mat'],'result')
% end
%%
fileList = dir('*.mat');

FR_summary = [];
noise_over_summary = [];
isolation_summary = [];
multi_summary = [];
Labels = [];
CV = [];
AMP = [];

for file=1:numel(fileList)
load(fileList(file).name);
FR_summary=[FR_summary;result.Avg_FR];
noise_over_summary = [noise_over_summary;result.noise_overlap];
isolation_summary=[isolation_summary;result.isolation];
multi_summary=[multi_summary result.multi];
CV=[CV result.peakCV];
Labels =[Labels; result.selectedClusters];
AMP = [AMP; cellfun(@max,result.P2P)];



end
CV=CV';
multi_summary=multi_summary';