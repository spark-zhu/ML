clear all
close all
DIR=dir('*.rhd');
for i=1:numel(DIR)
read_Intan_RHD2000_file_special(DIR(i).name);

if i==1
    amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);

data=amplifier_data;
% data2=stim_data;
else
        amplifier_data = amplifier_data - repmat(median(amplifier_data),size(amplifier_data,1),1);

    data=[data amplifier_data];
%     data2=[data2 stim_data];
end


end

y=[amplifier_channels.native_order];
selection = (y<48 & y >15);
amplifier_channels=amplifier_channels(selection);
data=data(selection,:);

%%

nexFileName=DIR(1).name;
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\Combine_To_Nex\nexFile_temp.mat')
nexFile.freq=frequency_parameters.amplifier_sample_rate;
nexFile.tend=size(data,2)/frequency_parameters.amplifier_sample_rate;
contvarTemp=nexFile.contvars;
contvarTemp.ADFrequency=nexFile.freq;
contvar=cell(size(data,1),1);
for i=1:size(data,1)
    contvarTemp.name=amplifier_channels(i).native_channel_name;
    contvarTemp.data=data(i,:)'/1000;
    contvar{i,1}=contvarTemp;
end
nexFile.contvars=contvar;
%specify correct path ,which is the nex folder in the upper folder. 
currentpath = cd('..');
parentpath = pwd();
cd(currentpath);
[result] = writeNexFile_special(nexFile, [parentpath '\nexFolder\' nexFileName '.nex']);

