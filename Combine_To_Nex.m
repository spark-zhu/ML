% the full package would need a 
% combine_Intan_RHD2000_files_V4.m,
% read_Intan_RHD2000_file.m
% nexFile_temp.mat
% writeNexFile.m 
% to run successfully. 

z=dir('*.rhd');
nexFileName=z(1).name(1:end-11);
combine_Intan_RHD2000_files_V4
read_Intan_RHD2000_file
load('C:\Users\Spark\Box Sync\ML_code\nexFile_temp.mat')
nexFile.freq=frequency_parameters.amplifier_sample_rate;
nexFile.tend=size(amplifier_data,2)/frequency_parameters.amplifier_sample_rate;
contvarTemp=nexFile.contvars;
contvarTemp.ADFrequency=nexFile.freq;
contvar=cell(size(amplifier_data,1),1);
for i=1:size(amplifier_data,1)
    contvarTemp.name=[amplifier_channels(i).native_channel_name '_' amplifier_channels(i).custom_channel_name];
    contvarTemp.data=amplifier_data(i,:)';
    contvar{i,1}=contvarTemp;
end
nexFile.contvars=contvar;
%specify correct path ,which is the nex folder in the upper folder. 
currentpath = cd('..');
parentpath = pwd();
cd(currentpath);
[result] = writeNexFile(nexFile, [parentpath '\nexFolder\' nexFileName '.nex']);




