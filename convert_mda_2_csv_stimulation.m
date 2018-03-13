Fs=20E3;
mdaStrut = dir('**/firings.mda');
mdaName =[mdaStrut(1).folder '\firings.mda'];
% [file, path, ~] = ...
% uigetfile('*.kwik', 'Select an KWIK', 'MultiSelect', 'off');
% KwikName=[path,file];%dir('*.kwik');
A=readmda(mdaName,0);

cluster=A(3,:);
MyTimes=A(2,:);



csvOut(:,1)=cluster';
csvOut(:,2)=1;
csvOut(:,3)=MyTimes'/Fs;
csvOut(:,7:46)=0;
%%

csvwrite('stimulation.csv',csvOut)