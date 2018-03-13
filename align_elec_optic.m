%[b,a]=butter(4,[1044/10E3,1048/10E3],'stop');
%fildata1=filtfilt(b,a,amplifier_data(1,:));
clear all
close all
read_Intan_RHD2000_file
amplifier_data=amplifier_data(1,:);% decrease the total memory consumed.
[pks,locs] = findpeaks(amplifier_data, 'MinPeakHeight',2800);
tp=(locs-5)/20E3;

store=zeros(24,20);
for i=1:24
store(i,:)=amplifier_data(1,locs(i)-10:locs(i)+9);
store(i,:)=store(i,:)-store(i,1);
end

x=tp;

%%
load('time.mat')
dif=diff(realTime);
tKick=find(dif>0.15);
realKO=realTime(tKick+1);

y=realKO;
%%
p = polyfit(x,y,1)
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal

%% we ultimately need to transform some time pts in the intan output excel file...
... with p(1)*x + b
clearvars -except p



% Investigation finds a perfect linear relationship between the laser off
% timepoints and the next laser peak time points. 
