%% import Intan Data.
read_Intan_RHD2000_file
d = fdesign.highpass('Fst,Fp,Ast,Ap',0.015,0.025,60,1);
Hd = design(d,'equiripple');
ch_num=17;
for i=1:ch_num
    i/ch_num
amplifier_data(i,:)=filter(Hd,amplifier_data(i,:));
end
%% import Spike Timing Data.
warning off;
cd 'C:\Users\Spark\Downloads\recording 2016-06-14\recording 2016-06-14\12 mins'
m=importdata('CombinedFiles.xls');
%% Group 1

h=figure('Visible','off');
plot(t_amplifier,amplifier_data(3,:))
hold on;
plot(t_amplifier,amplifier_data(1,:)-100,'r')
hold on;
plot(t_amplifier,amplifier_data(2,:)-200,'g')
hold on;
plot(t_amplifier,amplifier_data(4,:)-300,'k')
hold on;
xlabel('t(s)')


legend('5','3','4','6')

clear store
store=zeros(4,601,sum(m(:,1)==3));

for i=1:size(m,1)
    [i size(m,1)]
if m(i,1)==3 % look into the third channel
    t_start=(m(i,3)+8e-4);
scatter(t_start,100,5,'ro');
hold on;
axis([t_start-5e-3 t_start+25e-3 -400 150])
offsetIndex=find(m(:,1)==3);
offsetIndex = offsetIndex(1)-1;
store(1,:,i-offsetIndex)= amplifier_data(3,t_start*20e3-100:t_start*20e3+500);
store(2,:,i-offsetIndex)= amplifier_data(1,t_start*20e3-100:t_start*20e3+500);
store(3,:,i-offsetIndex)= amplifier_data(2,t_start*20e3-100:t_start*20e3+500);
store(4,:,i-offsetIndex)= amplifier_data(4,t_start*20e3-100:t_start*20e3+500);
%print(h,'-dpdf',[num2str(i) '.pdf'])
end
end

names=dir('*.pdf');
nameList=cell(size(names));
for i=1:numel(nameList)
nameList{i}=names(i,1).name;
end


MeanPlot=mean(store,3);
t_new=((1:601)-100)/20e3;

figure
plot(t_new,MeanPlot(3,:),'LineWidth',2)
hold on;
plot(t_new,MeanPlot(2,:)-70,'r','LineWidth',2)
hold on;
plot(t_new,MeanPlot(1,:)-140,'g','LineWidth',2)
hold on;
plot(t_new,MeanPlot(4,:)-210,'k','LineWidth',2)
hold on;
xlabel('Time(s)')
ylabel('Voltage(uV)')
legend('Ch4','Ch3','Ch5','Ch6')
title(['averaged over ' num2str(size(store,3)) 'occurances of spikes in Channel 5'])
axis([0.005 0.015 -260 50])

% axis([0 2*pi -100.5 100.5])
% print(h,'-dpdf','Figure2.pdf')
delete(nameList{:})
%% Group 2 
h=figure('Visible','off');
plot(t_amplifier,amplifier_data(6,:))
hold on;
plot(t_amplifier,amplifier_data(8,:)-100,'r')
hold on;
plot(t_amplifier,amplifier_data(7,:)-200,'g')
hold on;
xlabel('t(s)')
legend('10','13','11')
clear store
store=zeros(3,601,sum(m(:,1)==6));
for i=1:size(m,1)
    [i size(m,1)]
if m(i,1)==6 % look into the third channel
    t_start=(m(i,3)+8e-4);
scatter(t_start,100,5,'ro');
hold on;
axis([t_start-5e-3 t_start+25e-3 -400 150])
offsetIndex=find(m(:,1)==6);
offsetIndex = offsetIndex(1)-1;
store(1,:,i-offsetIndex)= amplifier_data(6,t_start*20e3-100:t_start*20e3+500);
store(2,:,i-offsetIndex)= amplifier_data(8,t_start*20e3-100:t_start*20e3+500);
store(3,:,i-offsetIndex)= amplifier_data(7,t_start*20e3-100:t_start*20e3+500);
% print(h,'-dpdf',[num2str(i) '.pdf'])
end
end

names=dir('*.pdf');
nameList=cell(size(names));
for i=1:numel(nameList)
nameList{i}=names(i,1).name;
end

MeanPlot=mean(store,3);
t_new=((1:601)-100)/20e3;

figure
plot(t_new,MeanPlot(3,:),'LineWidth',2)
hold on;
plot(t_new,MeanPlot(1,:)-70,'r','LineWidth',2)
hold on;
plot(t_new,MeanPlot(2,:)-140,'g','LineWidth',2)
hold on;
axis([0.005 0.015 -200 50])
xlabel('Time(s)')
ylabel('Voltage(uV)')
legend('Ch11','Ch10','Ch13')
title(['averaged over ' num2str(size(store,3)) 'occurances of spikes in Channel 10'])


% axis([0 2*pi -100.5 100.5])
% print(h,'-dpdf','Figure2.pdf')
%%
append_pdfs('Group_2.pdf',nameList{:})
delete(nameList{:})
%% Group3
h=figure('Visible','off');
plot(t_amplifier,amplifier_data(14,:))
hold on;
plot(t_amplifier,amplifier_data(15,:)-100,'r')
hold on;
plot(t_amplifier,amplifier_data(13,:)-200,'g')
hold on;
plot(t_amplifier,amplifier_data(16,:)-300,'k')
hold on;
xlabel('t(s)')
legend('59','60','58','61')

clear store
store=zeros(4,601,sum(m(:,1)==14));

for i=1:size(m,1)
    [i size(m,1)]
if m(i,1)==14 % look into the third channel
    t_start=(m(i,3)+8e-4);
scatter(t_start,100,5,'ro');
hold on;
axis([t_start-5e-3 t_start+25e-3 -400 150])
offsetIndex=find(m(:,1)==14);
offsetIndex = offsetIndex(1)-1;
store(1,:,i-offsetIndex)= amplifier_data(14,t_start*20e3-100:t_start*20e3+500);
store(2,:,i-offsetIndex)= amplifier_data(15,t_start*20e3-100:t_start*20e3+500);
store(3,:,i-offsetIndex)= amplifier_data(13,t_start*20e3-100:t_start*20e3+500);
store(4,:,i-offsetIndex)= amplifier_data(16,t_start*20e3-100:t_start*20e3+500);
% print(h,'-dpdf',[num2str(i) '.pdf'])
end
end

names=dir('*.pdf');
nameList=cell(size(names));
for i=1:numel(nameList)
nameList{i}=names(i,1).name;
end

MeanPlot=mean(store,3);
t_new=((1:601)-100)/20e3;

figure
plot(t_new,MeanPlot(1,:),'LineWidth',2)
hold on;
plot(t_new,MeanPlot(2,:)-70,'r','LineWidth',2)
hold on;
plot(t_new,MeanPlot(3,:)-140,'g','LineWidth',2)
hold on;
plot(t_new,MeanPlot(4,:)-210,'k','LineWidth',2)
hold on;
axis([0.005 0.015 -260 50])
xlabel('Time(s)')
ylabel('Voltage(uV)')
legend('Ch59','Ch60','Ch58','Ch61')
title(['averaged over ' num2str(size(store,3)) 'occurances of spikes in Channel 59'])
