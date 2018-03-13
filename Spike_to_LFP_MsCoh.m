% From Previous Defined TimePeriod For spike sorting, We extract the
% corresponding LFP,  Let's do a 12 mins averaged MScoherence for 1 to
% 100Hz, pairwise channel according to the channel distribution of that day
read_Intan_RHD2000_file_special
dend=floor(size(amplifier_data,2)/20000)*20000;
amplifier_data=amplifier_data(:,1:dend);

% data=amplifier_data(1,:);
[b,a]=butter(4,300/10000,'high');
[b1,a1]=butter(4,100/10000,'low');
 fprintf(1,'Filtering progress: ')
for i=1:size(amplifier_data,1)
    fprintf(1,'%4.2f',i/size(amplifier_data,1));   

amplifier_data_LFP(i,:) = filtfilt(b1,a1,amplifier_data(i,:));
amplifier_data(i,:) = filtfilt(b,a,amplifier_data(i,:)); 
 fprintf(1,'\b\b\b\b')
end
fprintf('\n')
% med = median(amplifier_data);
% 
% data=med;
% data=amplifier_data(25,:);
%  fil=filtfilt(b,a,data);
 fprintf(1,'feature calculating progress: ')
 for i=1:dend/20000
 fprintf(1,'%4.2f',i/dend*20000);   
 temp=amplifier_data(:,(i-1)*20e3+1:i*20e3);
 feat1(i)=median(mean((temp').^2));
 feat2(i)=max(max(abs(temp)));
 fprintf(1,'\b\b\b\b')
 end
fprintf('\n')
 
%  feat1=mean(newfil.^2);
%  feat1(feat1>400)=400;
%  
%  newfil=reshape(data,20000,numel(data)/20000);
%  feat2=max(abs(newfil));
%  feat2(feat2>600)=600;
 
 mod_feat1=repmat(feat1,20000,1);
 align_feat1=reshape(mod_feat1,dend,1);
 
 mod_feat2=repmat(feat2,20000,1);
 align_feat2=reshape(mod_feat2,dend,1);
 %% Key Part Aligning recording segment used for sorting with that of LFP analysis
load('param.mat')
data=amplifier_data_LFP(:,logical(align_feat1<param.thre_Median_MS)&logical(align_feat2<param.thre_Max_abs));
dur=6;
data=data(:,floor(size(data,2)/2)-dur*60*20e3+1:size(data,2)/2+dur*60*20e3);
%% Pairwise Ch mscohere. 
dF1=1;
dF2=100;
CohStore=XCMC_equal_Fre(data,20E3,hamming(20E3),dF1,dF2);
for ch=size(CohStore,1):-1:2
    for ch2=ch-1:-1:1
    CohStore{ch,ch2}=CohStore{ch2,ch};
    end
end
save('COH','CohStore')
%% Load Recording Specification
mkdir('LFP')
cd('LFP')
currentpath=pwd;
s1=strfind(currentpath,'recording');
if isempty(s1)
    folder_name = uigetdir;
    cd(folder_name)
    currentpath=pwd;
    s1=strfind(currentpath,'recording');
end
s2=strfind(currentpath,'\');
if numel(s2)>find(s2==(s1-1))
mainFolderPath = currentpath(1:s2(find(s2==(s1-1))+1)-1);
else
mainFolderPath = currentpath;
end


   if exist([mainFolderPath '\Dat_V_Map.mat'], 'file')
   load([mainFolderPath '\Dat_V_Map.mat'])
   else 
      read_Intan_RHD2000_file
      x(:,2)=[amplifier_channels.custom_order]';
      x(:,1)=1:size(x,1);
      Dat_V_Map=x;
      save([mainFolderPath '\Dat_V_Map.mat'],'Dat_V_Map')
   end
load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
Ch_Map=Ch_Map_new;

for i=1:size(Ch_Map,1)
for j=1:size(Ch_Map,2)
    if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j), 1))
        Ch_Map(i,j)=0;
    end
end
end

Ch_Map_2 = Ch_Map; 
for i=1:size(Ch_Map,1)
for j=1:size(Ch_Map,2)
if isempty(find(Dat_V_Map(:,2)==Ch_Map(i,j)))
Ch_Map_2(i,j)=0;
else 
Ch_Map_2(i,j)=find(Dat_V_Map(:,2)==Ch_Map(i,j));
end
end
end

%% Plotting 
for CurCh=1:size(data,1)
h=figure;
set(h,'PaperUnits','centimeters','PaperPositionMode','Auto')

Ox=0.6;
Oy=0.6;
gap=0.6;
Gw=2;
Gh=3;
numGx = 8;
numGy = 4;

Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
[X,Y] = meshgrid(Gx,Gy);

t=dF1:dF2;
for row=1:4
    for col=1:8

        ax=axes('Units','centimeters','position', [X(row,col) Y(row,col) Gw Gh]);
        set(ax,'FontSize',8)
        set(ax,'FontWeight','bold')

        if Ch_Map(row,col)~=0
            temp=find(Dat_V_Map(:,2)==Ch_Map(row,col));
            if CurCh~=temp
            plot(t,CohStore{CurCh,temp},'LineWidth',2)
            else 
                 plot(t,0.5*ones(size(t)),'r','LineWidth',10);
            end
            
        end
        axis([1 100 0 1])
        set(gca,'yTick',0:0.2:1)
    end
end
set(gcf,'Position',get(0,'Screensize'))
set(h,'PaperPosition',[0 0 40 27])
print(h,'-dpdf',[num2str(CurCh) '.pdf'],'-loose')
end
%%
DIR=dir('*.pdf');
append_pdfs('result_CoH.pdf',DIR.name)
delete(DIR.name)