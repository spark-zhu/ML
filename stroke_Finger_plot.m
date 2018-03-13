windowLen = 0.1*20E3; % window length in sample for Firing rate calculation
v= 1/windowLen*20E3*ones(1,windowLen);
oF = 100; % downsample fold. 
ignoredUnitNumber = [12 15 17 20 21 22 24 33 40 42 19]; % noise unit eliminated. 
load('10182017_tracking.mat')
    load('Finger_map.mat');
Ch_Map= Ch_Map_new-15;
PeakCh   = cellfun(@(x)  find(x==max(x)),result.P2P) ;



% validUnit = setdiff(result.list,ignoredUnitNumber);


Global_Time_Store = cell(4,8); % stored all the time stamps organized per channel  
Global_Unit_Store = cell(4,8); % store all the raw unit number(cluster label) per channel. 

for row=1:4
    for col=1:8
        unitOfThatCh = result.list(PeakCh==Ch_Map(row,col));
        validUnit = setdiff(unitOfThatCh,ignoredUnitNumber);

    if ~isempty(validUnit)
    Global_Unit_Store{row,col} = validUnit;
    TimeStore = cell(1,numel(validUnit));    
    for unit=1:numel(validUnit)
    
    
sample1 = result.time{result.list==validUnit(unit)}; % major mistake corrected here.  




binVersion=zeros(size(1:sample1(end)));
binVersion(sample1)=1;
w = conv(binVersion,v,'same'); 
w = downsample(w,oF);

     TimeStore{unit}=w;       
            
            
    end
        Global_Time_Store{row,col}=TimeStore;
    end
    
    end
end

%%

    h=figure;
        

    stepSize = 400;
    for finger=1:4
    subplot(1,4,finger)

    AllCh_record=Global_Time_Store(finger,:);
    for ele=1:8
        
     if isempty(AllCh_record{ele})
     plot([0 4000],[0 0]+(ele-1)*stepSize,'lineWidth',3)
     hold on
     end
     
     if numel(AllCh_record{ele})==1
         timeLine = AllCh_record{ele}{1};
     plot((1:numel(timeLine))./(20E3/oF),timeLine+(ele-1)*stepSize,'lineWidth',1)
     hold on
     end
     
     
     if numel(AllCh_record{ele})==2
         timeLine = AllCh_record{ele}{1};
     plot((1:numel(timeLine))./(20E3/oF),timeLine+(ele-1)*stepSize,'lineWidth',1)
     hold on
     timeLine = AllCh_record{ele}{2};
     plot((1:numel(timeLine))./(20E3/oF),timeLine+(ele-0.5)*stepSize,'lineWidth',1)
     hold on
     end
           
     
    
    end
    axis([2000 4100 0 stepSize*(ele+2)])
    end
% title(num2str(range))
set(gcf, 'Position', get(0,'Screensize'));
% saveas(h,num2str(range))
