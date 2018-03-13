function CohStore=XCMC_equal_Fre(data,winLen,win,desiredF1,desiredF2)
% assume winLen=Fs; desiredF!=0;
final=zeros(1,desiredF2);
sec=size(data,2)/winLen;
window=repmat(win,1,sec);
Klist=cell(1,desiredF2);
for F=1:numel(Klist) % frequency increment of 1 
k=0:(winLen-1);
k1=k*-1i*2*pi*F/winLen;
Klist{F}=exp(k1);
end

SigReshapeList = cell(size(data,1),1);

for ch=1:numel(SigReshapeList)
    listx=reshape(data(ch,:)',winLen,sec);
    listx=window.*listx;
    SigReshapeList{ch}=listx;
end

CohStore = cell(size(data,1));
for ch1=1:size(data,1)-1
    for ch2=ch1+1:size(data,1)
        [ch1 ch2]
        
        
        
      for F1=desiredF1:desiredF2  % through all croos comparison 
        F2=F1;
            
b=sum(abs(Klist{F1}*SigReshapeList{ch1}).^2)/sec;
c=sum(abs(Klist{F2}*SigReshapeList{ch2}).^2)/sec;

a=abs(sum((Klist{F1}*SigReshapeList{ch1}).*conj(Klist{F2}*SigReshapeList{ch2}))/sec).^2;

final(F1)=a/b/c;
        
     end  
           
    CohStore{ch1,ch2}=final;
    end
end
   
    
