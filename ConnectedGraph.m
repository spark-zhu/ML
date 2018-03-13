function finalAnswer = ConnectedGraph (M)
finalAnswer = '(0,0)';
for i=1:size(M,1)
    for j=1:size(M,2)
        if M(i,j)~=0
        a=[];
        b=[];
        c=[];
        d=[];
        e=[];
        
        if j+1<=size(M,2)
        a=M(i,j+1);
        end
        
        if i+1<=size(M,1)
        b=M(i+1,j);
        end
        
        if j-1>=1
        e=M(i,j-1);
        end
        
        if ~isempty(a)&& ~isempty(b)
            
        c=M(i+1,j+1);
            
        end
        
        if ~isempty(e)&& ~isempty(b)
          
            d=M(i+1,j-1);
    
        end
        
        if ~isempty(a)&&a~=0
          finalAnswer=strcat(finalAnswer,[', (' num2str(M(i,j)-1) ', ' num2str(a-1) ') ']);
        end
        if ~isempty(c)&&c~=0
          finalAnswer=strcat(finalAnswer,[', (' num2str(M(i,j)-1) ', ' num2str(c-1) ') ']);
        end
        if ~isempty(b)&&b~=0
          finalAnswer=strcat(finalAnswer,[', (' num2str(M(i,j)-1) ', ' num2str(b-1) ') ']);
        end
        if ~isempty(d)&&d~=0
          finalAnswer=strcat(finalAnswer,[', (' num2str(M(i,j)-1) ', ' num2str(d-1) ') ']);
        end
        end 
    end
end
end