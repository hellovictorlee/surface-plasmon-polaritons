function out = ccc_eparray(eparray,ycen,xcen,r)

for j=1:2*r+1
    
    for i=1:2*r+1
        
        if(j-r-1)^2+(i-r-1)^2 <= r^2;
        
        eparray(j+ycen-r,i+xcen-r)=11.575*8.85e-12;
        
        end
        
    end
    
end

out=eparray;