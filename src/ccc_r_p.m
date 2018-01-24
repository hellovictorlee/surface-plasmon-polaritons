function out = ccc_r_p(r_p,ycen,xcen,r)

for j=1:2*r+1
    
    for i=1:2*r+1
        
        if(j-r-1)^2+(i-r-1)^2 <= r^2;
        
        r_p(j+ycen-r,i+xcen-r)=2*pi*6.46e12;
        
        end
        
    end
    
end

out=r_p;