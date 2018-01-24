function out = ccc_omega_p(omega_p,ycen,xcen,r)

for j=1:2*r+1
    
    for i=1:2*r+1
        
        if(j-r-1)^2+(i-r-1)^2 <= r^2;
        
        omega_p(j+ycen-r,i+xcen-r)=2*pi*2.183e15;
        
        end
        
    end
    
end

out=omega_p;