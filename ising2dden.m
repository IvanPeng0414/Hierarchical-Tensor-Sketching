function [E] = ising2dden(J1,J2,d,n,x)
%x:matrix(d^2,1), need to reshape
E   = 0;
tmp = -1:2/(n-1):1;
x   = tmp(x); %-1 or 1
x = reshape(x,[d,d]);


for j = 1 : d
    for i = 1 : d
        if i==d-1
            E = E - J1*x(i,j)*x(i+1,j);
        elseif i<=d-2
            E = E - J1*x(i,j)*x(i+1,j) - J2*x(i,j)*x(i+2,j);
        else
            E = E - J1*x(d,j)*x(1,j);
        end
        if j == d-1
            E = E - J1*x(i,j)*x(i,j+1);
        elseif j <=d-2
            E = E - J1*x(i,j)*x(i,j+1) - J2*x(i,j)*x(i,j+2);
        else
            E = E - J1*x(i,d)*x(i,1);
        end
            
    end
    
end


end