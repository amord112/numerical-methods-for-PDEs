function A = getImplicitMatrix(Ngx,ia,ib,dt,c,dx,eta)

A = sparse(Ngx,Ngx); 
r = (dt^2)*(c^2)/(dx^2); 

for i = ia+1:ib-1
    A(i,i-1) = -eta*r;
    A(i,i  ) = 1 + 2*eta*r;
    A(i,i+1) = -eta*r;
end

A(ia,ia) = 1; 
A(ib,ib) = 1; 

A = decomposition(A,'LU'); 

end