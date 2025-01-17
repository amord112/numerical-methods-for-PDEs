function b = getRHS(un,unm1,J1,J2,I1,I2,Ngx,Ngy,dt,fn,eta,dpmx,dpmy,c)

b = zeros(Ngx,Ngy);
r = c^2*eta*dt^2;

b(J1,I2) = 2.*un(J1,I2) - unm1(J1,I2) + dpmy(un,J1,I2).*(c^2*dt^2 - 2*r) + r*dpmy(unm1,J1,I2);
b(I1,I2) = b(I1,I2) + dpmx(un,I1,I2).*(c^2*dt^2 - 2*r) + r*dpmx(unm1,I1,I2) + fn(I1,I2);


end

