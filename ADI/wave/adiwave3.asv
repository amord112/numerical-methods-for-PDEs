clear; clf; fontSize = 15;

%%% Solving the 2D wave equation
%%% Solve u_tt = c^2*(u_xx + u_yy) + f(x,y,t)
%%% using the ADI form of the scheme 
%%% dpmt u_i^n = c^2*(dpmx + dpmy)u_i^{n+1}

tFinal = 1;
c = 1;
ax = 0; bx = 1;
ay = 0; by = 1;
kx = 3*pi;
ky = 2*pi;
plotOption = 1;
ms = 'true';
numResolutions = 5;

%% Define a solution for testing
if strcmp(ms,'true') %homogeneous, zero BC's
    kt = c*sqrt((kx^2) + (ky^2)); 
    ue = @(x,y,t) sin(kx*x).*sin(ky*y).*sin(kt*t);
    ut = @(x,y,t) kt*sin(kx*x).*sin(ky*y).*cos(kt*t);
    f  = @(x,y,t) 0;
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
    u0  = @(x,y) ue(x,y,0);
    u1  = @(x,y) ut(x,y,0);
elseif strcmp(ms,'true1') %homogeneous, nonzero Dirichlet BC's 
    kt = c*sqrt((kx^2) + (ky^2)); 
    ue = @(x,y,t) sin(kx*x).*sin(ky*y).*sin(kt*t);
    ut = @(x,y,t) kt*sin(kx*x).*sin(ky*y).*cos(kt*t);
    f  = @(x,y,t) 0;
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
    u0  = @(x,y) ue(x,y,0);
    u1  = @(x,y) ut(x,y,0);
elseif strcmp(ms,'trig') 
    kx = 3.2*pi; ky = 1.4*pi; kt = 2.3*pi; 
    ue = @(x,y,t) cos(kx*x).*cos(ky*y).*cos(kt*t);
    ut = @(x,y,t) -kt*cos(kx*x).*cos(ky*y).*sin(kt*t);
    utt = @(x,y,t) -kt*kt*cos(kx*x).*cos(ky*y).*cos(kt*t);
    uxx = @(x,y,t) -kx*kx*cos(kx*x).*cos(ky*y).*cos(kt*t);
    uyy = @(x,y,t) -ky*ky*cos(kx*x).*cos(ky*y).*cos(kt*t);
    f  = @(x,y,t) utt(x,y,t) - (c^2)*(uxx(x,y,t) + uyy(x,y,t));
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
    u0  = @(x,y) ue(x,y,0);
    u1  = @(x,y) ut(x,y,0);
elseif strcmp(ms,'poly')
    c2 = 0; d2 = 1; d1 = 1; d0 = 1;
    ue  = @(x,y,t) (x.^2 + x + 1).*(d2*y.^2 + d1*y + d0).*(c2*t^2 + t + 1);
    ut  = @(x,y,t) (x.^2 + x + 1).*(d2*y.^2 + d1*y + d0).*(2*c2*t + 1);
    utt  = @(x,y,t) (x.^2 + x + 1).*(d2*y.^2 + d1*y + d0).*(2*c2);
    uxx = @(x,y,t) 2*(d2*y.^2 + d1*y + d0).*(c2*t^2 + t + 1);
    uyy = @(x,y,t) 2*d2*(x.^2 + x + 1).*(c2*t^2 + t + 1);
    u0  = @(x,y) ue(x,y,0);
    u1  = @(x,y) ut(x,y,0);
    f   = @(x,y,t) utt(x,y,t) - (c^2)*(uxx(x,y,t) + uyy(x,y,t));
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
else
    fprintf('Incorrect input for ms\n');
end

%% Grid refinement study
err = zeros(numResolutions,1);

for m = 1:numResolutions
    Nx = 10*(2^m);
    Ny = Nx;
    dx = (bx-ax)/Nx;
    dy = (by-ay)/Ny;

    iax = 1; ibx = iax + Nx;
    iay = 1; iby = iay + Ny;
    i1x = iax + 1; i1y = iay + 1; % first interior pt in x and y
    i2x = ibx - 1; i2y = iby - 1;

    Ngx = ibx; Ngy = iby; % total number of grid points in x and y
    Ngx1 = Ngx - 2; Ngy1 = Ngy - 2;

    I1 = i1x:i2x; I2 = i1y:i2y; 
    J1 = iax:ibx; J2 = iay:iby;

    x = zeros(Ngx,Ngy); y = x;

    for i = 1:Ngx
        for j = 1:Ngy
            x(i,j) = ax + (i - iax)*dx;
            y(i,j) = ay + (j - iay)*dy;
        end
    end

    dpmx = @(un,I,J) (1/(dx^2))*(un(I + 1,J) - 2*un(I,J) + un(I- 1,J));
    dpmy = @(un,I,J) (1/(dy^2))*(un(I,J+1) - 2*un(I,J) + un(I,J-1));

    % set the time-step
    dt = dx; % typically set dt according to CFL condition
    Nt = round(tFinal / dt);
    dt = tFinal / Nt;

    dMx = form1Dmatrix(Ngx, dx, dt, c, iax, ibx, I1);
    dMy = form1Dmatrix(Ngy, dy, dt, c, iay, iby, I2);

    ry = (c*dt/dy)^2; 
    
    % allocate space for solution
    unm1 = zeros(Ngx,Ngy);
    unp1 = zeros(Ngx,Ngy);
    vnp1 = zeros(Ngx,Ngy);
    rhs  = zeros(Ngx,Ngy);

    unm1 = u0(x,y);
    un = dt*u1(x,y) + unm1;

    fn = zeros(Ngx,Ngy);
    rhs = zeros(Ngx,Ngy);

    %% Start the time-stepping
    t = 0;

    for n = 2:Nt
        t = n*dt;

        I = I1; J = J2;
        fn(I,J) = f(x(I,J), y(I,J), t);

        rhs(I,J) = 2*un(I,J) - unm1(I,J) + (dt^2)*fn(I,J);
        rhs(iax,:) = gax(y(iax,:),t) - ry*(gax((y(iax,:) + dy),t) - 2*gax(y(iax,:),t) + gax((y(iax,:) - dy),t));
        rhs(ibx,:) = gbx(y(ibx,:),t) - ry*(gbx((y(ibx,:) + dy),t) - 2*gbx(y(ibx,:),t) + gbx((y(ibx,:) - dy),t));

        vnp1 = dMx\rhs;

        I = J1; J = I2;
        rhs(I,J) = vnp1(I,J);
        rhs(:,iay) = gay(x(:,iay),t);
        rhs(:,iby) = gby(x(:,iby),t);

        unp1 = (dMy\rhs')';

        ix=iax; un(ix,J2)=gax(y(ix,J2),t);
        ix=ibx; un(ix,J2)=gbx(y(ix,J2),t); 

        unm1 = un; 
        un = unp1;
    end

    uexact = ue(x,y,t);
    err(m) = max(max(abs(uexact - un)));

    fprintf('%s: t = %1.2f, Nx = %3d, Nt = %3d, dt = %1.2e, Max-err = %1.3e', ms, t, Nx, Nt, dt, err(m));
    if m > 1
        rate = log(err(m-1) / err(m)) / log(2);
        fprintf(' rate = %1.2f\n', rate);
    else
        fprintf('\n');
    end

    if plotOption == 1
        figure(1)
        surf(x, y, un); colorbar;
        title(sprintf('2D Wave Equation %s t=%9.3e (Nx=%d)', ms, t, Nx))
        xlabel('x'); ylabel('y'); zlabel('u'); set(gca, 'FontSize', fontSize); grid on;

        figure(2)
        surf(x, y, uexact - un); colorbar;
        title(sprintf('2D Wave Equation Error %s t=%9.3e (Nx=%d)', ms, t, Nx))
        xlabel('x'); ylabel('y'); zlabel('u'); set(gca, 'FontSize', fontSize); grid on;
    end
end

function dMz = form1Dmatrix(Ngz, dz, dt, c, iaz, ibz, I)

    Mz = sparse(Ngz,Ngz);
    rz = (c*dt/dz)^2;

    for i = I
         Mz(i,i-1) = -rz;
         Mz(i,i  ) = 1 + 2*rz;
         Mz(i,i+1) = -rz;
    end

    i = iaz; Mz(i,i) = 1; 
    i = ibz; Mz(i,i) = 1;

    dMz = decomposition(Mz,'lu');

end