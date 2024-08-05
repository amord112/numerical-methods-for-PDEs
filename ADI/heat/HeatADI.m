clear; clf; fontSize = 15;

%%% Solving the 2D heat equation
%%% Solve u_t = kappa*(u_xx + u_yy) + f(x,y,t)
%%% using BE-CD2 scheme

kappa = 0.05; % diffusion coefficient
tFinal = .5;
ax = 0; bx = 1;
ay = 0; by = 1;
kx = 3*pi;
ky = 2*pi;
plotOption = 1;
ms = 'true';
numResolutions = 3;

%% Define a solution for testing
if(strcmp(ms,'true'))
    ue = @(x,y,t) sin(kx*x).*sin(ky*y).*exp(-kappa*((kx^2) + (ky^2))*t);
    u0 = @(x,y) ue(x,y,0);
    f  = @(x,y,t) 0;
    gax = @(y,t) 0;
    gbx = @(y,t) 0;
    gay = @(x,t) 0;
    gby = @(x,t) 0;
elseif(strcmp(ms,'poly'))
    c2 = 0;
    ue  = @(x,y,t) (x.^2 + x + 1).*(y.^2 + y + 1).*(c2*t^2 + t + 1);
    ut  = @(x,y,t) (x.^2 + x + 1).*(y.^2 + y + 1).*(2*c2*t + 1);
    uxx = @(x,y,t) 2*(y.^2 + y + 1).*(c2*t^2 + t + 1);
    uyy = @(x,y,t) 2*(x.^2 + x + 1).*(c2*t^2 + t + 1);
    u0  = @(x,y) ue(x,y,0);
    f   = @(x,y,t) ut(x,y,t) - kappa*(uxx(x,y,t) + uyy(x,y,t));
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
else
    fprintf('incorrect input for ms\n');
end

%% Grid refinement study
err = zeros(numResolutions,1);

for m = 1:numResolutions
    Nx = 10*(2^m);
    Ny = Nx;
    dx = (bx - ax)/Nx;
    dy = (by - ay)/Ny;

    iax = 1; ibx = iax + Nx;
    iay = 1; iby = iay + Ny;
    Ngx = ibx; Ngy  = iby; % total number of grid points in x and y
    Ngx1 = Ngx - 2; Ngy1 = Ngy - 2;
    I1 = iax + 1:ibx - 1; J1 = iax:ibx;
    I2 = iay + 1:iby - 1; J2 = iay:iby;

    x = zeros(Ngx,Ngy); y = x;

    for i = 1:Ngx
        for j = 1:Ngy
            x(i,j) = ax + (i - iax)*dx;
            y(i,j) = ay + (j - iay)*dy;
        end
    end

    % set the time-step
    dt = dx;
    Nt = round(tFinal/dt);
    dt = tFinal/Nt;

    % form the lhs matrix 
    A = sparse(Ngx1,Ngx1);
    A1 = sparse(Ngy,Ngy);

    for i = iax+1:ibx-3
        A(i,i-1) = 1;
        A(i,i  ) = -2;
        A(i,i+1) = 1;
    end

     for i = iay+1:iby-1
        A1(i,i-1) = 1;
        A1(i,i  ) = -2;
        A1(i,i+1) = 1;
    end

    rx = dt*kappa/(dx^2);
    ry = dt*kappa/(dy^2);
    
    Mx = sparse(Ngx1, Ngx1);
    Mx = eye(Ngx1,Ngx1)-rx.*A;
    dMx = decomposition(Mx,'LU');


    My = sparse(Ngy1, Ngy1);
    My = eye(Ngy1,Ngy1)-ry.*A;
    dMy = decomposition(My,'LU');
    
    My1 = sparse(Ngy, Ngy);
    My1 = eye(Ngy,Ngy)-ry.*A1;
    
    % allocate space for solution
    unp1 = zeros(Ngx,Ngy);
    vnp1 = zeros(Ngx,Ngy);
    un = u0(x,y);

    %% Start the time-stepping
    t = 0;
    rhs = zeros(Ngx-2,Ngy-2);
    for n = 1:Nt

        t = n*dt;
        fn = f(x(I1,I2),y(I1,I2),t); 
        rhs = un(I1,I2) + dt*fn; 
     
        vnp1(I1,I2) = dMx\rhs;
   

        vnp1(iax,:) = 0;%My1.*gax(y(iax,:),t);
        vnp1(ibx,:) = 0;%My1.*gbx(y(ibx,:),t);


        unp1 = dMy\vnp1(I1,I2);

        un = unp1;

        un(iax,:) = gax(y(iax,:),t);
        un(ibx,:) = gbx(y(ibx,:),t);
        un(:,iay) = gay(x(:,iay),t);
        un(:,iby) = gby(x(:,iby),t);

    end

    uexact = ue(x,y,t);
    err(m) = max(max(abs(uexact - un)));

     fprintf('%s: t = %1.2f, Nx = %3d, Nt = %3d, dt = %1.2e, Max-err = %1.3e', ms, t, Nx, Nt, dt, err(m));
    if(m>1)
        rate = log(err(m-1)/err(m))/log(2);
        fprintf(' rate = %1.2f\n', rate);
    else
        fprintf('\n');
    end

    if(plotOption == 1)
        figure(1)
        surf(x,y,un); colorbar;
        title(sprintf('2D Heat Equation %s t=%9.3e (Nx=%d)', ms, t, Nx))
        xlabel('x'); ylabel('y'); zlabel('u');set(gca, 'FontSize', fontSize); grid on;

        figure(2)
        surf(x,y,uexact - un); colorbar;
        title(sprintf('2D Heat Equation Error %s t=%9.3e (Nx=%d)', ms, t, Nx))
        xlabel('x'); ylabel('y'); zlabel('u');set(gca, 'FontSize', fontSize); grid on;
    end
end