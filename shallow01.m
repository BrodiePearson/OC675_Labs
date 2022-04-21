%  OC 675 -- Laboratory exercises one and two:  Linear shallow water model.
%
%  March 2022
%  Jim Lerczak
%
%  For lab one, we solve the one-dimensional advection equation.  For both
%  the initial conditions, I will assume period boundary conditions.

clear all
more off

%  First some preliminaries:
U = 1 ;                %  Background advection speed.
dx = 1 ;               %  Spatial grid spacing.
dtvals = [1 0.5 0.2] ; %  Time steps to use in this lab.
Nx = 200 ;             %  Number of spatial grid points ;

x = dx*(0:1:Nx-1) ;

%  set up the top-hat initial condition
% u0 = 0*x ;
% nn = find((x>=5).*(x<10)) ;
% u0(nn) = 1 ;

%  Gaussian initial condition
u0 = exp(-0.5*((x - 60)/2).^2) ;


%%  Uncentered in space and time
figure(1)
clf

for idt = 1:3 ;
    subplot(3,1,idt) ;
    T = 0 ;
    un = u0 ;
    dt = dtvals(idt) ;
    Nt = 250/dt ;
    for n = 1:Nt ;
        T = T + dt ;
        un = unp1_uncent01(un,dt,dx,U,0) ;
        if abs((T-1))<0.5*dt ;
            plot(x,un,'k','linewidth',1) ;
            hold on
        end
        if abs((T-20))<0.5*dt ;
            plot(x,un,'b','linewidth',1) ;
        end
        if abs((T-100))<0.5*dt ;
            plot(x,un,'g','linewidth',1) ;
        end
    end
end

for ii = 1:3
    subplot(3,1,ii) ;
    axis([0 max(x) -0.25 1.25]) ;
    title(['dt = ' num2str(dtvals(ii))])
end

%%  Centered space and time
figure(2)
clf

for idt = 1:3 ;
    subplot(3,1,idt) ;
    T = 0 ;
    dt = dtvals(idt) ;
    Nt = 250/dt ;
    %  use the uncentered time, uncentered space to get initial un and unm1
    unm1 = u0 ;
    un = unp1_uncent01(u0,dt,dx,U,0) ;
    T = T + dt ;
    for n = 1:Nt ;
        T = T + dt ;
        unp1 = unp1_cent01(unm1,un,dt,dx,U,0) ;
        unm1 = un ;
        un = unp1 ;
        if abs((T-1))<0.5*dt ;
            plot(x,un,'k','linewidth',1) ;
            hold on
        end
        if abs((T-20))<0.5*dt ;
            plot(x,un,'b','linewidth',1) ;
            hold on
        end
        if abs((T-100))<0.5*dt ;
            plot(x,un,'g','linewidth',1) ;
        end
    end
end

for ii = 1:3
    subplot(3,1,ii) ;
    axis([0 max(x) -0.25 1.25]) ;
    title(['dt = ' num2str(dtvals(ii))])
end

%%  Uncentered time, centered space
figure(3)
clf

for idt = 1:3 ;
    subplot(3,1,idt) ;
    T = 0 ;
    un = u0 ;
    dt = dtvals(idt) ;
    Nt = 250/dt ;
    for n = 1:Nt ;
        T = T + dt ;
        un = unp1_uncent02(un,dt,dx,U,0) ;
        if abs((T-1))<0.5*dt ;
            plot(x,un,'k','linewidth',1) ;
            hold on
        end
        if abs((T-20))<0.5*dt ;
            plot(x,un,'b','linewidth',1) ;
        end
        if abs((T-100))<0.5*dt ;
            plot(x,un,'g','linewidth',1) ;
        end
    end
end

for ii = 1:3
    subplot(3,1,ii) ;
    axis([0 max(x) -0.25 1.25]) ;
    title(['dt = ' num2str(dtvals(ii))])
end

%% Boundary condition challenge

figure(4)
clf

idt = 1 ;
T = 0 ;
dt = dtvals(idt) ;
Nt = 250/dt ;

    %  use the uncentered time, uncentered space to get initial un and unm1
unm1c = u0 ;
unc = unp1_uncent01(u0,dt,dx,U,1) ;
unm1p = u0 ;
unp = unp1_uncent01(u0,dt,dx,U,1) ;

T = T + dt ;
for n = 1:Nt ;
    T = T + dt ;
    unp1c = unp1_cent01(unm1c,unc,dt,dx,U,0) ;
    unm1c = unc ;
    unc = unp1c ;
    unp1p = unp1_cent01(unm1p,unp,dt,dx,U,1) ;
    unm1p = unp ;
    unp = unp1p ;

    figure(4)
    clf
    subplot(2,1,1) ;
    plot(x,unc,'k','linewidth',1) ;
    axis([0 max(x) -0.25 1.25]) ;

    title(['T = ' num2str(T) ', clamped'])
    subplot(2,1,2) ;
    plot(x,unp,'k','linewidth',1) ;
    title(['T = ' num2str(T) ', periodic'])
    axis([0 max(x) -0.25 1.25]) ;
    pause(0.01) ;

end