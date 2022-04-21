clear all
more off

NL = 0 ;
T = 400 ;
%  Linear dt's:  4, 0.5, 0.1
dt = 0.1 ;
Nt = round(T/dt) ;
%  Nt = 1000 ;

%  domain preliminaries:
H = 5 ;
Lx = 6000 ;
dx = 20 ;
xu = 0:dx:Lx ;
xh = dx/2:dx:Lx+dx/2 ;
nmax = length(xu) ;

g = 9.81 ;                %  gravity
nu = 0 ;                 %  eddy viscosity
%  nu = 100 ;                 %  eddy viscosity

htheory = NaN*xu ;
utheory = NaN*xu ;

%  stencils for periodic BCs.
jneg = [nmax 1:nmax-1] ;
jpos = [2:nmax 1] ;

%  initial conditions
up = 0*xu ;
u = up ;
hp = 0*xh ;
h = hp ;
%  linear problem variables
if NL == 0
    disp(['linear simulation']) ;
    disp(' ') ;
    %  bump parameters
    hmin = 0.25 ;
    mu = 100 ;
    x0 = Lx/2 ;
    %  advection velocity (3 or 10)
    U = 3 ;

    c = sqrt(g*H) ;
    Fr = U/c ;
    
    %  Upstream Froude number
    disp(['Froude # = ' num2str(Fr)])
    disp([' ']) ;

    disp(['dx/(c + U) = ' num2str(dx/(c + U)) ' s'])
    disp(['dt = ' num2str(dt) ' s'])
    disp(['CFL:  (c + U)dt/dx = ' num2str((c+U)*dt/dx)])

    ylim = [-0.5 0.5] ;

    hsu = hmin*exp(-0.5*((xu - x0)/mu).^2) ; 
    hsh = hmin*exp(-0.5*((xh - x0)/mu).^2) ;
    hsxu = -((xu - x0)/mu^2).*hsu ;
    hsxh = -((xh - x0)/mu^2).*hsh ;

    htheory = (Fr^2)/(Fr^2 - 1)*hsu ;
    utheory = U*H./(htheory + H - hsu) - U ;
   
else
    disp(['nonlinear simulation']) ;
    disp(' ') ;
    %  bump parameters
    hmin = 0.9*H ;
    mu = 100 ;
    x0 = Lx/2 ;
    U = 1 ;

    c = sqrt(g*H) ;
    %  Upstream Froude number
    disp(['upstream Froude = ' num2str(U/c)])
    %  Initial Froude number at bump
    disp(['initial Froude number at bump = ' num2str( (U*H/(H-hmin)) / sqrt(g*(H-hmin)))]) ;
    disp([' ']) ;

    disp(['dx/(c + U) = ' num2str(dx/(c + U)) ' s'])
    disp(['dt = ' num2str(dt) ' s'])
    disp(['CFL:  (c + U)dt/dx = ' num2str((c+U)*dt/dx)])

    ylim = [-5 5] ;

    hsu = hmin*exp(-0.5*((xu - x0)/mu).^2) ; 
    hsh = hmin*exp(-0.5*((xh - x0)/mu).^2) ;
    hsxu = -((xu - x0)/mu^2).*hsu ;
    hsxh = -((xh - x0)/mu^2).*hsh ;
   
end

figure(1)
clf

tm = 0 ;
nplt = round(1/dt) ;
if nplt < 1
    nplt = 1 ;
end
for nt = 1:Nt ;

    if NL == 0
        uRHS = (1/dx)*(-(U/2)*(u(jpos) - u(jneg)) - g*(h - h(jneg))) + (nu/dx^2)*(u(jpos) - 2*u + u(jneg)) ;
        hRHS = (1/dx)*(-H*(u(jpos) - u) - (U/2)*(h(jpos) - h(jneg)) + U*hsxh*dx) + (nu/dx^2)*(h(jpos) - 2*h + h(jneg)) ;
    else
        arg = 0.5*(u + U).^2 ;
        uRHS = (1/dx)*(-(1/2)*(arg(jpos) - arg(jneg)) - g*(h - h(jneg))) + (nu/dx^2)*(u(jpos) - 2*u + u(jneg)) ;
        uh = U + 0.5*(u + u(jpos)) ;
        arg = uh.*(h + H - hsh) ;
        hRHS = -(1/(2*dx))*(arg(jpos) - arg(jneg)) + (nu/dx^2)*(h(jpos) - 2*h + h(jneg)) ;
    end

    tm = tm + dt ;
    up = u + dt*uRHS ;
    hp = h + dt*hRHS ;

    if mod(nt,nplt)==0
        figure(1)
        clf
        plot(xu,htheory,'b--') ;
        hold on
        plot(xu,utheory,'g--') ;
        plot(xu,up,'k') ;
        plot(xh,hp,'r') ;
        axis([0 Lx ylim])
        title(['t = ' num2str(tm) ' s']) ;
        pause(0.001) ;
    end

    u = up ;
    h = hp ;

end


return


uinit = 9;
hinit = 400.
delt = 2.
%
imax = 400;
u(1:imax) = uinit;
h(1:imax) = hinit;
um = u;
up = u;
htop = 200;
amnt = 2000;
gstr = 9.81*(20./1000);
alpha = 1.0;
nend = 1000;
dx = 200.0;
xx = [1:imax]*dx;
% delt = 0.2*dx/uinit
% nmax = (imax*dx)/(uint*delt)
nmax = 2000;
dhsdx(1:imax) = 0.0;
hter(1:imax) = 0.0;
for i=1:imax
  dhsdx(i) = -2.*(i-imax/2)*dx*htop*amnt*amnt/(amnt^2+((i-imax/2)*dx)^2)^2;
  hter(i) = htop*amnt*amnt/(amnt^2+((i-imax/2)*dx)^2);
end
% h = h-hter
dhsdx
pause(0.2);
hm = h;
hp = h;
ip = [2:imax 1];
in = [imax 1:imax-1];
% 
for nn=1:nmax
 ub = (u(in)+u)*0.5;

   udhdx = 0.5*((u+abs(u)).*(h-h(in))+...
                (u(ip)-abs(u(ip))).*(h(ip)-h))/dx;
%  udhdx = uinit.*(h-h(in))/dx;
%
hp = h-delt*(udhdx+hinit*(u(ip)-u)/dx);
%
%  ududx = 0.5*((u([imax 1:imax-1])+abs(u([imax 1:imax-1]))).*...
%              (u(2:imax 1]-u([2:imax 1]))+...
%              (ubar([2:imax 1])-abs(ubar([2:imax 1]))).*(u([2:imax 1])-u))/dx;

 ududx = 0.5*((ub+abs(ub)).*(u-u(in))+(ub(ip)-abs(ub(ip))).*(u(ip)-u))/dx;
%  ududx = 0.5*(ub(ip).^2-ub.^2)/dx;
%  ududx = uinit*(u-u(in))/dx;

 up = u - delt* (ududx + gstr*((h-h(in))/dx+dhsdx));
%
%  up = up + delt*100*(u([2:imax 1])+u([imax 1:imax-1])-2.*...
%                  u)/(dx*dx);
%  hp = hp + delt*100*(h([2:imax 1])+h([imax 1:imax-1])-2.*...
%                  h)/(dx*dx);
% u = u + 0.1*(up-2*u+um);
% h = h + 0.1*(hp-2*h+hm);
% um = u;
 u = up;
 h = hp;
 plot(xx(1:imax),h(1:imax)+hter(1:imax),'-k',xx(1:imax),hter(1:imax),'--k');
%  hold;
% plot(xx,u*10.,'--r')
%  hold;
%  pause;
%  plot(xx(1:imax),h(1:imax)+hter);
% axis([0 400 -1.5 1.5);
xlabel('Distance (m)','fontsize',14);
ylabel('h field','fontsize',14);
sum(u.*u)
pause(0.01);
% pause(0.15);
% plot(xx(1:imax),p(1:imax));
end

