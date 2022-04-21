%  lab03.m
%  Jim Lerczak
%  OC675
%
%  14 April 2022

clear all
more off

NL = 0 ;
T = 400 ;
%  Linear dt's:  4, 0.5, 0.1
dt = 0.25 ;
Nt = round(T/dt) ;
Nt = 2000 ;

U = 3 ;
V = 0 ;

g = 9.81 ;                %  gravity
K = 25 ;                  %  eddy viscosity
%  K = 100 ;                 %  eddy viscosity

%  domain preliminaries:
H = 5 ;
Lx = 6000 ;
Ly = Lx/10 ;
dx = 20 ;
dy = dx ;

hmtn = 0.25 ;
wmtn = 100 ;
x0 = Lx/2 ;
y0 = Ly/2 ;

xu = 0:dx:Lx ;
yu = ((0+dy/2):dy:Ly+dy/2)' ;

xv = (0:dx/2):dx:Lx+dx/2 ;
yv = (0:dy:Ly)' ;

xh = xv ;
yh = yu ;

%  stencils for periodic BCs.
Nu = length(xu) ;
Mu = length(yu) ;
imu = [Nu 1:Nu-1] ;
ipu = [2:Nu 1] ;
jmu = [Mu 1:Mu-1] ;
jpu = [2:Mu 1] ;
Nv = length(xv) ;
Mv = length(yv) ;
imv = [Nv 1:Nv-1] ;
ipv = [2:Nv 1] ;
jmv = [Mv 1:Mv-1] ;
jpv = [2:Mv 1] ;
Nh = length(xh) ;
Mh = length(yh) ;
imh = [Nh 1:Nh-1] ;
iph = [2:Nh 1] ;
jmh = [Mh 1:Mh-1] ;
jph = [2:Mh 1] ;

hsu = hmtn*exp(-0.5*((yu - y0)/wmtn).^2)*exp(-0.5*((xu - x0)/wmtn).^2) ;
hsv = hmtn*exp(-0.5*((yv - y0)/wmtn).^2)*exp(-0.5*((xv - x0)/wmtn).^2) ;
hsh = hmtn*exp(-0.5*((yh - y0)/wmtn).^2)*exp(-0.5*((xh - x0)/wmtn).^2) ;
% hsu = hmtn*ones(size(yu))*exp(-0.5*((xu - x0)/wmtn).^2) ;
% hsv = hmtn*ones(size(yv))*exp(-0.5*((xv - x0)/wmtn).^2) ;
% hsh = hmtn*ones(size(yh))*exp(-0.5*((xh - x0)/wmtn).^2) ;

%  initial conditions
u = zeros(Mu,Nu) ;
up = u ;
v = zeros(Mv,Nv) ;
vp = v ;
h = zeros(Mh,Nh) ;
hp = h ;


figure(1)
clf

tm = 0 ;
nplt = round(1/dt) ;
if nplt < 1
    nplt = 1 ;
end
for nt = 1:Nt ;

    if NL == 0
        uRHS =    -U*(u(:,ipu)-u(:,imu))/(2*dx) ...
               -   V*(u(jpu,:) - u(jmu,:))/(2*dy) ...
               -   g*(h - h(:,imh))/dx ;
               +   K*(u(:,ipu) - 2*u + u(:,imu))/(dx^2) ...
               +   K*(u(jpu,:) - 2*u + u(jmu,:))/(dy^2) ;

        vRHS =    -U*(v(:,ipv)-v(:,imv))/(2*dx) ...
               -   V*(v(jpv,:) - v(jmv,:))/(2*dy) ...
               -   g*(h - h(jmh,:))/dx ;
               +   K*(v(:,ipv) - 2*v + v(:,imv))/(dx^2) ...
               +   K*(v(jpv,:) - 2*v + v(jmv,:))/(dy^2) ;

        hRHS =    -U*((h(:,iph) - hsh(:,iph)) - (h(:,imh) - hsh(:,imh)))/(2*dx) ...
               -   V*((h(jph,:) - hsh(jph,:)) - (h(jmh,:) - hsh(jmh,:)))/(2*dy) ...
               -   H*((u(:,ipu) - u)/dx + (v(jpv,:) - v)/dy) ...
               +   K*(h(:,iph) - 2*h + h(:,imh))/(dx^2) ...
               +   K*(h(jph,:) - 2*h + h(jmh,:))/(dy^2) ;

    else
        disp(['blah']) ;
    end

    tm = tm + dt ;
    up = u + dt*uRHS ;
    vp = v + dt*vRHS ;
    hp = h + dt*hRHS ;

    if mod(nt,nplt)==0
        figure(1)
        clf
        pcolor(xh,yh,h)
        axis('equal')
        caxis([-.02 .02]) ;
        colormap redblue
        shading flat
        axis([0 Lx 0 Ly])
        title(['t = ' num2str(tm) ' s']) ;
        pause(0.001) ;
    end

    u = up ;
    v = vp ;
    h = hp ;

end


return
