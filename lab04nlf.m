%  lab04nlf.m
%  Jim Lerczak
%  OC675
%
%  14 April 2022

clear all
more off

%colormap for pcolor
cmap1 = flipud(cbrewer('div','RdBu',100,'linear'));

T = 400 ;
%  Linear dt's:  4, 0.5, 0.1
dt = 0.2 ;
Nt = round(T/dt) ;

g = 9.81 ;                   %  gravity
K = 5 ;                     %  eddy viscosity
%  K = 100 ;                 %  eddy viscosity

f = 5e-2 ;

%  domain preliminaries:
H = 5 ;
Lx = 6000 ;
Ly = Lx ;
dx = 20 ;
dy = dx ;

%  parameters for bottom bathymetry
hmtn = 0 ;
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

%%  set-up bottom bathymetry
hsu = hmtn*exp(-0.5*((yu - y0)/wmtn).^2)*exp(-0.5*((xu - x0)/wmtn).^2) ;
hsv = hmtn*exp(-0.5*((yv - y0)/wmtn).^2)*exp(-0.5*((xv - x0)/wmtn).^2) ;
hsh = hmtn*exp(-0.5*((yh - y0)/wmtn).^2)*exp(-0.5*((xh - x0)/wmtn).^2) ;

%%  set-up initial conditions
prob = 6 ;  

if prob == 1 ;  %  linear geostrophic adjustment
    NL = 0 ;
    U = 0 ;
    V = 0 ;
    u = zeros(Mu,Nu) ;
    up = u ;
    v = zeros(Mv,Nv) ;
    vp = v ;
    h0 = 4 ;
    wh = 100 ;
    x0 = Lx/2 ;
    y0 = Ly/2 ;
    h = h0*exp(-0.5*((yh - y0)/wh).^2)*exp(-0.5*((xh - x0)/wh).^2)  ;
    hp = h ;

elseif prob == 2 ;  %  nonlinear geostrophic adjustment
    NL = 1 ;
    U = 0 ;
    V = 0 ;
    u = zeros(Mu,Nu) ;
    up = u ;
    v = zeros(Mv,Nv) ;
    vp = v ;
    h0 = 4 ;
    wh = 100 ;
    x0 = Lx/2 ;
    y0 = Ly/2 ;
    h = h0*exp(-0.5*((yh - y0)/wh).^2)*exp(-0.5*((xh - x0)/wh).^2)  ;
    hp = h ;

elseif prob == 3 ;  %  balanced initial flow
    NL = 0 ;
    U = 0 ;
    V = 0 ;
    h0 = 2 ;
    wh = 100 ;
    x0 = Lx/2 ;
    y0 = Ly/2 ;
    u =  (g/f)*h0*(((yu - y0)/wh^2).*exp(-0.5*((yu - y0)/wh).^2))*exp(-0.5*((xu - x0)/wh).^2)  ;
    up = u ;
    v = -(g/f)*h0*exp(-0.5*((yv - y0)/wh).^2)*(((xv - x0)/wh^2).*exp(-0.5*((xv - x0)/wh).^2))  ;
    vp = v ;
    h = h0*exp(-0.5*((yh - y0)/wh).^2)*exp(-0.5*((xh - x0)/wh).^2)  ;
    hp = h ;

elseif prob == 4 ;  %  zonal jet (see notes Barotropic.pdf; not exactly sure where they are from).

    lR = sqrt(g*H)/f ;  % deformation radius
    T0 = 1/f ;          % inverse Coriolis parameter
    w = 0.5 ;           % nondimensional width
    lc = 1/(1 - exp(-2*w)) ;
    kc = sqrt(lc^2 - 1) ;
    U0 = 1 ;
    omi = kc*U0*exp(-2*lc*w) ;
    disp(' ') ;
    disp('Problem 4') ;
    disp(' ') ;
    disp(['Rossby deformation radius:          ' num2str(lR) ' m']) ;
    disp(['1/f:                                ' num2str(1/f) ' s']) ;
    w = w*lR ;
    disp(['Dimensional jet width:              ' num2str(w) 'm']) ;
    U0 = U0*lR/T0 ;
    disp(['Dimensional U0:                     ' num2str(U0) ' m/s']) ;
    disp(['Dimensional cross-wave scale:       ' num2str((2*pi/lc)*lR) ' m']) ;
    disp(['Dimensional instability wavelength: ' num2str((2*pi/kc)*lR) ' m']) ;
    disp(['Dimensional growth time scale:      ' num2str((2*pi/omi)*T0) ' s']) ;

    NL = 1 ;

    u = NaN*zeros(Mu,Nu) ;
    nn = find(abs(yu-y0)>w) ;
    u(nn,:) = (-U0*sign((yu(nn)-y0)).*exp(-(abs((yu(nn)-y0)) - w)/lR))*ones(size(xu)) ;
    nn = find(abs(yu-y0)<=w) ;
    u(nn,:) = (-U0*sinh((yu(nn)-y0)/lR)/sinh(w/lR))*ones(size(xu)) ;
%    u = u + 0.01*randn(Mu,Nu) ;
    up = u ;

    v = zeros(Mv,Nv) ;
    vp = v ;

    h = NaN*zeros(Mh,Nh) ;
    nn = find(abs(yh-y0)>w) ;
    h(nn,:) = (-U0*exp(-(abs((yh(nn)-y0)) - w)/lR))*ones(size(xh)) ;
    nn = find(abs(yh-y0)<=w) ;
    h(nn,:) = (U0*((cosh((yh(nn)-y0)/lR) - cosh(w/lR))/sinh(w/lR) - 1))*ones(size(xu)) ;
    hp = h ;

elseif prob == 5 ;  %  linear geostrophic adjustment zonal flow
    NL = 0 ;
    U = 0 ;
    V = 0 ;
    u = zeros(Mu,Nu) ;
    up = u ;
    v = zeros(Mv,Nv) ;
    vp = v ;
    h0 = 4 ;
    wh = 100 ;
    x0 = Lx/2 ;
    y0 = Ly/2 ;
    h = h0*exp(-0.5*((yh - y0)/wh).^2)*ones(size(xh)) ;  ;
    hp = h ;

elseif prob == 6 ;  %  nonlinear geostrophic adjustment zonal flow
    NL = 1 ;
    U = 0 ;
    V = 0 ;
    u = zeros(Mu,Nu) ;
    up = u ;
    v = zeros(Mv,Nv) ;
    vp = v ;
    h0 = 4 ;
    wh = 100 ;
    x0 = Lx/2 ;
    y0 = Ly/2 ;
    h = h0*exp(-0.5*((yh - y0)/wh).^2)*ones(size(xh)) ;  ;
    hp = h ;

end

figure(prob)
clf
plot_size(1,2,16,10) ;


tm = 0 ;
nplt = round(1/dt) ;
if nplt < 1
    nplt = 1 ;
end
for nt = 1:Nt ;

    %  This linear set up allows for infinitesimal bathymetry and advection
    %  by a spatially constant mean flow.
    if NL == 0
        vv = (v(jpv,:) + v)/2 ;
        vv = (vv + vv(:,imv))/2 ;
        uRHS =    -U*(u(:,ipu)-u(:,imu))/(2*dx) ...
               -   V*(u(jpu,:) - u(jmu,:))/(2*dy) ...
               +   f*(V + vv) ...
               -   g*(h - h(:,imh))/dx  ...
               +   K*(u(:,ipu) - 2*u + u(:,imu))/(dx^2) ...
               +   K*(u(jpu,:) - 2*u + u(jmu,:))/(dy^2) ;

        uu = (u(:,ipu) + u)/2 ;
        uu = (uu + uu(jmu,:))/2 ;
        vRHS =    -U*(v(:,ipv)-v(:,imv))/(2*dx) ...
               -   V*(v(jpv,:) - v(jmv,:))/(2*dy) ...
               -   f*(U + uu) ...
               -   g*(h - h(jmh,:))/dx ...
               +   K*(v(:,ipv) - 2*v + v(:,imv))/(dx^2) ...
               +   K*(v(jpv,:) - 2*v + v(jmv,:))/(dy^2) ;

        hRHS =    -U*((h(:,iph) - hsh(:,iph)) - (h(:,imh) - hsh(:,imh)))/(2*dx) ...
               -   V*((h(jph,:) - hsh(jph,:)) - (h(jmh,:) - hsh(jmh,:)))/(2*dy) ...
               -   H*((u(:,ipu) - u)/dx + (v(jpv,:) - v)/dy) ...
               +   K*(h(:,iph) - 2*h + h(:,imh))/(dx^2) ...
               +   K*(h(jph,:) - 2*h + h(jmh,:))/(dy^2) ;

    %  This is fully nonlinear
    else
        %  get v values on the u grid.
        vv = (v(jpv,:) + v)/2 ;
        vv = (vv + vv(:,imv))/2 ;
        uRHS =    -u.*(u(:,ipu)-u(:,imu))/(2*dx) ...
               -   vv.*(u(jpu,:) - u(jmu,:))/(2*dy) ...
               +   f*vv ...
               -   g*(h - h(:,imh))/dx  ...
               +   K*(u(:,ipu) - 2*u + u(:,imu))/(dx^2) ...
               +   K*(u(jpu,:) - 2*u + u(jmu,:))/(dy^2) ;

        %  get u values on the v grid
        uu = (u(:,ipu) + u)/2 ;
        uu = (uu + uu(jmu,:))/2 ;
        vRHS =    -uu.*(v(:,ipv)-v(:,imv))/(2*dx) ...
               -   v.*(v(jpv,:) - v(jmv,:))/(2*dy) ...
               -   f*uu ...
               -   g*(h - h(jmh,:))/dx ...
               +   K*(v(:,ipv) - 2*v + v(:,imv))/(dx^2) ...
               +   K*(v(jpv,:) - 2*v + v(jmv,:))/(dy^2) ;

        %  get u and v values on the h grid.
        uu = (u(:,ipu) + u)/2 ;
        vv = (v(jpv,:) + v)/2 ;
        argu = (h + H - hsh).*uu ;
        argv = (h + H - hsh).*vv ;
        hRHS =    -(argu(:,iph) - argu(:,imh))/(2*dx) ...
               -   (argv(jph,:) - argv(jmh,:))/(2*dy) ...
               +   K*(h(:,iph) - 2*h + h(:,imh))/(dx^2) ...
               +   K*(h(jph,:) - 2*h + h(jmh,:))/(dy^2) ;
    end

    tm = tm + dt ;
    up = u + dt*uRHS ;
    vp = v + dt*vRHS ;
    hp = h + dt*hRHS ;

    if mod(nt,nplt)==0
        figure(prob)
        clf
        tiledlayout(1,4) ;
        nexttile([1 3]) ;
        pcolor(xh,yh,h)
        axis('equal')
        caxis([-.1 .1]) ;
        colormap(cmap1)
        shading flat
        axis([0 Lx 0 Ly])
        title(['t = ' num2str(tm) ' s']) ;
        nexttile
        plot(h(:,round(Mh/2)),yh,'k','linewidth',2)
        hold on
        plot([0 0],[0 Ly],'k--')
        axis([-1 1 0 Ly])
        pause(0.001) ;
    end

    u = up ;
    v = vp ;
    h = hp ;

end


return
