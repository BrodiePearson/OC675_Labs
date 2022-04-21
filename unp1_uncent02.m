%  unp1_uncent02
%  Jim Lerczak
%  27 March
%
%  bc = 0:  clamped
%  bc = 1:  periodic

function [unp1] = unp1_uncent02(un,dt,dx,U,bc) ;

    Nx = length(un) ;
    jj = 2:Nx-1 ;
    unp1 = un ;
    unp1(jj) = un(jj) - (dt*U/(2*dx))*(un(jj+1) - un(jj-1)) ;

    
    if bc == 0      %  clamped
        unp1(1) = un(1) - (dt*U/(2*dx))*(un(2) - 0) ;
        unp1(Nx) = un(Nx) - (dt*U/(2*dx))*(0 - un(Nx-1)) ;
    elseif bc == 1  %  periodic 
        unp1(1) = un(1) - (dt*U/(2*dx))*(un(2) - un(Nx)) ;
        unp1(Nx) = un(Nx) - (dt*U/(2*dx))*(un(1) - un(Nx-1)) ;
    end

  
end