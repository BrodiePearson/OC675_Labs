%  unp1_cent01
%  Jim Lerczak
%  27 March
%
%  bc = 0:  clamped
%  bc = 1:  periodic

function [unp1] = unp1_cent01(unm1,un,dt,dx,U,bc) ;

    Nx = length(un) ;
    jj = 2:Nx-1 ;
    unp1 = un ;
    unp1(jj) = unm1(jj) - (dt*U/dx)*(un(jj+1) - un(jj-1)) ;

    
    if bc == 0      %  clamped
        unp1(1) = unm1(1) - (dt*U/dx)*(un(2) - 0) ;
        unp1(Nx) = unm1(Nx) - (dt*U/dx)*(0 - un(Nx-1)) ;
    elseif bc == 1  %  periodic 
        unp1(1) = unm1(1) - (dt*U/dx)*(un(2) - un(Nx)) ;
        unp1(Nx) = unm1(Nx) - (dt*U/dx)*(un(1) - un(Nx-1)) ;
    end

  
end