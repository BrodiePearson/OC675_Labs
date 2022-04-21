%  unp1_uncent01
%  Jim Lerczak
%  27 March
%
%  bc = 0:  clamped
%  bc = 1:  periodic

function [unp1] = unp1_uncent01(un,dt,dx,U,bc) ;

    Nx = length(un) ;
    jj = 2:Nx ;
    unp1 = un ;
    unp1(jj) = un(jj) - (dt*U/dx)*(un(jj) - un(jj-1)) ;

    
    if bc == 0      %  clamped
        unp1(1) = un(1) - (dt*U/dx)*(un(1) - 0) ;
    elseif bc == 1  %  periodic 
        unp1(1) = un(1) - (dt*U/dx)*(un(1) - un(Nx)) ;
    end

  
end