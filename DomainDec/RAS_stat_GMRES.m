% RAS method used as a PRECONDITIONER. 
% Then it will be used by GMRES method

function [uw1] = RAS_stat_GMRES(uw0,b,nsub_x,nsub_t,A,Alocal,Rlocal,Rtlocal,iter_max,tol)
nsub=nsub_x*nsub_t;
incr=tol+1;
niter=0;
    while(max(incr) > tol && niter<iter_max)      
        uw1= uw0+ precondAction(b-A*uw0,Alocal,Rlocal,Rtlocal,nsub);
        incr=uw1-uw0; %it is not rigorous to check the increment!
                      % anyway it is not a problem since the loop run only
                      % for one iteration (iter_max=1)
        uw0=uw1;
        niter=niter+1;
    end
end