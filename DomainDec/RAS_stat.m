% RAS method used as a SOLVER. 
% Stationary Richardson method with the RAS preconditioner

function [uw1,relres2P_vec,relresinfP_vec,relresinf_vec,niter,solves,solves_vec,mat_res_pipe] = RAS_stat(uw0,b,nsub_x,nsub_t,A,Alocal,Rlocal,Rtlocal,iter_max,tol,tol_pipe,it_wait,pipe)
% incr=tol+1;
solves_vec=0;
mat_res_pipe=[];   % all residuals, it will be used for postprocessing
nsub=nsub_x*nsub_t;
res=tol+1;
niter=0;
solves=0;
relres2P_vec=[];    % relative residual preconditioned, L2 
relresinfP_vec=[];  % relative residual preconditioned, Linf
relresinf_vec=[];   % relative residual, Linf

bb=precondAction(b,Alocal,Rlocal,Rtlocal,nsub); 
Pb2=norm(bb,2);
Pbinf=norm(bb,inf);
binf=norm(b,inf);

if pipe== 1
    relres2P_vec=0;
    relresinfP_vec=0;
    solves_vec = 0;
    sx=1;
    zone=2;
    it_waited=0;
    uw1=uw0;
    while( res > tol && niter<iter_max)  
        v=b-A*uw1;        
        mat_res_pipe=[mat_res_pipe,v];
        [z,sx,zone,it_waited,solves]=precondAction_pipe(v,Alocal,Rlocal,Rtlocal,nsub_x,nsub_t,sx,tol_pipe,it_wait,it_waited,zone,solves);
        res=norm(v,inf);
        niter=niter+1;
        relresinf_vec(niter)=res;
        solves_vec(niter)=solves;
        uw0=uw1;
        uw1= uw0+ z;
               
    end    
else    
    z=precondAction(b-A*uw0,Alocal,Rlocal,Rtlocal,nsub);
    while( res > tol && niter<iter_max)                 
        uw1= uw0+ z;
        v=b-A*uw1;
        z=precondAction(v,Alocal,Rlocal,Rtlocal,nsub);
        res=norm(z/Pb2,2);
        niter=niter+1;
        relres2P_vec(niter)=res;
        relresinfP_vec(niter)=norm(z/Pbinf,inf);
        relresinf_vec(niter)=norm(v/binf,inf);
        uw0=uw1;
        
    end
    solves=niter*nsub;

   
end
end