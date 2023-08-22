% run domain decomposition methods: RAS, RAS+GMRES, RAS PIPELINED 
% return perf: it is a struct with all residuals, iterations, time, solved
% subdomains of the methods. All the performances

function [perf,DataDD]=DD_run(A,b,Alocal,Rlocal,Rtlocal,DataDD,meths,it_max,it_max_pipe,it_wait,tol,tol_sx)                        

%initialize struct
perf.it_ras=0;
perf.time_ras=0;
perf.relres2P_vec_ras=0;
perf.relresinfP_vec_ras=0;
perf.relresinf_vec_ras=0;
perf.solved_dom_ras=0;  % solved subdomains

perf.it_gmres=0;
perf.time_gmres=0;
perf.relres2P_vec_gmres=0;
perf.solved_dom_gmres=0;

perf.it_pipe=0;
perf.time_pipe=0;
perf.resinf_vec_pipe=0;
perf.solved_dom_pipe=0;

if any(meths==1) %RAS
    uw0=zeros(2*DataDD.NT*DataDD.Nx*DataDD.nln,1);     
    tic  
    % action of preconditioner on residual
    [~,relres2P_vec,relresinfP_vec,relresinf_vec,niter,solves,~,~] = RAS_stat(uw0,b,DataDD.nsub_x,DataDD.nsub_t,A,Alocal,Rlocal,Rtlocal,it_max,tol,0,0,0);
    perf.time_ras=toc;
    perf.it_ras=niter;
    perf.relres2P_vec_ras=relres2P_vec;
    perf.relresinfP_vec_ras=relresinfP_vec;
    perf.relresinf_vec_ras=relresinf_vec;
    perf.solved_dom_ras=solves;
   
end 
    
if any(meths==2)    %GMRES 
    Afun=@(v) v-RAS_stat_GMRES(v,b*0,DataDD.nsub_x,DataDD.nsub_t,A,Alocal,Rlocal,Rtlocal,1,tol);    
    bb=RAS_stat_GMRES(0*b,b,DataDD.nsub_x,DataDD.nsub_t,A,Alocal,Rlocal,Rtlocal,1,tol);    
    tic
    [~,~,~,it,resvec]=gmres(Afun,bb,[],tol,it_max); 
    perf.time_gmres=toc;
    perf.it_gmres=it(2);     
    perf.relres2P_vec_gmres=resvec/norm(bb,2); 
    perf.solved_dom_gmres=DataDD.nsub_x*DataDD.nsub_t*perf.it_gmres;

end   
if any(meths==3)  %RAS PIPELINE      
    uw0=15*ones(2*DataDD.NT*DataDD.Nx*DataDD.nln,1);    
    tic         
    % action of preconditioner on residual
    [~,~,~,resinf_vec,niter,solves,solves_vec,mat_res] = RAS_stat(uw0,b,DataDD.nsub_x,DataDD.nsub_t,A,Alocal,Rlocal,Rtlocal,it_max_pipe,tol,tol_sx,it_wait,1);
    perf.time_pipe=toc; 
    perf.it_pipe=niter;
    perf.resinf_vec_pipe=resinf_vec;   
    perf.solved_dom_pipe=solves;
    perf.solves_vec = solves_vec;
    perf.mat_res=mat_res;  % all residuals, for postprocessing
   
    
end
end











