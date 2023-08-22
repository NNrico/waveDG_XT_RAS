
XT_DG for the 1D wave equation, cartesian mesh only
Domain decomposition methods: RAS, GMRES with RAS preconditioner, RAS PIPELINE


from
A space-time discontinuous Galerkin method for the wave equation
Authors: Francesco Songia, Enrico Zardi
mailto: francesco.songia@mail.polimi.it, enrico.zardi@mail.polimi.it 
Advisors: Ilario Mazzieri, Gabriele Ciaramella
Project in Numerical Analysis for Partial Differential Equations
Academic year 2021-22


In this file:
 - code main structure
 - 'tests' folder description
 - result structures descriptions


For the DG discretization we have used as a base structure a given library for solving 
acoustic poroelastic problems with polygonal meshes. We have kept only the acoustic 
structure but then we have modified all the matrices computation since we exploit a 
space time strategy.
Those was provided by our advisors. 
We have created from scratch the domain decomposition library (DomainDec) with all the 
methods. 

-------------------------------------------------------------------------------------
Code main structure

main.m
   DG discretization (XT_DG_run.m)
	
   domain decomposition (run_domainDec_withA.m)
	
   residual evolution in dd methods

-------------------------------------------------------------------------------------
'tests' folder

In the folder 'tests' we stored the results of domain decomposition methods (if 
saveinfo = 1 in main.m). 
In this folder it is created a new folder named with the creation time, the mesh
used (es mesh_0105_10020 stands for a [0,1]x[0,5] domain (space x time) with 100 and 
20 finite elements for space and time respectively), the formulation used {IP,IPH} 
and the DG stability coefficient mu.

k = number of domain decomposition choice exploited in the test
Inside this new folder there are:
 - k .png images with the decompsed domain
 - k .mat files named as the mesh (backup of a single decomposition choice). Here there 
   is the backup struct with all the test information
 - test_result (.mat and .xls) with a compact description of all results

--------------------------------------------------------------------------------------
Result structures description

### test_results
tests_result=[info,i,dd,result_ras,result_gmres,result_pipe] (is a vector)
     info: [X,T,NX,NT,prob,mu,alpha,form,it_max,it_max_pipe,tol,tol_sx,it_wait]
     i: index of test with i=1,..,k, k = how many decomposition choice
     dd: [nsub_x,nsub_t,m,n,ot,ox] 
     result_ras: [it_ras, time_ras, relres2P_vec_ras(end),relresinfP_vec_ras(end),
                  relresinf_vec_ras(end),solved_dom_ras]        
     result_gmres: [it_gmres,time_gmres,relres2P_vec_gmres(end),solved_dom_gmres]   
     result_pipe: [it_pipe,time_pipe,resinf_vec_pipe(end),solved_dom_pipe]

X,T: domain [0,X]x[0,T]
NX, NT: finite elements
prob: number of input data problem (InputData\DataTest.m)
mu, alpha: DG stability coeff, mu = alpha * fem_degree/h, h = X/NX
form: IP (0) or IPH (1) formulation
it_max: max number of iterations for RAS, GMRES
it_max_pipe: max number of iterations for RAS PIPELINE
tol: tolerance for convergence
tol_sx: only for RAS PIPELINE. Tolerance to reach for move the left edge of the 
        subdomains window
it_wait: only for RAS PIPELINE. Wait it_wait iterations before update the right edge  
         of the subdomains window
i: index of test with i=1,..,k, k = how many decomposition choice
nsub_x,nsub_t: number of subdomains
m, n: each subdomain is composed by n*m finite elements (space*time) 
ot, ox: number of finite elements that form the overlap. They are an average if  
        there isn't a perfect decomposition
result_***: contains iterations used, time, final residuals and solved subdomains.
            relres2P stands for relative residual L2 norm, P stands for 
            preconditioned residual.

### backup
     info: the same as above,
     dd: the same as above,
     perf: performance of the methods. It contains the same result as above but  
           there are all the residual vectors (all the history).

--------------------------------------------------------------------------------------

   
	
