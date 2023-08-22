% #######################################################################

% XT-DG WAVE_1D + DOM_DEC METHODS (RAS, RAS+GMRES, RAS PIPELINED)

% #######################################################################

% -------------------MAIN----------------------------- 
% sections:
    % addpath
    % assembly dg matrix
    % domain decomposition
    % residual evolution dd methods
    
% run one section at a time
% ----------------------------------------------------

%% addpath
currentFolder = pwd;

pathInput          = fullfile(currentFolder,'InputData');
pathAssembly       = fullfile(currentFolder,'Assembly');
pathErrors         = fullfile(currentFolder,'Errors');
pathMeshGeneration = fullfile(currentFolder,'MeshErrorAnalysis');
pathFEspace        = fullfile(currentFolder,'FEspace');
pathPostProcessing = fullfile(currentFolder,'PostProcessing');
pathDomainDec      = fullfile(currentFolder,'DomainDec'); 

addpath(genpath(pathInput));
addpath(genpath(pathAssembly));
addpath(genpath(pathErrors));
addpath(genpath(pathMeshGeneration));
addpath(genpath(pathFEspace));
addpath(genpath(pathPostProcessing));
addpath(genpath(pathDomainDec));

%% discontinuos Galerkin - obtain matrices A,b
% to set:
    % problem name (test)
    % domain [0,X] x [0,T]
    % number of finite elements (NT, NX)
    % damping (damp)
    % meshname
    % formulation 0: IP, 1: IPH 
    % dg stability coeff (mu)

test='Test11';   %from InputData %Test11 is the one of the report
Data=DataTest(test);
Data.Degree=2;

% mesh in MeshErrorAnalysis 
mesh='ProvaMONO_0105_20100_el.mat';  % domain [0,1]x[0,5], 20 elem in space, 100 in time
% mesh='ProvaMONO_0101_20_el.mat';
Data.X=1;
Data.T=5;
Data.NT =100; 
Data.NX = 20;
Data.damp=0; % DAMPING

formulation=1; %0: IP, 1: IPH 
mu=1000; %mettere anche c %1000
alpha=mu*(Data.X/Data.NX)/(Data.Degree*Data.Degree); %do not touch 
plot_sol=1;  %if 1 the solution is plotted 

Data.meshfile=mesh;
[A,b,femregion,Solutions, ERROR, INFO,bkData] = XT_DG_run(Data,alpha,formulation,plot_sol);
% (if plot_sol = 1). If there is not an exact solution in Data, then the analytical plot doesn't make sense

%%
% create name of DG discretized problem
if formulation ==0
    f='IP';
else
    f='IPH';
end
if Data.damp>0
    name=strcat('mesh_0',num2str(Data.X),'0',num2str(Data.T),'_',num2str(Data.NX),num2str(Data.NT),'_',f,'_mu_',num2str(mu),'_damp');
else
    name=strcat('mesh_0',num2str(Data.X),'0',num2str(Data.T),'_',num2str(Data.NX),num2str(Data.NT),'_',f,'_mu_',num2str(mu));
end

%--------------------------------------------------------------------------
% --------possibility to save dg matrices in a separate folder ------------

% %% save A and info from the test above
% pathname=fullfile(pwd,'folder_to_create\');

% matfile=fullfile(pathname, name);
% save(matfile,'A','b','Data','femregion','formulation','mesh','Solutions');

% %% load A (if previously saved)
% load(name)
% % back to previous folder
%--------------------------------------------------------------------------

%% domain decomposition

% -------- SET DOMAIN AND DG MATRIX PROPERTIES
% Ainfo is a struct that stores the data of the problem discretized with DG
Ainfo.name=name; 
Ainfo.form=0;   % to set
Ainfo.mu=400;   % to set
Ainfo.prob=11;  % must be a number, for example in the first section I have 'Test11', so here set 11
Ainfo.X=1;      % to set
Ainfo.T=5;      % to set
Ainfo.NX=20;    % to set
Ainfo.NT=100;   % to set


% -------- SET DOMAIN DECOMPOSITION DATA FOR TEST
% (all to set)

%   methods:
% meths= [1 2 3]              if I want to exploit all methods (1 RAS, 2 GMRES, 3 PIPELINE)
% meths= [1 0 0] or [0 2 3]   If I want to exploit (RAS) or (GMRES and PIPE), and all other combinations

%   decomposition:
% m: number of TIME finite elements of a subdomain
% n: number of SPACE finite elements of a subdomain
% - once I have set m and n, it remains two degrees of fredoom for defining
% the decomposition: the overlap length and the number of subdomains. We
% have to fix one of those and consequently the other one is adjusted. 
%     -fix overlap: overlap_sub = 0 and set overlap in time and space ot,ox;
%     or
%     -fix subdomains: overlap_sub = 1 and nsub_x, nsub_t.            
% ------> IT IS BETTER TO SET THE NUMBER OF SUBDOMAINS (nsub_x, nsub_t) !!!

% - m,n and (ot,ox) or (nsub_x, nsub_t) must be of the same length: length=k with k = how many 
% decomposition choices I want to exploit

newtests.meths=[1 0 0]; %[0 2 0] [0 2 3]
newtests.m=[6];%[6,15];   %[15,6,6];
newtests.n=[12];%[5,8];   %[8,8,4]; 
overlap_nsub=1;    %0 overlap - 1 nsubx,t
if overlap_nsub==0
    newtests.ot=[2,10];
    newtests.ox=[12,10];
    % requires a perfect division if overlap is set. It is better to choose
    % the number of subdomain (overlap_nsub=1)
else
    newtests.nsub_x=[2];%[5,3];%[2,2,5]; 
    newtests.nsub_t=[20];%[20,7];%[8,20,20];  %doesn't work for only one subt
end

% m, n and number of subs must be correctly chosen
[ok] = check_dd_choice(newtests,overlap_nsub,Ainfo.NX,Ainfo.NT);
if ok==0
    return
end

newtests.it_max=10; % max iterations for RAS and GMRES
newtests.tol=1e-10; 
% PIPELINE info
newtests.it_max_pipe=300; % max iterations
newtests.tol_sx=1e-10; % tolerance to be reached for move left edge of moving window S
newtests.it_wait=3; % how many iterations I have to wait until move the right edge of moving window S (aggressive version)

% final call for domain decomposition test
saveinfo=0; %{0,1} if 1 it creates a folder 'tests' and inside another folder with all test results
[tests_result,all_tests_info] = run_domainDec_withA(A,b,femregion,Data,newtests,Ainfo,overlap_nsub,saveinfo);



% -------RESULTS STRUCTURES DESCRIPTION------------------------------------
% tests_result=[info,i,dd,result_ras,result_gmres,result_pipe] (is a vector)
%     info: [X,T,NX,NT,prob,mu,alpha,form,it_max,it_max_pipe,tol,tol_sx,it_wait]
%     i: index of test with i=1,..,k, k=how many decomposition choice
%     dd: [nsub_x,nsub_t,m,n,ot,ox] where ot, ox are an average if there isn't a perfect decomposition
%     result_ras: [it_ras, time_ras, relres2P_vec_ras(end),relresinfP_vec_ras(end),relresinf_vec_ras(end),solved_dom_ras]        
%     result_gmres: [it_gmres,time_gmres,relres2P_vec_gmres(end),solved_dom_gmres]   
%     result_pipe: [it_pipe,time_pipe,resinf_vec_pipe(end),solved_dom_pipe]
%         relres2P stands for relative residual L2 norm, P stands for preconditioned residual

% all_tests_info is a struct that contains k backups
% backup
%     info: the same as above,
%     dd: the same as above,
%     perf: performance of the methods. contains the same result as above but there are all 
%         the residual vectors.
%
% ------------------------------------------------------------------------



%% dd postprocessing - residual evolution dd methods

%load performance from all_tests_info
test_to_analyze=1; %from 1 to length(m)
backup=all_tests_info{test_to_analyze};

pipe=1; %0 for RAS, RAS+GMRES. 1 for RAS PIPELINE
plot_residual_evolution(backup,pipe,femregion) %femregion needed. If you have saved the dg matrices 
                                               %there is also the femregion

% For the plot of RAS PIPELINE it plots the residual in the domain
% to see the moving subdomain window but only for few iterations. Those
% iterations can be modified inside plot_residual_evolution (not general..) 


%

%% plot solution ----------------------------------------------
% load meshinfo
n=20*100; % nx nt
coord=zeros(n,2);
for i=1:length(femregion.coords_element)
    vals=femregion.coords_element{i};
    xc=mean(vals(:,1));
    yc=mean(vals(:,2));
    coord(i,1)=xc;
    coord(i,2)=yc;
end
% u=vec(1:end/2);
u=Solutions.phi_h;
m=reshape(full(u),6,n);
mean_u=sum(m)/6;

figure('Position', [200 200 700 250])
scatter(coord(:,2),coord(:,1),50,mean_u','filled')
title('displacement u')
xlabel('time')
ylabel('space')
colorbar()

%%
% Solutions.phi_h=vec(1:end/2);
% Solutions.dot_phi_h=vec(end/2+1:end);

[GUp,GWp,Gphi,GUe] = plot_solution_dis(bkData,femregion,Solutions,0.1); %t=0.1
scatter_plot(GUp,GWp,Gphi,GUe,0,bkData); 