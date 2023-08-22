% run domain decomposition methods with DG matrix A
% input:
    % A, b, femregion, Dati: from dg discretization
    % newtests: tests parameters. There are all domain decomposition choices
    % Ainfo: other struct to describe domain and DG matrix properties
    % overlap_nsub: {0,1} if the decomposition is set based on overlap or number of subdomains
    % saveinfo: {0,1} if results must be saved in tests folder


function [test_result, all_tests_info]=run_domainDec_withA(A,b,femregion,Dati,newtests,Ainfo,overlap_nsub,saveinfo)

all_tests_info={};

name=Ainfo.name;

if saveinfo==1
    c=clock;
    dir=strcat(num2str(c(1)),num2str(c(2)),num2str(c(3)),num2str(c(4)),num2str(c(5)),name);
    mkdir('tests',dir);
    pathname=fullfile(pwd,strcat('tests\',num2str(c(1)),num2str(c(2)),num2str(c(3)),num2str(c(4)),num2str(c(5)),name));
end

% SET DG AND DECOMPOSITION DATA -----------------------------
form=Ainfo.form;
if strcmp(form,'IP')
    form=0;
elseif strcmp(form,'IPH') 
    form=1;
end

mu=Ainfo.mu;
prob=Ainfo.prob; 
X=Ainfo.X;
T=Ainfo.T;
NX=Ainfo.NX;
NT=Ainfo.NT;
alpha=mu*(X/NX)/(Dati.Degree*Dati.Degree); 

meths=newtests.meths;
m=newtests.m;
n=newtests.n;
if overlap_nsub==0
    ot=newtests.ot;
    ox=newtests.ox;
else
    nsub_x=newtests.nsub_x;
    nsub_t=newtests.nsub_t;
end

it_max=newtests.it_max;
it_max_pipe=newtests.it_max_pipe;
tol=newtests.tol;
tol_sx=newtests.tol_sx;
it_wait=newtests.it_wait;

info=[X,T,NX,NT,prob,mu,alpha,form,it_max,it_max_pipe,tol,tol_sx,it_wait];
DataDD.theta = 0.5;  % weight to add in restriction and prolungation matrices (Rtilde) for overlapped elements
DataDD.Nx=NX;
DataDD.NT=NT;
DataDD.X=X;
DataDD.T=T;

test_result=zeros(length(m),length(info)+6+7+4+4);

for i=1:length(m)
    DataDD.m =m(i); 
    DataDD.n =n(i); 
    if overlap_nsub==0
        DataDD.ot =ot(i);
        DataDD.ox =ox(i);
        DataDD.nsub_x=1+(NX-DataDD.n)/(DataDD.n-DataDD.ox);
        DataDD.nsub_t=1+(NT-DataDD.m)/(DataDD.m-DataDD.ot);
        if check_prefect_division(DataDD)==0 % check if there is a perfect division 
            return
        end
    else
        DataDD.nsub_x=nsub_x(i);
        DataDD.nsub_t=nsub_t(i);
    end     
    
    DataDD.nln=femregion.nln; 
    
    disp(strcat('Domain decomposition choice ',num2str(i)))
    disp(strcat('nsubX: ',num2str(DataDD.nsub_x),', nsubT: ',num2str(DataDD.nsub_t)))
    disp(strcat('m: ',num2str(DataDD.m),', n: ',num2str(DataDD.n)))
    disp('Computing ...');
    % DECOMPOSE THE DOMAIN -----------------------------------------------    
    [DataDD]=  createDomDec(DataDD); 
    plotdd(DataDD,femregion)   % plot decomposed domain
    % subs are ordered in the spatial direction
    
    if saveinfo==1
        mat=fullfile(pathname, strcat(name(1:end-4),'_subt',num2str(DataDD.nsub_t),'_subx',num2str(DataDD.nsub_x),'.png'));
        saveas(gcf,mat);
    end    
    
    % PROLUNGATION AND RESTRICTION MATRICES -------------------------------
    [Alocal,Rlocal,Rtlocal]=local_matrices(A,DataDD); 
     
    % TEST RUN ------------------------------------------------------------
    [perf,DataDD]=DD_run(A,b,Alocal,Rlocal,Rtlocal,DataDD,meths,it_max,it_max_pipe,it_wait,tol,tol_sx); %%%%%%%%%%%%
    
    % STORE ALL RESULTS --------------------------------------------------
    if overlap_nsub==0
        dd=[DataDD.nsub_x,DataDD.nsub_t,DataDD.m,DataDD.n,DataDD.ot,DataDD.ox];
    else
        p=length(DataDD.list_ot_forw);
        dd=[DataDD.nsub_x,DataDD.nsub_t,DataDD.m,DataDD.n,mean(DataDD.list_ot_forw(1:round(p/2))),mean(DataDD.list_ox_forw(1:round(p/2)))];
    end
    result_ras=[perf.it_ras, perf.time_ras, perf.relres2P_vec_ras(end),perf.relresinfP_vec_ras(end),perf.relresinf_vec_ras(end),perf.solved_dom_ras];        
    result_gmres=[perf.it_gmres,perf.time_gmres,perf.relres2P_vec_gmres(end),perf.solved_dom_gmres];     
    result_pipe=[perf.it_pipe,perf.time_pipe,perf.resinf_vec_pipe(end),perf.solved_dom_pipe];           
    row_data = [info,i,dd,result_ras,result_gmres,result_pipe];
    test_result(i,1:length(row_data))=row_data;
    
    backup.info=info;
    backup.dd=dd;
    backup.perf=perf;
    
    all_tests_info{i}=backup;
    

    if saveinfo==1
        if i<10
            matfile = fullfile(pathname, strcat(name(1:end-4),'_',num2str(0),num2str(i)));
        else
            matfile = fullfile(pathname, strcat(name(1:end-4),'_',num2str(i)));
        end
        save(matfile,'backup') 
    end
    disp('Done')
    disp('------------------------------------------------------')
   
end

if saveinfo==1
    matfile=fullfile(pathname, 'test_results');
    save(matfile,'test_result')
    xlswrite(matfile,test_result)
    disp('result saved in tests folder')
end

% close all

end


