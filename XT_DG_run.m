% DG discretization
% IP and IPh formulation

function [A,b,femregion,Solutions, ERROR, INFO,Data] = XT_DG_run(Data,alpha,formulation,plot_sol)

Data.fem = Data.Degree;
Data.alp=alpha;
% set up quadrature nodes
Data.nqn = 2*Data.fem + 1;


%% load meshfile
load(Data.meshfile);
disp('Uploading mesh...');
disp('------------------------------------------------------')


Data.domain = Dati.domain;
region.id = E_tag;

%% checking tags 
% keep tag A that comes from older implentation

for i = 1 : length(E_tag)   
    
    for j = 1 : length(Data.tag_ac)
        if E_tag(i) == Data.tag_ac(j)
            region.tag(i,1) = 'A';
        end
    end

end


%% computing neighbourung elements
[neighbour] = neighbours_new(region,Data);


%% plot mesh
plot_poly_mesh(Data,region,neighbour)
hold on;
for i = 1 : size(region.coord,1)
    text(region.coord(i,1),region.coord(i,2), num2str(i));
end


%%  creation of the finite element space
[femregion] = create_dof_new(Data,region);

%%  MATRIX ASSEMBLY

disp('Matrices computation ... ');

if formulation==1
    [MatricesAc]  = matrix2D_iph(femregion,neighbour,Data);
    disp('IPH formulation')
elseif formulation==0
    [MatricesAc]  = matrix2D_ip(femregion,neighbour,Data);
    disp('IP formulation')
else
    disp('formulation error')
end


disp('Done')
disp('------------------------------------------------------')


%% RIGHT HAND SIDE

disp('Computing rhs ... ');
if formulation==1
    [fx_a] = evaluate_f_acu_iph(neighbour,femregion,Data);
elseif formulation==0
    [fx_a] = evaluate_f_acu_ip(neighbour,femregion,Data);
else
    disp('formulation error')
end
disp('Done')
disp('------------------------------------------------------')


%% Solve
uh_wh=MatricesAc.block\fx_a;
m=size(uh_wh)/2;
wh=uh_wh((m+1):end);
uh=uh_wh(1:m);
Solutions.phi_h=uh;
Solutions.dot_phi_h=wh;


%% scatter plot of the solution 
if plot_sol==1

    [GUp,GWp,Gphi,GUe] = plot_solution_dis(Data,femregion,Solutions,0.1); %t=0.1
    scatter_plot(GUp,GWp,Gphi,GUe,0,Data); 
end


%% Compute Errors
% missing part not studied in this project

ERROR=-1;
INFO=-1;


%% return
A=MatricesAc.block;
b=fx_a;

end
