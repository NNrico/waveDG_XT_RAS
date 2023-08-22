% #######################################################################

% CREATE MATRICES FOR APSC PROJECT COURSE

% #######################################################################

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
mu=1000;
alpha=mu*(Data.X/Data.NX)/(Data.Degree*Data.Degree); %do not touch 
plot_sol=1;  %if 1 the solution is plotted 

Data.meshfile=mesh;
[A,b,femregion,Solutions, ERROR, INFO,bkData] = XT_DG_run(Data,alpha,formulation,plot_sol);
% (if plot_sol = 1). If there is not an exact solution in Data, then the analytical plot doesn't make sense

%%
% save matrices
[i,j,v]=find(A);
A_txt=[i j v];
writematrix(A_txt, "A_1_5.txt");

[i,j,v]=find(b);
b_txt=[i j v];
writematrix(b_txt,"b_1_5.txt");


%%
% save coords for plot
n=Data.NT *Data.NX;
coord=zeros(n,2);
for i=1:length(femregion.coords_element)
    vals=femregion.coords_element{i};
    xc=mean(vals(:,1));
    yc=mean(vals(:,2));
    coord(i,1)=xc;
    coord(i,2)=yc;
end
writematrix(coord, "coord_1_5.txt");