function [A, B, C] = assembling_matrices(Data,Mat_Poro,Mat_Ac,Mat_El,...
                                Mat_PoroAc,Mat_PoroEl,Mat_ElAc)  
                            
                            
[~,np1] = size(Mat_Poro.M_P_rho);
[~,np2] = size(Mat_Poro.M_P_rhow);
[~,nac] = size(Mat_Ac.M_A);
[~,nel] = size(Mat_El.M_P_rho);
 
[r,s] = size(Mat_PoroAc.C1_P);
[~,m] = size(Mat_ElAc.C1_A);



A = [Mat_Poro.M_P_rho,   Mat_Poro.M_P_rhof,  sparse(np1,nac),  sparse(np1,nel);
     Mat_Poro.M_P_rhof,  Mat_Poro.M_P_rhow,  sparse(np2,nac),  sparse(np2,nel);
     sparse(nac,np1),    sparse(nac,np2),    Mat_Ac.M_A,       sparse(nac,nel);
     sparse(nel,np1),    sparse(nel,np2),    sparse(nel,nac),  Mat_El.M_P_rho;];
 
 
% C = [Mat_Poro.A_E + Mat_Poro.A_P_beta2, Mat_Poro.A_P_beta, sparse(np1,nac),  Mat_PoroEl.C1_P;    
%      Mat_Poro.A_P_beta,                 Mat_Poro.A_P,      sparse(np2,nac),  sparse(np2,nel);
%      sparse(nac,np1),                   sparse(nac,np2),   Mat_Ac.A_A,       sparse(nac,nel);
%      Mat_PoroEl.C1_E,                   sparse(nel,np2),   sparse(nel,nac),  Mat_El.A_E + Mat_El.Ddis];
%+Mat_Poro.BT_beta_E

C = [Mat_Poro.A_E + Mat_Poro.A_P_beta2, Mat_Poro.A_P_beta,    sparse(np1,nac),  Mat_PoroEl.C_P;    
     Mat_Poro.A_P_beta,                 Mat_Poro.A_P,         sparse(np2,nac),  sparse(np2,nel); 
     sparse(nac,np1),                   sparse(nac,np2),      Mat_Ac.A_A,       sparse(nac,nel);
     Mat_PoroEl.C_E,                    sparse(nel,np2),      sparse(nel,nac),  Mat_El.A_E + Mat_El.Ddis+Mat_El.Rdis];

  % [~,q] = size(Mat_ElAc.C1_A);
 
if (Data.tau ~=0)
    
  D = Mat_Poro.M_P_eta_kper + (1-Data.tau)/Data.tau*Mat_Poro.D_P;
  
%   B = [sparse(r,r),        sparse(r,r),          Mat_PoroAc.C1_P,  sparse(r,m);
%       sparse(r,r),         D,                    Mat_PoroAc.C2_P;  sparse(r,m);
%       Mat_PoroAc.C1_A,     Mat_PoroAc.C2_A,      sparse(s,s),      Mat_ElAc.C1_A;
%       sparse(m,r),         sparse(m,r)           Mat_ElAc.C1_E     Mat_El.Dvel];

  B = [Mat_Poro.ABC_UU,     Mat_Poro.ABC_UW,      Mat_PoroAc.C1_P,  sparse(r,m);
       Mat_Poro.ABC_WU,     Mat_Poro.ABC_WW + D,  Mat_PoroAc.C2_P;  sparse(r,m);
       Mat_PoroAc.C1_A,     Mat_PoroAc.C2_A,      sparse(s,s),      Mat_ElAc.C1_A;
       sparse(m,r),         sparse(m,r)           Mat_ElAc.C1_E     Mat_El.Dvel+Mat_El.Svel];
else
    
  D = Mat_Poro.M_P_eta_kper;
    
%   B = [sparse(r,r),       sparse(r,r),  Mat_PoroAc.C1_P,  sparse(r,m);
%        sparse(r,r),       D,            sparse(r,s),      sparse(r,m);
%        Mat_PoroAc.C1_A,   sparse(s,r),  sparse(s,s),      Mat_ElAc.C1_A;
%        sparse(m,r),       sparse(m,r)   Mat_ElAc.C1_E     Mat_El.Dvel];
       
  B = [Mat_Poro.ABC_UU,   Mat_Poro.ABC_UW,       Mat_PoroAc.C1_P,  sparse(r,m);
       Mat_Poro.ABC_WU,   Mat_Poro.ABC_WW + D,   sparse(r,s),      sparse(r,m);
       Mat_PoroAc.C1_A,   sparse(s,r),           sparse(s,s),      Mat_ElAc.C1_A;
       sparse(m,r),       sparse(m,r)            Mat_ElAc.C1_E     Mat_El.Dvel+Mat_El.Svel];
        
 end





% if (Data.problem_type == 1)
%     
%     
%     
%     C = [Mat_.A_E + Mat_.A_P_beta2, Mat_.A_P_beta, sparse(r,s);
%         Mat_.A_P_beta,                  Mat_.A_P,      sparse(r,s);
%         sparse(s,r),                        sparse(s,r),       Mat_.A_A];
%     
% else
%     
%     A = [Mat_.M_P_rho,   Mat_.M_P_rhof
%         Mat_.M_P_rhof,  Mat_.M_P_rhow ];
%     
%     B = [sparse(r,r),  sparse(r,r),         ;
%         sparse(r,r),  Mat_.M_P_eta_kper];
%     
%     C = [Mat_.A_E + Mat_.A_P_beta2,  Mat_.A_P_beta;
%         Mat_.A_P_beta,                  Mat_.A_P];
%     
% end