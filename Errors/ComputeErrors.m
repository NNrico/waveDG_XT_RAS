function [ERROR,INFO] = ComputeErrors(Solutions,Data,time,...
                                      femregion,neighbour,MatricesPoro,...
                                      MatricesAc,MatricesEl)

[up0,wp0,phi0,ue0] = evaluate_solution_new(Data,femregion);

                                  
%% compute errors new version - matrix form
% compute exact solution U at final time
t = time;
up_ex   = up0   * eval(Data.exact_upt);
wp_ex   = wp0   * eval(Data.exact_wpt);
phi_ex  = phi0  * eval(Data.exact_phit);
ue_ex   = ue0  * eval(Data.exact_uet);


% compute \dot{U} at final time
dot_up_ex   = up0   * eval(Data.exact_dupt);
dot_wp_ex   = wp0   * eval(Data.exact_dwpt);
dot_phi_ex  = phi0  * eval(Data.exact_dphit);
dot_ue_ex   = ue0   * eval(Data.exact_duet);

% switch from nodal to modal representation of the solution coefficients
up_ex    = MatricesPoro.MPrjP\up_ex;
wp_ex    = MatricesPoro.MPrjP\wp_ex;
phi_ex   = MatricesAc.MPrjA\phi_ex;
ue_ex    = MatricesEl.MPrjP\ue_ex;

dot_up_ex   = MatricesPoro.MPrjP\dot_up_ex;
dot_wp_ex   = MatricesPoro.MPrjP\dot_wp_ex;
dot_phi_ex  = MatricesAc.MPrjA\dot_phi_ex;
dot_ue_ex   = MatricesEl.MPrjP\dot_ue_ex;


% Definition of errors - differences between modal coeff.
error_up       = up_ex  - Solutions.up_h;
error_wp       = wp_ex  - Solutions.wp_h;
error_phi      = phi_ex - Solutions.phi_h;
error_ue       = ue_ex  - Solutions.ue_h;
%error_bu_w    = beta * error_u + error_w;

error_dot_up   = dot_up_ex  - Solutions.dot_up_h;
error_dot_wp   = dot_wp_ex  - Solutions.dot_wp_h;
error_dot_phi  = dot_phi_ex - Solutions.dot_phi_h;
error_dot_ue   = dot_ue_ex  - Solutions.dot_ue_h;

%% DG errors
% Calculation of dg errors as in the paper
% Check the definition of the matrices

error_dGep   = error_up'    * MatricesPoro.DGe * error_up;
error_dGp    = error_up'    * MatricesPoro.DGp_beta * error_up + error_wp' * MatricesPoro.DGp * error_wp;
error_dGa    = error_phi'   * MatricesAc.DGa * error_phi;
error_dGp_w  = error_wp'    * MatricesPoro.DGp * error_wp;
error_dGe    = error_ue'    * MatricesEl.DGe * error_ue;

if isempty(error_dGep); error_dGep = 0; end
if isempty(error_dGp);  error_dGp = 0;  end
if isempty(error_dGa);  error_dGa = 0;  end
if isempty(error_dGe);  error_dGe = 0;  end
    
error_dG     = error_dGep + error_dGp + error_dGa + error_dGe;

%% L2 velocity errors
% Calculation of L^2 velocity errors
% Check the definition of the matrices

error_L2_dot_up   = error_dot_up'   * (    MatricesPoro.M_P_rho     ) * error_dot_up;
error_L2_dot_wp   = error_dot_wp'   * (    MatricesPoro.M_P_rhow    ) * error_dot_wp;
error_L2_dot_uwp  = error_dot_wp'   * (2 * MatricesPoro.M_P_rhof    ) * error_dot_up;
error_L2_dot_phi  = error_dot_phi'  * (    MatricesAc.M_A           ) * error_dot_phi;
error_L2_dot_ue   = error_dot_ue'   * (    MatricesEl.M_P_rho       ) * error_dot_ue;

if isempty(error_L2_dot_up);  error_L2_dot_up = 0;  end
if isempty(error_L2_dot_wp);  error_L2_dot_wp = 0;  end
if isempty(error_L2_dot_uwp); error_L2_dot_uwp = 0; end
if isempty(error_L2_dot_phi); error_L2_dot_phi = 0; end
if isempty(error_L2_dot_ue);  error_L2_dot_ue = 0;  end


error_L2_vel     = error_L2_dot_up + error_L2_dot_wp + error_L2_dot_uwp ...
                    + error_L2_dot_phi + error_L2_dot_ue;

%% L2 displacement errors                
                
error_L2_up   = error_up'   * MatricesPoro.MPrjP * error_up;
error_L2_wp   = error_wp'   * MatricesPoro.MPrjP * error_wp;
error_L2_phi  = error_phi'  * MatricesAc.MPrjA   * error_phi;
error_L2_ue   = error_ue'   * MatricesEl.MPrjP   * error_ue;

if isempty(error_L2_up);  error_L2_up = 0;  end
if isempty(error_L2_wp);  error_L2_wp = 0;  end
if isempty(error_L2_phi); error_L2_phi = 0; end
if isempty(error_L2_ue);  error_L2_ue = 0;  end

error_L2     = error_L2_up + error_L2_wp + error_L2_phi + error_L2_ue;

%% interface errors for poroacoustics
error_B_wp = error_wp' * (MatricesPoro.M_P_eta_kper) * error_wp;
if (Data.tau ~= 0 && Data.tau ~= 1)
    error_B_interface = compute_Berror_interface(Data,femregion,neighbour,Solutions.wp_h,t,Data.tau);
else
    error_B_interface = 0;
end
error_B   = error_B_wp + error_B_interface;
if isempty(error_B); error_B = 0; end

%% Total errors 
% Energy error
error_Energy = sqrt(error_L2_vel + error_dG + error_B);
error_L2_v   = sqrt(error_L2_vel);
error_L2_d   = sqrt(error_L2);

% compute poroelastic pressure
u_h1 = Solutions.up_h(1:femregion.ndof_p);
u_h2 = Solutions.up_h(femregion.ndof_p+1:2*femregion.ndof_p);
w_h1 = Solutions.wp_h(1:femregion.ndof_p);
w_h2 = Solutions.wp_h(femregion.ndof_p+1:2*femregion.ndof_p);

if isempty(u_h1)
    p_h = 0;
else
    p_h = MatricesPoro.MPrjP_1\(MatricesPoro.P1_beta*u_h1 + MatricesPoro.P2_beta*u_h2 + MatricesPoro.P1*w_h1 + MatricesPoro.P2*w_h2);
end
% if exact poroelastic pressure = 0
error_L2_p_h = p_h' * MatricesPoro.MPrjP_1 * p_h;
error_L2_pressure = sqrt(error_L2_p_h + error_L2_dot_phi);



%% Savings

INFO.Nel = femregion.ne;
INFO.h   = Data.h;
INFO.tau = Data.tau;
INFO.p   = Data.fem;

ERROR.error_L2_p_h      = error_L2_p_h;
ERROR.error_L2_pressure = error_L2_pressure;
ERROR.error_L2_v        = error_L2_v;
ERROR.error_L2_d        = error_L2_d;
ERROR.error_Energy      = error_Energy;
ERROR.error_B           = error_B;
ERROR.error_B_interface = error_B_interface;
ERROR.error_L2_vel      = error_L2_vel;
ERROR.error_L2_dot_phi  = error_L2_dot_phi;
ERROR.error_L2_dot_uwp  = error_L2_dot_uwp;
ERROR.error_L2_dot_wp   = error_L2_dot_wp;
ERROR.error_L2_dot_up   = error_L2_dot_up;
ERROR.error_L2_dot_ue   = error_L2_dot_ue;
ERROR.error_dG          = error_dG;
ERROR.error_dGp_w       = error_dGp_w;
ERROR.error_dGa         = error_dGa;
ERROR.error_dGp         = error_dGp;
ERROR.error_dGep        = error_dGep;
ERROR.error_dGe         = error_dGe;
ERROR.error_L2_phi      = error_L2_phi;
ERROR.error_L2_wp       = error_L2_wp;
ERROR.error_L2_up       = error_L2_up;
ERROR.error_L2_ue       = error_L2_ue;


