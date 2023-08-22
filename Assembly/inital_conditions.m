function [SolutionsOld] = inital_conditions(Data,femregion,time,...
                                   MatricesPoro,MatricesAc,MatricesEl)


[up0,wp0,phi0,ue0] = evaluate_solution_new(Data,femregion);

t = time;

up0t   = up0   * eval(Data.exact_upt);
wp0t   = wp0   * eval(Data.exact_wpt);
phi0t  = phi0  * eval(Data.exact_phit);
ue0t   = ue0   * eval(Data.exact_uet);

dup0t   = up0   * eval(Data.exact_dupt);
dwp0t   = wp0   * eval(Data.exact_dwpt);
dphi0t  = phi0  * eval(Data.exact_dphit);
due0t   = ue0   * eval(Data.exact_duet);

up0t   = MatricesPoro.MPrjP\up0t;
wp0t   = MatricesPoro.MPrjP\wp0t;
phi0t  = MatricesAc.MPrjA  \phi0t;
ue0t   = MatricesEl.MPrjP  \ue0t;

dup0t   = MatricesPoro.MPrjP\dup0t;
dwp0t   = MatricesPoro.MPrjP\dwp0t;
dphi0t  = MatricesAc.MPrjA  \dphi0t;
due0t   = MatricesEl.MPrjP  \due0t;

SolutionsOld.up0   = up0t;
SolutionsOld.wp0   = wp0t;
SolutionsOld.phi0  = phi0t;
SolutionsOld.ue0   = ue0t;
SolutionsOld.dup0  = dup0t;
SolutionsOld.dwp0  = dwp0t;
SolutionsOld.dphi0 = dphi0t;
SolutionsOld.due0  = due0t;
