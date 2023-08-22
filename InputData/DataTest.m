function [Data]=DataTest(test)
Data.name = 'Convergence_Ac';


Data.tag_ac      = [1];
Data.tag_ac_bc   = [2 3 4 5]; % Numbers
Data.lab_ac_bc   = ['DDDD']; % do not modify, it was useful for older implementation


if strcmp(test,'Test1')
    Data.rho_a = [1];
    Data.c     = [1];
    Data.b=[0;1];
    Data.source_phi  = '(pi*pi+1)*exp(-y).*sin(pi*x)';
    Data.source_phit =  '0.*x';
    % exact solution --> used to compute the initial conditions

    Data.exact_phi  =  'exp(-y).*sin(pi*x)'; 
    Data.exact_phit =  '0.*x'; 
    Data.exact_dphit = '0.*x';

    % exact gradient --> used for the error analysis
    Data.exact_dphi_x =  'pi*exp(-y).*cos(pi*x)';
    Data.exact_dphi_y =  '-exp(-y).*sin(pi*x)';
    
elseif strcmp(test,'Test2')   
    Data.rho_a = [1];
    Data.c     = [1];
    Data.b=[0;1];
    Data.source_phi  = '- 3*exp(x.^2).*cos(y) - 4*x.^2.*exp(x.^2).*cos(y)';
    Data.source_phit =  '0.*x';
    %exact solution --> used to compute the initial conditions

    Data.exact_phi  =  'exp(x.^2).*cos(y)';
    Data.exact_phit =  '0.*x'; 
    Data.exact_dphit = '0.*x';

    %exact gradient --> used for the error analysis
    Data.exact_dphi_x =  '2*x.*exp(x.^2).*cos(y)';
    Data.exact_dphi_y =  '-exp(x.^2).*sin(y)';
    
    
elseif strcmp(test,'Test10') 
    Data.rho_a = [1];
    Data.c     = [1];
    Data.b=[0;1];
    Data.source_phi  = 'pi^2*sin(10*y).*cos(pi*x) - 100*sin(10*y).*cos(pi*x) +  10*(10*cos(10*y).*cos(pi*x))'; %damping
    Data.source_phit =  '0.*x';
    % exact solution --> used to compute the initial conditions

    Data.exact_phi  =  'sin(10*y).*cos(pi*x)';
    Data.exact_phit =  '0.*x'; 
    Data.exact_dphit = '0.*x';

    % exact gradient --> used for the error analysis
    Data.exact_dphi_x =  '-pi*sin(10*y).*sin(pi*x)';
    Data.exact_dphi_y =  '10*cos(10*y).*cos(pi*x)';
    Data.exact_dphi_y0 =  '10*cos(pi*x)';
    
elseif strcmp(test,'Test11')
    Data.rho_a = [1];
    Data.c     = [1];
    Data.b=[0;1];
    Data.source_phi  = '0*x.*y'; 
    Data.source_phit =  '0.*x';
    % exact solution --> used to compute the initial conditions

    Data.exact_phi  =  '(x.^2).*(1-x).*(sin(pi*x).^2)';
    Data.exact_phit =  '0.*x'; 
    Data.exact_dphit = '0.*x';

    % exact gradient --> used for the error analysis
    Data.exact_dphi_x =  '0*x.*y';
    Data.exact_dphi_y =  '0*x.*y';
    Data.exact_dphi_y0 =  '0*x.*y';
    
elseif strcmp(test,'Test12') 
    Data.rho_a = [1];
    Data.c     = [1];
    Data.b=[0;1];
    Data.source_phi  = '(pi^2+1)*cos(pi*x).*exp(-y)';
    Data.source_phit =  '0.*x';
    % exact solution --> used to compute the initial conditions

    Data.exact_phi  =  'exp(-y).*cos(pi*x)';
    Data.exact_phi0  =  'cos(pi*x)';
    Data.exact_phit =  '0.*x'; 
    Data.exact_dphit = '0.*x';

    % exact gradient --> used for the error analysis
    Data.exact_dphi_x =  '-pi*exp(-y).*sin(pi*x)';
    Data.exact_dphi_y =  '-exp(-y).*cos(pi*x)';
    Data.exact_dphi_y0 =  '-cos(pi*x)';

end


