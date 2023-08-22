function [E_L2, E_SEMI_H1] = C_error_L2_H1(femregion, uh, Dati)
%% [E_L2, E_SEMI_H1]=C_error_L2_H1(femregion, uh, Dati)
%==========================================================================
% Compute L2 and semi-H1 errors
%==========================================================================
%    called in C_compute_errors.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(nfod,1))  solution vector
%
%    OUTPUT:
%          E_L2        : (real) L2 error 
%          E_SEMI_H1   : (real) H1 error - seminorm


nln = femregion.nln;
ne = femregion.ne;



% shape functions
%[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
%[nodes_1D, w_1D] = C_quadrature(Dati);
[nodes_1D,w_1D] = xwlgl(nln);


% evaluation of shape bases 
%[dphiq,Grad] = C_evalshape(basis,nodes_1D);
%[Grad_lgl_std] = derlgl(nodes_1D,nln);
[dphiq,Grad] = basis_and_der_at_lgl(nodes_1D,nln);

E_SEMI_H1_LOC = zeros(ne,1);
E_L2_LOC = zeros(ne,1);
%%%%%%%%%%%%%%%%% CICLO SUGLI ELEMENTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for ie=1:ne
         
    i = femregion.connectivity(1:nln,ie);   
    
    
    % Determinant of the Jacobian and physical coordinates 
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(i,:), nodes_1D);
   
        
    local_uh=uh(i);
    
    % Solution and gradient at physical nodes
	x = pphys_1D;
	t = Dati.T;
    local_exact = eval(Dati.exact_sol)';
    local_grad_exact = eval(Dati.grad_exact)';
       
    % Approximated solution and gradient at quadrature nodes
    local_grad_aprox = zeros(length(w_1D),1);
    local_aprox = zeros(1,length(w_1D));
    for k = 1:length(w_1D)
        for s=1:nln
             local_aprox(k)=local_aprox(k) + dphiq(1,k,s).*local_uh(s);
             local_grad_aprox(k,:) = local_grad_aprox(k,:) + Grad(k,:,s).*local_uh(s);
        end
    end
    
    % Evaluation integrals for errors
    for k=1:length(w_1D)
        Jdet = BJ;               % determinant
        dx = abs(Jdet).*w_1D(k); % wheight
        Binv = 1./BJ;            % inverse
 
        pointwise_diff(k,:)=(local_grad_exact(k,:))-(local_grad_aprox(k,:)*Binv);
          
            
        E_SEMI_H1_LOC(ie)= E_SEMI_H1_LOC(ie) + (pointwise_diff(k,:) * transpose(pointwise_diff(k,:))).*dx;
        E_L2_LOC(ie)=E_L2_LOC(ie) + ((local_aprox(k)-local_exact(k)).^2).*dx;
    end
end

% OUTPUT - h1 and l2 error
E_SEMI_H1=sqrt(sum(E_SEMI_H1_LOC));
E_L2=sqrt(sum(E_L2_LOC));