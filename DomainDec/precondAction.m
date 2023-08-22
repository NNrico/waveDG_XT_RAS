% action of the preconditioner on x (the residual)
% preconditioner: sum(Rtilde_i * A_i * R_i')

function [z]=precondAction(x,Alocal,Rlocal,Rtlocal,nsub)  
    z=0;
    for k=1:nsub        
        xk=Rlocal{k}*x;    
        uk=Alocal{k}\xk;
        zk=transpose(Rtlocal{k})*uk;
        z=z+zk;
    end
    
end


%issymmetric(A)
% isposdef = all(eig(A) > 0)