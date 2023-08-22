% compute prolungations and restriction matrices Rlocal
% Rtlocal: (Rtilde) there are weights in the overlapped elements

function [Alocal,Rlocal,Rtlocal]=local_matrices(A,DataDD)
    nsub=DataDD.nsub_t*DataDD.nsub_x;
    Z=zeros(DataDD.m*DataDD.n*DataDD.nln, DataDD.NT*DataDD.Nx*DataDD.nln);
   
    Alocal=cell(1,nsub);
    Rlocal=cell(1,nsub);
    Rtlocal=cell(1,nsub);
    
    % subdomain k are ordered in the spatial direction
    for k=1:nsub
        [Rk,Rkt] =  createRk(DataDD,k);
        Rk=[Rk Z;Z Rk];
        Rkt=[Rkt Z;Z Rkt];            
        Ak=Rk*A*transpose(Rk);
        Alocal{k}=Ak;
        Rlocal{k}=Rk;
        Rtlocal{k}=Rkt;
    end



end