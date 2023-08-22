function [Rk,Rkt] =  createRk(DataDD,k)

    % restriction and prolungation matrices for subdomain k. 
    % in Rkt (Rtilde) there are weights (DataDD.theta) in the overlapped elements
    
    m = DataDD.m;
    n = DataDD.n;
    NT = DataDD.NT;
    Nx = DataDD.Nx;
    nln = DataDD.nln;
    theta = DataDD.theta;
    elem = DataDD.list_sub(k);

    ot_forw = DataDD.list_ot_forw(k);
    ox_forw = DataDD.list_ox_forw(k);
    ot_back = DataDD.list_ot_back(k);
    ox_back = DataDD.list_ox_back(k);

    indexrow = 1:nln*m*n;
    indexcol = zeros(1,nln*m*n);
    value = zeros(1,nln*m*n);
    for k =0:n-1
        indexcol(k*nln*m+1:nln*m*(k+1)) = k*NT*nln+1+(elem-1)*nln:(k*NT+m)*nln+(elem-1)*nln; 
    end

    for k =0:ox_back-1
    value(k*nln*m+1:nln*m*(k+1)) = [ones(1,ot_back*nln)*(1-theta)/2 ones(1,(m-ot_forw-ot_back)*nln)/2 ones(1,ot_forw*nln)*theta/2];
    end
    for k=ox_back:n-1-ox_forw
       value(k*nln*m+1:nln*m*(k+1)) = [ones(1,ot_back*nln)*(1-theta) ones(1,(m-ot_forw-ot_back)*nln) ones(1,ot_forw*nln)*theta];
    end 
    for k =n-ox_forw:n-1
    value(k*nln*m+1:nln*m*(k+1)) = [ones(1,ot_back*nln)*(1-theta)/2 ones(1,(m-ot_forw-ot_back)*nln)/2 ones(1,ot_forw*nln)*theta/2];
    end
        
    Rkt = sparse(indexrow,indexcol,value,nln*m*n,nln*NT*Nx);
    Rk =sparse(indexrow,indexcol,ones(1,nln*m*n),nln*m*n,nln*NT*Nx);

end