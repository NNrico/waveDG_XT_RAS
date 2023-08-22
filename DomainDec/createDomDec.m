%  creation of the subdomain decomposition

function [DataDD] =  createDomDec(DataDD)    

    m=DataDD.m;
    n = DataDD.n;
    n_sub_x = DataDD.nsub_x;
    n_sub_t = DataDD.nsub_t;
    NT = DataDD.NT;
    Nx = DataDD.Nx;

    % compute average ot (temporal overlap) and average ox (spacial overlap) and also the rest     
    if n_sub_t > 1 
    ot_mean = fix(m- (NT -m)/(n_sub_t-1));
    resto_ot = (m-ot_mean)*(n_sub_t-1)- (NT -m);
    else
        ot_mean = 0;
        resto_ot = 0;
    end

    if n_sub_x > 1 
    ox_mean = fix(n- (Nx -n)/(n_sub_x-1));
    resto_ox = (n-ox_mean)*(n_sub_x-1)- (Nx -n);
    else
        ox_mean = 0;
        resto_ox = 0;
    end

    % the following lists collects the overlaps for all the subdomains
    
    %list_ot_forw
    a = ones(1,resto_ot)*ot_mean+1;
    b = ones(1,n_sub_t-resto_ot-1)*ot_mean;
    list_ot_forw = repmat([a b 0],1,n_sub_x);

    %list_ox_forw
    a = ones(resto_ox,n_sub_t)*ox_mean+1;
    b = ones(n_sub_x-resto_ox-1,n_sub_t)*ox_mean;
    c = zeros(1,n_sub_t);
    list_ox_forw = reshape([a;b;c]',1,n_sub_t*n_sub_x);


    %list_ot_back (from list_ot_forw)
    a = ones(1,resto_ot)*ot_mean+1;
    b = ones(1,n_sub_t-resto_ot-1)*ot_mean;
    list_ot_back = repmat([0 a b],1,n_sub_x);

    %list_ox_back (from list_ox_frow)
    a = ones(resto_ox,n_sub_t)*ox_mean+1;
    b = ones(n_sub_x-resto_ox-1,n_sub_t)*ox_mean;
    c = zeros(1,n_sub_t);
    list_ox_back = reshape([c;a;b]',1,n_sub_t*n_sub_x);

    list_elem = zeros(n_sub_x,n_sub_t);
    m_list_ot_forw = reshape(list_ot_forw,n_sub_t,n_sub_x)';
    m_list_ox_forw = reshape(list_ox_forw,n_sub_t,n_sub_x)';


    for i = 1:n_sub_x
        if i == 1
            list_elem(i,1) = 1;
        else
            j=n_sub_t;
            list_elem(i,1) = list_elem(i-1,1) + (n-m_list_ox_forw(i-1,j))*NT;
        end
        for j = 2:n_sub_t
            list_elem(i,j) = list_elem(i,j-1)+m-m_list_ot_forw(i,j-1);
        end
    end

    list_elem = reshape(list_elem',1,n_sub_t*n_sub_x);
    % list_elem = uint32(list_elem);

    DataDD.list_sub = list_elem;
    DataDD.list_ot_forw = list_ot_forw;
    DataDD.list_ot_back = list_ot_back;
    DataDD.list_ox_forw = list_ox_forw;
    DataDD.list_ox_back = list_ox_back;
    
end



