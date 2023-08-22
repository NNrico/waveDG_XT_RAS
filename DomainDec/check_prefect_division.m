
% Check if there is a perfect division. This means if there is a costant
% overlap for all subdomains.

function [ok]=check_prefect_division(DataDD)
ok=1;
[isComp_x,isComp_t]= compatible_division(DataDD); 
if(isComp_x && isComp_t)
    disp('perfect division')
%     disp('compatible division with grids that are multiple of the subdomain')
else
    disp('error! Not perfect division')
%     disp('NO compatible division with grids that are multiple of the subdomain')
    ok=0;
end
disp('------------------------------------------------------');
end

function [isComp_x,isComp_t]= compatible_division(DataDD)
    isComp_x=false;
    isComp_t=false;
    if(DataDD.m>0 && DataDD.m>DataDD.ot && DataDD.n>0 && DataDD.n>DataDD.ox && DataDD.ox>=0 && DataDD.ot>=0 && DataDD.ot<=m/2 && DataDD.ox<=n/2)
        last_t=DataDD.NT;
        last_x=DataDD.Nx;
        isComp_x=isCompatible(last_x,DataDD.n,DataDD.ox);
        isComp_t=isCompatible(last_t,DataDD.m,DataDD.ot);
    end
end

function [isComp]=isCompatible(last,m,ov)
    curr=m-ov;
    isComp=false;
    while(~isComp && curr+m<=last)
        if(curr+m ==last)
            isComp=true;
        else
            curr=curr+m-ov;            
        end
    end
end
