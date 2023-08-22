% sx: number of starting subdomain in the subdomain window
% check if the corresponding subdomains have reached a convergence

function [sx]=check_sx(x,Rlocal,subx,subt,sx,tol)
if sx > subt 
    disp('err in the definiton of left edge of subdomain window');
    return
end

k_to_check=[sx : subt : subt*(subx-1)+1];
fail=0;

for i=k_to_check    
    err=norm(Rlocal{i}*x,inf);
    if(err > tol)
        fail=1;
        break
    end
end

if fail==0
    sx=sx+1;
end

% to generalize left edge could be composed by more time subdomains

end