% action of the preconditioner on x (the residual)
% preconditioner: sum(Rtilde_i * A_i * R_i') but the sum is over only few
%                 subdomains inside the subdomain window 

% two steps:
% - first define the subdomains in the subdomain window: check for
% convergence the left edge and update the right edge after it_wait
% iterations (aggressive version)
% - compute the action of the preconditioner on x



function [z,sx,zone,it_waited,solves]=precondAction_pipe(x,Alocal,Rlocal,Rtlocal,subx,subt,sx,tol,it_wait,it_waited,zone,solves)   
    % sx: number of starting subdomain
    % zone: subdomain window width
    
    z=0; 
    [sx1]=check_sx(x,Rlocal,subx,subt,sx,tol); %check if we can update left edge of subdomain window
    f=0; %flag to control update sx
    if sx1 >sx
        sx=sx1;
        zone=zone-1;
        f=1;
    end

    
    % isend: flag to check if we are at the end
    if zone+sx==subt+1
        isend=1;
    else
        isend = 0;
    end

    if ~isend
     if zone==0
        zone = zone + 1;
     end
    end

    %update dx 
    if ~isend
        if it_waited>=it_wait %I have waited enough, now update right edge of subdomain window
            zone=zone+1;        
            it_waited=0;
        elseif f==0 && it_waited<it_wait
            it_waited=it_waited+1;          
        elseif f==1 && it_waited<it_wait
            it_waited=0;
            zone=zone+1;
        end
    end

    if zone<=-1
        disp('zone is negative')
        return
    end
       
    %nsubt > = sx+2 hypotesis

    dx=sx+zone;
    sub_zone=set_zone(subx,subt,sx,dx);  % set subdomains inside the subdomain window
    solves=solves+length(sub_zone);
    
    % action of the preconditioner. now the preconditioner sees only
    % sub_zone subdomains
    for k=sub_zone        
        xk=Rlocal{k}*x;
        uk=Alocal{k}\xk;
        zk=transpose(Rtlocal{k})*uk;
        z=z+zk;
    end


end