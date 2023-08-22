function [ok] = check_dd_choice(test,overlap_nsub,nx,nt)
ok=1;    
if overlap_nsub==0
    if length(test.ot)~=length(test.ox) || length(test.m)~=length(test.n) || length(test.m)~=length(test.ot) 
        ok=0;
    end
    if ok==1
        for i=1:max([length(test.m),length(test.n),length(test.ot)])       
            if test.ot(i) > test.m(i) || test.ox(i) > test.n(i)
                ok=0;
            end
        end
    end
else
    if length(test.nsub_x)~=length(test.nsub_t) || length(test.m)~=length(test.n) || length(test.m)~=length(test.nsub_t) 
        ok=0;
    end
    if ok==1
        for i=1:max([length(test.m),length(test.n),length(test.nsub_x)])
            if nt > test.m(i)*test.nsub_t(i)+2 || nx > test.n(i)*test.nsub_x(i)+2
                ok=0;
            end
        end
    end
end

if ok==0
    disp('error in dd choice')

end