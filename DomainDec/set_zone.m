 % set subdomains inside the subdomain window (zone)
function sub_zone = set_zone(subx,subt,sx,dx)

matrix_sub = reshape(1:subx*subt,subt,subx);
zonematrix= matrix_sub(sx:(dx-1),:);
sub_zone = reshape(zonematrix,1,[]);

end