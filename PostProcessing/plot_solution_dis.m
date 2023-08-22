function [GUp,GWp,Gphi,GUe] = plot_solution_dis(Dati,femregion,Solutions,time)

t = time;

nln=femregion.nln;

% 1D and 2D quadrature nodes and weights
[~, ~, nodes_2D, w_2D]=quadrature(Dati.nqn);

% G = [ x | y | uh(x,y) | vh(x,y) | uex(x,y) | vex(x,y) ];
GUp  = zeros(1,6);
GWp  = zeros(1,6);
% Gphi = zeros(1,4);
Gphi = zeros(1,6);
GUe  = zeros(1,6);

kindp = 0;
kinda = 0;
kinde = 0;
for ie=1:femregion.ne % loop over elements
    %     ie
%     id_ie = femregion.id(ie);
    tag_ie = femregion.tag(ie);
    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    coords_elem=femregion.coords_element{ie};
    
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    
    

        
    if tag_ie == 'A'
        
        index_a = index ;
        local_phi_h = Solutions.phi_h(index_a);
        local_w_h = Solutions.dot_phi_h(index_a);
        
        for iTria = 1:size(Tria,1)
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
            
            %                 for k=1:length(w_2D)
            
            lw_2d = length(w_2D);   % loop over quadrature nodes
            x = pphys_2D(1:lw_2d,1);
            y = pphys_2D(1:lw_2d,2);
            
%             local_exact_phi=eval(Dati.exact_phi).*eval(Dati.exact_phit);
            local_exact_phi=eval(Dati.exact_phi);
            local_exact_w=eval(Dati.exact_dphi_y);
            local_aprox_phi=0;
            local_aprox_w=0;
            
            %for s=1:nln  % reconstruct the discrete solution at the quadrature nodes
            %    local_aprox_phi = local_aprox_phi + dphiq(k,s).*local_phi_h(s);
            %end
            local_aprox_phi = local_aprox_phi + dphiq(1:lw_2d,1:nln)*local_phi_h(1:nln);
            local_aprox_w = local_aprox_w + dphiq(1:lw_2d,1:nln)*local_w_h(1:nln);
            
            
            Gphi(kinda + 1 : kinda + lw_2d,:) = [x, y, local_aprox_phi, local_exact_phi,local_aprox_w, local_exact_w];
            kinda = kinda + lw_2d;
            
            
            %                 end
        end
       
    end
    
    
    
end
end


%
%
