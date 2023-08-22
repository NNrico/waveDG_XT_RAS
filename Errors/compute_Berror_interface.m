function [Ew_interface] = compute_Berror_interface(Dati,femregion,neighbour,w_h,time,tau)

% initialization
Ew_interface = 0;

% 1D and 2D quadrature nodes and weights
[nodes_1D, w_1D, ~, ~] = quadrature(Dati.nqn);
nln = femregion.nln;
nqn_1D=length(w_1D);


for ie=1:femregion.ne % loop over elements
    
%     id_ie = femregion.id(ie);
    tag_ie = femregion.tag(ie);

    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    coords_elem=femregion.coords_element{ie};
    
    [normals,meshsize]=get_normals_meshsize_faces(coords_elem);
    
    neigh_ie = neighbour.neigh{ie};
    %     neighedges_ie = neighbour.neighedges{ie};
    
    if tag_ie == 'P'
        
        for iedg = 1: neighbour.nedges(ie) % loop over faces
            
            
            if (neigh_ie(iedg) >0 && femregion.tag(neigh_ie(iedg)) == 'A')
                
                
                local_wh1 = w_h(index);
                local_wh2 = w_h(index+femregion.ndof_e);
                
                if iedg<neighbour.nedges(ie)
                    p1 = coords_elem(iedg,:)'; p2 = coords_elem(iedg+1,:)';
                    mean = 0.5*(p1+p2);
                    vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
                    v_hat = [-vy;vx];
                    p3 = mean+v_hat;
                    v = [p1';p2';p3'];
                else
                    p1 = coords_elem(iedg,:)'; p2 = coords_elem(1,:)';
                    mean = 0.5*(p1+p2);
                    vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
                    v_hat = [-vy;vx];
                    p3 = mean+v_hat;
                    v = [p1';p2';p3'];
                end
                
                [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);
                [B_edge,~] = evalshape2D(femregion, ie, pphys_1D);
                
                
%                 for k=1:nqn_1D   % loop over 1D quadrature nodes
                    
                    ds = meshsize(iedg)*w_1D(1:nqn_1D);
                    x = pphys_1D(1:nqn_1D,1);
                    y = pphys_1D(1:nqn_1D,2);
                    t = time;
                    
                    Bedge = B_edge(1:nqn_1D,:);
                    
                    local_exact_w1 = eval(Dati.exact_wp1).*eval(Dati.exact_wpt)';
                    local_exact_w2 = eval(Dati.exact_wp2).*eval(Dati.exact_wpt)';
                    
%                     local_aprox_w1 = 0;
%                     local_aprox_w2 = 0;
                    
%                     for s = 1 : nln  % reconstruct the discrete solution and at the quadrature nodes
%                         
%                         local_aprox_w1 = local_aprox_w1 + Bedge(:,s).*local_wh1(s);
%                         local_aprox_w2 = local_aprox_w2 + Bedge(:,s).*local_wh2(s);
%                         
%                     end

                     local_aprox_w1 = Bedge*local_wh1;
                     local_aprox_w2 = Bedge*local_wh2;


                    Ew_interface = Ew_interface + ds * ((1-tau)/tau * ...
                        ((local_aprox_w1 - local_exact_w1) * normals(1,iedg) ...
                        + (local_aprox_w2 - local_exact_w2) * normals(2,iedg)).^2);
                    
%                 end
                
            end
            
        end
        
    end
    
end

