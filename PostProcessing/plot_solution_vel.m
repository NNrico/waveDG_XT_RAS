function [dot_GUp,dot_GWp,dot_Gphi,dot_GUe] = plot_solution_vel(Dati,femregion,Solutions,time)

t = time;

nln=femregion.nln;

% 1D and 2D quadrature nodes and weights
[~, ~, nodes_2D, w_2D]=quadrature(Dati.nqn);

% G = [ x | y | uh(x,y) | vh(x,y) | uex(x,y) | vex(x,y) ];
GUp  = zeros(1,6);
GWp  = zeros(1,6);
Gphi = zeros(1,4);
GUe  = zeros(1,6);

kindp = 0;
kinda = 0;
kinde = 0;
for ie=1:femregion.ne % loop over elements
    %     ie
%     id_ie = femregion.id(ie);
    tag_ie = femregion.tag(ie);
    
%     BBox_ie = femregion.BBox(ie,:);
    %     if( reciever(1) >= BBox_ie(1) && reciever(1) <= BBox_ie(2) ...
    %             &&  reciever(2) >= BBox_ie(3) && reciever(2) <= BBox_ie(4))
    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    coords_elem=femregion.coords_element{ie};
    
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    
    
    if tag_ie == 'P' %id_ie == 1
        
        local_uh1 = Solutions.dot_up_h(index);
        local_uh2 = Solutions.dot_up_h(index+femregion.ndof_p);
        local_wh1 = Solutions.dot_wp_h(index);
        local_wh2 = Solutions.dot_wp_h(index+femregion.ndof_p);
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
            
            lw_2d = length(w_2D);
            x=pphys_2D(1:lw_2d,1);
            y=pphys_2D(1:lw_2d,2);
            
            local_exact_u1 = eval(Dati.exact_up1)*eval(Dati.exact_upt);
            local_exact_u2 = eval(Dati.exact_up2)*eval(Dati.exact_upt);
            local_exact_w1 = eval(Dati.exact_wp1)*eval(Dati.exact_wpt);
            local_exact_w2 = eval(Dati.exact_wp2)*eval(Dati.exact_wpt);
            
            local_aprox_u1 = dphiq(1:lw_2d,1:nln)*local_uh1(1:nln);
            local_aprox_u2 = dphiq(1:lw_2d,1:nln)*local_uh2(1:nln);
            local_aprox_w1 = dphiq(1:lw_2d,1:nln)*local_wh1(1:nln);
            local_aprox_w2 = dphiq(1:lw_2d,1:nln)*local_wh2(1:nln);
            
            GUp(kindp + 1 : kindp + lw_2d,:) = [x, y, local_aprox_u1, local_aprox_u2, local_exact_u1, local_exact_u2];
            GWp(kindp + 1 : kindp + lw_2d,:) = [x, y, local_aprox_w1, local_aprox_w2, local_exact_w1, local_exact_w2];
            
            kindp = kindp + lw_2d;
            %                 end
            
        end
        %             for iTria = 1:size(Tria,1)
        %                 v1 = coords_elem(Tria(iTria,1),:);
        %                 v2 = coords_elem(Tria(iTria,2),:);
        %                 v3 = coords_elem(Tria(iTria,3),:);
        %                 [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        %                 Jdet = det(BJ);
        %                 [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
        %
        %
        %                 for k=1:length(w_2D) % loop over quadrature nodes
        %                     x=pphys_2D(k,1);
        %                     y=pphys_2D(k,2);
        %
        %                     local_exact_u1 = eval(Dati.exact_up1)'*eval(Dati.exact_upt)';
        %                     local_exact_u2 = eval(Dati.exact_up2)'*eval(Dati.exact_upt)';
        %                     local_exact_w1 = eval(Dati.exact_wp1)'*eval(Dati.exact_wpt)';
        %                     local_exact_w2 = eval(Dati.exact_wp2)'*eval(Dati.exact_wpt)';
        %
        %                     local_aprox_u1 = 0;
        %                     local_aprox_u2 = 0;
        %                     local_aprox_w1 = 0;
        %                     local_aprox_w2 = 0;
        %
        %                     for s=1:nln  % reconstruct the discrete solution at the quadrature nodes
        %                         local_aprox_u1 = local_aprox_u1 + dphiq(k,s).*local_uh1(s);
        %                         local_aprox_u2 = local_aprox_u2 + dphiq(k,s).*local_uh2(s);
        %                         local_aprox_w1 = local_aprox_w1 + dphiq(k,s).*local_wh1(s);
        %                         local_aprox_w2 = local_aprox_w2 + dphiq(k,s).*local_wh2(s);
        %                     end
        %
        %                     GUp(kindp,:) = [x, y, local_aprox_u1, local_aprox_u2, local_exact_u1, local_exact_u2];
        %                     GWp(kindp,:) = [x, y, local_aprox_w1, local_aprox_w2, local_exact_w1, local_exact_w2];
        %
        %                     kindp = kindp +1 ;
        %                 end
        
        
    elseif tag_ie == 'A'
        
        index_a = index - femregion.ndof_p;
        local_phi_h = Solutions.dot_phi_h(index_a);
        
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
            
            local_exact_phi=eval(Dati.exact_phi).*eval(Dati.exact_phit);
            local_aprox_phi=0;
            
            %for s=1:nln  % reconstruct the discrete solution at the quadrature nodes
            %    local_aprox_phi = local_aprox_phi + dphiq(k,s).*local_phi_h(s);
            %end
            local_aprox_phi = local_aprox_phi + dphiq(1:lw_2d,1:nln)*local_phi_h(1:nln);
            
            
            Gphi(kinda + 1 : kinda + lw_2d,:) = [x, y, local_aprox_phi, local_exact_phi];
            kinda = kinda + lw_2d;
            
            
            %                 end
        end
        %             for iTria = 1:size(Tria,1)
        %                 v1 = coords_elem(Tria(iTria,1),:);
        %                 v2 = coords_elem(Tria(iTria,2),:);
        %                 v3 = coords_elem(Tria(iTria,3),:);
        %                 [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
        %
        %                 [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
        %
        %                 for k=1:length(w_2D) % loop over quadrature nodes
        %                     x=pphys_2D(k,1);
        %                     y=pphys_2D(k,2);
        %                     local_exact_phi = eval(Dati.exact_phi)'*eval(Dati.exact_phit)';
        %
        %                     local_aprox_phi = 0;
        %                     local_aprox_phi = local_aprox_phi + dphiq(k,1:nln)*local_phi_h(1:nln);
        %
        %                     Gphi(kinda,:) = [x, y, local_aprox_phi, local_exact_phi];
        %                     kinda = kinda+1;
        %
        %
        %                 end
        %             end
        
    elseif tag_ie == 'E' %id_ie == 1
        
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        
        local_uh1 = Solutions.dot_ue_h(index_e);
        local_uh2 = Solutions.dot_ue_h(index_e+femregion.ndof_e);
        
        
        for iTria = 1:size(Tria,1)
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
            
            
            %                 for k=1:length(w_2D) % loop over quadrature nodes
            lw_2d = length(w_2D);
            x     = pphys_2D(1:lw_2d,1);
            y     = pphys_2D(1:lw_2d,2);
            
            local_exact_u1 = eval(Dati.exact_ue1).*eval(Dati.exact_uet)';
            local_exact_u2 = eval(Dati.exact_ue2).*eval(Dati.exact_uet)';
            
%             local_aprox_u1 = 0;
%             local_aprox_u2 = 0;
%             
%             for s=1:nln  % reconstruct the discrete solution at the quadrature nodes
%                 local_aprox_u1 = local_aprox_u1 + dphiq(k,s).*local_uh1(s);
%                 local_aprox_u2 = local_aprox_u2 + dphiq(k,s).*local_uh2(s);
%             end

            local_aprox_u1 = dphiq(1:lw_2d,1:nln)*local_uh1(1:nln);
            local_aprox_u2 = dphiq(1:lw_2d,1:nln)*local_uh2(1:nln);

            GUe(kinde + 1 : kinde + lw_2d,:) = [x, y, local_aprox_u1, local_aprox_u2, local_exact_u1, local_exact_u2];           
%             GUe(kinde,:) = [x, y, local_aprox_u1, local_aprox_u2, local_exact_u1, local_exact_u2];
            
            kinde = kinde + lw_2d ;
            %                 end
            
        end
    end
    
    
    
end


%
%
