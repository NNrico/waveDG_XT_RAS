function [Pressure]= compute_pressure(Dati,femregion,Solutions)

nln = femregion.nln;

% 1D and 2D quadrature nodes and weights
[~, ~, nodes_2D, w_2D]=quadrature(Dati.nqn);

ph_poro = zeros(1,3);
ph_acu  = zeros(1,3);
ph_el   = zeros(1,3);


kindp = 0;
kinda = 0;
kinde = 0;

id_shift = max(Dati.tag_poro);
if(isempty(id_shift)); id_shift = 0; end

id_shift_e = max([Dati.tag_ac,Dati.tag_poro]);
if(isempty(id_shift_e)); id_shift_e = 0; end

prog = 0;
counter_el_p = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

for ie=1:femregion.ne % loop over elements
    
    id_ie = femregion.id(ie);
    tag_ie = femregion.tag(ie);
    
    counter_el_p = counter_el_p + 1;
    prog = ( 100*(counter_el_p/femregion.ne) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    %BBox_ie = femregion.BBox(ie,:);
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    coords_elem=femregion.coords_element{ie};
    
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    
    
    if strcmp(tag_ie,'P')
        
        
        local_uh1 = Solutions.up_h(index);
        local_uh2 = Solutions.up_h(index+femregion.ndof_p);
        local_wh1 = Solutions.wp_h(index);
        local_wh2 = Solutions.wp_h(index+femregion.ndof_p);
        
        for iTria = 1:size(Tria,1)
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [~,Grad]= evalshape2D(femregion, ie, pphys_2D);
            
            lw_2d = length(w_2D);
            x  = pphys_2D(1:lw_2d,1);
            y  = pphys_2D(1:lw_2d,2);
            beta = Dati.beta(id_ie);
            m = Dati.m(id_ie);
            
            grad_x = Grad(1:lw_2d,1,:);
            grad_y = Grad(1:lw_2d,2,:);
            
            local_aprox_du1_x = 0; local_aprox_du2_y = 0;
            local_aprox_dw1_x = 0; local_aprox_dw2_y = 0;
            
            for i = 1 : nln
                local_aprox_du1_x = local_aprox_du1_x + grad_x(:,:,i)*local_uh1(i);
                local_aprox_du2_y = local_aprox_du2_y + grad_y(:,:,i)*local_uh2(i);
                local_aprox_dw1_x = local_aprox_dw1_x + grad_x(:,:,i)*local_wh1(i);
                local_aprox_dw2_y = local_aprox_dw2_y + grad_y(:,:,i)*local_wh2(i);
            end
            
            
            local_u_div = local_aprox_du1_x + local_aprox_du2_y;
            local_w_div = local_aprox_dw1_x + local_aprox_dw2_y;
            
            ph_poro(kindp + 1 : kindp + lw_2d,:) = [x, y, - m * (beta * local_u_div + local_w_div)];
            kindp = kindp + lw_2d;
            
            
        end
        
    elseif strcmp(tag_ie,'A')
        
        index_a = index - femregion.ndof_p;
        rho_a = Dati.rho_a(id_ie-id_shift);
        
        local_dot_phi_h = Solutions.dot_phi_h(index_a);
        
        for iTria = 1:size(Tria,1)
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [dphiq,~]= evalshape2D(femregion, ie, pphys_2D);
            
            lw_2d = length(w_2D);   % loop over quadrature nodes
            x = pphys_2D(1:lw_2d,1);
            y = pphys_2D(1:lw_2d,2);
            
            local_aprox_phi=0;
            
            local_aprox_phi = local_aprox_phi + dphiq(1:lw_2d,1:nln)*local_dot_phi_h(1:nln)*rho_a;
            
            
            ph_acu(kinda + 1 : kinda + lw_2d,:) = [x, y, local_aprox_phi];
            kinda = kinda + lw_2d;
            
        end
        
    elseif strcmp(tag_ie,'E')
        
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        
        local_uh1 = Solutions.ue_h(index_e);
        local_uh2 = Solutions.ue_h(index_e+femregion.ndof_e);
        
        
        for iTria = 1:size(Tria,1)
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [~, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            
            [~,Grad]= evalshape2D(femregion, ie, pphys_2D);
            
            lw_2d = length(w_2D);
            x  = pphys_2D(1:lw_2d,1);
            y  = pphys_2D(1:lw_2d,2);
            
            lambda = Dati.lam_el(id_ie-id_shift_e);
            
            grad_x = Grad(1:lw_2d,1,:);
            grad_y = Grad(1:lw_2d,2,:);
            local_aprox_du1_x = 0; local_aprox_du2_y = 0;
            
            for i = 1 : nln
                local_aprox_du1_x = local_aprox_du1_x + grad_x(:,:,i)*local_uh1(i);
                local_aprox_du2_y = local_aprox_du2_y + grad_y(:,:,i)*local_uh2(i);
            end            
            local_u_div = local_aprox_du1_x + local_aprox_du2_y;
            
            ph_el(kinde + 1 : kinde + lw_2d,:) = [x, y, -lambda * (local_u_div)];
            kinde = kinde + lw_2d;
            
            
        end
        
    end
    
end


Pressure.poro = ph_poro;
Pressure.acu  = ph_acu;
Pressure.el   = ph_el;




