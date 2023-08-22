function [up,wp,phi_a,ue] = evaluate_solution_new(Data,femregion)

[~, ~, nodes_2D, w_2D]=quadrature(Data.nqn);


up1 = sparse(femregion.ndof_p,1);
up2 = sparse(femregion.ndof_p,1);
wp1 = sparse(femregion.ndof_p,1);
wp2 = sparse(femregion.ndof_p,1);

phi_a = sparse(femregion.ndof_a,1);

ue1 = sparse(femregion.ndof_e,1);
ue2 = sparse(femregion.ndof_e,1);


index_shift=0;

prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);
counter_el = 0;


for ie=1:femregion.ne
    
    counter_el = counter_el + 1;
    prog = ( 100*(counter_el/femregion.ne) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    
    
    tag_ie = femregion.tag(ie);
    id_ie  = femregion.id(ie);
    
    index         = (ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift   = index_element(end);
    
    coords_elem = femregion.coords_element{ie};
    
    edges    = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io       = Tria_Del.inOutStatus();
    Tria     = Tria_Del.Triangulation(io==1,:);
    
    
    if tag_ie == 'P'
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            Jdet = det(BJ);
            
            % Funzioni di Base
            [dphiq,~] = evalshape2D(femregion, ie, pphys_2D);
            
            %             for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            lw_2d = length(w_2D);
            
            dx = w_2D(1:lw_2d)*Jdet;
            x  = pphys_2D(1:lw_2d,1);
            y  = pphys_2D(1:lw_2d,2);
            
            phi = dphiq(1:lw_2d,:);
            
            U1 = eval(Data.exact_up1);
            U2 = eval(Data.exact_up2);
            W1 = eval(Data.exact_wp1);
            W2 = eval(Data.exact_wp2);
            
            for i = 1:femregion.nln % loop over scalar shape functions
                
                up1(index(i)) = up1(index(i)) + dx * (U1 .* phi(:,i));
                up2(index(i)) = up2(index(i)) + dx * (U2 .* phi(:,i));
                wp1(index(i)) = wp1(index(i)) + dx * (W1 .* phi(:,i));
                wp2(index(i)) = wp2(index(i)) + dx * (W2 .* phi(:,i));
                
                
            end
            %             end
        end
        
    elseif tag_ie == 'A'
        
        index_a = index - femregion.ndof_p;
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            Jdet=det(BJ);
            
            % Funzioni di Base
            [dphiq,~] = evalshape2D(femregion, ie, pphys_2D);
            
            %             for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            lw_2d = length(w_2D);
            
            dx = w_2D(1:lw_2d)*Jdet;
            x  = pphys_2D(1:lw_2d,1);
            y  = pphys_2D(1:lw_2d,2);
            
            phi = dphiq(1:lw_2d,:);
            
            Phi_ex = eval(Data.exact_phi);
            
            for i=1:femregion.nln
                
                phi_a(index_a(i)) = phi_a(index_a(i)) + dx * (Phi_ex .* phi(:,i));
                
            end
            %             end
        end
        
    elseif tag_ie == 'E'
        
        index_e = index - femregion.ndof_p - femregion.ndof_a;
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            Jdet=det(BJ);
            
            % Funzioni di Base
            [dphiq,~] = evalshape2D(femregion, ie, pphys_2D);
            
            %             for k=1:length(w_2D) % loop over 2D quadrature nodes
            
            lw_2d = length(w_2D);
            
            dx = w_2D(1:lw_2d)*Jdet;
            x  = pphys_2D(1:lw_2d,1);
            y  = pphys_2D(1:lw_2d,2);
            
            phi = dphiq(1:lw_2d,:);
            
            U1 = eval(Data.exact_ue1);
            U2 = eval(Data.exact_ue2);
            
            for i=1:femregion.nln % loop over scalar shape functions
                
                ue1(index_e(i)) = ue1(index_e(i)) + dx * (U1 .* phi(:,i));
                ue2(index_e(i)) = ue2(index_e(i)) + dx * (U2 .* phi(:,i));
                
            end
            
            %             end
        end
        
    end
end
up = [up1; up2];
wp = [wp1; wp2];
ue = [ue1; ue2];

fprintf('\n');

