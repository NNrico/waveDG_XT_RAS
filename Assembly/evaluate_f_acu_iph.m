% force term for ip formulation

function [f_a] = evaluate_f_acu_iph(neighbour,femregion,Dati)

[nodes_1D, w_1D, nodes_2D, w_2D]=quadrature(Dati.nqn);
nqn_1D=length(w_1D);

penalty_scaled2 = Dati.alp * (2)^2/(Dati.X/Dati.NX);

f_a = sparse(femregion.ndof_a,1);
f1_a = sparse(femregion.ndof_a,1);

index_shift=0;
id_shift =0;

for ie=1:femregion.ne % loop over elements
    
    id_ie = femregion.id(ie);
    tag_ie = femregion.tag(ie);
    
    index=(ie-1)*femregion.nln*ones(femregion.nln,1) + [1:femregion.nln]';
    index_element = index_shift + [1:1:femregion.nedges(ie)]';
    index_shift = index_element(end);
    
    neigh_ie = neighbour.neigh{ie};
    neighedges_ie = neighbour.neighedges{ie};
    coords_elem = femregion.coords_element{ie};
    
    [normals,meshsize]=get_normals_meshsize_faces(coords_elem);
    edges = [[1:femregion.nedges(ie)]' [2:femregion.nedges(ie) 1]'];
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    
    
    if tag_ie == 'A'
        
        
        index_a = index ;
        
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D);
            Jdet=det(BJ);                       % determinant
            
            % Funzioni di Base
            [dphiq,~] = evalshape2D(femregion, ie, pphys_2D);
            
            %             for k=1:length(w_2D) % loop over 2D quadrature nodes
            lw_2D = length(w_2D);
            
            dx = w_2D(1:lw_2D)*Jdet;
            c=Dati.c;
            x  = pphys_2D(1:lw_2D,1);
            y  = pphys_2D(1:lw_2D,2);
            
            phi = dphiq(1:lw_2D,:);
            
           
            F_A = eval(Dati.source_phi);
            
            
            for i=1:femregion.nln % loop over scalar shape functions
                
                f_a(index_a(i)) = f_a(index_a(i)) + dx * (F_A .* phi(:,i));
                
            end
        end
        %     end
        
        for iedg=1:neighbour.nedges(ie) % loop over faces
            
            neigedge=neighedges_ie(iedg);    % index of neighbour edge
            
           
            
            if iedg<neighbour.nedges(ie)
                p1 = coords_elem(iedg,:)';
                p2 = coords_elem(iedg+1,:)';
                mean = 0.5*(p1+p2);
                vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
                v_hat = [-vy;vx];
                p3 = mean+v_hat;
                v = [p1';p2';p3'];
            else
                p1 = coords_elem(iedg,:)';
                p2 = coords_elem(1,:)';
                mean = 0.5*(p1+p2);
                vx = p2(1)-p2(1); vy = p2(2)-p2(2); v =[vx;vy];
                v_hat = [-vy;vx];
                p3 = mean+v_hat;
                v = [p1';p2';p3'];
            end
            
            [pphys_1D] = get_jacobian_physical_points_faces(v, nodes_1D);
            [B_edge,G_edge] = evalshape2D(femregion, ie, pphys_1D);
                      
            
            ds = meshsize(iedg)*w_1D(1:nqn_1D);
            
            Bedge   = B_edge(1:nqn_1D,:);
            Gedge_x = G_edge(1:nqn_1D,1,:);
            Gedge_y = G_edge(1:nqn_1D,2,:);
            
            x = pphys_1D(1:nqn_1D,1);
            y = pphys_1D(1:nqn_1D,2);
                     
            
            rho_a = Dati.rho_a(id_ie-id_shift);
            
            aa = 0.5  * normals(1,iedg);
            bb = 0.5  * normals(2,iedg);
            
            if  neigh_ie(iedg) == -1 % boundary conditions
                
                for i=1:femregion.nln % loop over scalar shape functions
                    gdphi = eval(Dati.exact_phi);
                    wd= eval(Dati.exact_dphi_y);%                     
                    if (normals(2,iedg)>-10e-3 && normals(2,iedg)<10e-3) 
                        f_a(index_a(i)) = f_a(index_a(i)) + ds * (penalty_scaled2/2 .* Bedge(:,i) .* gdphi ...
                                                                        - (2*(c*c*aa*Gedge_x(:,:,i)).*gdphi )); %S,I 
                       
                    end
                    if(normals(1,iedg)>-10e-3 && normals(1,iedg)<10e-3 && normals(2,iedg)==-1) % to generalize to get low boundary                
                        f_a(index_a(i)) = f_a(index_a(i)) +ds * (  wd.*Bedge(:,i)); 
                        f1_a(index_a(i)) = f1_a(index_a(i)) +ds * (  gdphi.*Bedge(:,i));  
                  
                    end
                   
                end
            end
                        
          
        end
        
        
    end
end
f_a=[f_a;f1_a];





