%--------------------------------------------------------------------
% PURPOSE:
%
% Matrix assembly with IPH formulation
%
%
%--------------------------------------------------------------------

function [MatricesAc]= matrix2D_iph(femregion,neighbour,Dati)
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Dati.nqn);
nqn_1D = length(w_1D);
magn1=0;
magn2=0;
penalty_scaled2 = Dati.alp * femregion.fem^2/(Dati.X/Dati.NX);


%-------------ACOUSTIC -------------------------------
M=sparse(femregion.ndof_a,femregion.ndof_a);  % \int_{\Omega} w*v dx
B=sparse(femregion.ndof_a,femregion.ndof_a);  % \int_{\Omega} c2*ux*vx dx
A=sparse(femregion.ndof_a,femregion.ndof_a);  % \int_{\Omega} wt*v dx

I=sparse(femregion.ndof_a,femregion.ndof_a);  % \int_{vertical} {c2*vx} . [u] ds
S=sparse(femregion.ndof_a,femregion.ndof_a);  % space stability \int_{vertical} gamma/2*[u]*[v] ds
S1=sparse(femregion.ndof_a,femregion.ndof_a); % space stability \int_{vertical} (1/2gamma)*[ux]*[vx] ds

W=sparse(femregion.ndof_a,femregion.ndof_a);  % time stability  \int_{low horizontal} [utemp]*v ds

index_shift=0;
ne_a = femregion.ne_a;

id_shift = 0;

counter_el_a = 0;

prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);

%% matrix computations
for ie = 1 : femregion.ne % loop over elements
    
    
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
    Tria_Del = DelaunayTri(coords_elem(:,1),coords_elem(:,2), edges);  %Creates a Delaunay triangulation from a set of points
    io = Tria_Del.inOutStatus();
    Tria = Tria_Del.Triangulation(io==1,:);
    
    
    %%%%%%%%%%%%%%%%%%%
    
    % CALCOLO MATRICI 
    if tag_ie == 'A'
        
        counter_el_a = counter_el_a + 1;
   
        prog = ( 100*(counter_el_a/femregion.ne_a) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
        
        index_a = index ;
        
        for iTria = 1:size(Tria,1)
            
            v1 = coords_elem(Tria(iTria,1),:);
            v2 = coords_elem(Tria(iTria,2),:);
            v3 = coords_elem(Tria(iTria,3),:);
            
            [BJ, ~, pphys_2D] = get_jacobian_physical_points([v1;v2;v3], nodes_2D); %non ricevo BJinv
            Jdet=det(BJ);                       % determinant
            
            % Funzioni di Base
            [dphiq,Grad] = evalshape2D(femregion, ie, pphys_2D);
            
            lw_2D = length(w_2D);                           
            
            dx = w_2D(1:lw_2D)*Jdet;
            x = pphys_2D(1:lw_2D,1);
            y = pphys_2D(1:lw_2D,2);
            
            phi = dphiq(1:lw_2D,:);
            grad_x = Grad(1:lw_2D,1,:);
            grad_y = Grad(1:lw_2D,2,:);
            
            rho_a = Dati.rho_a(id_ie-id_shift);
            c = Dati.c(id_ie-id_shift);
            
            for i=1:femregion.nln % loop over scalar shape functions
        
                for j=1:femregion.nln % loop over scalar shape functions
                                        
                    A_loc=dx *(grad_y(:,:,j).*phi(:,i));  
                    B_loc=dx * (c*c*(grad_x(:,:,j) .* grad_x(:,:,i))); 
                    M_loc=dx * ( phi(:,j) .* phi(:,i) ); 
                
                    A(index_a(i),index_a(j)) = A(index_a(i),index_a(j)) +A_loc;
                    B(index_a(i),index_a(j)) = B(index_a(i),index_a(j)) +B_loc;
                    M(index_a(i),index_a(j)) = M(index_a(i),index_a(j)) +M_loc;
                end
            end
            
        end
        
        IN = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
        SN = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
        S1N = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
        WN = zeros(femregion.nln, femregion.nln, neighbour.nedges(ie));
       
        for iedg=1:neighbour.nedges(ie) % loop over faces
            
            neigedge=neighedges_ie(iedg);    % index of neighbour edge
            id_el_neigh = neigh_ie(iedg);
            if id_el_neigh > 0
                id_neigh = femregion.id(id_el_neigh);
                tag_neigh = femregion.tag(id_el_neigh);
            else
                id_neigh = 0;
                tag_neigh = 'NaN';
            end

            
            
         
            
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
            
            if neigh_ie(iedg) ~= -1
                [B_edge_neigh,G_edge_neigh] = evalshape2D(femregion, neigh_ie(iedg), pphys_1D);
             
            end
                    
            
            ds = meshsize(iedg)*w_1D(1:nqn_1D);
            x = pphys_1D(1:nqn_1D,1);
            y  =pphys_1D(1:nqn_1D,2);
            
            Bedge = B_edge(1:nqn_1D,:);
            Gedge_x = G_edge(1:nqn_1D,1,:);
            Gedge_y = G_edge(1:nqn_1D,2,:);              
            
            aa = 0.5 * normals(1,iedg);            
                        
            if (neigh_ie(iedg) == -1 || femregion.tag(neigh_ie(iedg)) == 'A') %femregion.id(neigh_ie(iedg)) == 2)
                if(normals(2,iedg)>-10e-3 && normals(2,iedg)<10e-3) %vertical edge
                    % to generalize for other meshes
                    
                    for i=1:femregion.nln % loop over scalar shape functions
                        for j=1:femregion.nln % loop over scalar shape functions                                              
                            S(index_a(i),index_a(j)) = S(index_a(i),index_a(j)) + ds * ((penalty_scaled2/2) .* Bedge(:,j) .* Bedge(:,i));
                            S1(index_a(i),index_a(j)) = S1(index_a(i),index_a(j)) + ds * ((1/(2*penalty_scaled2)) .* Gedge_x(:,:,j) .* Gedge_x(:,:,i));
                            magn1=magn1+S(index_a(i),index_a(j)) ;
                            magn2=magn2+S1(index_a(i),index_a(j)) ;
                            
                            if neigh_ie(iedg) >0 % internal faces                         
                                Bedgeneigh = B_edge_neigh(1:nqn_1D,:);
                                Gedgeneigh_x=G_edge_neigh(1:nqn_1D,1,:);
%                                 
                                I(index_a(i),index_a(j)) = I(index_a(i),index_a(j)) + ds * ( c*c*aa*Gedge_x(:,:,i).*Bedge(:,j) ...
                                                                                                +c*c*aa*Gedge_x(:,:,j).*Bedge(:,i));
                                IN(i,j,iedg) = IN(i,j,iedg) + ds*c*c*aa *(Gedgeneigh_x(:,:,j).*Bedge(:,i)  ...   
                                                                                      -Gedge_x(:,:,i).*Bedgeneigh(:,j));                                                  
                                                                                                                                                                                                                            %                                 
                                 SN(i,j,iedg) = SN(i,j,iedg) - ds * ((penalty_scaled2/2) .* Bedge(:,i) .* Bedgeneigh(:,j));
                                 S1N(i,j,iedg) = S1N(i,j,iedg) - ds * ((1/(2*penalty_scaled2)) .* Gedgeneigh_x(:,:,j) .* Gedge_x(:,:,i));
                                 magn1=magn1+SN(i,j,iedg) ;
                                 magn2=magn2+S1N(i,j,iedg) ;
                                
                            elseif neigh_ie(iedg) == -1 % boundary faces
                                I(index_a(i),index_a(j)) = I(index_a(i),index_a(j)) + ds *2*( c*c*aa*Gedge_x(:,:,i).*Bedge(:,j) ... 
                                                                                                   +c*c*aa*Gedge_x(:,:,j).*Bedge(:,i) );
%                                                                                                   
                            end
                            
                        end
                    end
                elseif (normals(1,iedg)>-10e-3 && normals(1,iedg)<10e-3 &&normals(2,iedg)==-1 ) % low horizontal edges
                    
                    
                    for i=1:femregion.nln % loop over scalar shape functions
                        for j=1:femregion.nln % loop over scalar shape functions
                                W(index_a(i),index_a(j)) = W(index_a(i),index_a(j)) + ds * ( Bedge(:,j) .* Bedge(:,i));                  
                                
                                if neigh_ie(iedg) >0                   
                                    Bedgeneigh = B_edge_neigh(1:nqn_1D,:);
                                    WN(i,j,iedg) = WN(i,j,iedg) -ds * ( Bedge(:,i) .* Bedgeneigh(:,j));                         
                                end                      
                        end
                    end                      
                    
                    
                    
                end
                
                

                
    
                
            end
            
        end
        
        %%% ASSEMBLAGGIO 
        [I] = assemble_neigh_A(I, index_a, neigh_ie, IN, femregion.nln, neighbour.nedges(ie), ne_a); % assemble the neighbours local matrices
        [W] = assemble_neigh_A(W, index_a, neigh_ie, WN, femregion.nln, neighbour.nedges(ie), ne_a); % assemble the neighbours local matrices
        [S] = assemble_neigh_A(S, index_a, neigh_ie, SN, femregion.nln, neighbour.nedges(ie), ne_a); % assemble the neighbours local
        [S1] = assemble_neigh_A(S1, index_a, neigh_ie, S1N, femregion.nln, neighbour.nedges(ie), ne_a); % assemble the neighbours local
    
        
    end
    
end
% disp('\nmagn IPH mu/2: ');
% disp(magn1);
% disp('\nmagn IPH 1/2mu: ');
% disp(magn2);

%% COSTRUZIONE MATRICI


MatricesAc = struct('A', A,...
    'M', M,...
    'B', B,...
    'I', I,...
    'W', W,...
    'S', S,...
    'block', [(B-I+S-S1) (A+W+Dati.damp*M);(A+W) -M],... 
    'MPrjA',-1,...
    'DGa',-1); 


fprintf('\n');

