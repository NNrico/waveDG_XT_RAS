function [errors]=ComputeErrors_new(Dati,femregion,solutions)
[nodes_1D, w_1D, nodes_2D, w_2D] = quadrature(Dati.nqn);
nqn_1D = length(w_1D);


W=sparse(femregion.ndof_a,femregion.ndof_a);  % \int_{orizz_bassi} [utemp]*v ds

index_shift=0;
ne_a = femregion.ne_a;
ne_p = femregion.ne_p;

id_shift = max(Dati.tag_poro);
if(isempty(id_shift)); id_shift = 0; end

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
    
    % CALCOLO MATRICI ACUSTICHE
    if tag_ie == 'A'
        
        counter_el_a = counter_el_a + 1;       
        
        index_a = index - femregion.ndof_p;
        
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
            
            uex=eval(Data.exact_phi);
            
            
            for i=1:femregion.nln % loop over scalar shape functions
        
                for j=1:femregion.nln % loop over scalar shape functions
                    
                    %sembra non venga usato BJinv, non c'e rifereimento?
                    A_loc=dx *(grad_y(:,:,j).*phi(:,i));  %w_t * v
                    B_loc=dx * (c*c*(grad_x(:,:,j) .* grad_x(:,:,i)));  %c2*u_x*v_x
                    M_loc=dx * ( phi(:,j) .* phi(:,i) ); %w*v
                
                    A(index_a(i),index_a(j)) = A(index_a(i),index_a(j)) +A_loc;
                    B(index_a(i),index_a(j)) = B(index_a(i),index_a(j)) +B_loc;
                    M(index_a(i),index_a(j)) = M(index_a(i),index_a(j)) +M_loc;
                end
            end
            %             end
        end
        
        
                    
                    
    end
                
                

                
                %                 end
                
          
            
        
        %%% ASSEMBLAGGIO IT ACUSTICHE
        [I] = assemble_neigh_A(I, index_a, neigh_ie, IN, femregion.nln, neighbour.nedges(ie), ne_p, ne_p + ne_a); % assemble the neighbours local matrices
        [W] = assemble_neigh_A(W, index_a, neigh_ie, WN, femregion.nln, neighbour.nedges(ie), ne_p, ne_p + ne_a); % assemble the neighbours local matrices
        [S] = assemble_neigh_A(S, index_a, neigh_ie, SN, femregion.nln, neighbour.nedges(ie), ne_p, ne_p + ne_a); % assemble the neighbours local
        %    end
        
   
    
end


%% COSTRUZIONE MATRICI


MatricesAc = struct('A', A,...
    'M', M,...
    'B', B,...
    'I', I,...
    'W', W,...
    'S', S,...
    'block', [(B-I+S) (A+W+10*M);(A+W) -M],... %
    'MPrjA',S,...
    'DGa',S); %assemblo a blocchi fuori, le ultime due metto S a caso. poi capire cosa servono
                % c'è ancora da fare f, perchè non farla qua?

%[(B -transpose(I) -I+W+S+S3) (-transpose(A)+W1);(transpose(A)-W1-transpose(W1)+S) (M+S2)]      
                
% [(B -transpose(I) -I +transpose(W)+W+S) (-transpose(A));(transpose(A)-W1+S) M]                
% [(B -transpose(I) -I +S) (-transpose(A)+W1);A M]


fprintf('\n');


end