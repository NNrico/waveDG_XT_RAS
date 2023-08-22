function [neighbor]= neighbours_new(region,Data)

ne =region.ne;
connectivity=region.connectivity;
neigh = cell(1,ne);
neighedges = cell(1,ne);

for i=1:ne
    neigh{i}=-ones(size(connectivity{i}));
    neighedges{i}=-ones(size(connectivity{i}));
end

for i=1:(ne-1)
    edges =[];
    n_edges = length(connectivity{i});
    
    for vertices = 1:n_edges
        v(vertices)=connectivity{i}(vertices);
    end
    
    for e = 1:n_edges-1
        edges(e,:)=[v(e) v(e+1)];
    end
    edges(n_edges,:) = [v(n_edges) v(1)];

    for j=(i+1):ne
        edgesn =[];
        n_edgesn = length(connectivity{j});
        for verticesn = 1:n_edgesn
            vn(verticesn)=connectivity{j}(verticesn);
        end
        
        for e = 1:n_edgesn-1
            edgesn(e,:)=[vn(e) vn(e+1)];
        end
        edgesn(n_edgesn,:) = [vn(n_edgesn) vn(1)];

        for s = 1:size(edges,1)
            for t = 1:size(edgesn,1)
                if (edges(s,1) == edgesn(t,2) && edges(s,2) == edgesn(t,1))
                    neigh{i}(s)=j;
                    neigh{j}(t)=i;
                    neighedges{i}(s)=t;
                    neighedges{j}(t)=s;
                end
            end
        end
        
    end
    
    
end


% Tag_boundary(Data.tag_poro_bc) = Data.lab_poro_bc;
Tag_boundary(Data.tag_ac_bc)   = Data.lab_ac_bc;
% Tag_boundary(Data.tag_el_bc)   = Data.lab_el_bc;

counter = 0;
for i = 1 : ne
%     i
    for j = 1 : size(neigh{i},2)
%         disp([i,j,neigh{i}(j)]);
        if (neigh{i}(j) == -1)
            counter = counter +1;
            n_edges = length(connectivity{i});
    
            for vertices = 1:n_edges
                v(vertices)=connectivity{i}(vertices);
            end
    
            for e = 1:n_edges-1
                edges(e,:)=[v(e) v(e+1)];
            end
            edges(n_edges,:) = [v(n_edges) v(1)];
            
%             disp([j,edges(j,:)]);
            for k = 1 : size(region.connectivity_bc,1) 
                if ((region.connectivity_bc(k,1) == edges(j,1) && (region.connectivity_bc(k,2) == edges(j,2))) || ...
                    (region.connectivity_bc(k,2) == edges(j,1) && (region.connectivity_bc(k,1) == edges(j,2))) )
                 
                    tag_bc = region.bc_tag(k);
                    lab_bc = Tag_boundary(tag_bc);
                    switch lab_bc
                        case('D')
                            neigh{i}(j) = -1;
                        case('N')
                            neigh{i}(j) = -2;
                        case('A')
                            neigh{i}(j) = -3;
                        otherwise
                            disp('Bc not known!')
                    end
                end
            end
        end
    end
end
                            
 
                
    


%===================================================================================
% COSTRUZIONE STRUTTURA NEIGHBOUR
%===================================================================================
neighbor.nedges = region.nedges;
neighbor.neigh = neigh;
neighbor.neighedges = neighedges;