function [] = scatter_plot_phis(GUp,GWp,GUe,Data)


if(~isempty(Data.tag_poro))
    figure(200);
    tri = delaunay(GUp(:,1),GUp(:,2));
    hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,3));
    title('numerical u_{p,1}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside; hold on;
    
    
    figure(300);
    tri = delaunay(GUp(:,1),GUp(:,2));
    hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,4));
    title('numerical u_{p,2}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside; hold on;
    

    figure(250);
    tri = delaunay(GWp(:,1),GWp(:,2));
    hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,3));
    title('numerical w_{p,1}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside; hold on;
    
    
    figure(350);
    tri = delaunay(GWp(:,1),GWp(:,2));
    hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,4));
    title('numerical w_{p,2}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside; hold on;
    
    
    
    
end


if(~isempty(Data.tag_el))
    figure(200);
    tri = delaunay(GUe(:,1),GUe(:,2));
    hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,3));
    title('numerical u_{e,1}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;

    figure(300);
    tri = delaunay(GUe(:,1),GUe(:,2));
    hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,4));
    title('numerical u_{e,2}')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
       
    
end