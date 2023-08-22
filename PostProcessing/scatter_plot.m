function [] = scatter_plot(GUp,GWp,Gphi,GUe,~,Data) %no pressure ~


if(~isempty(Data.tag_ac))
    figure;
    tri = delaunay(Gphi(:,1),Gphi(:,2));
    subplot(1,2,1)
    hfig = trisurf(tri, Gphi(:,1), Gphi(:,2), Gphi(:,3));
    title('numerical \phi')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
    subplot(1,2,2)
    hfig = trisurf(tri, Gphi(:,1), Gphi(:,2), Gphi(:,4));
    title('analytical')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
    
    figure;
    tri = delaunay(Gphi(:,1),Gphi(:,2));
    subplot(1,2,1)
    hfig = trisurf(tri, Gphi(:,1), Gphi(:,2), Gphi(:,5));
    title('numerical velocity')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
    subplot(1,2,2)
    hfig = trisurf(tri, Gphi(:,1), Gphi(:,2), Gphi(:,6));
    title('analytical')
    axis vis3d;
    lfig = light('Position',[-50 -15 29]);
    lighting phong; shading interp; colorbar EastOutside;
    
    
   
%     figure(500)
%     hfig = trisurf(tri, Pressure.acu(:,1), Pressure.acu(:,2), Pressure.acu(:,3));
%     title('pressure')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     hold on;    
    
    
end


% if(~isempty(Data.tag_poro))
%     figure;
%     tri = delaunay(GUp(:,1),GUp(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,3));
%     title('numerical u_{p,1}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,5));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     
%     
%     figure;
%     tri = delaunay(GUp(:,1),GUp(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,4));
%     title('numerical u_{p,2}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GUp(:,1), GUp(:,2), GUp(:,6));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     
%     figure;
%     tri = delaunay(GWp(:,1),GWp(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,3));
%     title('numerical w_{p,1}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,5));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     
%     
%     figure;
%     tri = delaunay(GWp(:,1),GWp(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,4));
%     title('numerical w_{p,2}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GWp(:,1), GWp(:,2), GWp(:,6));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     
%     
% %     figure(500)
% %     hfig = trisurf(tri, Pressure.poro(:,1), Pressure.poro(:,2), Pressure.poro(:,3));
% %     title('pressure')
% %     axis vis3d;
% %     lfig = light('Position',[-50 -15 29]);
% %     lighting phong; shading interp; colorbar EastOutside;
% %     hold on;    
% 
%     
%     
%     
% end
% 
% 
% if(~isempty(Data.tag_el))
%     figure;
%     tri = delaunay(GUe(:,1),GUe(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,3));
%     title('numerical u_{e,1}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,5));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
% 
%     figure;
%     tri = delaunay(GUe(:,1),GUe(:,2));
%     subplot(1,2,1)
%     hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,4));
%     title('numerical u_{e,2}')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;
%     subplot(1,2,2)
%     hfig = trisurf(tri, GUe(:,1), GUe(:,2), GUe(:,6));
%     title('analytical')
%     axis vis3d;
%     lfig = light('Position',[-50 -15 29]);
%     lighting phong; shading interp; colorbar EastOutside;  
%     
% %     
% %     figure(500)
% %     hfig = trisurf(tri, Pressure.el(:,1), Pressure.el(:,2), Pressure.el(:,3));
% %     title('pressure')
% %     axis vis3d;
% %     lfig = light('Position',[-50 -15 29]);
% %     lighting phong; shading interp; colorbar EastOutside;
% %     hold on;    
% 

    
    
end