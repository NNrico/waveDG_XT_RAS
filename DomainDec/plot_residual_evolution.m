% plot residual evolution for RAS and GMRES (pipe=0) or for RAS PIPELINE
% (pipe=1). For the plot of RAS PIPELINE it plot the residual in the domain
% to see the moving subdomain window but only for few iterations. Those
% iterations can be modified in the code.

function [] = plot_residual_evolution(backup,pipe,femregion)

if pipe==0 %RAS,GMRES
 
    figure('Position', [0 0  700 700])
    sgtitle(strcat('Residuals evolution over iterates - mu: ',num2str(backup.info(6))))
    subplot(2,1,1)
    semilogy([1:backup.perf.it_ras],backup.perf.relres2P_vec_ras,'ro-',[1:backup.perf.it_ras],backup.perf.relresinfP_vec_ras,'bo-');
    legend('2','inf')
    xlabel('iterates')
    ylabel('relative preconditoned residual')
    title('XT-RAS')
    subplot(2,1,2)
    semilogy([1:backup.perf.it_gmres],backup.perf.relres2P_vec_gmres(2:end),'ro-');
    legend('2')
    xlabel('iterates')
    ylabel('relative preconditoned residual')
    title('GMRES') 

    % % if you want to save the image
    % mat=fullfile(pwd, strcat('resRAS_GMRES',name(1:end-4),'_subt',num2str(backup.dd(2)),'_subx',num2str(backup.dd(1)),'.png'));
    % saveas(gcf,mat);
else
 
    it_pipe=backup.perf.it_pipe;
    iter=[2:round(it_pipe/5):it_pipe]; %never take it1
    n=backup.info(3)*backup.info(4);  %nx*nt (number of finite elements)
    coord=zeros(n,2);
    for i=1:length(femregion.coords_element)
        vals=femregion.coords_element{i};
        xc=mean(vals(:,1));
        yc=mean(vals(:,2));
        coord(i,1)=xc;
        coord(i,2)=yc;
    end

    if length(iter)==1
        fig=figure('Position', [200 200 700 210]);
    else 
        fig=figure;
    end
    k=1;
    for i= iter
        subplot(length(iter),1,k)
        k=k+1;
        r=backup.perf.mat_res(1:end/2,i);
        m=reshape(full(r),6,n);
        mean_r=sum(m)/6;
        scatter(coord(:,2),coord(:,1),50,mean_r,'filled')
        title(strcat('XT-RAS pipelined iterates=',num2str(i),' #solves= ',num2str(backup.perf.solves_vec(i))))
        colorbar()
        caxis([-10 10])
        sgtitle('Residual evolution in PIPE')
    end
    
end


end