function plot_cos_time_double(time_cosm, time_cosm_exp, time_cosm_mp, time_cosm_mp_s,...
                 time_cosm_tay, time_cosm_pol, perm, Tcolors,...
                 Tmarkers, msize, lw, print_yticks, print_legend,...
                 legend_perm, png_out, tikz_out, y_lim, n)
% This function generates the plots of the execution time of the algorithms
% for computing the cosine in double precision. 
%
%               %%Description of some of the inputs%% 
% PERM: a vector sorting the errors in decreasing order of TIME_COSM_ADV
% Y_LIM: upper range of the y-axis in the plot
    figure
    clf;
    n_mats = size(time_cosm_mp,2);     
    mats = 1:1:n_mats;
  
    hold on;
    
    semilogy(mats, time_cosm_mp_s(perm),...
                Tmarkers{6}, 'Color', Tcolors(6,:), 'MarkerSize', msize, 'Linewidth', lw);
    semilogy(mats, time_cosm_mp(perm),...
                Tmarkers{5}, 'Color', Tcolors(5,:), 'MarkerSize', msize, 'Linewidth', lw);
    semilogy(mats, time_cosm(perm),...
                Tmarkers{3}, 'Color', Tcolors(3,:), 'MarkerSize', msize, 'Linewidth', lw);
    semilogy(mats, time_cosm_exp(perm),...
                Tmarkers{4}, 'Color', Tcolors(4,:), 'MarkerSize', msize, 'Linewidth', lw);
    semilogy(mats, time_cosm_tay(perm),...
                Tmarkers{2}, 'Color', Tcolors(2,:), 'MarkerSize', msize, 'Linewidth', lw);
    semilogy(mats, time_cosm_pol(perm),...
                Tmarkers{1}, 'Color', Tcolors(1,:), 'MarkerSize', msize, 'Linewidth', lw);
%     set(gca,'XTick',0)
 
    if ~print_yticks
        set(gca,'YTickLabels',{})
    end
    
    if(print_legend)
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(legend_perm)));
        methods_names = {'cosm\_mp\_s', 'cosm\_mp', 'cosm', 'cosm\_exp', 'cosm\_tay', 'cosm\_pol'};
        legend(methods_names{legend_perm}, 'interpreter', 'latex', 'Location', 'NW', 'FontSize', 11);
        set(gca,'linewidth',1.2)
        set(gca,'fontsize',12)
    end
    xlim([min(mats)-1, max(mats)+1]);
    ylim([0 y_lim])
    
    xlabel_name = sprintf('n=%d', n);
    
    xlabel(xlabel_name,'interpreter','latex','FontWeight','normal','fontsize',18)
%     ylabel('Execution time','fontsize',18)
    box on     
    
    filename_suffix = 'time_double';
        
    % save .png
    if png_out
        filename = sprintf('../figs/cos_%s_n%d', filename_suffix, n);
        saveas(gcf, filename, 'png');
    end

    % tikz
    if tikz_out
        filename = sprintf('../figs/cos_%s_n%d', filename_suffix, n);
        matlab2tikz(filename, 'showInfo', false);
    end

end