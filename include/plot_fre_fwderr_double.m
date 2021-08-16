function plot_fre_fwderr_double(condest1u, perm, fwderr_fre_blk,...
                            fwderr_fre_exp, fwderr_fre_mp, fwderr_fre_mp_s,...
                            prec, ymin, ymax, Tcolors, color_cond, Tmarkers,...
                            ls_cond, msize, lw, lw_cond, print_yticks, print_legend,...
                            legend_perm, png_out, tikz_out)
% This function generates the plots of the errors of the computed Frechet
% derivative in double precision. 
%
%               %%Description of some of the inputs%% 
% CONDESTLU: a vector storing condition number of the matrices * unit roundoff
% PERM: a vector sorting the errors in decreasing order of CONDESTLU
% PREC: digits of precision
% [YMIN, YMAX]: range of the errors in the plots

    fwderr_fre_blk(fwderr_fre_blk > ymax) = ymax;
    fwderr_fre_exp(fwderr_fre_exp > ymax) = ymax;
    fwderr_fre_mp(fwderr_fre_mp > ymax) = ymax;
    fwderr_fre_mp_s(fwderr_fre_mp_s > ymax) = ymax;
    
    fwderr_fre_blk(fwderr_fre_blk < ymin) = ymin;
    fwderr_fre_exp(fwderr_fre_exp < ymin) = ymin;
    fwderr_fre_mp(fwderr_fre_mp < ymin) = ymin;
    fwderr_fre_mp_s(fwderr_fre_mp_s < ymin) = ymin;
    
    % forward error
    condest1u(isnan(condest1u)) = Inf;
    figure
    clf;
    mats = 1:1:length(condest1u);
    hold on;

    mp_semilogy(mats, condest1u(perm),...
                [0 200 ymin ymax], 0, ls_cond, 'Color', color_cond, 'MarkerSize', msize,'Linewidth', lw_cond);
    mp_semilogy(mats, fwderr_fre_mp(perm),...
                [0 200 ymin ymax], 0, Tmarkers{1}, 'Color', Tcolors(1,:), 'MarkerSize', msize, 'Linewidth', lw);
    mp_semilogy(mats, fwderr_fre_blk(perm),...
                [0 200 ymin ymax], 0, Tmarkers{2}, 'Color', Tcolors(2,:), 'MarkerSize', msize, 'Linewidth', lw);  
    mp_semilogy(mats, fwderr_fre_mp_s(perm),...
                [0 200 ymin ymax], 0, Tmarkers{3}, 'Color', Tcolors(3,:), 'MarkerSize', msize, 'Linewidth', lw);
    mp_semilogy(mats, fwderr_fre_exp(perm),...
                [0 200 ymin ymax], 0, Tmarkers{4}, 'Color', Tcolors(4,:), 'MarkerSize', msize, 'Linewidth', lw);
    
%     set(gca,'XTick',0)
    if ~print_yticks
        set(gca,'YTickLabels',{})
    end
    mp.Digits(prec);

    if(print_legend)
        local_perm = [1, legend_perm+1];
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(local_perm)));
        methods_names = {'$\kappa_{L}(A,E)u$', 'cosm\_fre\_mp', 'cosm\_fre\_blk',...
            'cosm\_fre\_mp\_s', 'cosm\_fre\_exp'};
        legend(methods_names{local_perm}, 'interpreter', 'latex',...
               'Location', 'NE', 'FontSize', 11);
        set(gca,'linewidth',1.2)
        set(gca,'fontsize',12)
    end
    box on

    % save .png
    if png_out
        filename = sprintf('../figs/fre_fwderr_double');
        saveas(gcf, filename, 'png');
    end

    % save .tikz
    if tikz_out
        filename = sprintf('../figs/fre_fwderr_double.tikz');
        matlab2tikz(filename, 'showInfo', false,...
                    'extraTikzpictureOptions',...
                    {'trim axis left', 'trim axis right'});
    end

end