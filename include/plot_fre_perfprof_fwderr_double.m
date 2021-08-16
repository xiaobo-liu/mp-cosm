function plot_fre_perfprof_fwderr_double(fwderr_fre_blk, fwderr_fre_exp,...
                          fwderr_fre_mp, fwderr_fre_mp_s, prec, Tcolors,...
                          Tstyles, print_legend, legend_perm, png_out,...
                          tikz_out, theta_lim)
% This function generates the performance profile of the errors in the 
% computed Frechet derivative in double precision. 
%
%               %%Description of some of the inputs%% 
% PREC: digits of precision
% THETA_LIM: upper range of x-axis in the performance profile 

    mp.Digits(prec);
    perfprof_norm = @(x)(max(x, x*(1-5e-2) + 5e-2*mp('eps')/2));

    % performance profiles
    figure
    clf
    T = [perfprof_norm(fwderr_fre_mp);...
         perfprof_norm(fwderr_fre_blk);...
         perfprof_norm(fwderr_fre_mp_s);...
         perfprof_norm(fwderr_fre_exp)]';
    perfprof(T, theta_lim, Tcolors, Tstyles, 1.0);
    axis([1 theta_lim 0 1])
    if(print_legend)
        local_perm = legend_perm;
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(local_perm)));
        methods_names = {'cosm\_fre\_mp', 'cosm\_fre\_blk', 'cosm\_fre\_mp\_s', 'cosm\_fre\_exp'};
        lgd = legend(methods_names{local_perm}, 'interpreter', 'latex',...
               'Location', 'SE', 'Orientation', 'vertical', 'FontSize', 11);
        lgd.NumColumns = 1;
        set(gca,'linewidth',1.2)
        set(gca,'fontsize',12)
    end
    xlabel('$\theta$','interpreter','latex','FontWeight','normal','fontsize',18)

    % save .png
    if png_out
        filename = sprintf('../figs/fre_perfprof_fwderr_double');
        saveas(gcf, filename, 'png');
    end

    % tikz
    if tikz_out
        filename = sprintf('../figs/fre_perfprof_fwderr_double.tikz');
        matlab2tikz(filename, 'showInfo', false);
    end
end