function plot_cos_perfprof_time_mp(fwderr_cosm_adv,fwderr_cosm_exp,...
                          fwderr_cosm_mp, fwderr_cosm_mp_s,...
                          prec, Tcolors, Tstyles, print_legend, legend_perm,...
                          png_out, tikz_out, theta_lim, n)
% This function generates the performance profile of the execution time of 
% the algorithms for computing the cosine in multiprecision. 
%
%               %%Description of some of the inputs%% 
% THETA_LIM: upper range of x-axis in the performance profile

    mp.Digits(prec);
    % rescale values smaller than u
    perfprof_norm = @(x)(max(x, x*(1-5e-2) + 5e-2*mp('eps')/2));

    % performance profiles
    figure
    clf
    T = [perfprof_norm(fwderr_cosm_mp);...
         perfprof_norm(fwderr_cosm_exp);...
         perfprof_norm(fwderr_cosm_mp_s);...
         perfprof_norm(fwderr_cosm_adv)]';
    perfprof(T, theta_lim, Tcolors, Tstyles, 1.0);
    axis([1 theta_lim 0 1])
    if(print_legend)
        local_perm = legend_perm;
        ghH = get(gca, 'Children');
        ghH = flipud(ghH);
        set(gca, 'Children', flipud(ghH(local_perm)));
        methods_names = {'cosm\_mp', 'cosm\_exp', 'cosm\_mp\_s', 'cosm\_adv'};
        lgd = legend(methods_names{local_perm}, 'interpreter', 'latex',...
               'Location', 'SE', 'Orientation', 'vertical', 'FontSize', 11);
        lgd.NumColumns = 1;
        set(gca,'linewidth',1.2)
        set(gca,'fontsize',12)
    end
    xlabel('$\theta$','interpreter','latex','FontWeight','normal','fontsize',18)

    filename_suffix = 'time_mp';
        
    % save .png
    if png_out
        filename = sprintf('../figs/cos_perfprof_%s_%d_n%d', filename_suffix, prec, n);
        saveas(gcf, filename, 'png');
    end

    % tikz
    if tikz_out
        filename = sprintf('../figs/cos_perfprof_%s_%d_n%d', filename_suffix, prec, n);
        matlab2tikz(filename, 'showInfo', false);
    end
end