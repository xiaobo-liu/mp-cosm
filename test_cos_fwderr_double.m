% this code tests the error of the computed cosine in double precision

format compact

addpath('data');
addpath('external');
addpath('include');
warning('off');

%% the following codes should be commented when run all the tests at one time

initialize_tests

reproduce_results = false; % reproduce the test results?


%% set output parameters and run the tests

theta_lim = 25; % the range of x-axis in the performance profile generated

png_out = true;

tikz_out = false;

ymin = 1e-18;
ymax = 1e-2;

compute_condest = true;

ids_min = 1;
[~, ids_max] = testmats();
ids_max = ids_max * 2; % size doubled by multiplying the imaginary unit i

digits_work = 16;
digits_ref = 10 * digits_work;


if reproduce_results

    n_mats = ids_max - ids_min + 1;

    fwderr_cosm = zeros(1, n_mats);
    fwderr_cosm_adv = zeros(1, n_mats);
    fwderr_cosm_exp = zeros(1, n_mats);
    fwderr_cosm_mp = zeros(1, n_mats);
    fwderr_cosm_mp_s = zeros(1, n_mats);
    fwderr_cosm_tay = zeros(1, n_mats);
    fwderr_cosm_pol = zeros(1, n_mats);

    condest1u = zeros(1, n_mats);

    rng(1);

    main_loop = tic; % record the time consumption
    
    for k = ids_min : ids_max

        fprintf('\n* Matrix id: %d\n', k);

        
        if k <= n_mats/2
            i_mats = k;
            A = testmats(k,n);
        elseif k <= n_mats
            i_mats = k - n_mats/2;
            A = testmats(i_mats,n)*1i;
        end
        
        mp.Digits(digits_work);
        
        fprintf('Processing cosm...\n');
        X_cosm = cosmahr(A);
        
        fprintf('Processing cosmadv...\n');
        X_cosm_adv = cosm(mp(A)); % cosm from Advanpix tool box takes 'mp' input
        
        fprintf('Processing cosmexp...\n');
        X_cosm_exp = cosmexp_double(A);
        
        fprintf('Processing cosmmp...\n');
        X_cosm_mp = cosm_double(A,'algorithm','transfree');

        fprintf('Processing cosmmp_s...\n');
        X_cosm_mp_s = cosm_double(A,'algorithm','realschur');
        
        fprintf('Processing cosmtay...\n');
        X_cosm_tay = cosmtay(A);
        
        fprintf('Processing cosmpol...\n');
        X_cosm_pol = cosmpol(A);
        
        mp.Digits(digits_ref);
        fprintf('Computing the reference solution...\n');
        exact = cosm_mp(mp(A),'algorithm','transfree');
        exactnorm1 = norm(exact, 1);
        
        fwderr_cosm(k) = double(norm(X_cosm - exact, 1) / exactnorm1);
        fwderr_cosm_adv(k) = double(norm(X_cosm_adv - exact, 1) / exactnorm1);
        fwderr_cosm_exp(k) = double(norm(X_cosm_exp - exact, 1) / exactnorm1);
        fwderr_cosm_mp(k) = double(norm(X_cosm_mp - exact, 1) / exactnorm1);
        fwderr_cosm_mp_s(k) = double(norm(X_cosm_mp_s - exact, 1) / exactnorm1);
        fwderr_cosm_tay(k) = double(norm(X_cosm_tay - exact, 1) / exactnorm1);
        fwderr_cosm_pol(k) = double(norm(X_cosm_pol - exact, 1) / exactnorm1);
        
        fprintf('Processing condest...\n');
        if (compute_condest)
            condest1u(k) = funm_condest1(double(A),@cosmmct) * eps/2;
        end

    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    
    dataname = sprintf('data/cos_fwderr_double.mat');
    save(dataname, 'condest1u', 'fwderr_cosm', 'fwderr_cosm_adv','fwderr_cosm_exp',...
            'fwderr_cosm_mp', 'fwderr_cosm_mp_s', 'fwderr_cosm_tay', 'fwderr_cosm_pol');
else
    dataname = sprintf('data/cos_fwderr_double.mat');
    load(dataname)
end


%% Plots
print_legend = true;
print_yticks = true;

Tcolors = [color_cosm_tay; color_cosm_pol; color_cosm_mp; color_cosm_exp;...
           color_cosm; color_cosm_adv; color_cosm_mp_s];
Tstyles = {ls_cosm_tay; ls_cosm_pol; ls_cosm_mp; ls_cosm_exp;...
           ls_cosm; ls_cosm_adv; ls_cosm_mp_s};
Tmarkers = {marker_cosm_tay; marker_cosm_pol; marker_cosm_mp; marker_cosm_exp;...
           marker_cosm; marker_cosm_adv; marker_cosm_mp_s};
numalg = size(Tmarkers,1);
Tstyles_perfprof = cell(numalg,1);
for k=1:numalg,  Tstyles_perfprof(k) = append(Tstyles(k),Tmarkers(k)); end
        
legend_perm = [1:numalg];

[~, perm] = sort(condest1u, 'descend'); % sort the matrices according to condu

% Plot performance profile and error

plot_cos_perfprof_fwderr_double(fwderr_cosm, fwderr_cosm_adv, fwderr_cosm_exp,... 
    fwderr_cosm_mp, fwderr_cosm_mp_s, fwderr_cosm_tay, fwderr_cosm_pol,...
    digits_work, Tcolors, Tstyles_perfprof, print_legend, legend_perm,...
    png_out, tikz_out, theta_lim);

plot_cos_fwderr_double(condest1u, perm, fwderr_cosm, fwderr_cosm_adv,...
    fwderr_cosm_exp, fwderr_cosm_mp, fwderr_cosm_mp_s, fwderr_cosm_tay,...
    fwderr_cosm_pol, digits_work, ymin, ymax, Tcolors, color_cond, Tmarkers,...
    ls_cond, msize, lw, lw_cond, print_yticks, print_legend, legend_perm, png_out, tikz_out);