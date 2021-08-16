% this code tests the error in the computed Frechet derivative in 
% double precision

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

    fwderr_fre_blk = zeros(1, n_mats);
    fwderr_fre_exp = zeros(1, n_mats);
    fwderr_fre_mp = zeros(1, n_mats);
    fwderr_fre_mp_s = zeros(1, n_mats);
   

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
        
        E = randn(size(A,1), class(A));
        E = E / norm(E, 1);
        
        mp.Digits(digits_work);
        
        fprintf('Processing fre_blk...\n');
        X_fre_blk = cosm_frechet_blk(A, E);
        
        fprintf('Processing fre_exp...\n');
        X_fre_exp = cosm_frechet_exp(A, E);
        
        fprintf('Processing fre_mp...\n');
        [~, X_fre_mp] = cosm_frechet_mp(A, E, 'algorithm','transfree');

        fprintf('Processing fre_mp_s...\n');
        [~, X_fre_mp_s] = cosm_frechet_mp(A, E, 'algorithm','realschur');
        
        fprintf('Computing the reference solution...\n');
        [~, exact] = cosm_frechet_mp(A, E, 'algorithm', 'transfree', ...
            'precision', digits_ref, 'maxdegree', 2500);
        exactnorm1 = norm(exact, 1);
        
        fwderr_fre_blk(k)  = double(norm(X_fre_blk - exact, 1) / exactnorm1);
        fwderr_fre_exp(k)  = double(norm(X_fre_exp - exact, 1) / exactnorm1);
        fwderr_fre_mp(k)   = double(norm(X_fre_mp - exact, 1) / exactnorm1);
        fwderr_fre_mp_s(k) = double(norm(X_fre_mp_s - exact, 1) / exactnorm1);
        
        fprintf('Processing condest...\n');
        if (compute_condest)
            condest1u(k) = fd_condest(@cosm_frechet_mp, A, E) * eps/2;            
        end

    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    
    dataname = sprintf('data/fre_fwderr_double.mat');
    save(dataname, 'condest1u', 'fwderr_fre_blk', 'fwderr_fre_exp',...
            'fwderr_fre_mp', 'fwderr_fre_mp_s');
else
    dataname = sprintf('data/fre_fwderr_double.mat');
    load(dataname)
end


%% Plots
print_legend = true;
print_yticks = true;

Tcolors = [color_fre_mp; color_fre_blk; color_fre_mp_s; color_fre_exp];
Tstyles = {ls_fre_mp; ls_fre_blk; ls_fre_mp_s; ls_fre_exp};
Tmarkers = {marker_fre_mp; marker_fre_blk; marker_fre_mp_s; marker_fre_exp};
numalg = size(Tmarkers,1);
Tstyles_perfprof = cell(numalg,1);
for k=1:numalg,  Tstyles_perfprof(k) = append(Tstyles(k),Tmarkers(k)); end
        
legend_perm = [1:numalg];

[~, perm] = sort(condest1u, 'descend'); % sort the matrices according to condu

% Plot performance profile and error

plot_fre_perfprof_fwderr_double(fwderr_fre_blk, fwderr_fre_exp,... 
    fwderr_fre_mp, fwderr_fre_mp_s, digits_work, Tcolors, Tstyles_perfprof,...
    print_legend, legend_perm, png_out, tikz_out, theta_lim);

plot_fre_fwderr_double(condest1u, perm, fwderr_fre_blk, fwderr_fre_exp,...
    fwderr_fre_mp, fwderr_fre_mp_s, digits_work, ymin, ymax, Tcolors, ...
    color_cond, Tmarkers, ls_cond, msize, lw, lw_cond, print_yticks, print_legend,...
    legend_perm, png_out, tikz_out);