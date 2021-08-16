% this code tests the execution times for computing the cosine in multiprecision

format compact

addpath('data');
addpath('external');
addpath('include');
warning('off');

%% the following commented when run all the tests at one time

initialize_tests

% n = 100; % change to reset size of the test matrices

reproduce_results = false; % reproduce the test results


%% set output parameters and run the tests

% precisions
prec_vec = [256];

theta_lim = 10; % the range of x-axis in the performance profile generated

% the range of y-axis in the time plots
if n==16,  y_lim = 1; end
if n==100, y_lim = 100;  end

png_out = true;
tikz_out = false;

ids_min = 1;
[~, ids_max] = testmats_var();

n_mats = ids_max - ids_min + 1;

n_prec = length(prec_vec);


if reproduce_results
    
    % initializations
    time_cosm_adv = zeros(n_prec, n_mats);
    time_cosm_exp = zeros(n_prec, n_mats);
    time_cosm_mp = zeros(n_prec, n_mats);
    time_cosm_mp_s = zeros(n_prec, n_mats);

    rng(1);

    main_loop = tic; % record the time consumption
    
    for j = 1:n_prec
        
        prec_work = prec_vec(j); % digits of the working precision
 
        mp.Digits(prec_work);
        epsilon = mp('eps');
        
        for i = ids_min : ids_max
            fprintf('\n* Matrix id: %d\n', i);
            
            A = testmats_var(i,n);
            
            A = mp(A);
        
            fprintf('Processing cosm_adv...\n');
            [X_cosm_adv, time] = smart_timer(@cosm, A);
            time_cosm_adv(j, i) = time;
        
            fprintf('Processing cosm_exp...\n');
            [X_cosm_exp, time] = smart_timer(@cosmexp_mp, A);
            time_cosm_exp(j, i) = time;
        
            fprintf('Processing cosm_mp...\n');
            [X_cosm_mp, time] = smart_timer(@cosm_mp, A,'algorithm','transfree');
            time_cosm_mp(j, i) = time;
     
            fprintf('Processing cosm_mp_s...\n');
            [X_cosm_mp_s, time] = smart_timer(@cosm_mp, A,'algorithm','complexschur');
            time_cosm_mp_s(j, i) = time;
            
        end
        dataname = sprintf('data/cos_time_mp_%d_n%d.mat', prec_vec(j), n);
        save(dataname, 'time_cosm_adv', 'time_cosm_exp', 'time_cosm_mp', 'time_cosm_mp_s');
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);    
end

%% Plots
print_legend = true;
print_yticks = true;

Tcolors = [color_cosm_mp; color_cosm_exp; color_cosm_mp_s; color_cosm_adv];
Tstyles = {ls_cosm_mp; ls_cosm_exp; ls_cosm_mp_s; ls_cosm_adv};
Tmarkers = {marker_cosm_mp; marker_cosm_exp; marker_cosm_mp_s; marker_cosm_adv};
numalg = size(Tmarkers,1);
    
T_styles_appended = cell(numalg,1);
for i=1:numalg,  T_styles_appended(i) = append(Tstyles(i),Tmarkers(i)); end
        
legend_perm = [1:numalg];
    
for j = 1:n_prec
    dataname = sprintf('data/cos_time_mp_%d_n%d.mat', prec_vec(j), n);
    load(dataname)
    [~, perm] = sort(time_cosm_adv(j,:), 'ascend'); % sorting in the time plots by the times of cosm_adv
    
    time_cosmadv_id1 = time_cosm_adv(j,:);
    time_cosm_exp_id1 = time_cosm_exp(j,:);
    time_cosm_mp_id1 = time_cosm_mp(j,:);
    time_cosm_mp_s_id1 = time_cosm_mp_s(j,:);

    % Plot performance profile
    plot_cos_perfprof_time_mp(time_cosmadv_id1, time_cosm_exp_id1, time_cosm_mp_id1,...
                    time_cosm_mp_s_id1, prec_vec(j), Tcolors, T_styles_appended,...
                    print_legend, legend_perm, png_out, tikz_out, theta_lim, n);
    % Plot the execution time 
    plot_cos_time_mp(time_cosmadv_id1, time_cosm_exp_id1, time_cosm_mp_id1,...
                    time_cosm_mp_s_id1, perm, Tcolors, Tmarkers, msize,...
                            lw, print_yticks, print_legend, legend_perm,...
                            png_out, tikz_out, y_lim, n, prec_vec(j))
end