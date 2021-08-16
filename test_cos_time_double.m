% this code tests in double precision the execution times for computing the
% cosine the algorithms on test set V.

format compact

addpath('external');
addpath('include');
addpath('data');
warning('off');


%% the following commented when run all the tests at one time

initialize_tests

% n = 100; % change to reset the size of the test matrices

reproduce_results = false; % reproduce the test results


%% set output parameters and run the tests

theta_lim = 30; % the range of x-axis in the performance profile generated

% the range of y-axis in the time plots
if n==16,  y_lim = 0.02; end
if n==100, y_lim = 0.06;  end

png_out = true;

tikz_out = false;

ids_min = 1;
[~, ids_max] = testmats_var();

digits_work = 16; 

if reproduce_results

    n_mats = ids_max - ids_min + 1;

    time_cosm = zeros(1, n_mats);
    time_cosm_exp = zeros(1, n_mats);
    time_cosm_mp = zeros(1, n_mats);
    time_cosm_mp_s = zeros(1, n_mats);
    time_cosm_tay = zeros(1, n_mats);
    time_cosm_pol = zeros(1, n_mats);
    time_cosm_adv = zeros(1, n_mats);

    rng(1);

    main_loop = tic; % record the time consumption
    
    for i = ids_min : ids_max

        fprintf('\n* Matrix id: %d\n', i);

        A = testmats_var(i,n);

        mp.Digits(digits_work);
        
        fprintf('Processing cosmadv...\n');
        [X_cosm_adv, time] = smart_timer(@cosm, mp(A));
        time_cosm_adv(i) = time;
        
        fprintf('Processing cosm...\n');
        [X_cosm, time] = smart_timer(@cosmahr, A);
        time_cosm(i) = time;
         
        fprintf('Processing cosmexp...\n');
        [X_cosm_exp, time] = smart_timer(@cosmexp_double, A);
        time_cosm_exp(i) = time;
        
        fprintf('Processing cosmmp...\n');
        [X_cosm_mp, time] = smart_timer(@cosm_double, A,'algorithm','transfree');
        time_cosm_mp(i) = time;
        
        fprintf('Processing cosmmp_s...\n');
        [X_cosm_mp_s, time] = smart_timer(@cosm_double, A,'algorithm','complexschur');
        time_cosm_mp_s(i) = time;
        
        fprintf('Processing cosmtay...\n');
        [X_cosm_tay, time] = smart_timer(@cosmtay, A);
        time_cosm_tay(i) = time;
        
        fprintf('Processing cosmpol...\n');
        [X_cosm_pol, time] = smart_timer(@cosmpol, A);
        time_cosm_pol(i) = time;
        
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    
    dataname = sprintf('data/cos_time_double_n%d.mat', n);
    save(dataname, 'time_cosm', 'time_cosm_exp', 'time_cosm_mp', 'time_cosm_mp_s',...
         'time_cosm_tay', 'time_cosm_pol', 'time_cosm_adv');
else % loading existing data
    dataname = sprintf('data/cos_time_double_n%d.mat', n);
    load(dataname)
end

%% Plots
print_legend = true;
print_yticks = true;

Tcolors = [color_cosm_pol; color_cosm_tay; color_cosm; color_cosm_exp;...
           color_cosm_mp; color_cosm_mp_s; color_cosm_adv];
Tstyles = {ls_cosm_pol; ls_cosm_tay; ls_cosm; ls_cosm_exp;...
           ls_cosm_mp; ls_cosm_mp_s; ls_cosm_adv};
Tmarkers = {marker_cosm_pol; marker_cosm_tay; marker_cosm; marker_cosm_exp;...
            marker_cosm_mp; marker_cosm_mp_s; marker_cosm_adv};
        
numalg = size(Tmarkers,1);
T_styles_appended = cell(numalg,1);
for i=1:numalg,  T_styles_appended(i) = append(Tstyles(i),Tmarkers(i)); end

legend_perm = (1:numalg);
legend_perm2 = (1:numalg-1); % time of cosm_adv is not plotted (too slow)

[~, perm] = sort(time_cosm_mp_s(), 'ascend');

% Plot performance profile
plot_cos_perfprof_time_double(time_cosm, time_cosm_exp, time_cosm_mp, time_cosm_mp_s,...
                 time_cosm_tay, time_cosm_pol, time_cosm_adv, Tcolors, T_styles_appended,... 
                 print_legend, legend_perm, png_out, tikz_out, theta_lim, n);
             
% Plot the execution time 
plot_cos_time_double(time_cosm, time_cosm_exp, time_cosm_mp, time_cosm_mp_s,...
                 time_cosm_tay, time_cosm_pol, perm, Tcolors,...
                 Tmarkers, msize, lw, print_yticks, print_legend,...
                 legend_perm2, png_out, tikz_out, y_lim, n);