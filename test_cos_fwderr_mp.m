% this code tests the error of the computed cosine in multiprecision

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

% precisions
prec_vec = [256, 1024];
ymins = [1e-258, mp('1e-1026')];
ymaxs = [1e-242, mp('1e-1010')];

% a quicker test
% prec_vec = 64; ymins = 1e-66; ymaxs = 1e-50;
% another one
% prec_vec = 256; ymins = 1e-258; ymaxs = 1e-242;


compute_condest = true;

ids_min = 1;
[~, ids_max] = testmats();
ids_max = ids_max * 2; % size doubled by multiplying the imaginary unit i

n_mats = ids_max - ids_min + 1;

n_prec = length(prec_vec);

if reproduce_results
    
    % digits used only in initializations
    mp.Digits(max(prec_vec));
    
    fwderr_cosm_adv = zeros(n_prec, n_mats, 'mp');
    fwderr_cosm_exp = zeros(n_prec, n_mats, 'mp');
    fwderr_cosm_mp = zeros(n_prec, n_mats, 'mp');
    fwderr_cosm_mp_s = zeros(n_prec, n_mats, 'mp');
    
    condest1 = zeros(1, ids_max, 'mp');
    condest1u = zeros(n_prec, n_mats, 'mp');
   
    rng(1);
    
    main_loop = tic; % record the time consumption
    
    for j = 1:n_prec
        
        digits_work = prec_vec(j); % digits of the working precision
        fprintf('Digits of working precision is %.4d...\n', digits_work);
        
        digits_ref = floor(digits_work * 2);
        
        mp.Digits(digits_work);
        epsilon = mp('eps');
        
        for k = ids_min : ids_max
            fprintf('\n* Matrix id: %d\n', k);
            
            if k <= n_mats/2
                i_mats = k;
                A = mp(testmats(i_mats,n));
            elseif k <= n_mats
                i_mats = k - n_mats/2;
                A = mp(testmats(i_mats,n)*1i);
            end
       
            fprintf('Processing cosm_adv...\n');
            X_cosm_adv = cosm(A);
        
            fprintf('Processing cosm_exp...\n');
            X_cosm_exp = cosmexp_mp(A);
        
            fprintf('Processing cosm_mp...\n');
            X_cosmtrf = cosm_mp(A,'algorithm','transfree');
     
            fprintf('Processing cosm_mp_s...\n');
            X_cosm_mp_s = cosm_mp(A,'algorithm','realschur');
            
            % Compute reference solution.
            mp.Digits(digits_ref);
            fprintf('Computing the reference solution...\n');
            exact = cosm_mp(mp(A),'algorithm','transfree','maxdegree',2500);
            exactnorm1 = norm(exact, 1);        
    
            fwderr_cosm_adv(j, k) = norm(X_cosm_adv - exact, 1) / exactnorm1;
            fwderr_cosm_exp(j, k) = norm(X_cosm_exp - exact, 1) / exactnorm1;
            fwderr_cosm_mp(j, k) = norm(X_cosmtrf - exact, 1) / exactnorm1;
            fwderr_cosm_mp_s(j, k) = norm(X_cosm_mp_s - exact, 1) / exactnorm1;
            
             if (j == 1) % compute cond only once and multiply it by different eps 
                 old_digits = mp.Digits();
                 mp.Digits(34); % storing the values of cond in mp
                 fprintf('Processing condest...\n');
                 condest1(k) = mp(funm_condest1(double(A), @cosmmct));
                 mp.Digits(old_digits);
             end
             
             mp.Digits(digits_work);
             condest1u(j, k) = condest1(k) * epsilon;
            
        end
        dataname = sprintf('data/cos_fwderr_mp_%d.mat', prec_vec(j));
        save(dataname, 'condest1u', 'fwderr_cosm_adv', 'fwderr_cosm_exp',...
        'fwderr_cosm_mp', 'fwderr_cosm_mp_s');
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60); 
end

%% Plots
print_legend = true;
print_yticks = true;
Tcolors = [color_cosm_mp; color_cosm_exp; color_cosm_adv; color_cosm_mp_s];
Tstyles = {ls_cosm_mp; ls_cosm_exp; ls_cosm_adv; ls_cosm_mp_s};
Tmarkers = {marker_cosm_mp; marker_cosm_exp; marker_cosm_adv; marker_cosm_mp_s};

numalg = size(Tmarkers,1);
Tstyles_perfprof = cell(numalg,1);
for k=1:numalg,  Tstyles_perfprof(k) = append(Tstyles(k),Tmarkers(k)); end
        
legend_perm = [1:numalg];
    
for j = 1:n_prec
    dataname = sprintf('data/cos_fwderr_mp_%d.mat', prec_vec(j));
    load(dataname)
    [~, perm] = sort(condest1u(j,:), 'descend');
   
    % Plot performance profile    
    plot_cos_perfprof_fwderr_mp(fwderr_cosm_adv(j,:), fwderr_cosm_exp(j,:),...
        fwderr_cosm_mp(j,:), fwderr_cosm_mp_s(j,:), prec_vec(j), Tcolors,... 
        Tstyles_perfprof, print_legend, legend_perm, png_out, tikz_out, theta_lim);

    plot_cos_fwderr_mp(condest1u(j,:), perm, fwderr_cosm_adv(j,:),...
        fwderr_cosm_exp(j,:), fwderr_cosm_mp(j,:), fwderr_cosm_mp_s(j,:),... 
        prec_vec(j), ymins(j), ymaxs(j), Tcolors, color_cond, Tmarkers,... 
        ls_cond, msize, lw, lw_cond, print_yticks, print_legend, legend_perm, png_out, tikz_out);
                
end