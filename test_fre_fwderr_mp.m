% this code tests the error in the computed Frechet derivative in 
% multiprecision

format compact

addpath('data');
addpath('external');
addpath('include');
warning('off');

%% the following codes should be commented when run all the tests at one time

initialize_tests

produce_results = false; % reproduce the test results?


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

compute_condest = true;

ids_min = 1;
[~, ids_max] = testmats();
ids_max = ids_max * 2; % size doubled by multiplying the imaginary unit i

n_mats = ids_max - ids_min + 1;

n_prec = length(prec_vec);


if produce_results
    
    % digits used only in initializations
    mp.Digits(max(prec_vec));
    
    fwderr_fre_blk = zeros(n_prec, n_mats, 'mp');
    fwderr_fre_mp = zeros(n_prec, n_mats, 'mp');
    fwderr_fre_mp_s = zeros(n_prec, n_mats, 'mp');
    
    condest1 = zeros(1, ids_max, 'mp');
    condest1u = zeros(n_prec, n_mats, 'mp');
   
    rng(1);
    
    main_loop = tic; % record the time consumption
    
    for j = 1:n_prec
        
        prec_work = prec_vec(j); % digits of the working precision
        fprintf('Digits of working precision is %.4d...\n', prec_work);
        
        prec_ref = floor(prec_work * 2);
        
        mp.Digits(prec_work);
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
            
            E = randn(size(A,1), class(A));
            E = E / norm(E, 1);
       
            fprintf('Processing fre_blk...\n');
            X_fre_blk = cosm_frechet_blk(A, E);
        
            fprintf('Processing fre_mp...\n');
            [~, X_fre_mp] = cosm_frechet_mp(A, E, 'algorithm','transfree');
     
            fprintf('Processing fre_mp_s...\n');
            [~, X_fre_mp_s] = cosm_frechet_mp(A, E, 'algorithm','realschur');
            
            % Compute reference solution.
            mp.Digits(prec_ref);
            fprintf('Computing the reference solution...\n');
            [~, exact] = cosm_frechet_mp(A, E, 'algorithm', 'transfree',...
                'precision', prec_ref, 'maxdegree', 2500);
            exactnorm1 = norm(exact, 1);        
    
            fwderr_fre_blk(j, k) = norm(X_fre_blk - exact, 1) / exactnorm1;
            fwderr_fre_mp(j, k) = norm(X_fre_mp - exact, 1) / exactnorm1;
            fwderr_fre_mp_s(j, k) = norm(X_fre_mp_s - exact, 1) / exactnorm1;
            
             if (j == 1) % compute cond only once and multiply it by different eps 
                 old_digits = mp.Digits();
                 mp.Digits(34); % storing the values of cond in mp
                 fprintf('Processing condest...\n');
                 condest1(k) = fd_condest(@cosm_frechet_mp, A, E);
                 mp.Digits(old_digits);
             end
             
             mp.Digits(prec_work);
             condest1u(j, k) = condest1(k) * epsilon;
            
        end
 
        dataname = sprintf('data/fre_fwderr_mp_%d.mat', prec_vec(j));
        save(dataname, 'condest1u', 'fwderr_fre_blk', 'fwderr_fre_mp',...
        'fwderr_fre_mp_s');
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60); 
end

%% Plots
print_legend = true;
print_yticks = true;
Tcolors = [color_fre_mp; color_fre_blk; color_fre_mp_s];
Tstyles = {ls_fre_mp; ls_fre_blk; ls_fre_mp_s};
Tmarkers = {marker_fre_mp; marker_fre_blk; marker_fre_mp_s};

numalg = size(Tmarkers,1);
Tstyles_perfprof = cell(numalg,1);
for k=1:numalg,  Tstyles_perfprof(k) = append(Tstyles(k),Tmarkers(k)); end
        
legend_perm = [1:numalg];
    
for j = 1:n_prec
    dataname = sprintf('data/fre_fwderr_mp_%d.mat', prec_vec(j));
    load(dataname)
    [~, perm] = sort(condest1u(j,:), 'descend');
   
    % Plot performance profile    
    plot_fre_perfprof_fwderr_mp(fwderr_fre_blk(j,:), fwderr_fre_mp(j,:),...
        fwderr_fre_mp_s(j,:), prec_vec(j), Tcolors, Tstyles_perfprof,...
        print_legend, legend_perm, png_out, tikz_out, theta_lim);

    plot_fre_fwderr_mp(condest1u(j,:), perm, fwderr_fre_blk(j,:),...
        fwderr_fre_mp(j,:), fwderr_fre_mp_s(j,:), prec_vec(j), ...
        ymins(j), ymaxs(j), Tcolors, color_cond, Tmarkers, ls_cond,...
        msize, lw, lw_cond, print_yticks, print_legend, legend_perm, png_out, tikz_out);
                
end