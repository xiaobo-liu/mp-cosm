function [A, n_mats] = testmats_var(k, n)
%TESTMATS_VAR   104 test matrices from TESTMATS whose size is variable.
%Matrices causing overflow problems are excluded. Matrices formed by these
%matrices multipying the imaginary unit are also tested. In total, 104
%matrices are used in the experiments for comparing the execution time.

    n_mats = 104;
    if nargin < 1
        A = 104; % total number of matrices
        return;
    end
    if k <= n_mats/2
        A = full(testmats_var_r(k, n));
    elseif k <= n_mats
        A = full(testmats_var_r(k - n_mats/2,n))*1i;
    end
    
    function [A, n_mats] = testmats_var_r(k, n)
    n_mats = 52;  
    switch k
    case 01, sA='core/gfpp'; A=anymatrix(sA,n); 
    case 02, sA='core/hessfull01'; A=anymatrix(sA,n);
    case 03, sA='core/nilpotri'; A=anymatrix(sA,n); %sparse
    case 04, sA='core/rschur'; A=anymatrix(sA,n);
    case 05, sA='core/soules'; A=anymatrix(sA,n);
    case 06, sA='core/vand'; A=anymatrix(sA,n);
    case 07, sA='gallery/chebvand'; A=anymatrix(sA,n);
    case 08, sA='gallery/chow'; A=anymatrix(sA,n);
    case 09, sA='gallery/clement'; A=anymatrix(sA,n);
    case 10, sA='gallery/condex'; A=anymatrix(sA,n,3); %lower triangular
    case 11, sA='gallery/cycol'; A=anymatrix(sA,n);
    case 12, sA='gallery/dorr'; A=anymatrix(sA,n); %sparse
    case 13, sA='gallery/dramadah'; A=anymatrix(sA,n,1);
    case 14, sA='gallery/dramadah'; A=anymatrix(sA,n,2);
    case 15, sA='gallery/dramadah'; A=anymatrix(sA,n,3);
    case 16, sA='gallery/forsythe'; A=anymatrix(sA,n);
    case 17, sA='gallery/frank'; A=anymatrix(sA,n);
    case 18, sA='gallery/gearmat'; A=anymatrix(sA,n);
    case 19, sA='gallery/grcar'; A=anymatrix(sA,n);
    case 20, sA='gallery/hanowa'; A=anymatrix(sA,n);
    case 21, sA='gallery/jordbloc'; A=anymatrix(sA,n);
    case 22, sA='gallery/kahan'; A=anymatrix(sA,n);
    case 23, sA='gallery/leslie'; A=anymatrix(sA,n);
    case 24, sA='gallery/lesp'; A=anymatrix(sA,n);
    case 25, sA='gallery/lotkin'; A=anymatrix(sA,n);
    case 26, sA='gallery/normaldata'; A=anymatrix(sA,n,10);
    case 27, sA='gallery/orthog'; A=anymatrix(sA,n,4);
    case 28, sA='gallery/orthog'; A=anymatrix(sA,n,-2);
    case 29, sA='gallery/parter'; A=anymatrix(sA,n);
    case 30, sA='gallery/randcolu'; A=anymatrix(sA,n);
    case 31, sA='gallery/randhess'; A=anymatrix(sA,n);
    case 32, sA='gallery/rando'; A=anymatrix(sA,n,1);
    case 33, sA='gallery/rando'; A=anymatrix(sA,n,2);
    case 34, sA='gallery/rando'; A=anymatrix(sA,n,3);
    case 35, sA='gallery/rando'; A=10*triu(anymatrix(sA,n),1); % nilpotent, triangular.
    case 36, sA='gallery/randsvd'; A=anymatrix(sA,n,1);
    case 37, sA='gallery/randsvd'; A=anymatrix(sA,n,2);
    case 38, sA='gallery/randsvd'; A=anymatrix(sA,n,3);
    case 39, sA='gallery/randsvd'; A=anymatrix(sA,n,4);
    case 40, sA='gallery/randsvd'; A=anymatrix(sA,n,5);
    case 41, sA='gallery/redheff'; A=anymatrix(sA,n);
    case 42, sA='gallery/riemann'; A=anymatrix(sA,n);
    case 43, sA='gallery/smoke'; A=anymatrix(sA,n);
    case 44, sA='gallery/smoke'; A=anymatrix(sA,n,1);
    case 45, sA='gallery/toeppen'; A=anymatrix(sA,n); %sparse
    case 46, sA='gallery/triw'; A=anymatrix(sA,n,-1);
    case 47, sA='gallery/triw'; A=anymatrix(sA,n,-2);
    case 48, sA='gallery/uniformdata'; A=anymatrix(sA,n,1e3);
    case 49, A = zeros(n); A(n+1:n+1:n^2) = 1:n-1; % log of Cholesky factor of Pascal matrix. See \cite{edst04}.
    case 50, A = anymatrix('gallery/triw',n,1);  m = (n-1)/2; % \cite[p.~9, Ex II]{pang85}
             A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1); % Fasi's interpretation for arbitrary n.
             for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end % N = 31 corresponds to the matrix in above ref.
    case 51, A = anymatrix('gallery/triw',n,1,1); % \cite[p.~10, Ex III]{pang85}
             A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2);
    case 52, alpha = 1; beta = 1;  % % \cite{kuda10}, no values are given in the paper.
             A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
             A(1,2) = beta; A(n,n-1) = beta; 
    end
    end
end