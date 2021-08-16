function [A, n_mats] = testmats(k, n)
%TESTMATS_AM   99 test matrices from ANYMATRIX matrix collection (including
%   matrices from the matrices-mp-cosm collection) that are nonnormal, 
%   square and not causing overflow problems in computing cos(A).

    n_mats = 99;
    if nargin < 1
        A = n_mats; % total number of matrices
        return;
    end
    if nargin < 2
        n = 16;
    end

    switch k
    case 01, sA='core/edelman27'; A=anymatrix(sA); % 27 by 27
    case 02, sA='core/gfpp'; A=anymatrix(sA,n); 
    case 03, sA='core/hessfull01'; A=anymatrix(sA,n);
    case 04, sA='core/nilpotri'; A=anymatrix(sA,n); % sparse
    case 05, sA='core/rschur'; A=anymatrix(sA,n);
    case 06, sA='core/soules'; A=anymatrix(sA,n);
    case 07, sA='core/vand'; A=anymatrix(sA,n);
    case 08, sA='gallery/binomial'; A=anymatrix(sA,n);
    case 09, sA='gallery/chebspec'; A=anymatrix(sA,n);
    case 10, sA='gallery/chebvand'; A=anymatrix(sA,n);
    case 11, sA='gallery/chow'; A=anymatrix(sA,n);
    case 12, sA='gallery/clement'; A=anymatrix(sA,n);
    case 13, sA='gallery/condex'; A=anymatrix(sA,4,1); % 4 by 4
    case 14, sA='gallery/condex'; A=anymatrix(sA,3,2); % 3 by 3
    case 15, sA='gallery/condex'; A=anymatrix(sA,n,3); % lower triangular
    case 16, sA='gallery/cycol'; A=anymatrix(sA,n);
    case 17, sA='gallery/dorr'; A=anymatrix(sA,n); % sparse
    case 18, sA='gallery/dramadah'; A=anymatrix(sA,n,1);
    case 19, sA='gallery/dramadah'; A=anymatrix(sA,n,2);
    case 20, sA='gallery/dramadah'; A=anymatrix(sA,n,3);
    case 21, sA='gallery/forsythe'; A=anymatrix(sA,n);
    case 22, sA='gallery/forsythe'; A=anymatrix(sA,10,1e-10,0); % \cite[Test 4]{ward77}.
    case 23, sA='gallery/frank'; A=anymatrix(sA,n);
    case 24, sA='gallery/gearmat'; A=anymatrix(sA,n);
    case 25, sA='gallery/grcar'; A=anymatrix(sA,n);
    case 26, sA='gallery/hanowa'; A=anymatrix(sA,n);
    case 27, sA='gallery/invhess'; A=anymatrix(sA,n);
    case 28, sA='gallery/invol'; A = anymatrix(sA,8)*8*pi; % \cite[exmp.~3]{hism03}
    case 29, sA='gallery/invol'; A = anymatrix(sA,13); A = triu(schur(A,'complex'),1); % nilpotent, triangular
    case 30, sA='gallery/jordbloc'; A=anymatrix(sA,n);
    case 31, sA='gallery/kahan'; A=anymatrix(sA,n);
    case 32, sA='gallery/leslie'; A=anymatrix(sA,n);
    case 33, sA='gallery/lesp'; A=anymatrix(sA,n);
    case 34, sA='gallery/lotkin'; A=anymatrix(sA,n);
    case 35, sA='gallery/neumann'; A=anymatrix(sA,ceil(sqrt(n))^2); % sparse, 25 by 25 for n=20;
    case 36, sA='gallery/normaldata'; A=anymatrix(sA,n,10);
    case 37, sA='gallery/orthog'; A=anymatrix(sA,n,4);
    case 38, sA='gallery/orthog'; A=anymatrix(sA,n,-2);
    case 39, sA='gallery/parter'; A=anymatrix(sA,n);
    case 40, sA='gallery/randcolu'; A=anymatrix(sA,n);
    case 41, sA='gallery/randhess'; A=anymatrix(sA,n);
    case 42, sA='gallery/rando'; A=anymatrix(sA,n,1);
    case 43, sA='gallery/rando'; A=anymatrix(sA,n,2);
    case 44, sA='gallery/rando'; A=anymatrix(sA,n,3);
    case 45, sA='gallery/rando'; A=10*triu(anymatrix(sA,n),1); % nilpotent, triangular.
    case 46, sA='gallery/randsvd'; A=anymatrix(sA,n,1);
    case 47, sA='gallery/randsvd'; A=anymatrix(sA,n,2);
    case 48, sA='gallery/randsvd'; A=anymatrix(sA,n,3);
    case 49, sA='gallery/randsvd'; A=anymatrix(sA,n,4);
    case 50, sA='gallery/randsvd'; A=anymatrix(sA,n,5);
    case 51, sA='gallery/randsvd'; A=anymatrix(sA,8,1e14); % \cite{fahi19}
    case 52, sA='gallery/redheff'; A=anymatrix(sA,n);
    case 53, sA='gallery/riemann'; A=anymatrix(sA,n);
    case 54, sA='gallery/sampling'; A=anymatrix(sA,n);
    case 55, sA='gallery/smoke'; A=anymatrix(sA,n);
    case 56, sA='gallery/smoke'; A=anymatrix(sA,n,1);
    case 57, sA='gallery/toeppen'; A=anymatrix(sA,n); % sparse
    case 58, sA='gallery/triw'; A=anymatrix(sA,n,-1);
    case 59, sA='gallery/triw'; A=anymatrix(sA,n,-2);
    case 60, sA='gallery/uniformdata'; A=anymatrix(sA,n,1e3);
    case 61, sA='gallery/wilk'; A=anymatrix(sA,3); % 3 by 3
    case 62, sA='gallery/wilk'; A=anymatrix(sA,4); % 4 by 3
    case 63, sA='matlab/pascal'; A=anymatrix(sA,n,1);
    case 64, sA='matlab/pascal'; A=anymatrix(sA,n,2);
    case 65, sA='nessie/spl0708b'; A=anymatrix(sA); % sparse, 41 by 41
    % matrices from matrices-mp-cosm collection
    case 66, sA='mpcosm/almo2r1'; A=anymatrix(sA);
    case 67, sA='mpcosm/almo2r2'; A=anymatrix(sA);
    case 68, sA='mpcosm/almo4'; A=anymatrix(sA);
    case 69, sA='mpcosm/dahi03'; A=anymatrix(sA);
    case 70, sA='mpcosm/dipa00'; A=anymatrix(sA);
    case 71, sA='mpcosm/edst04'; A=anymatrix(sA);
    case 72, sA='mpcosm/eigt7'; A=anymatrix(sA);
    case 73, sA='mpcosm/fahi19'; A=anymatrix(sA);
    case 74, sA='mpcosm/fasi7'; A=anymatrix(sA);
    case 75, sA='mpcosm/hism03'; A=anymatrix(sA);
    case 76, sA='mpcosm/jemc05r1'; A=anymatrix(sA);
    case 77, sA='mpcosm/jemc05r2'; A=anymatrix(sA);
    case 78, sA='mpcosm/kase99'; A=anymatrix(sA);
    case 79, sA='mpcosm/kela89'; A=anymatrix(sA);
    case 80, sA='mpcosm/kela98r1'; A=anymatrix(sA);
    case 81, sA='mpcosm/kela98r2'; A=anymatrix(sA);
    case 82, sA='mpcosm/kela98r3'; A=anymatrix(sA);
    case 83, sA='mpcosm/kela98r4'; A=anymatrix(sA);
    case 84, sA='mpcosm/kuda10'; A=anymatrix(sA);
    case 85, sA='mpcosm/lara17r1'; A=anymatrix(sA);
    case 86, sA='mpcosm/lara17r2'; A=anymatrix(sA);
    case 87, sA='mpcosm/lara17r3'; A=anymatrix(sA);
    case 88, sA='mpcosm/lara17r4'; A=anymatrix(sA);
    case 89, sA='mpcosm/lara17r5'; A=anymatrix(sA);
    case 90, sA='mpcosm/lara17r6'; A=anymatrix(sA);
    case 91, sA='mpcosm/mopa03r1'; A=anymatrix(sA);
    case 92, sA='mpcosm/mopa03r2'; A=anymatrix(sA);
    case 93, sA='mpcosm/naha95'; A=anymatrix(sA);
    case 94, sA='mpcosm/pang85r1'; A=anymatrix(sA);
    case 95, sA='mpcosm/pang85r2'; A=anymatrix(sA);
    case 96, sA='mpcosm/pang85r3'; A=anymatrix(sA);
    case 97, sA='mpcosm/trem05'; A=anymatrix(sA);
    case 98, sA='mpcosm/ward77r1'; A=anymatrix(sA);
    case 99, sA='mpcosm/ward77r2'; A=anymatrix(sA);          
    end
    A = full(A);
end