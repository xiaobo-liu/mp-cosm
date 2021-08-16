function [varargout] = cosm_mp(A, varargin)
%COSM_MP  Multiprecision algorithm for the matrix cosine.
%
%   This code requires the Advanpix Multiprecision Computing Toolbox (see
%   www.advanpix.com).
%
%   [C, s, m] = cosm_mp(A) computes the matrix cosine C of the square
%   matrix A. The output parameters s and m are the scaling parameter and 
%   the degree of the approximant.

%   [...] = cosm_mp(...,'precision',DIGITS) specifies the number of digits to be
%   used in the computation. Default is mp.Digits() if A is of class mp. The
%   computation is performed in single or double arithmetic if A is of class
%   single or double, respectively, and the precision is not specified.
%
%   [...] = cosm_mp(...,'epsilon',EPSILON) specifies the tolerance to be
%   used in evaluating the Taylor approximants. Default is machine epsilon 
%   of the precision of A if A is of class 'single' or 'double', mp.eps() 
%   if A is of class 'mp'.
%
%   [...] = cosm_mp(...,'maxdegree',MAXDEGREE) specifies the maximum degree of
%   the Taylor approximants. Default is 200.
%
%   [...] = cosm_mp(...,'preprocessing',PREPRCS) determines whether 
%   preprocessing is performed on the input matrix A: PREPRCS = True/False. 
%   Preprocessing is turned off by default.
%
%   [...] = cosm_mp(...,'algorithm', KIND), specifies whether the algorithm
%   should be performed on the full matrix (KIND='transfree'), on the upper
%   triangular factor of the complex Schur decomposition of the matrix (KIND =
%   'complexschur'), or on the upper quasi-triangular factor of the real Schur
%   decomposition of the matrix. If the input matrix has complex entries,
%   'complexschur' is used instead of 'realschur'. Default is 'transfree'.
%
%   Reference: 
%   [1] M. Fasi, mpexpm, https://github.com/mfasi/mpexpm/blob/master/include/expm_mp.m.

% Copyright (c) 2017-2021, Massimiliano Fasi, Xiaobo Liu
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
% EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Input Validation.
    p = inputParser;
    addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
    addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
    addParameter(p, 'maxdegree', 500, @(x)(x == round(x) && x > 0));
    addParameter(p, 'preprocessing', 0, @(x) islogical(x)); 
    addParameter(p, 'algorithm', 'transfree',...
                 @(x)(ischar(x) && strcmp(x, 'transfree')...
                      || strcmp(x, 'realschur') || strcmp(x, 'complexschur')));
  
    parse(p,A,varargin{:});

    A = p.Results.A;
    digits = p.Results.precision;
    epsilon = p.Results.epsilon;
    maxdegree = p.Results.maxdegree;
    pre_prcs = p.Results.preprocessing;
    
    [n1, n] = size(A);
    if n1 ~= n, error('The matrix ``A'' must be square.'); end

    if nargout > 3
        error('This function returns at most three values.');
    end
    
%% Determine whether the computation will be single, double or mp.
    current_digits = mp.Digits();
    if ~isempty(digits) % working precision specified
        mp.Digits(digits);
        A = mp(A);
    else
        if isa(A, 'double')
            digits = 16;
        elseif isa(A, 'single')
            digits = 8;
        else
            digits = mp.Digits();
        end
        mp.Digits(digits);
    end
    if isempty(epsilon)
        epsilon = myeps(class(A));
    end
    curr_epsilon = epsilon;
    
%% Use diagonalization if A is hermitian.
    if ishermitian(A)
        [V, d] = eig(A,'vector'); 
        if nargout >= 0
            C = V*(cos(d).*V');
            [varargout{1}] = (C + C') / 2; 
        end
        if nargout > 1
            [varargout{2}] = 0;
        end
        if nargout > 2
            [varargout{3}] = 0;
        end
        return 
    end
        
%% Form complex Schur form if required.
    compute_schur = false;
    recompute_diag_blocks = false;
    switch p.Results.algorithm
      case 'transfree'
        if matlab.internal.math.isschur(A)
            recompute_diag_blocks = true;
        else
            recompute_diag_blocks = false;
        end
      case 'realschur'
        if ~matlab.internal.math.isschur(A)
            compute_schur = true;
            schur_string = 'real';
        end
      case 'complexschur'
        if ~istriu(A)
            compute_schur = true;
            schur_string = 'complex';
        end
    end
    
    if compute_schur 
        % if A is complex the complex Schur form is returned and in this
        % case the second input `schur_string' is ignored by MATLAB.
        [Q, A] = schur(A, schur_string);
    end

%% preprocessing the input matrix if asked to.
    if pre_prcs
        mp_pi = mp('pi');
        mu_prcs = trace(A) / (n*mp_pi);
        q = [0, floor(mu_prcs), ceil(mu_prcs)];
        normi = zeros(1,3);
        for qi=1:3, normi(qi) = norm(A(1:n+1:n^2) - mp_pi*q(qi), 1); end
        [~, index] = min(normi);
        A(1:n+1:n^2) = A(1:n+1:n^2) - mp_pi*q(index); 
        % balancing
        [D, M] = balance(A);
        balancing = false;
        if norm(M, 1) < norm(A, 1) 
            A = M;
            balancing = true;
        end
    end

    I = eye(n,'mp');
    B = A * A; % working with A^2;
    B_double = double(B);
    
    if ~any(B_double(:)) % if B is a zero matrix
        if nargout > 0
            [varargout{1}] = I;
        end
        if nargout > 1
            [varargout{2}] = 0;
        end
        if nargout > 2
            [varargout{3}] = 0;
        end
        return
    end
    
    alphas = zeros(1,maxdegree);
    alphas(1) = norm(B_double, 1);
    
    degrees = opt_degrees(maxdegree); % vector of optimal degrees for current approximant.
    maxdegree = degrees(end);        % max degree
    
    i_deg = 1; % i_deg: index of current degree - we start from 1
    m = degrees(i_deg);           % current degree
    s = 0;
    Bsq_powers = zeros(n,n,ceil(sqrt(maxdegree))+1,'mp');
    Bsq_powers(:,:,1) = I;
    Bsq_powers(:,:,2) = B;
    Bsq_powers_lgth = 2; % current length of sequence B is 2 (B_0 = I, B_1 = B);
    Bsq_powers_double = double(Bsq_powers); % form double precision powers of B.
    
    
    abs_err = false; % use relative error bound
    if ~abs_err, compute_normcosm = true;  end
    
    % optional: compute the sum in t^{ch}_{2m} in the error bound using some extra precision
    ext_prec = true; % slight improves the accuracy, and brings no difference to speed. 
    prec_factor = 1.2; % the factor that determines the extra precision
    
    degree_found = false;
    delta_pre = Inf;

    % Prepare the factorials required in evaluating the bound - the 
    % initialization "factorials = factorial(2*mp(0:maxdegree));  
    % factorials_double = double(factorials);" is not used since its
    % execution time can be significant when maxdegree is large.
    % Instead, we initialize as following and compute more factorials 
    % only when needed.
    
    factorials = zeros(1,maxdegree+1,'mp');
    factorials_double = zeros(1,maxdegree+1); % powers of B in double precision for evaluating norm of cos(A/2^s)
    factorials(1:m+1) = factorial(2*mp(0:m));
    factorials_double(1:m+1) = double(factorials(1:m+1));
    
    tempnormcosm1 = Inf; % 'previous' approximation to the norm of cos(A/2^s) used in the function eval_bound
    [delta_nxt, Bsq_powers, tempnormcosm1] = eval_bound(Bsq_powers, alpha(s, m), m, ext_prec);
    if ~abs_err, curr_epsilon = epsilon * tempnormcosm1; end
    if delta_nxt < curr_epsilon, degree_found = true; end  % Awad (< instead of <=)
    
    order_k = 3;
    while ~degree_found && m < maxdegree 
         if abs(delta_pre) <= abs(delta_nxt)^order_k   % Awad
            s = s + 1;
            compute_normcosm = true; % s increased, recompute the 1-norm of cos(A/2^s)
        else
            i_deg = i_deg + 1; % m_i
            m_pre = m; % store the previous degree m
            m = degrees(i_deg);
            % update the factorials required in evaluating the bound
            factorials(m_pre+2:m+1) = factorial(2*mp(m_pre+1:m));
            factorials_double(m_pre+2:m+1) = double(factorials(m_pre+2:m+1));
        end
        delta_pre = delta_nxt;
        [delta_nxt, Bsq_powers, tempnormcosm1] = eval_bound(Bsq_powers, alpha(s, m), m, ext_prec);
        if ~abs_err, curr_epsilon = epsilon * tempnormcosm1; end
        if delta_nxt <= curr_epsilon, degree_found = true; end
    end
    
    % s and m found, now we compute the Taylor polynomial of matrix cosine
    C = PSeval_cos(Bsq_powers, s, m);
    if recompute_diag_blocks
            C = recompute_diagonals(A/(2^s), C);
    end
    
    %% revocering phase.
    for r = 1:s
            C = 2 * C * C - I;
            if recompute_diag_blocks
                C = recompute_diagonals(2^(r-s)*A, C);
            end
    end
    
    %% postprocessing
    if pre_prcs
        C = (-1)^q(index) * C;
        if balancing
            C = D * C / D;
        end
    end
    
    if compute_schur
        C = Q * C * Q';
    end
    if isreal(A)
        C = real(C);
    end

%% Prepare output.
    if nargout >= 0
        [varargout{1}] = C;
    end
    if nargout > 1
        [varargout{2}] = s;
    end
    if nargout > 2
        [varargout{3}] = m;
    end
    
 %% Restore mp working precision.
    mp.Digits(current_digits);
    
    
    
%% -----------------------       SUBFUNCTIONS       -----------------------

    function [err_bnd, Bsq_powers, tempnormcosm] = eval_bound(Bsq_powers, alpha, m, ext_prec)
    %EVAL_BOUND   Error bound and current approx. of the norm of cos(A/2^s).
    %
    %   EVAL_BOUND(BSQ_POWERS, ALPHA, M, EXT_PREC) computes an upper bound 
    %   ERR_BND on the truncation error of the Taylor approximant of 
    %   degree M to COS(A/2^s) and an approximation TEMPNORMCOSM to its 
    %   norm using available powers of B stored in BSQ_POWERS. EXT_PREC 
    %   is a Boolean variable that determines whether we compute the 
    %   scalar coefficients in the error bound using some extra precision. 
    %   ALPHA represents the value alpha^*_m(B/4^s).
    
        nu = floor(sqrt(m));
        for i = (Bsq_powers_lgth+1):nu+1
            Bsq_powers(:,:,i) = Bsq_powers(:,:,i-1) * Bsq_powers(:,:,2);
        end
        % update the double precision powers
        Bsq_powers_double(:,:,Bsq_powers_lgth+1:nu+1) = double(Bsq_powers(:,:,Bsq_powers_lgth+1:nu+1));
        
        Bsq_pwlg_incrd = false; % Bsq_powers_lgth increased = false
        Bsq_pwlg_pre = Bsq_powers_lgth; 
        Bsq_powers_lgth = nu + 1;
        if Bsq_powers_lgth>Bsq_pwlg_pre, Bsq_pwlg_incrd = true; end
        
        digits_old = mp.Digits();
        if ext_prec, mp.Digits(prec_factor*digits_old); end
        alpha = sqrt(mp(alpha));
        t_ch =  sum(alpha.^(2*(0:m))./factorials(1:m+1)); 
        err_bnd = abs(cosh(alpha) - t_ch);
        mp.Digits(digits_old);
        if compute_normcosm % estimate and update 1-norm of cos(A/2^s) using potentially lower precision (double)
            reshaped_vector = reshape(4.^(-s*(0:Bsq_powers_lgth-1))./...
                factorials_double(1:Bsq_powers_lgth),[1,1,Bsq_powers_lgth]);
            approx = sum(bsxfun(@times,...
                                Bsq_powers_double(:,:,1:Bsq_powers_lgth),... % reshaped_matrix,...
                                reshaped_vector),3);
            tempnormcosm = norm(approx, 1);
            % for the same s turn off the estimation if the value of tempnormcosm is
            % converging (when Bsq_powers_lgth increased but diff is small)!
            if Bsq_pwlg_incrd && abs(tempnormcosm1 - tempnormcosm) / ...
                    abs(tempnormcosm) < 0.1
                compute_normcosm = false;
            end
        else
            tempnormcosm = tempnormcosm1;
        end
    end

    function [est,mv] = normBd(B,d)
    %NORMAM   Estimate of 1-norm of power of matrix.
    %   NORMBD(B,d) estimates norm(B^d,1).
    %
    %   If B has nonnegative elements the estimate is exact.
    %   [EST,MV] = NORMBD(B,d) returns the estimate EST and the number MV of
    %   matrix-vector products computed involving B or B^*.

    %   Reference: 
    %   [1] A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 
    %   31(3):970-989, 2009.
    %   [2] M. Fasi and N. J. Higham, A Arbitrary Precision Scaling and 
    %   Squaring Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 
    %   40(4):1233-1256, 2019
    
        nc = 1; % number of columns in the iteration matrix in NORMEST1
        n1 = size(B,1);
        if isequal(B,abs(B))
            e = ones(n1,1);
            for c=1:d         % for positive matrices only
                e = B'*e;
            end
            est = norm(e,inf);
            mv = d;
        else
            % Use powers of B in Bsq_powers in order to efficienlty compute
            % B^d * z.
            % mult_index stores indeces of powers of B for multiplications to be done later
            mult_index = zeros(1, Bsq_powers_lgth);  
            d_dec = d;
            lgth_dec = min(Bsq_powers_lgth, d+1); 
            while d_dec > 0
                mult_index(lgth_dec) = floor(d_dec / (lgth_dec-1));
                d_dec = mod(d_dec, lgth_dec-1);
                lgth_dec = min(lgth_dec - 1, d_dec+1);
            end
        [est,~,~,it] = normest1(@afun_power,nc);
        mv = it(2)*nc*d;
        end

        function prod = afun_power(flag,Z)
        %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

            if isequal(flag,'dim')
                prod = n1;
            elseif isequal(flag,'real')
                prod = isreal(B);
            else
                [n2,~] = size(Z);
                if n2 ~= n1
                    error('Dimension mismatch')
                end
                if isequal(flag,'notransp')
                    % do the matrix-vector multiplication according to mult_index
                    for i = find(mult_index(2:end))
                        for ii = 1:mult_index(i+1)  
                            Z = double(Bsq_powers(:,:,i+1)) * Z; % Bsq_powers starts from the identity
                        end
                    end
                    if mult_index(1) ~= 0
                        for ii = 1:mult_index(1)
                            Z = B * Z;
                        end
                    end
                elseif isequal(flag,'transp')
                    for i = find(mult_index(2:end))
                        for ii = 1:mult_index(i+1)
                            Z = double(Bsq_powers(:,:,i+1)') * Z;
                        end
                    end
                    if mult_index(1) ~= 0
                        for ii = 1:mult_index(1)
                            Z = B' * Z;
                        end
                    end
                end
                prod = Z;
            end
        end
    end

    function y = alpha(s, m)
    %ALPHA    Estimate 4^{-s}*max(||B^d||^(1/d), ||X^(d+1)||^(1/(d+1))).
    %   ALPHA(S,M) estimates 4^{-s}*max(||B^d||^(1/d), ||B^(d+1)||^(1/(d+1))) 
    %   in the 1-norm, where d is the largest integer such that d(d-1) < m+2.
        
        d_star = floor((1 + sqrt(4 * m + 5)) / 2);

        % Experimental strategy
        if alphas(d_star+1) == 0 % The second term is not known, we try the bound.
            if alphas(d_star) == 0
                alphas(d_star) = normBd(B_double, d_star)^(1/d_star);
            end
            known = find(alphas ~= 0);
            low = min(known);
            high = max(known);
            bin_counter = 0;
            found_upper_bound = false;
            while low < high
                if (low + high == d_star+1)
                    if (alphas(d_star) > alphas(low) * alphas(high))
                        found_upper_bound = true;
                        break;
                    end
                end
                if bin_counter
                    low = low + 1;
                else
                    high = high - 1;
                end
                bin_counter = mod(bin_counter + 1, 2);
            end
            if found_upper_bound
                %             fprintf('.\n');
                y = alphas(d_star) / 4^s;
                %             i = j;
                return
            else
                assert(alphas(d_star+1) == 0)
                alphas(d_star+1) = normBd(B_double, d_star+1)^(1/(d_star+1));
            end
        end
        [y,~] = max(alphas(d_star:d_star+1));
        y = y / 4^s;
    end

    function C = PSeval_cos(Bsq_powers, s, m)
    %EVAL_POLY_PS   modified Paterson-Stockmeyer method for matrix polynomials.
    %
    %   EVAL_POLY_PS(Bsq_powers, S, M) evaluates p^c_m(4^{-s}B) using the
    %   Paterson-Stockmeyer method.

        ss = floor(sqrt(m));
        if ss == 0
            C = I;
            return
        end
        rr = floor(m/ss);

        scal = 4^s;

        C = Bsq_powers(:,:,ss+1) * (-1) / scal / ((2*m) *(2*m-1)) ;
        for j = ss-1:-1:1
            C = (C + Bsq_powers(:,:,j+1)) * (-1) / scal / ...
                ( (2*(m-ss+j)) *(2*(m-ss+j)-1) );
        end
        
        C = Bsq_powers(:,:,ss+1) * C;
        for k = rr - 2:-1:1
            for j = ss:-1:1
                C = (C + Bsq_powers(:,:,j+1)) * (-1) / scal /...
                    ( (2*(k*ss+j)) *(2*(k*ss+j)-1) );
            end
            C = Bsq_powers(:,:,ss+1) * C;
        end
        
        for j = ss:-1:1
            C = (C + Bsq_powers(:,:,j+1)) * (-1) / scal / ((2*j) *(2*j-1));
        end
        C(1:n + 1:n^2) = C(1:n + 1:n^2) + 1; % add the identity matrix
    end

    function degs = opt_degrees(mmax)
    %OPT_DEGREES    Vector DEGS stores optimal degrees of Taylor
    %approximants to the cosine.
        degs = [1,    2,    4,    6,    9,   12,   16,   20,   25,...
                30,   36,   42,   49,   56,   64,   72,   81,   90,  100,...
                110,  121,  132,  144,  156,  169,  182,  196,  210,  225,...
                240,  256,  272,  289,  306,  324,  342,  361,  380,  400,...
                420,  441,  462,  484,  506,  529,  552,  576,  600,  625,...
                650,  676,  702,  729,  756,  784,  812,  841,  870,  900,...
                930,  961,  992, 1024, 1056, 1089, 1122, 1156, 1190, 1225,...
                1260, 1296, 1332, 1369, 1406, 1444, 1482, 1521, 1560, 1600,...
                1640, 1681, 1722, 1764, 1806, 1849, 1892, 1936, 1980, 2025,...
                2070, 2116, 2162, 2209, 2256, 2304, 2352, 2401, 2450, 2500];
        degs = degs(1:find(degs <= mmax, 1, 'last'));
    end

    function F = recompute_diagonals(T, F)
    %Recomputation of diagonal blocks.
        i = 1;
        while i <= size(T, 1)
            if (n == i + 1) || ((i <= n - 2) && (T(i+2,i+1) < epsilon))
                % start of 2-by-2 block
                if T(i+1,i) < epsilon % triangular block
                    F(i:i+1,i:i+1) = cosm2by2_tri(T(i:i+1,i:i+1));
                else % full block
                    F(i:i+1,i:i+1) = cosm2by2_full(T(i:i+1,i:i+1));
                end
                i = i + 2; % skip next element if it is 2-by-2 block
            else % start of 1-by-1 block, no superdiagonal
                F(i,i) = cos(T(i,i));
                i = i + 1;
            end
        end
    end

    function F = cosm2by2_tri(T)
    %COSM2BY2_TRI	Trigonometric cosine of 2-by-2 upper triangular matrix.
    %   COSM2BY2_TRI(T) computes the cosine of the 2-by-2 upper
    %   triangular matrix T.

    % Equation (3.1) and (3.3) of
    % A. H. Al-Mohy, N. J. Higham, and S. D. Relton. New Algorithms for
    % Computing the Matrix Sine and Cosine Separately or Simultaneously. 
    % SIAM J. Sci. Comput. 37(1):A456-A487, 2015. 
    
        diagT = diag(T);
        cosT = cos(diagT);
        F = diag(cosT);
        
        if abs(diagT(1)-diagT(2))<epsilon % diagT(1) = diagT(2)
            cos_divdiff = - sin(diagT(1));
        else
            cos_divdiff = - sin( (diagT(1)+diagT(2))/2 ) * ...
                sin( (diagT(1)-diagT(2))/2 ) / ( (diagT(1)-diagT(2))/2 );
        end
        F(1,2) = T(1,2) * cos_divdiff;
    end

    function F = cosm2by2_full(X)
    %COSM2BY2_FULL   Trigonometric cosine of 2-by-2 full matrix.
    %   COSM2BY2_FULL(T) computes the cosine of the 2-by-2 full matrix X.

    % Equation (3.6) of
    % A. H. Al-Mohy, N. J. Higham, and S. D. Relton. New Algorithms for
    % Computing the Matrix Sine and Cosine Separately or Simultaneously. 
    % SIAM J. Sci. Comput. 37(1):A456-A487, 2015. 

        a = X(1,1);
        b = X(1,2);
        c = X(2,1);
        % X(2,2) = a;
        theta = sqrt(-b*c);
        
        F(1,1) = cos(a) * cosh(theta);
        F(2,1) = - c * sin(a) * sinh(theta) / theta;
        F(1,2) = - b * sin(a) * sinh(theta) / theta;
        F(2,2) = F(1,1);

    end

    function e = myeps(curr_class)
    %COMP_EPS	Machine epsilon.
    %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
        if(strcmp(curr_class, 'mp'))
            e = mp('eps');
        else
            e = eps(curr_class)/2;
        end
    end

end