function [condrel] = fd_condest(fun_fd, A, E, varargin)
%FD_CONDEST Estimates condition number of Frechet derivative Lf(A,E)

%   Reference: Nicholas J. Higham and Samuel D. Relton, Estimating the 
%   Condition Number of the Frechet Derivative of a Matrix Function, 
%   SIAM J. Sci. Comput., 36(6):C617-C634, 2014

n = length(A);

%Get the function value, FD and rel. condition number
[F, L, cn] = feval(fun_fd, A, E, varargin{:});
% [F, L, cn] = fun_fd(A,E);
nrmL = norm(L, 1);
nrmE = norm(E, 1);
nrmA = norm(A, 1);

%Make the condition number absolute
cn = cn * norm(F, 1)/nrmA;

%Scaling factor to get correct rel. condition number
s = nrmA / nrmE;

%Perform the power method on the large matrix using normest1
maxpert = normest1(@second_fd, 2, []);
    
%Convert this into the rel. condition number for the FD Lf(A,E)
condrel = (cn + s * maxpert) * nrmE / nrmL;

%Nested Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [vecL2] = second_fd(flag, vecdA)
        t = size(vecdA, 2);
        %When called by normest1 need to use flags
        skip = false;
        if strcmp(flag, 'real')
            vecL2 = isreal(A) && isreal(E); skip = true;
        elseif strcmp(flag, 'dim')
            vecL2 = n^2; skip = true;
        elseif strcmp(flag ,'transp')
            D = zeros(n,n,t); %The second FD
            dA = zeros(n,n,t);
            for k = 1:t
                dA(:,:,k) = unvec(vecdA(:,k));
                dA(:,:,k) = dA(:,:,k)'; %Transpose
            end
        elseif strcmp(flag, 'notransp')
            D = zeros(n,n,t); %The derivatives
            dA = zeros(n,n,t);
            for k = 1:t
                dA(:,:,k) = unvec(vecdA(:,k));
            end
        end
        %If we only need to return some data from the flags then skip
        if ~skip
            for k = 1:t
                %Compute the second FD
    %             temp = fun_fd([A, dA(:,:,k);zeros(n),A], blkdiag(E,E));
%                 dA(:,:,k)
                temp = feval(fun_fd,[A, dA(:,:,k);zeros(n),A],blkdiag(E,E),varargin{:});
                D(:,:,k) = temp(1:n,n+1:end);
            end
            if strcmp(flag, 'transp')
                for k = 1:t
                    D(:,:,k) = D(:,:,k)';
                end
            end
            vecL2 = zeros(n^2, t);
            for k = 1:t
                vecL2(:,k) = vec(D(:,:,k));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    function V = unvec(vec)
    %UNVEC Unvectorizes the column vector vec
    V = zeros(n);
        for b = 0:n-1
            V(:,b+1) = vec(1+b*n:(b+1)*n);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    function v = vec(mat)
    v = mat(:);
    end

    end