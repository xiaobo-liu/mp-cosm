function L = cosm_frechet_blk(A, E)
% This function compute the Frechet derivative L_cos(A,E) by applying 
% cosm_mp to a 2 by 2 block matrix.


n = size(A, 1);
O = zeros(n,class(A));

Z = [A E; O A];
F = cosm_mp(Z);

L = F(1:n,n+1:2*n);

end