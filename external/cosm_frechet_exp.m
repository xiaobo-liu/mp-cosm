function L = cosm_frechet_exp(A, E)
% This function compute the Frechet derivative L_cos(A,E) by exploiting
% L_cos(A,E) = i/2 * (L_exp(iA,E) - L_exp(-iA,E)) and applying the algorithm
% EXPM_FRECHET_PADE to computing the Frechet derivative of the matrix
% exponential. This algorithm is intended for double precision only.

L = (expm_frechet_pade(1i*A,E) - expm_frechet_pade(-1i*A,E)) * 1i / 2;

end