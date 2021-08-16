function [C, S] = cosmexp_mp(A)
%COS_VEF   Compute the matrix sine and cosine in arbitrary precision via 
%   the matrix exponential, where the exponential is computed using the
%   arbitrary precision scaling and squareing algorithm expm_mp developed
%   in
%
%   M. Fasi and N. J. Higham, A Arbitrary Precision Scaling and 
%   Squaring Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 
%   40(4):1233-1256, 2019


X = expm_mp(1i*A);
if isreal(A)
    C = real(X);
    S = imag(X);
else
    invX = expm_mp(-1i*A);
    C = (X + invX) / 2;
    S = (X - invX) / (2*1i);
end