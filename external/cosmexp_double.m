function [C, S] = cosmexp_double(A)
%COSMEXP_DOUBLE   Compute the matrix sine and cosine in double precision  
%   via the matrix exponential, where the exponential is computed using 
%   the built-in EXPM function.

X = expm(1i*A);
if isreal(A)
    C = real(X);
    S = imag(X);
else
    invX = expm(-1i*A);
    C = (X + invX) / 2;
    S = (X - invX) / (2*1i);
end