# Arbitrary precision algorithms for computing the matrix cosine and its Fréchet derivative. 

Function `include/cosm_mp.m` computes the matrix cosine of a square matrix in arbitrary precision. Function `include/cosm_double.m` is the counterpart that works in double precision and uses no `mp` computation. Function `include/cosm_frechet_mp.m` computes the matrix cosine and its Fréchet derivative simoutaneously in arbitrary precision. 

Details on the underlying algorithms can be found in the technical report:

Awad H. Al-Mohy, N. J. Higham and X. Liu. Arbitrary Precision Algorithms for Computing the Matrix Cosine and its Fréchet Derivative, MIMS EPrint 2021.x, 2021.

All codes used for generating the data in the above report are included in this repository.

## Dependencies

The code in this repository may require the Advanpix Multiprecision Computing
Toolbox for MATLAB (www.advanpix.com).

## License

See `license.txt` for licensing information.
