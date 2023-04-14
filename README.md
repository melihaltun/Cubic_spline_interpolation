# Cubic_spline_interpolation
A cubic spline interpolation implementation that returns cubic polynomial coefficients for each segment of inputs.

This code takes an x and a corresponding y sequence, each with length n, and returns cubic polynomial coefficients for each segment. Therefore the returned c0[], c1[], c2[] and c3[] arrays all have size n-1.

The additonal function that evaulates these coefficients and returns an interpolated y2 sequence for a given x2 sequence will be added later. 

Cholesky_Decomposition repository is required to perform matrix inversion and some other matrix operations:
https://github.com/melihaltun/Cholesky_Decomposition
