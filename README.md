# Cubic_spline_interpolation
A cubic spline interpolation implementation that returns cubic polynomial coefficients for each segment of inputs.
Once polynomial coefficients are found, a secondary function can be called which outputs spline interpolated values for given inputs.

This code takes an x and a corresponding y sequence, each with length n, and returns cubic polynomial coefficients for each segment. Therefore the returned c0[], c1[], c2[] and c3[] arrays all have size n-1.

Once the coefficient arrays are calculated evaluate_polynomial function will calculate interpolated y2 values for an x2 series, which should cover the same range as the original x, but the number of elements can be different. 

Cholesky_Decomposition repository is required to perform matrix inversion and some other matrix operations:
https://github.com/melihaltun/Cholesky_Decomposition
