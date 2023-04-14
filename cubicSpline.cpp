/**
* @file cubicSpline.cpp
* @author Melih Altun @2015
**/

#include "cubicSpline.h"

/* Finds polynomial coefficients for cubic spline interpolation
parameters: (outputs) c0, c1, c2, c3 poly coefficients for each spline interval
			(inputs) x coordinates, y coordinates, number of points */
void cubicSpline(double c0[], double c1[], double c2[], double c3[], double x[], double y[], int n)
{
	int i;
	double *diff_x, *diff_y, *A, *A_inv, *m, *b;
	diff_x = new double[n - 1];
	diff_y = new double[n - 1];
	A = new double[(n - 2)*(n - 2)];
	A_inv = new double[(n - 2)*(n - 2)];
	b = new double[n - 2];
	m = new double[n];

	memset(A, 0, (n - 2)*(n - 2)*sizeof(double));

	//dx, dy
	for (i = 0; i < n - 1; i++) {
		diff_x[i] = x[i + 1] - x[i];
		diff_y[i] = (y[i + 1] - y[i]) / diff_x[i];
	}

	//construct tri-diagonal matrix of basis functions and difference of 1st derivative of y
	for (i = 0; i < n - 2; i++) {
		if (i>0)
			A[lin_index(i, i - 1, n-2)] = diff_x[i];
		A[lin_index(i, i, n-2)] = 2 * (diff_x[i] + diff_x[i + 1]);
		if (i < n - 3)
			A[lin_index(i, i + 1, n-2)] = diff_x[i + 1];
		b[i] = 6 * (diff_y[i + 1] - diff_y[i]);
	}

	cholesky_inverse(A_inv, A, n - 2);
	multiply_matrix_with_vector(m, A_inv, b, n - 2, n - 2);  // m = A \ b

	//second derivatives with zero boudary conditions
	for (i = n - 2; i>0; i--)
		m[i] = m[i - 1];
	m[0] = m[n - 1] = 0;  //shift m and zero pad on both sides

	//calculate 0th, 1st, 2nd and 3rd poly coeffs. for each interval
	for (i = 0; i < n-1; i++) {
		c0[i] = y[i];
		c1[i] = diff_y[i] - diff_x[i] * (2 * m[i] + m[i + 1]) / 6;
		c2[i] = m[i] / 2;
		c3[i] = (m[i + 1] - m[i]) / (6 * diff_x[i]);
	}

	//clean up
	delete[] diff_x;
	delete[] diff_y;
	delete[] A;
	delete[] A_inv;
	delete[] m;
	delete[] b;
}
