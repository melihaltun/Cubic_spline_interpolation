/**
* @file cubicSpline.cpp
* @author Melih Altun @2015
**/

#include "cubicSpline.h"

/* Finds polynomial coefficients for cubic spline interpolation
parameters: (outputs) c0, c1, c2, c3 poly coefficients for each spline interval
			(inputs) x coordinates, y coordinates, number of points */
void cubicSpline(float c0[], float c1[], float c2[], float c3[], float x[], float y[], int n)
{
	int i;
	float *diff_x, *diff_y, *A, *A_inv, *m, *b;
	diff_x = new float[n - 1];
	diff_y = new float[n - 1];
	A = new float[(n - 2)*(n - 2)];
	A_inv = new float[(n - 2)*(n - 2)];
	b = new float[n - 2];
	m = new float[n];

	memset(A, 0, (n - 2)*(n - 2)*sizeof(float));

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



// Evaluate the cubic polynomial for a single segment using the segment's coefficients and a given x value
float evaluate_segment(float x, float coeffs_segment[]) {
	float y = coeffs_segment[3] * pow(x, 3) + coeffs_segment[2] * pow(x, 2) + coeffs_segment[1] * x + coeffs_segment[0];
	return y;
}


// finds the x segment where calculations will take place
int getSegmentIndex(float x1[], int n1, float x)
{
	for (int i = n1 - 1; i >= 0; i--) {
		if (x >= x1[i])
			return i;
	}
	return 0;
}

// Evaluate the cubic polynomial using all of the segment coefficients and a given x value
// parameters: (output) y2: interpolated y value;
//             (inputs) x2: new x values for interpolation, n2: number of x2 and y2 elements, x1: x coordinates before interpolation, 
//             (inputs) c0, c1, c2, c3: poly coefficients for each spline interval
void evaluate_polynomial(float y2[], float x2[], int n2, float x1[], int n1, float c0[], float c1[], float c2[], float c3[])
{
	float x;
	int n_segments, segment_index, start_index, end_index;
	float segment_coeffs[4];
	n_segments = n1 - 1;
	for (int i = 0; i < n2; i++) {
		x = x2[i];
		segment_index = getSegmentIndex(x1, n1, x);

		if (segment_index == n_segments)
			segment_index = n_segments - 1;

		segment_coeffs[0] = c0[segment_index];
		segment_coeffs[1] = c1[segment_index];
		segment_coeffs[2] = c2[segment_index];
		segment_coeffs[3] = c3[segment_index];
	
		y2[i] = evaluate_segment(x-x1[segment_index], segment_coeffs);
	}
}
