/**
* @file testCubicSpline.h
* @author Melih Altun @2015-2023
**/

#include "cubicSpline.h"

int main()
{
	float x1[] = { 0, 1, 2, 3, 4 };
	float y1[] = { 4, -1, 2, 5, 1 };
	float x2[] = { 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4 };
	float y2[9], c0[4], c1[4], c2[4], c3[4];   // coeff sizes are size(x1) - 1

	cubicSpline(c0, c1, c2, c3, x1, y1, 5);

	evaluate_polynomial(y2, x2, 9, x1, 5, c0, c1, c2, c3);

	return 0;
}

