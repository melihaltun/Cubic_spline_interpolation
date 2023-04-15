/**
* @file cubicSpline.h
* @author Melih Altun @2015-2023
**/

#include "matrixOperations.h"
#include "Cholesky.h"

void cubicSpline(float c0[], float c1[], float c2[], float c3[], float x[], float y[], int n);
void evaluate_polynomial(float y2[], float x2[], int n2, float x1[], int n1, float c0[], float c1[], float c2[], float c3[]);
