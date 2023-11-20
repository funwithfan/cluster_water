#include <stdio.h>
#include <stdbool.h>

void histogram(int *array, int size, int *bins);
void initializeIntArray(int size, int *array, int value);
void initializeDoubleArray(int size, double *array, double value);
void initializeBoolArray(int size, bool *array, bool value);
void initializeIntMatrix(int rows, int cols, int **array, int value);
void initializeDouble2DMatrix(int rows, int cols, double **array, double value);
int findMaxValueIndex(int *array, int size);
double meanDoubleArray(double *array, int length);
double maxDouble2DMatrix(int rows, int cols, double **array);
double minDouble2DMatrix(int rows, int cols, double **array);
void sumRowsInt2DMatrix(int rows, int cols, int **array, int *sumRows);
void doubleArraysDotProduct(int size, double *array1, double *array2, double *product);
void intArraysMutiplyConstant(int size, int *array, double constant, double *product);
void doubleArraysDotDivision(int size, double *array1, double *array2, double *product);