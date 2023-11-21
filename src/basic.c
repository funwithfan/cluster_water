#include "basic.h"

/** Count frequency **/
void histogram(int *array, int size, int *count) {
    int i, value;
    //initializeIntArray(size, count, 0);

    // Count occurrences of values in array and update bins
    for ( i = 0; i < size; i++) {
        value = array[i];
        count[value]++;
    }
}

/** Initialize int array **/
void initializeIntArray(int size, int *array, int value) {
    int i;
    for (i = 0; i < size; i++) {
        array[i] = value;
    }
}

/** Initialize double array **/
void initializeDoubleArray(int size, double *array, double value) {
    int i;
    for (i = 0; i < size; i++) {
        array[i] = value;
    }
}

/** Initialize bool array **/
void initializeBoolArray(int size, bool *array, bool value) {
    int i;
    for (i = 0; i < size; i++) {
        array[i] = value;
    } 
}

/** Initialize int matrix **/
void initializeIntMatrix(int rows, int cols, int **array, int value) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            array[i][j] = value;
        }
    }
}

/** Initialize double matrix **/
void initializeDouble2DMatrix(int rows, int cols, double **array, double value) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            array[i][j] = value;
        }
    }
}

// maxmum value in a double matrix
double maxDouble2DMatrix(int rows, int cols, double **array) {
    int i, j;
    double max = array[0][0];

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (array[i][j] > max) {
                max = array[i][j];
            }
        }
    }
    return max;
}

// minmum value in a double matrix
double minDouble2DMatrix(int rows, int cols, double **array) {
    int i, j;
    double min = array[0][0];

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (array[i][j] < min) {
                min = array[i][j];
            }
        }
    }
    return min;
}


/** Index of max value in an int array **/
int findMaxValueIndex(int *array, int size) {
    int maxIndex = 0;
    int i;
    for (i = 0; i < size; i++) {
        if (array[i] > array[maxIndex]) {
            maxIndex = i;
        }
    }
    return maxIndex;
}

/** Average value of the array **/
double meanDoubleArray(double *array, int length) {
    int i;
    double mean = 0;
    for (i = 0; i < length; i++) {
        mean = mean + array[i];
    }
    mean = mean / length;
    return mean;
}

// Sum of rows in a 2D int matrix
void sumRowsInt2DMatrix(int rows, int cols, int **array, int *sumRows) {
    int i, j;
    for(i = 0;i < rows; i++) {
        sumRows[i] = 0;
        for(j = 0; j < cols; j++) {
            sumRows[i] += array[i][j];
        }
    }
}

// Dot product of double arrays
void doubleArraysDotProduct(int size, double *array1, double *array2, double *product) {
    int i;
    for (i = 0; i < size; i++) {
        product[i] = array1[i] * array2[i];
    }
}

// Product of a int array and a double constant
void intArraysMutiplyConstant(int size, int *array, double constant, double *product) {
    int i;
    for (i = 0; i < size; i++) {
        product[i] = array[i] * constant;
    }
}

// pointwise division of two double arrays
void doubleArraysDotDivision(int size, double *array1, double *array2, double *product) {
    int i;
    for (i = 0; i < size; i++) {
        product[i] = array1[i] / array2[i];
    }
}

// Write 2D double matrix to file
void saveDouble2DMatrix(char *filename, int rows, int cols, double **array) {
    int i, j;
    FILE *fpw = fopen(filename, "w");
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
            fprintf(fpw, "%f ", array[i][j]);
        }
        fprintf(fpw, "\n");
    }
    fclose(fpw);
    printf("Matrix saved to %s\n", filename);
}

// Write 2D int matrix to file
void saveInt2DMatrix(char *filename, int rows, int cols, int **array) {
    int i, j;
    FILE *fpw = fopen(filename, "w");
    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
            fprintf(fpw, "%d ", array[i][j]);
        }
        fprintf(fpw, "\n");
    }
    fclose(fpw);
    printf("Matrix saved to %s\n", filename);
}

// Write 1D int array to file
void saveInt1DArray(char *filename, int length, int *array) {
    int i;
    FILE *fpw = fopen(filename, "w");
    for(i = 0; i < length; i++) {
        fprintf(fpw, "%d\n", array[i]);
    }
    fclose(fpw);
    printf("Array saved to %s\n", filename);
}