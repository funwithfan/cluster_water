#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "header.h"
#include "basic.h"

void readAtomInformation(char *filename, int *numMol, double *boxLength, struct MOL *mol);
double trueDistance(double distance, double box_length);
void findDistanceMatrix(double **r, struct MOL *mol, int mol1_id, int mol2_id, double boxLength, bool *isWithinCutoff);
double intermolecularPE(double **r);
void calculatePE_perMol(struct MOL *mol, int numMol, double **peMatrix, double boxLength);
void calculateKE_perMol(struct MOL *mol, int numMol);
void calculateCOMvelcoity_perMol(struct MOL *mol, int numMol);
double relativeKE(struct MOL *mol, int mol1_id, int mol2_id);
void getAdjacencyMatrix_Hill(struct MOL *mol, int **adjacencyMatrix, double **peMatrix, int numMol);
void getAdjacencyMatrix_distance(struct MOL *mol, int **adjacencyMatrix, int numMol, double boxLength, double cutoff);
void energy_perCluster_size(struct MOL *mol, int numClusters, int *clusterSizes, int **clusters, double *PE_perCluster_size, double *KE_perCluster_size);
void inner_energy_perCluster_size(struct MOL *mol, double **peMatrix, int numClusters, int *clusterSizes, int **clusters, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated);
void save_result(char *filename, int max_size, int numFrame, int *size_count_accumulated, int *degree_count_accumulated, double *PE_s_accumulated, double *KE_s_accumulated);
