#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "readMDinformation.h"
#include "header.h"
#include "basic.h"
#include "graph.h"

/** initialize gobal vars **/
// force field parameters
double mass_mol = 18.0154;
double rc = 10;
double mass[3] = {15.9994,1.008,1.008};
double eps = 0.1553;
double sig = 3.166;
double charge[3] = {-0.8476,0.4238,0.4238};

// unit conversion
double NA = 6.02e23;
double kcal_mol2j = 6.9477e-21; 
double A2m = 1e-10;
double e2c = 1.60217663e-19;
double a_fs2m_s = 1e5;
double g_mol2kg = 1.6611e-27;

int main(int argc, char ** argv) {
    int N = 8544;
    struct MOL mol[N];
    int i;
    
    double **peMatrix = (double **)malloc(N * sizeof(double *));
    int **adjacencyMatrix = (int **)malloc(N * sizeof(int *));
    int **clusters = (int **)malloc(N * sizeof(int *));
    for (i = 0; i < N; ++i) {
        peMatrix[i] = (double *)malloc(N * sizeof(double));
        adjacencyMatrix[i] = (int *)malloc(N * sizeof(int));
        clusters[i] = (int *)malloc(N * sizeof(int));
    }
    int numClusters;
    int clusterSizes[N];
    int nodesDegree[N];

    FILE *fpr_position;
    FILE *fpr_velocity;
    char resultDir[FILENAME_MAX];
    char simulationDir[FILENAME_MAX];
    char filename_position[FILENAME_MAX];
    char filename_velocity[FILENAME_MAX];
    char filename_result[FILENAME_MAX];

    int p ,t, numFrame, save_every;

    if(argc == 7) {
        p = atoi(argv[1]);
        t = atoi(argv[2]);
        numFrame = atoi(argv[3]);
        save_every = atoi(argv[4]);
        strcpy(resultDir, argv[5]);
        strcpy(simulationDir, argv[6]);
    }
    else {
        // comment out these two lines if debugging locally
        printf("Not enough input variables. Quiting...\n");
        return 0; 

        printf("Local debugging:\n");
        p = 240;
        t = 300;
        numFrame = 10;
        save_every = 5;
        sprintf(resultDir, "/Users/jingcunfan/Documents/data_cp_network/Hill_SPCE");
        sprintf(simulationDir, "/Users/jingcunfan/Documents/MD_simulations/sc_water/cluster");
        
    }
    printf("P = %d atm\n", p);
    printf("T = %d K\n", t);
    printf("Number of frames = %d\n", numFrame);
    printf("Saving frequency: every %d frames\n", save_every);
    printf("Result file directory: %s\n", resultDir);
    printf("Simulation file directory: %s\n", simulationDir);

    // Create result directory
    char cmd[100];
    sprintf(resultDir, "%s/%datm", resultDir, p);
    sprintf(cmd, "mkdir -p %s", resultDir);
    printf("%s\n", cmd);
    system(cmd);
    sprintf(filename_result, "%s/%datm_%dK.txt", resultDir, p, t);
    
    // Open simulation dump file
    sprintf(filename_position, "%s/SPCE_%datm/%dK/position.xmol", simulationDir, p, t);
    sprintf(filename_velocity, "%s/SPCE_%datm/%dK/velocity.xmol", simulationDir, p, t);
    fpr_position = fopen(filename_position,"r");
    fpr_velocity = fopen(filename_velocity,"r");

    int size_count_accumulated[N], degree_count_accumulated[N];
    initializeIntArray(N, size_count_accumulated, 0);
    initializeIntArray(N, degree_count_accumulated, 0);
    
    double PE_s_accumulated[N], KE_s_accumulated[N];
    initializeDoubleArray(N, PE_s_accumulated, 0);
    initializeDoubleArray(N, KE_s_accumulated, 0);

    // Traverse through each snapshot
    printf("-----------------------------------------------------\n");
    printf("Traverse through each frame:\n");
    int frame;
    for (frame = 0; frame < numFrame; frame++) {
        int numMol;
        double boxLength;

        // Load atom information  
        readAtomInformation(fpr_position, fpr_velocity, &numMol, &boxLength, mol);

        // Calculate molecule energies
        calculatePE_perMol(mol, numMol, peMatrix, boxLength);
        calculateKE_perMol(mol, numMol);
        calculateCOMvelcoity_perMol(mol, numMol);

        // Find clusters
        getAdjacencyMatrix(mol, adjacencyMatrix, peMatrix, numMol);
        findConnectedComponents(adjacencyMatrix, numMol, &numClusters, clusterSizes, clusters);
        sumRowsInt2DMatrix(numMol, numMol, adjacencyMatrix, nodesDegree);
        
        // Cluster size distribution and degree distribution
        histogram(clusterSizes, numClusters, size_count_accumulated);
        histogram(nodesDegree, numMol, degree_count_accumulated);

        // Energy per cluster vs. cluster size
        energy_perCluster_size(mol, numClusters, clusterSizes, clusters, PE_s_accumulated, KE_s_accumulated);

        printf("Finished with frame %d. There are %d clusters in this frame.\n", frame, numClusters);

        // Save every certain frames
        if((frame + 1) % save_every == 0){
            save_result(filename_result, N, frame + 1, size_count_accumulated, degree_count_accumulated, PE_s_accumulated, KE_s_accumulated);
            printf("Result averaged over all %d frames and saved to %s\n", frame + 1, filename_result);
        }
    }
    fclose(fpr_position);
    fclose(fpr_velocity);

    // Save result
    save_result(filename_result, N, numFrame, size_count_accumulated, degree_count_accumulated, PE_s_accumulated, KE_s_accumulated);
    printf("Result averaged over all %d frames and saved to %s\n", numFrame, filename_result);

    free(peMatrix);
    free(adjacencyMatrix);
    free(clusters);

    printf("All finished.\n");
    return 0;
}