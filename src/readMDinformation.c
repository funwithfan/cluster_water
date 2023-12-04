#include "readMDinformation.h"

void readAtomInformation(char *filename, int *numMol, double *boxLength, struct MOL *mol) {
    FILE *fpr;
    int i, j, index, mol_id;
    double lo, hi;
    char str[100];

    fpr = fopen(filename, "r");

    // read lines of dump file
    fgets(str,100,fpr); // TIMESTEP
    fgets(str,100,fpr); // timestep
    fgets(str,100,fpr); // ITEM: NUMBER OF ATOMS
    fgets(str,100,fpr); // number of atoms
    *numMol = atoi(str) / 3;
    fgets(str,100,fpr); // ITEM: BOX BOUNDS pp pp pp
    fscanf(fpr,"%s",str); // xlo
    lo = atof(str);
    fscanf(fpr,"%s",str); // xhi
    hi = atof(str);
    *boxLength = hi - lo;
    fgets(str,100,fpr); // skip line?
    fgets(str,100,fpr); // ylo yhi
    fgets(str,100,fpr); // zlo zhi
    fgets(str,100,fpr); // ITEM: ATOMS id mol type x y z ix iy iz vx vy vz element

    // read atom informations
    for(i = 0; i < *numMol * 3; i++) {
        fscanf(fpr,"%s",str); // id
        index = (atoi(str) - 1) % 3; // O H1 H2 for 0 1 2
        mol_id = ceil(atof(str)/3) - 1;
        fscanf(fpr, "%s", str); // type
        fscanf(fpr, "%s", str); // x
        mol[mol_id].position[index][0] = atof(str);
        fscanf(fpr, "%s", str); // y
        mol[mol_id].position[index][1] = atof(str);
        fscanf(fpr, "%s", str); // z
        mol[mol_id].position[index][2] = atof(str);
        fscanf(fpr, "%s", str); // vx
        mol[mol_id].velocity[index][0] = atof(str);
        fscanf(fpr, "%s", str); // vy
        mol[mol_id].velocity[index][1] = atof(str);
        fscanf(fpr, "%s", str); // vz
        mol[mol_id].velocity[index][2] = atof(str);
        fscanf(fpr, "%s", str); // element
    }
    fclose(fpr);
}


// Find the true distance considering pbc
double trueDistance(double distance, double box_length){
    if (distance > box_length/2) {
        distance -= box_length;
    }
    else if (distance < -box_length/2) {
        distance += box_length;
    }
    return distance;
}

// Calculate distance matrix between two molcules
void findDistanceMatrix(double **r, struct MOL *mol, int mol1_id, int mol2_id, double boxLength, bool *isWithinCutoff) {
    int i, j;
    double dx, dy, dz;
    *isWithinCutoff = true;

    for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
            dx = mol[mol1_id].position[i][0] - mol[mol2_id].position[j][0];
            dy = mol[mol1_id].position[i][1] - mol[mol2_id].position[j][1];
            dz = mol[mol1_id].position[i][2] - mol[mol2_id].position[j][2];
            dx = trueDistance(dx, boxLength);
            dy = trueDistance(dy, boxLength);
            dz = trueDistance(dz, boxLength);
            r[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
            if (r[i][j] > rc) {
                *isWithinCutoff = false;
                break;
            }
        }
        if (!*isWithinCutoff) {
            break;
        }
    }
} 


// Calculate inter-molecular potential energy between two molecules, unit: kcal/mol
double intermolecularPE(double **r) {
    double pe, vdw, coul;
    int i, j;

    // van der Waals
    vdw = 4 * eps * (pow(sig/r[0][0], 12) - pow(sig/r[0][0], 6));

    // Coulomb
    double Pi = 3.1415926;
    double eps_0 = 8.85418782e-12; // F/m, dielectric permittivity
    double C = 1 / (4 * Pi * eps_0); // energy-conversion constant
    int dielectric = 1 ; // dielectric constant
    coul = 0;
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            coul += C * charge[i] * charge[j] / (dielectric * r[i][j]);   
        }
    }

    pe = vdw + coul * e2c * e2c / A2m / kcal_mol2j;
    return pe;
}

// Calculate inter-molecular potential energy of each molecule, unit: kcal/mol
void calculatePE_perMol(struct MOL *mol, int numMol, double **peMatrix, double boxLength) {
    int i, mol1_id, mol2_id;
    double pe;
    double **r = (double **)malloc(3 * sizeof(double *));
    for (i = 0; i < 3; ++i) {
        r[i] = (double *)malloc(3 * sizeof(double));
    }
    bool isWithinCutoff;

    
    // Initialize PE of each molecule
    for (i = 0; i < numMol; i++) {
        mol[i].PE = 0;
    }
    
    // Traverse through each pair
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id + 1; mol2_id < numMol; mol2_id++) {
            findDistanceMatrix(r, mol, mol1_id, mol2_id, boxLength, &isWithinCutoff);
            if (isWithinCutoff) {
                pe = intermolecularPE(r);
            }
            else {
                pe = 0;
            }
            peMatrix[mol1_id][mol2_id] = pe;
            peMatrix[mol2_id][mol1_id] = pe;
            mol[mol1_id].PE += 0.5 * pe;
            mol[mol2_id].PE += 0.5 * pe;
        }
    }
    free(r);
}

// Calculate kinetic energy of each molecule, unit: kcal/mol
void calculateKE_perMol(struct MOL *mol, int numMol) {
    int molId, atomIdx, direction;
    double mv2;

    for (molId = 0; molId < numMol; molId++) {
        mv2 = 0;
        for (direction = 0; direction < 3; direction++) {
            for (atomIdx = 0; atomIdx < 3; atomIdx++) {
                mv2 += mass[atomIdx] * mol[molId].velocity[atomIdx][direction] * mol[molId].velocity[atomIdx][direction];
            } 
        }
        mol[molId].KE = 0.5 * mv2 * g_mol2kg * a_fs2m_s * a_fs2m_s / kcal_mol2j;
    }
}

// Calculate COM velocity each molecule
void calculateCOMvelcoity_perMol(struct MOL *mol, int numMol) {
    int molID, atomIdx, direction;

    for (molID = 0; molID < numMol; molID++) {
        for (direction = 0; direction < 3; direction++) {
            mol[molID].velocity_com[direction] = 0;
            for (atomIdx = 0; atomIdx < 3; atomIdx++) {
                mol[molID].velocity_com[direction] +=  mol[molID].velocity[atomIdx][direction] * mass[atomIdx];
            } 
            mol[molID].velocity_com[direction] /=  mass_mol;
        }
    }
}

// Calculate relative kinetic energy between two molecules, unit: kcal/mol
double relativeKE(struct MOL *mol, int mol1_id, int mol2_id){
    double v_relative_square, relative_ke;
    double v_relative[3];
    int i;
    v_relative_square = 0;

    for (i = 0; i < 3; i++){ // each direction
        v_relative[i] = mol[mol1_id].velocity_com[i] - mol[mol2_id].velocity_com[i];
        v_relative_square += v_relative[i] * v_relative[i];
    }
    relative_ke = 0.5 * mass_mol * v_relative_square * g_mol2kg * a_fs2m_s * a_fs2m_s / kcal_mol2j;
    return relative_ke;
}

// Calculate adjacency matrix accoding to Hill's energy criterion
void getAdjacencyMatrix_Hill(struct MOL *mol, int **adjacencyMatrix, double **peMatrix, int numMol) {
    int mol1_id, mol2_id;
    double ke, pe;
    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id + 1; mol2_id < numMol; mol2_id++) {
            adjacencyMatrix[mol1_id][mol2_id] = 0;
            adjacencyMatrix[mol2_id][mol1_id] = 0;
            pe = peMatrix[mol1_id][mol2_id];
            if (pe < 0) {
                ke = relativeKE(mol, mol1_id, mol2_id);
                if(ke + pe < 0) {
                    adjacencyMatrix[mol1_id][mol2_id] = 1;
                    adjacencyMatrix[mol2_id][mol1_id] = 1;
                }
            }
        }
    }
}

// Calculate adjacency matrix accoding to O-O distance cutoff
void getAdjacencyMatrix_distance(struct MOL *mol, int **adjacencyMatrix, int numMol, double boxLength, double cutoff) {
    int mol1_id, mol2_id;
    double dx, dy, dz, distance;

    for (mol1_id = 0; mol1_id < numMol; mol1_id++) {
        for (mol2_id = mol1_id + 1; mol2_id < numMol; mol2_id++) {
            adjacencyMatrix[mol1_id][mol2_id] = 0;
            adjacencyMatrix[mol2_id][mol1_id] = 0;

            dx = mol[mol1_id].position[0][0] - mol[mol2_id].position[0][0];
            dy = mol[mol1_id].position[0][1] - mol[mol2_id].position[0][1];
            dz = mol[mol1_id].position[0][2] - mol[mol2_id].position[0][2];
            dx = trueDistance(dx, boxLength);
            dy = trueDistance(dy, boxLength);
            dz = trueDistance(dz, boxLength);

            distance = sqrt(dx*dx + dy*dy + dz*dz);

            if (distance < cutoff) {
                adjacencyMatrix[mol1_id][mol2_id] = 1;
                adjacencyMatrix[mol2_id][mol1_id] = 1;
            }
        }
    }
}


// Energy per cluster as a function of cluster size
void energy_perCluster_size(struct MOL *mol, int numClusters, int *clusterSizes, int **clusters, double *PE_perCluster_size_accumulated, double *KE_perCluster_size_accumulated) {
    int i, cluster_id, size, molID;

    for (cluster_id = 0; cluster_id < numClusters; cluster_id++) {
        size = clusterSizes[cluster_id];
        for (i = 0; i < size; i++) {
            molID = clusters[cluster_id][i];
            mol[molID].clusterID = cluster_id;
            mol[molID].clusterSize = size;
            PE_perCluster_size_accumulated[size] += mol[molID].PE;
            KE_perCluster_size_accumulated[size] += mol[molID].KE;
        }
    }
}

void save_result(char *filename, int max_size, int numFrame, int *size_count_accumulated, int *degree_count_accumulated, double *PE_s_accumulated, double *KE_s_accumulated) {
    int size;
    FILE *fpw;
    fpw = fopen(filename, "w");
    double size_count_mean, degree_count_mean, PE_s_mean, KE_s_mean;

    // Write header
    fprintf(fpw, "number of frames = %d\n", numFrame);
    fprintf(fpw, "s\tN_c\tN_k\tPE(s)\tKE(s)\n");

    // Write data
    for (size = 0; size < max_size; size++) {
        size_count_mean = (1.0 * size_count_accumulated[size]) / (1.0 * numFrame);
        degree_count_mean = (1.0 * degree_count_accumulated[size]) / (1.0 * numFrame);
        PE_s_mean = PE_s_accumulated[size] / (1.0 * size_count_accumulated[size]);
        KE_s_mean = KE_s_accumulated[size] / (1.0 * size_count_accumulated[size]);

        fprintf(fpw, "%d\t%e\t%e\t%e\t%e\n", size, size_count_mean, degree_count_mean, PE_s_mean, KE_s_mean);
    }
    
    fclose(fpw);
}
