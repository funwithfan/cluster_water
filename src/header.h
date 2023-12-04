#ifndef __COMMON_H__
#define __COMMON_H__
struct MOL {
    double position[3][3];
    double velocity[3][3];
    //double com[3];
    double velocity_com[3];
    double PE;
    double KE;
    int clusterID;
    int clusterSize;
};

// unit conversion
double NA;
double kcal_mol2j; 
double A2m;
double e2c;
double a_fs2m_s;
double g_mol2kg;

// force field parameters
double mass_mol;
double rc;
double mass[3];
double eps;
double sig;
double charge[3];

#endif


