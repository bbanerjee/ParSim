#include <stdio.h>
#include <stdlib.h>


void diamm_calc_(int* nblk, int* ninsv, double* dt, double* ui, 
                 double* stress, double* d, double* sv, double* usm);

// MIG utilities stubs for standalone test
/*
void faterr_(char* mes1, char* mes2, int len_mes1, int len_mes2) {
    printf("FATAL ERROR: %.*s %.*s ", len_mes1, mes1, len_mes2, mes2);
    exit(1);
}
void bombed_(char* mes, int len_mes) {
    printf("BOMBED: %.*s ", len_mes, mes);
    exit(1);
}
void logmes_(char* mes, int len_mes) {
    printf("LOG: %.*s ", len_mes, mes);
}
*/

int main() {
    int nblk = 1;
    int ninsv = 26;
    double dt = 1.0e-6;
    double ui[47] = {0};
    double stress[6] = {10.0e6, -5.0e6, 2.0e6, 1.0e6, 0.0, 0.0};
    double d[6] = {1000.0, -500.0, 200.0, 100.0, 0.0, 0.0}; // strain rate
    double sv[26] = {0};
    double usm = 0.0;

    // Set some properties
    ui[0] = 100.0e9; // B0
    ui[3] = 75.0e9;  // G0
    ui[14] = 1000.0; // R0
    ui[15] = 294.0;  // T0
    sv[7] = 1000.0;  // Density
    sv[9] = 1.0;     // Jacobian

    printf("Running Fortran diamm_calc...  ");
    diamm_calc_(&nblk, &ninsv, &dt, ui, stress, d, sv, &usm);
    printf("Success. USM = %g ", usm);

    return 0;
}
