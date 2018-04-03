#ifndef _OROPRECIPLIB_H
#define MKL_ILP64   // Use 64-bit integers in MKL

#define _OROPRECIPLIB_H

#define IY 0
#define IX 1

#define EQC 0
#define EQS 1

#define NEQ 2

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <mkl.h>
#include <mkl_cluster_sparse_solver.h>


/* Thermodynamic parameters of the precipitation model */
struct cfg_thermodyn {
    double rho_ref;    // Density at reference level; kg m^-3
    double qsat_ref;   // Moisture saturation at reference level; kg kg^-1
    double Hm;         // Scale height for moisture content; km
};

/* Configuration parameters of the precipitation model */
struct cfg_options {
    int nx, ny;            // Resolution
    double dx, dy;         // Grid spacing
    double Lx, Ly;         // Domain size
    double U[2];           // Wind speed, y-, x-dir
    double tauc, tauf;     // Moisture conversion / precipitation time delays
    double Sbg;            // Background precipitation source
    int sink_at_downslope; // Allow negative precipitation at downslope side?
    int downsample;        // downsample topography by this value
    double z0;             // reference level (sealevel)
};


/* CSR (3 arr) storage container for the coefficient matrix,
 * for Intel MKL sparse solvers */
struct SCSR3_t {
    MKL_INT nval;       // num of stored values (len of values, columns)
    MKL_INT n, m;       // size of global matrix (n x m)
    double *values;
    MKL_INT *columns;
    MKL_INT *rowIndex;
    MKL_INT allocsize;  // allocated size for values, columns
    MKL_INT *perm;
    int metaonly;
    MKL_INT blksize;
};

/* Routines related to Intel MKL sparse solver and
 * matrix storage */
void free_SCSR3(struct SCSR3_t **scsr3);
void printformat_SCSR3(struct SCSR3_t *s);
void print_SCSR3(struct SCSR3_t *scsr3);
int setdiags_SCSR3(struct SCSR3_t *s, MKL_INT nd, MKL_INT *ld, MKL_INT halfwidth);
int create_SCSR3(struct SCSR3_t **ptrscsr3, MKL_INT n, MKL_INT m, MKL_INT blksize, int metaonly);
void shiftdata_SCSR3(struct SCSR3_t *scsr3, MKL_INT loc);
void shiftrows_SCSR3(struct SCSR3_t *scsr3, MKL_INT from);
void addval_SCSR3(struct SCSR3_t *scsr3, MKL_INT loc, MKL_INT icol, MKL_INT irow, double val);
int setval_SCSR3(struct SCSR3_t *scsr3, MKL_INT irow, MKL_INT icol, double val);
MKL_INT * get_iparm();
int solve_sparse_mkl(struct SCSR3_t *s, double *b, double *x, int cleanup);

/* Oroprecip related routines */
MKL_INT gidx(int i, int j, int eq, int nx);
void oroprecip(double *h, struct cfg_options *co, struct cfg_thermodyn *ct, double *p);


#endif
