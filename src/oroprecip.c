#include <oropreciplib.h>

MKL_INT gidx(int i, int j, int eq, int nx) {
    return (MKL_INT)NEQ*((MKL_INT)i*(MKL_INT)nx + (MKL_INT)j) + (MKL_INT)eq;
	//return NEQ * (i*nx + j) + eq;
}

void oroprecip(double *h, struct cfg_options *co, struct cfg_thermodyn *ct, double *p) {
	/*
	 * Solves numerically the orographic precipitation as described in
	 * Smith and Barstad 2004 (Journal of the American Meteorological Society)
	 *
	 *  - "No dynamics" version
	 *  - S decays as function of z (relevant for high topography),
	 *    w is vertically constant
	 *
	 * Input:
	 * 		h	Topography, size: co.nx * co.ny
	 * 		co	Configuration options
	 * 		ct	Thermodynamic options
	 *
	 * 	Output:
	 * 		p	Precipitation, size: co.nx * co.ny
	 *
	 * 	Input and output is in row-major format.
	 *
	 *
	 * 	struct cfg_options *co has following fields:
	 * 		nx, ny [int]	Num of grid points in x- ja y-dir
	 * 		dx, dy [double]	Grid spacing
	 * 		Lx, Ly [double]	Domain size
	 * 		                (make sure that (nx-1)*dx == Lx)
	 * 		U [double[2]]	Wind speed in y-, x-dir
	 * 		tauc [double]	Conversion delay 
	 * 		tauf [double]	Rainfall dealy
	 * 		Sbg [double]	Background moisture source value
	 * 		sink_at_downslope [int]
	 * 						0: Source term may become negative at downslope
	 * 						1: Source term always >= 0
	 *
	 * 	struct cfg_thermodyn *ct has following fields:
	 * 		rho_ref [double]	Density of air at reference level
	 * 		qsat_ref [double]	Saturation moisture content at reference level
	 * 		Hm [double]			Scale height for moisture content
	 *
	 *
	 * TODO:
	 *
	 * 	- Coefficient matrix need not to be build every time, make it static
	 */
	int iproc, nproc;

	struct SCSR3_t *A;
	double *x, *b;
	const int And = 6;
	MKL_INT Ald[6];

    int nx, ny;
	int N, dof;
    double dx, dy;
    double *qc, *qs, *S, *grh;
    double U[2];
	double tauc, tauf;
    double ddx, ddy;

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

	/* Shorthand notation */
	nx = co->nx; ny = co->ny;
	dx = co->dx; dy = co->dy;
	U[IX] = co->U[IX]; U[IY] = co->U[IY];
	tauc = co->tauc; tauf = co->tauf;

    /* Calculate derived config variables */
    N = nx*ny;
    dof = NEQ*N;
    ddx = 2.0 * dx;
    ddy = 2.0 * dy;

    /* Set up fields */
    qc  = (double *)malloc(sizeof(double)*N);
    qs  = (double *)malloc(sizeof(double)*N);
    S   = (double *)malloc(sizeof(double)*N);
    grh = (double *)malloc(sizeof(double)*2*N);

	assert (qs != NULL);
	assert (qc != NULL);
	assert (S != NULL);
	assert (grh != NULL);

    /* Set up matrices/vectors */
	A = NULL;
	if (iproc == 0) {
		create_SCSR3(&A, dof, dof, 5*dof, 0);
	} else {
		create_SCSR3(&A, dof, dof, 5*dof, 1);
	}
	Ald[0] = -NEQ*nx;
	Ald[1] = -NEQ;
	Ald[2] = -1; // EQC-EQS
	Ald[3] = 0;
	Ald[4] = NEQ;
	Ald[5] = NEQ*nx;
	setdiags_SCSR3(A, And, Ald, 3);
	b = calloc(sizeof(double), dof);
    x = calloc(sizeof(double), dof);

	assert (b != NULL);
	assert (x != NULL);

	/* Calculate grad h */
    for (int i = 1; i < ny-1; i++) {
        for (int j = 1; j < nx-1; j++) {
			grh[2*(i*nx + j) + IY] = (h[(i+1)*nx+j] - h[(i-1)*nx+j]) / ddy;
			grh[2*(i*nx + j) + IX] = (h[i*nx+j+1] - h[i*nx+j-1]) / ddx;
        }
		grh[2*(i*nx + 0   ) + IY] = (h[(i+1)*nx+0] - h[(i-1)*nx+0]) / ddy;
		grh[2*(i*nx + 0   ) + IX] = (h[i*nx+1] - h[i*nx+0]) / dx;               // gives some error...
		grh[2*(i*nx + nx-1) + IY] = (h[(i+1)*nx+nx-1] - h[(i-1)*nx+nx-1]) / ddy;
		grh[2*(i*nx + nx-1) + IX] = (h[i*nx+nx-1] - h[i*nx+nx-2]) / dx;         // ... this one, too
    }
	for (int j = 1; j < nx-1; j++) {
		grh[2*(0*nx + j) + IY] = (h[1*nx+j] - h[0*nx+j]) / dy;                  // gives some error ...
		grh[2*(0*nx + j) + IX] = (h[0*nx+j+1] - h[0*nx+j-1]) / ddx;
		grh[2*((ny-1)*nx + j) + IY] = (h[(ny-1)*nx+j] - h[(ny-2)*nx+j]) / dy;             // ... this one, too
		grh[2*((ny-1)*nx + j) + IX] = (h[(ny-1)*nx+j+1] - h[(ny-1)*nx+j-1]) / ddx;
	}
	grh[2*(0*nx + 0) + IY] = grh[2*(0*nx + 1) + IY];
	grh[2*(0*nx + 0) + IX] = grh[2*(1*nx + 0) + IY];
	grh[2*(0*nx + nx-1) + IY] = grh[2*(0*nx + nx-2) + IY];
	grh[2*(0*nx + nx-1) + IX] = grh[2*(1*nx + nx-1) + IX];
	grh[2*((ny-1)*nx + 0) + IY] = grh[2*((ny-1)*nx + 1) + IY];
	grh[2*((ny-1)*nx + 0) + IX] = grh[2*((ny-2)*nx + 0) + IX];
	grh[2*((ny-1)*nx + nx-1) + IY] = grh[2*((ny-1)*nx + nx-2) + IY];
	grh[2*((ny-1)*nx + nx-1) + IX] = grh[2*((ny-2)*nx + nx-1) + IY];

    /* Calculate water source based on terrain gradient and background source */
    for (int i = 1; i < ny-1; i++) {
        for (int j = 1; j < nx-1; j++) {
            S[i*nx + j] = (U[IY] * grh[2*(i*nx + j) + IY] + U[IX] * grh[2*(i*nx + j) + IX]);
            S[i*nx + j] *= ct->rho_ref * ct->qsat_ref * exp(-h[i*nx + j] / ct->Hm);
			S[i*nx + j] += co->Sbg;
			if (co->sink_at_downslope == 0 && S[i*nx +j] < 0.0) S[i*nx +j] = 0.0;
        }
    }


	/* Form coefficient matrix */

    /* Coefs, main grid points */
	for (int i = 0; i < dof; i++) b[i] = 0.0;

    for (int i = 1; i < ny-1; i++) {
        for (int j = 1; j < nx-1; j++) {
			setval_SCSR3(A, gidx(i, j, EQC, nx), gidx(i-1, j  , EQC, nx), -tauc * U[IY] / ddy);
			setval_SCSR3(A, gidx(i, j, EQC, nx), gidx(i  , j-1, EQC, nx), -tauc * U[IX] / ddx);
			setval_SCSR3(A, gidx(i, j, EQC, nx), gidx(i  , j  , EQC, nx), tauc * 1.0 / tauc);
			setval_SCSR3(A, gidx(i, j, EQC, nx), gidx(i+1, j  , EQC, nx), tauc * U[IY] / ddy);
			setval_SCSR3(A, gidx(i, j, EQC, nx), gidx(i  , j+1, EQC, nx), tauc * U[IX] / ddx);
            b[gidx(i, j, EQC, nx)] = tauc * S[i*nx + j];
            
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i-1, j  , EQS, nx), -tauc * U[IY] / ddy);
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i  , j-1, EQS, nx), -tauc * U[IX] / ddx);
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i  , j  , EQS, nx), tauc * 1.0 / tauf);
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i+1, j  , EQS, nx), tauc * U[IY] / ddy);
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i  , j+1, EQS, nx), tauc * U[IX] / ddx);
			setval_SCSR3(A, gidx(i, j, EQS, nx), gidx(i  , j  , EQC, nx), -tauc * 1.0 / tauc);
            b[gidx(i, j, EQS, nx)] = tauc * 0.0;
        }
    }

    /* B.C. for EQC
     * natural */
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < 2; j++) {
			setval_SCSR3(A, gidx(i, j*(nx-1), EQC, nx), gidx(i, j*(nx-1)  ,         EQC, nx),  1.0);
			setval_SCSR3(A, gidx(i, j*(nx-1), EQC, nx), gidx(i, j*(nx-1) + (1-2*j), EQC, nx), -1.0);
            b[gidx(i, j*(nx-1), EQC, nx)] = 0.0;
        }
    }

    for (int j = 0; j < nx; j++) {
        for (int i = 0; i < 2; i++) {
			setval_SCSR3(A, gidx(i*(ny-1), j, EQC, nx), gidx(i*(ny-1)          , j, EQC, nx),  1.0);
			setval_SCSR3(A, gidx(i*(ny-1), j, EQC, nx), gidx(i*(ny-1) + (1-2*i), j, EQC, nx), -1.0);
            b[gidx(i*(ny-1), j, EQC, nx)] = 0.0;
        }
    }

    /* B.C. for EQS
     * natural */
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < 2; j++) {
			setval_SCSR3(A, gidx(i, j*(nx-1), EQS, nx), gidx(i, j*(nx-1)  ,         EQS, nx),  1.0);
			setval_SCSR3(A, gidx(i, j*(nx-1), EQS, nx), gidx(i, j*(nx-1) + (1-2*j), EQS, nx), -1.0);
            b[gidx(i, j*(nx-1), EQS, nx)] = 0.0;
        }
    }

    for (int j = 0; j < nx; j++) {
        for (int i = 0; i < 2; i++) {
			setval_SCSR3(A, gidx(i*(ny-1), j, EQS, nx), gidx(i*(ny-1)          , j, EQS, nx),  1.0);
			setval_SCSR3(A, gidx(i*(ny-1), j, EQS, nx), gidx(i*(ny-1) + (1-2*i), j, EQS, nx), -1.0);
            b[gidx(i*(ny-1), j, EQS, nx)] = 0.0;
        }
    }


	/* Solve it! */

	solve_sparse_mkl(A, b, x, 0);

	/* Transform solution vector to qc and qs components */

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            qc[i*nx + j] = x[NEQ*(i*nx + j) + EQC];
            qs[i*nx + j] = x[NEQ*(i*nx + j) + EQS];
			p[i*nx + j] = qs[i*nx + j] / tauf;
        }
    }

	free(qs);
	free(qc);
	free(S);
	free(grh);
	free(b);
	free(x);
	free_SCSR3(&A);
}


#if OROPRECIP_STANDALONE == 1
int main(int argc, char **argv) {
	int iproc, nproc;
	FILE *fp;
	struct cfg_thermodyn ct;
	struct cfg_options co;
	double *h, *p;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

	/* ======================= */
    /* Read config             */
	/* ======================= */

	co.nx = 512;
	co.ny = 512;
	co.Lx = 250e3;
	co.Ly = 300e3;
	co.dx = co.Lx / (double)(co.nx-1);
	co.dy = co.Ly / (double)(co.ny-1);

    co.U[IY] = 0.0; 
    co.U[IX] = 15.0; 
	co.tauc = 1000;
	co.tauf = 1000;
	co.Sbg = 0.0;
	co.sink_at_downslope = 0;

	ct.rho_ref = 1.2;
	ct.qsat_ref = 8e-3;
	ct.Hm = 2500;

	/* ======================= */
	/* End of Read config      */
	/* ======================= */

    /* Form terrain */
	h = calloc(sizeof(double), co.nx*co.ny);
    for (int i = 0; i < co.ny; i++) {
        for (int j = 0; j < co.nx; j++) {
            double x, y, r;
			y = co.Ly * (double)i / (double)(co.ny-1) - 0.5*co.Ly;
			x = co.Lx * (double)j / (double)(co.nx-1) - 0.5*co.Lx;
			h[i*co.nx + j] = 500.0 * exp(-pow(x/35e3,2.0)-pow(y/35e3,2.0));
        }
    }

	p = calloc(sizeof(double), co.nx*co.ny);
	oroprecip(h, &co, &ct, p);

	if (iproc == 0) {
		fp = NULL;
		fp = fopen("outnumtopo.txt", "w");
		if (fp == NULL) {
			fprintf(stderr, "ERROR opening file for writing\n");
			exit(1);
		}
		for (int i = 0; i < co.ny; i++) {
			for (int j = 0; j < co.nx; j++) {
				fprintf(fp, "%d,%d,0,%e\n", i,j,h[i*co.nx + j]);
			}
		}
		fclose(fp);

		fp = NULL;
		fp = fopen("outnumprecip.txt", "w");
		if (fp == NULL) {
			fprintf(stderr, "ERROR opening file for writing\n");
			exit(1);
		}
		for (int i = 0; i < co.ny; i++) {
			for (int j = 0; j < co.nx; j++) {
				fprintf(fp, "%d,%d,0,%e\n", i,j,p[i*co.nx + j]);
			}
		}
		fclose(fp);
	}

	free(p);

	MPI_Finalize();
	return 0;
}
#endif


