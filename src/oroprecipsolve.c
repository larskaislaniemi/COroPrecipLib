#include <oropreciplib.h>

void free_SCSR3(struct SCSR3_t **scsr3) {
    struct SCSR3_t *s;
    assert(*scsr3 != NULL);

    s = *scsr3;
    if (s->values != NULL) free(s->values);
    if (s->columns != NULL) free(s->columns);
    if (s->rowIndex != NULL) free(s->rowIndex);
    if (s->perm != NULL) free(s->perm);
    free(*scsr3);
}

void printformat_SCSR3(struct SCSR3_t *s) {
    double *A;
    A = calloc(sizeof(double), s->n * s->m);

    for (MKL_INT i = 0; i < s->n; i++) {
        printf("%ld:\t%ld\t", i, s->rowIndex[i+1] - s->rowIndex[i]);
        for (MKL_INT j = s->rowIndex[i]; j < s->rowIndex[i+1]; j++) {
            printf("%ld ", s->columns[j]);
        }
        printf("\n");
    }

    printf("\n");
}

void print_SCSR3(struct SCSR3_t *scsr3) {
    printf("n, m: %ld, %ld; allocsize: %ld; nval: %ld\n", scsr3->n, scsr3->m, scsr3->allocsize, scsr3->nval);
    printf("Values: ");
    for (int i = 0; i < scsr3->nval; i++) printf("%f ", scsr3->values[i]);
    printf("\nColumns: ");
    for (int i = 0; i < scsr3->nval; i++) printf("%ld ", scsr3->columns[i]);
    printf("\nrowIndex: ");
    for (int i = 0; i <= scsr3->n; i++) printf("%ld ", scsr3->rowIndex[i]);
    printf("\n\n");
}

int setdiags_SCSR3(struct SCSR3_t *s, MKL_INT nd, MKL_INT *ld, MKL_INT halfwidth) {
    // fill matrix with given diagonals (zero values)
    // nd: total number of diagonals
    // ld: list of diagonals; len: nd
    // halfwidth: halfwidth of the band, excluding main diagonal (nl+1+nu)
    //      if unsymmetric width, then greater of the two
    //  nb! overwrites all existing values
    MKL_INT nval;
    MKL_INT addcount;

    if (s->metaonly) return 0;

    nval = 0;
    for (MKL_INT i = 0; i < nd; i++) {
        nval = nval + s->n - abs(ld[i]);
    }

    s->allocsize = nval;
    s->nval = nval;
    if (s->values != NULL) free(s->values);
    if (s->columns != NULL) free(s->columns);
    s->values = calloc(nval, sizeof(double));
    s->columns = calloc(nval, sizeof(MKL_INT));

    assert(s->values != NULL);
    assert(s->columns != NULL);

    addcount = 0;
    for (MKL_INT i = 0; i < s->n; i++) {
        s->rowIndex[i+1] = s->rowIndex[i];
        for (MKL_INT id = 0; id < nd; id++) {
            MKL_INT j = i + ld[id];
            if (j >= 0 && j < s->n) {
                addcount += 1;
                s->columns[s->rowIndex[i+1]] = j;
                s->rowIndex[i+1] += 1;
            }
        }
    }

    assert(addcount == nval);

    return 0;
}

int create_SCSR3(struct SCSR3_t **ptrscsr3, MKL_INT n, MKL_INT m, MKL_INT blksize, int metaonly) {
    /* Reserve memory for a matrix in CSR3 format, single precision;
     */
    struct SCSR3_t *scsr3;

    assert (*ptrscsr3 == NULL);
    assert (blksize > 0);
    assert (n > 0); 
    assert (m > 0);

    *ptrscsr3 = calloc(sizeof(struct SCSR3_t), 1);
    scsr3 = *ptrscsr3;

    scsr3->nval = 0;
    scsr3->n = n;
    scsr3->m = m;
    scsr3->metaonly = metaonly;
    scsr3->blksize = blksize;

    if (metaonly) {
        scsr3->values = NULL;
        scsr3->columns = scsr3->rowIndex = scsr3->perm = NULL;

        scsr3->allocsize = 0;
    } else {
        scsr3->values = calloc(sizeof(double), scsr3->blksize);
        scsr3->columns = calloc(sizeof(MKL_INT), scsr3->blksize);
        scsr3->rowIndex = calloc(sizeof(MKL_INT), scsr3->n+1);
        scsr3->perm = calloc(sizeof(MKL_INT), scsr3->n);

        assert(scsr3->values != NULL);
        assert(scsr3->columns != NULL);
        assert(scsr3->rowIndex != NULL);
        assert(scsr3->perm != NULL);

        scsr3->allocsize = scsr3->blksize;
    }

    return 0;
}

void shiftdata_SCSR3(struct SCSR3_t *scsr3, MKL_INT loc) {
    // add space for a new value in the values/columns array
    // by shifting existing values to the right. 
    // realloc more memory if needed.
    
    if (scsr3->metaonly) return;

    while (scsr3->allocsize <= scsr3->nval) {
        scsr3->values = realloc(scsr3->values, sizeof(double)*(scsr3->allocsize + scsr3->blksize));
        scsr3->columns = realloc(scsr3->columns, sizeof(MKL_INT)*(scsr3->allocsize + scsr3->blksize));

        assert(scsr3->values != NULL);
        assert(scsr3->columns != NULL);

        scsr3->allocsize += scsr3->blksize;
    }

    for (MKL_INT i = scsr3->nval-1; i >= loc; i--) {
        scsr3->values[i+1] = scsr3->values[i];
        scsr3->columns[i+1] = scsr3->columns[i];
    }

    //scsr3->values[loc] = 0.0;
    //scsr3->columns[loc] = 0;

    scsr3->nval += 1;
}

void shiftrows_SCSR3(struct SCSR3_t *scsr3, MKL_INT from) {
    // "shift" rows, i.e. add +1 to all row indexes after given row.
    // used after adding a value/column on a row.
    
    if (scsr3->metaonly) return;

    for (MKL_INT i = from; i < scsr3->n+1; i++) scsr3->rowIndex[i] += 1;
    assert (scsr3->rowIndex[scsr3->n] == scsr3->nval);
}

void addval_SCSR3(struct SCSR3_t *scsr3, MKL_INT loc, MKL_INT icol, MKL_INT irow, double val) {
    // add a value/column at given location and readjust row indices
    if (scsr3->metaonly) return;

    shiftdata_SCSR3(scsr3, loc);
    scsr3->columns[loc] = icol;
    scsr3->values[loc] = val;
    shiftrows_SCSR3(scsr3, irow+1);
}

int checkrows(struct SCSR3_t *s) {
    for (MKL_INT i = 1; i < s->n+1; i++) {
        if (s->rowIndex[i] < s->rowIndex[i-1]) return 1;
    }
    return 0;
}

int setval_SCSR3(struct SCSR3_t *scsr3, MKL_INT irow, MKL_INT icol, double val) {
    // set value of the array at given row/col
    int hasval;

    if (scsr3->metaonly) return 0;

    assert (irow >= 0 && irow < scsr3->n);
    assert (icol >= 0 && irow < scsr3->n);

    hasval = 0;

    if (scsr3->rowIndex[irow] == scsr3->rowIndex[irow+1]) {
        // this row has no values at the moment
        //printf("NO ROW\n");
        MKL_INT loc = scsr3->rowIndex[irow];
        addval_SCSR3(scsr3, loc, icol, irow, val);
        hasval = 1;
    } else {
        // row exists, check if column is already filled
        assert (scsr3->rowIndex[irow] < scsr3->rowIndex[irow+1]);  // TODO: FAILS if reallocs done in shiftdata (why?!)
        for (MKL_INT loc = scsr3->rowIndex[irow]; loc < scsr3->rowIndex[irow+1]; loc++) {
            if (scsr3->columns[loc] == icol) {
                // column already has a value, replace
                scsr3->values[loc] = val;
                hasval = 1;
                break;
            } else if (scsr3->columns[loc] > icol) {
                // column does not exist, a column to the right was found
                printf("NO COL\n");
                addval_SCSR3(scsr3, loc, icol, irow, val);
                hasval = 1;
                break;
            } else if (loc == scsr3->rowIndex[irow+1]-1) {
                printf("COL AT END\n");
                // column was not found before the end of the row
                addval_SCSR3(scsr3, loc+1, icol, irow, val);
                hasval = 1;
                break;
            } 
        }
    }

    assert(hasval == 1);

    return 0;
}

MKL_INT * get_iparm() {
    MKL_INT * iparm;
    
    iparm = calloc(sizeof(MKL_INT), 64);
    assert(iparm != NULL);

    iparm[0]  = 1;  // non-default values
    iparm[1]  = 2;  // METIS
    iparm[2]  = 1;  // n threads(?)
    iparm[4]  = 0;  // no user permutation
    iparm[5]  = 0;  // solution goes to x
    iparm[7]  = 2;  // max num of refinements
    iparm[9]  = 13; // pivoting perturbation default
    iparm[10] = 1;  // enable scaling
    iparm[12] = 1;  // enable matching
    iparm[17] = -1; // do  output num of nonzeros LU
    iparm[18] = -1; // do  output: Mflops of factorization
    iparm[26] = 1;  // do  check the sparse matrix representation for errors
    iparm[27] = 0;  // DOUBLE precision
    iparm[34] = 1;  // ZERO BASED indexing
    iparm[36] = 0;  // CSR matrix format (three array variation)
    iparm[39] = 0;  // data stored on master

    return iparm;
}

int solve_sparse_mkl(struct SCSR3_t *s, double *b, double *x, int cleanup) {
    /* State */
    static int firstrun = 1;
    static void *pt[64] = { 0 };

    /* MPI */
    int iproc, nproc;
    MPI_Comm comm;

    /* MKL */
    MKL_INT *iparm;
    const MKL_INT mtype = 11;
    const MKL_INT nrhs = 1;
    const MKL_INT msglvl = 1;
    const MKL_INT maxfct = 1;
    const MKL_INT mnum = 1;
    MKL_INT phase;
    MKL_INT ierr;

    assert (s != NULL);

    //mkl_verbose(1);
    mkl_set_num_threads(1);  // Use pure MPI
    iparm = get_iparm();

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    comm = MPI_Comm_c2f(MPI_COMM_WORLD);

    if (firstrun) {
        phase = 12;
        //cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &(s->n), s->values, s->rowIndex, s->columns,
        //  s->perm, &nrhs, iparm, &msglvl, b, x, &comm, &ierr);
        firstrun = 0;
    }

    if (cleanup) {
        phase = -1;
        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &(s->n), s->values, s->rowIndex, s->columns,
            s->perm, &nrhs, iparm, &msglvl, b, x, &comm, &ierr);
    } else {
        phase = 13;
        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &(s->n), s->values, s->rowIndex, s->columns,
            s->perm, &nrhs, iparm, &msglvl, b, x, &comm, &ierr);
    }

    return ierr;
}

