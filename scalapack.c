#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <f2c.h>
#include "blacs.h"
#include "scalapack.h"

/* a parameter that defines the global matrix size */
#define NTOT 100

/*
?* Generate a diagonal matrix whose eigenvalues are known.
?* The matrix has diagonal values of 1--NTOT
?* the off-diagonal values are 1E-16/sqrt(i*j), where i and j are the
?* global indices of the off-diagonal element
?* (expressly made non-zero to help avoid floating-point troubles)
?*
?* This function also distributes the matrix in a cyclic fashion
?*/
static int genmat(int ictxt, float **mat, int desc[]) {
? int nprow, npcol, myprow, mypcol;
? int i, j, ig, jg, blk, lld, mb, nb, irsrc, icsrc;
? float *tmp;

? /* get information about the BLACS grid */
? blacs_gridinfo(&ictxt, &nprow, &npcol, &myprow, &mypcol);

? /* allocate memory for the local part of the matrix */
? mb = NTOT / nprow + (NTOT%nprow<=myprow?0:1);
? nb = NTOT / npcol + (NTOT%npcol<=mypcol?0:1);
? lld = ((mb+3)/4)*4;
? *mat = (float*)malloc(sizeof(float)*lld*nb);
? if (*mat == (float*)0) {
??? i = fprintf(stderr, "Unable to allocate matrix on processor (%d, %d)\n",
?? ??? ?myprow, mypcol);
??? return -1;
? }

? /* set the descriptor for the matrix */
? j = NTOT; blk = 1; irsrc=0; icsrc=0;
? descinit(desc, &j, &j, &blk, &blk, &irsrc, &icsrc, &ictxt, &lld, &i);
? if (i != 0) return -2;

? /*
?? * assign values to the local part of the matrix according to the
?? * above distribution
?? */
? myprow++; mypcol++;
? for (j=0; j<nb; j++) {
??? jg = j*npcol + mypcol;
??? tmp = *mat + j*lld;
??? for (i=0; i<mb; i++) {
????? ig = i*nprow + myprow;
????? if (ig != jg) /* off-diagonal */
?? ?*tmp = 1.0e-16 / sqrt((double)(ig*jg));
????? else
?? ?*tmp = ig;
????? tmp++;
??? }
? }
? return 0;
} /* end of genmat */

int main(int argc, char *argv[]) {
? int i, mypnum, nproc, ictxt, nprow, npcol;
? int n, neig, ia, ja, il, iu, m, nz, iz, jz, lwork, liwork;
? int desc[10], *iwork, ifail[NTOT], *icluster;
? char porder='R', jobZ='V', range='I', uplo='U';
? float vl, vu, abstol, orfac, *mat, *Z, *work, *gap, eigval[NTOT];

? /* Initialize MPI */
? MPI_Init(&argc, &argv);

? /* setup BLACS grid */
? i = 0;
? blacs_setup(&mypnum, &nproc);
? blacs_get(&i, &i, &ictxt);

? nprow = sqrt((double) nproc);
? npcol = nproc / nprow;
? if (mypnum == 0) {
??? printf(" BLACS will setup %d X %d processor grid\n", nprow, npcol);
??? if (npcol*nprow != nproc)
????? printf(" BLACS grid (%d X %d) does not cover all processors (%d).\n",
?? ????? nprow, npcol, nproc);
??? fflush(NULL);
? }

? blacs_gridinit(&ictxt, &porder, &nprow, &npcol);

? /* generate the matrix */
? i = genmat(ictxt, &mat, desc);
? if (i != 0) blacs_abort(&ictxt, &i);

? /*
?? * parameters needed by pssyevx -- setup to compute the five smallest
?? * eigenvalues of the above generated symmetric matrix
?? * most array indices start with 1
?? */
? n = NTOT; ia = 1; ja = 1; iz = 1; jz = 1; il = 1; iu = 5; neig = -1;
? abstol = 0.0; orfac = 10.0;
? vl = 0.0; vu = 0.0;
? i = 4*n;
? liwork = (i>14)?i:14;
? i = 3*n + nproc + 1;
? liwork = ((liwork>i)? liwork : i) + 2 * n;
? icluster = (int*)malloc(sizeof(int) * 2*nprow*npcol);
? iwork = (int*)malloc(sizeof(int) * liwork);
? gap = (float*)malloc(sizeof(float) * nprow * npcol);
? Z = (float*)malloc(sizeof(float) * desc[8] * (1+n/npcol));
? lwork = 16*n;
? work = (float*)malloc(sizeof(float) * lwork);
? if (iwork==(int*)0 || icluster==(int*)0 || gap==(float*)0 ||
????? work==(float*)0) {
??? i = fprintf(stderr, "Unable to allocate workspace needed for pssyevx\n");
??? blacs_abort(&ictxt, &i);
? }

? /* call ScaLAPACK to compute the five smallest eigenvalues */
? i = 0;
? pssyevx(&jobZ, &range, &uplo, &n, mat, &ia, &ja, desc, &vl, &vu,
?? ?????? &il, &iu, &abstol, &neig, &nz, eigval, &orfac, Z, &iz, &jz,
?? ?????? desc,? work, &lwork, iwork, &liwork, ifail, icluster, gap,
?? ?????? &i, 1, 1, 1);

? if (i == -21 || i&4 == 4) { /* need more workspace */
??? lwork = -work[0];
??? free(work);
??? work = (float*)malloc(sizeof(float)*lwork);
??? if (work == (float*)0) {
????? i = fprintf(stderr, "Unable to reallocate work[%d]\n", lwork);
????? blacs_abort(&ictxt, &i);
??? }
??? pssyevx(&jobZ, &range, &uplo, &n, mat, &ia, &ja, desc, &vl, &vu,
?? ??? ?&il, &iu, &abstol, &neig, &nz, eigval, &orfac, Z, &iz,
?? ??? ?&jz, desc, work, &lwork, iwork, &liwork, ifail, icluster,
?? ??? ?gap, &i, 1, 1, 1);
? }

? /* output somethings from the first processor */
? if (mypnum == 0) {
??? iz = printf(" PSSYEV found %d eigenvalues (5 request), return code %d\n",
?? ??? ?neig, i);
??? if (neig > 0) {
????? for (i=0; i<neig; i++)
?? ?printf("E[%d] = %f, error=%e\n", i, eigval[i], eigval[i]-i-1);
??? }
??? fflush(NULL);
? }

? blacs_gridinfo(&ictxt, &nprow, &npcol, &ia, &ja);
? /* check the eigenvectors */
? if (mypnum == 0)
??? iz = printf(" PSSYEV found %d eigenvectors (5 request)\n", nz);
? if (nz > 0) {
??? for (i=0; i<nz; i++) {
????? /* Z(i,i) is located on processor (iz, jz) */
????? iz = i % nprow; jz = i % npcol;
????? if (iz==ia && jz==ja) {
?? ?vl = Z[i/nprow+desc[8]*(i/npcol)];
?? ?printf("Z(%d,%d) = %f, error = %e\n", i, i,
?? ??????? vl, (vl>0.0?vl-1.0:vl+1.0));
????? }
??? }
??? fflush(NULL);
? }

? /* cleaning up */
? free(work); free(gap); free(iwork); free(icluster); free(mat);
? i = 0;
? blacs_exit(&i); /* calls blacs_gridexit() and MPI_Finalize() */
? return i;
} /* end of main */
 