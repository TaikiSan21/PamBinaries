#include <R.h>

#ifdef USING_R
#include "Rinternals.h"
#include "R_ext/Lapack.h"

#include "sapa_lapack.h"

SEXP R_sapa_dpss(SEXP nsample, SEXP ntaper, SEXP bandwidth)
{
  int     Ntaper = INTEGER(ntaper)[0];
  int     N = INTEGER(nsample)[0];
  double  NW = REAL(bandwidth)[0];

  double *D;
  double *E;
  double *W;
  int    *IBLOCK;
  int    *ISPLIT;
  double *Z;
  double *WORK1;
  double *WORK2;
  int    *IWORK1;
  int    *IWORK2;
  int    *IFAIL;
  int     INFO;
  double  fac;
  double  fac2;
  int     i;
  int     j;
  double *pd_taper;

  SEXP    tapers;
  // SEXP    dimnames;
  double *tapervals;

  double  vl = 1.0;
  double  vu = 1.0;
  double  abstol = 1.0;
  int     LDZ = N;
  int     M1 = 1;
  int     M2 = Ntaper;
  int     il = N - Ntaper + 1;
  int     iu = N;
  int     nsplit = 1;

  /* allocate memory */

  D     = (double *) R_alloc(N, sizeof(double));
  E     = (double *) R_alloc(N-1, sizeof(double));
  W     = (double *) R_alloc(N, sizeof(double));
  Z     = (double *) R_alloc(LDZ * M2, sizeof(double));
  WORK1 = (double *) R_alloc(4*N, sizeof(double));
  WORK2 = (double *) R_alloc(5*N, sizeof(double));

  IBLOCK = (int *) R_alloc(N, sizeof(int));
  ISPLIT = (int *) R_alloc(N, sizeof(int));
  IWORK1 = (int *) R_alloc(3*N, sizeof(int));
  IWORK2 = (int *) R_alloc(N, sizeof(int));
  IFAIL  = (int *) R_alloc(M2, sizeof(int));

  /* initialize variables */

  fac = (2.0 * M_PI * NW)/ (double) N;

  /* D: generate N diagonal elements of symmetric tridiagonal matrix

     E: generate N-1 off-diagonal elements; note that actually N elements
     are generated here; only the first N-1 are used by dstebz,
     but the documentation for dstein indicates that it needs
     a vector with N elements, the first N-1 of which are the
     off-diagonal elements. */

  for (i = 0; i < N; i++){

    fac2 = (double) (N - 1 - 2 * i) / 2.0;
    D[i] = cos(fac) * fac2 * fac2;
    W[i] = 0.0;
    IBLOCK[i] = 0;
    ISPLIT[i] = 0;

    if (i < N - 1)
      E[i] = (double) ((i + 1) * (N - i - 1)) / 2.0;
  }

  /* DSTEBZ - compute the eigenvalues of a symmetric tridiagonal matrix */
  F77_CALL(dstebz)(
    (char *) "I",
    (char *) "B",
    &N,
    &vl,
    &vu,
    &il,
    &iu,
    &abstol,
    D,
    E,
    &M1,
    &nsplit,
    W,
    IBLOCK,
    ISPLIT,
    WORK1,
    IWORK1,
    &INFO);
  //  if (INFO != 0)
  //  error(_("error code %d from LAPACK routine '%s'"), info, "dstebz");

  /* DSTEIN - compute the eigenvectors of a real symmetric        */
  /* tridiagonal matrix T corresponding to specified eigenvalues, */
  /* using inverse iteration                                      */
  F77_CALL(dstein)(
    &N,
    D,
    E,
    &Ntaper,
    W,
    IBLOCK,
    ISPLIT,
    Z,
    &N,
    WORK2,
    IWORK2,
    IFAIL,
    &INFO);
  // if (INFO != 0)
  //  error(_("error code %d from LAPACK routine '%s'"), info, "dstein");

  PROTECT(tapers = allocMatrix(REALSXP, N, Ntaper));
  tapervals = REAL(tapers);
  pd_taper = Z;

  /* fill the taper matrix */
  for (i = 0; i < Ntaper; i++){

    for(j = 0; j < N; j++){

      *tapervals = *pd_taper;
      tapervals++;
      pd_taper++;
    }
  }

  // PROTECT(dimnames = allocVector(VECSXP, 2));
  // SET_VECTOR_ELT(dimnames, 1, install(""));
  // dimnamesgets(gradient, dimnames);

  UNPROTECT(1);
  return tapers;
}
#endif
