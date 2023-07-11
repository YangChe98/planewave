#include "phg.h"
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <phg/map.h>
#include <phg/matvec.h>
#include <phg/phg-mpi.h>
#include <phg/utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#define check phgPrintf("HELLO I AM HERE LINE=%d\n", __LINE__)
#define Pi 3.14159265359
#define h_bar 1
int min(int a, int b) {
  if (a < b) {
    return a;
  } else {
    return b;
  }
}
void matrix_expontenial(MAT *A, VEC *v, FLOAT deltat, int m, int mmax, int smax,
                        double tol, VEC **Av) {
  MAT *Adeltat;
  phgMatAXPBY(deltat, A, 0.0, &Adeltat);
  VEC *v1;
  phgMatVec(0.0, 1.0, Adeltat, v, 0.0, &v1);
  int j;
  for (j = 2; j < (m + 2); j++) {
    phgMatVec(0.0, 1.0, Adeltat, v1, 0.0, &v1);
  }
  int s;
  double nfac[mmax];
  nfac[0] = 1;
  for (j = 1; j < mmax; j++) {

    nfac[j] = nfac[j - 1] / (j + 1);
  }
  s = ceil(pow(phgVecNorm1(v1, 0, NULL) * nfac[m] / tol, 1. / (m + 1)));

  double p = m * s;
  double s1, p1;

  int f = 0;
  while (f == 0 && m < mmax) {
    m++;
    phgMatVec(0.0, 1.0, Adeltat, v1, 0.0, &v1);
    s1 = ceil(pow(phgVecNorm1(v1, 0, NULL) * nfac[m] / tol, 1. / (m + 1)));
    p1 = m * s1;
    if (p1 <= p) {
      p = p1;
      s = s1;
    } else {
      m = m - 1;
      f = 1;
    }
  }
  s = min(s, smax);
  p = m * s;

  VEC *w;
  phgVecCopy(v, &w);
  phgMatVec(0.0, 1.0, Adeltat, v, 0.0, &v1);
  for (j = 0; j < m; j++) {
    phgMatVec(0.0, 1.0, Adeltat, v1, 0.0, &v1);
    phgVecAXPBY(nfac[j] / pow(s, (j + 1)), v1, 1., &w);
  }
  MAT *Adeltats;
  phgMatAXPBY(1. / s, Adeltat, 0.0, &Adeltats);
  int i;
  VEC *v2;
  for (i = 0; i < (s - 1); i++) {
    phgVecCopy(w, &v2);
    for (j = 0; j < m; j++) {
      phgMatVec(0.0, 1.0, Adeltats, v2, 0.0, &v2);
      phgVecAXPBY(nfac[j], v2, 1.0, &w);
    }
  }
  phgVecCopy(w, Av);
}

int main(int argc, char *argv[]) {

  INT Nbasis = 100;

  phgInit(&argc, &argv);
  // phgPrintf("%e \n", phgNProcs);
  INT nprocs = phgNProcs;
  INT vec1procs = 0;
  if (Nbasis * Nbasis * Nbasis % nprocs == 0) {
    vec1procs = Nbasis * Nbasis * Nbasis / nprocs;
  } else {
    return 0;
  }
  phgPrintf("%d \n", vec1procs);
  MAP *vecmap =
      phgMapCreateSimpleMap(phgComm, vec1procs, Nbasis * Nbasis * Nbasis);

  VEC *cini = phgMapCreateVec(vecmap, 1);
  INT max_step = 1e3;
  FLOAT tol = 1e-2, alpha = 0.8;
  VEC *c = NULL;
  INT tN = 100;

  FLOAT x3value, x2value, x1value;

  MAT *Amat_ini = phgMapCreateMat(vecmap, vecmap);
  MAT *DiagMat = phgMapCreateMat(vecmap, vecmap);
  MAT *x1Mat_imag = phgMapCreateMat(vecmap, vecmap);
  MAT *x2Mat_imag = phgMapCreateMat(vecmap, vecmap);
  MAT *x3Mat_imag = phgMapCreateMat(vecmap, vecmap);

  phgMatDisassemble(DiagMat);
  phgMatDisassemble(x1Mat_imag);
  phgMatDisassemble(x2Mat_imag);
  phgMatDisassemble(x3Mat_imag);

  INT i, j, k, k1, k2, j1, j2, i1, i2, h;
  INT location, location2, location1;
  if (phgRank == 0) {
    for (i = 0; i < Nbasis; i++) {
      for (j = 0; j < Nbasis; j++) {
        for (k = 0; k < Nbasis; k++) {
          location = i * Nbasis * Nbasis + j * Nbasis + k;
          i1 = i - Nbasis / 2;
          j1 = j - Nbasis / 2;
          k1 = k - Nbasis / 2;
          phgMatAddGlobalEntry(DiagMat, location, location,
                               -(i1 ^ 2 + j1 ^ 2 + k1 ^ 2) / 2);
        }
      }
    }
    for (i = 0; i < Nbasis; i++) {
      for (j = 0; j < Nbasis; j++) {
        for (k1 = 0; k1 < Nbasis; k1++) {
          for (k2 = 0; k2 < Nbasis; k2++) {
            location1 = i * Nbasis * Nbasis + j * Nbasis + k1;
            location2 = i * Nbasis * Nbasis + j * Nbasis + k2;
            h = k1 + k2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }

            phgMatAddGlobalEntry(x3Mat_imag, location1, location2, x3value);
            phgMatAddGlobalEntry(x3Mat_imag, location2, location1, x3value);
          }
        }
      }
    }

    for (k = 0; k < Nbasis; k++) {
      for (j = 0; j < Nbasis; j++) {
        for (i1 = 0; i1 < Nbasis; i1++) {
          for (i2 = 0; i2 < Nbasis; i2++) {
            location1 = i1 * Nbasis * Nbasis + j * Nbasis + k;
            location2 = i2 * Nbasis * Nbasis + j * Nbasis + k;
            h = i1 + i2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }
            phgMatAddGlobalEntry(x1Mat_imag, location1, location2, x3value);
            phgMatAddGlobalEntry(x1Mat_imag, location2, location1, x3value);
          }
        }
      }
    }
    for (i = 0; i < Nbasis; i++) {
      for (k = 0; j < Nbasis; j++) {
        for (j1 = 0; j1 < Nbasis; j1++) {
          for (j2 = 0; j2 < Nbasis; j2++) {
            location1 = i * Nbasis * Nbasis + j1 * Nbasis + k;
            location2 = i * Nbasis * Nbasis + j2 * Nbasis + k;
            h = j1 + j2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }
            phgMatAddGlobalEntry(x2Mat_imag, location1, location2, x3value);
            phgMatAddGlobalEntry(x2Mat_imag, location2, location1, x3value);
          }
        }
      }
    }
  }
  phgMatAssemble(DiagMat);
  phgMatAssemble(x1Mat_imag);
  phgMatAssemble(x2Mat_imag);
  phgMatAssemble(x3Mat_imag);
  char *diag_mat_name = "diagmat";
  char *output_diag_mat_name = "diagmat.m";

  char *x1_mat_name = "x1mat";
  char *output_x1_mat_name = "x1mat.m";

  char *x2_mat_name = "x2mat";
  char *output_x2_mat_name = "x2mat.m";

  char *x3_mat_name = "x3mat";
  char *output_x3_mat_name = "x3mat.m";

  phgMatDumpMATLAB(DiagMat, diag_mat_name, output_diag_mat_name);
  phgMatDumpMATLAB(x1Mat_imag, x1_mat_name, output_x1_mat_name);
  phgMatDumpMATLAB(x2Mat_imag, x2_mat_name, output_x2_mat_name);
  phgMatDumpMATLAB(x3Mat_imag, x3_mat_name, output_x3_mat_name);
  // ode45_solver_parallel(vecmap, c, cini, Nbasis, max_step, tol, alpha, tN);
  phgFinalize();

  return 0;
}
