
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
void ode45_solver_parallel(MAP *vecmap, VEC *c, VEC *cini, INT Nbasis,
                           INT max_step, FLOAT tol, FLOAT alpha, INT tN) {

  VEC *c_pre = NULL;
  VEC *c_err = NULL;

  FLOAT c_err_norm2;
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
  /*
  MAT *Amat = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat1 = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat2 = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat3 = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat4 = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat5 = phgMapCreateMat(vecmap, vecmap);
  MAT *Amat6 = phgMapCreateMat(vecmap, vecmap);
*/
  MAT *Amat;
  MAT *Amat1;
  MAT *Amat2;
  MAT *Amat3;
  MAT *Amat4;
  MAT *Amat5;
  MAT *Amat6;
  phgVecDisassemble(cini);

  phgPrintf("%d \n", cini->nvec);
  // phgVecAddGlobalEntry(cini, 0, 1, 1.0);
  INT *I = phgAlloc(Nbasis * sizeof(INT));
  // FLOAT *Addvec = phgAlloc(vec1procs * sizeof(FLOAT));

  if (phgRank == 0) {
    phgVecAddGlobalEntry(cini, 0, 0, 1.0);
  }

  /*
  for (i = 0; i < vec1procs; i++) {
    *(I + i) = i;
    *(Addvec + i) = i * i;
    phgPrintf("%e,\n", *(Addvec + i));
  }*/
  // phgVecAddEntries(cini, 0, vec1procs, I, Addvec);
  phgVecAssemble(cini);

  // phgVecCopy(cini, &c_pre);
  char *vec_name = "vec1";
  char *outputname = "vec11.m";
  char *vec_name1 = "vec1_pre";
  char *outputname1 = "vec11_pre.m";
  char *mat_name = "mat1";
  char *output_mat_name = "mat11.m";
  char *wc_mat_name1 = "wcmat1";
  char *output_wc_mat_name1 = "wcmat11.m";

  char *wc_mat_name2 = "wcmat2";
  char *output_wc_mat_name2 = "wcmat22.m";

  char *wc_mat_name3 = "wcmat3";
  char *output_wc_mat_name3 = "wcmat33.m";

  char *wc_mat_name4 = "wcmat4";
  char *output_wc_mat_name4 = "wcmat44.m";

  char *wc_mat_name5 = "wcmat5";
  char *output_wc_mat_name5 = "wcmat55.m";

  char *wc_mat_name6 = "wcmat6";
  char *output_wc_mat_name6 = "wcmat66.m";

  FLOAT vecnorm2;

  phgVecNorm2(cini, 0, &vecnorm2);

  phgPrintf("%e,\n", vecnorm2);

  FLOAT *AddMat = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat1 = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat2 = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat3 = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat4 = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat5 = phgAlloc(Nbasis * sizeof(FLOAT));
  FLOAT *AddMat6 = phgAlloc(Nbasis * sizeof(FLOAT));
  // FLOAT *AddMat7 = phgAlloc(Nbasis * sizeof(FLOAT));

  FLOAT t_ini = 0;
  FLOAT deltat = 1e-3;

  FLOAT t_now, t1, t2, t3, t4, t5, t6;

  t_now = t_ini;
  /// INT j, k;
  // for (j = 0; j < Nbasis; j++) {
  //  *(I + j) = j;
  //}

  /* if (c1 == NULL) {
     phgPrintf("c1=NULL");
   }*/
  printf("%d\n", (int)phgRank);
  FLOAT omega = 1.0;
  /*for (j = 0; j < vec1procs; j++) {
    check;
    // I[j] = phgMapL2G(vecmap, j);
    check;
    printf("%d %d\n", phgMapL2G(vecmap, j), j);
  }
*/

  FLOAT tmin = deltat / (FLOAT)max_step;
  FLOAT t_delta_change, total_t_delta, total_err, normw, err;

  for (j = 0; j < Nbasis; j++) {
    *(I + j) = j;
  }
  // check;
  c_pre = phgVecCopy(cini, NULL);
  phgVecDumpMATLAB(c_pre, vec_name1, outputname1);
  // phgVecDumpMATLAB(c_pre, vec_name, outputname);
  //  check;
  VEC *c1 = NULL;
  VEC *c2 = NULL;
  VEC *c3 = NULL;
  VEC *c4 = NULL;
  VEC *c5 = NULL;
  VEC *c6 = NULL;
  VEC *c7 = NULL;
  VEC *wc1 = NULL;

  VEC *wc2 = NULL;
  VEC *wc3 = NULL;
  VEC *wc4 = NULL;
  VEC *wc5 = NULL;
  VEC *wc6 = NULL;
  INT t;
  for (t = 0; t < tN; t++) {
    t_now = t_ini + t * deltat;
    for (i = 0; i < max_step; i++) { // for (i = 0; i < tN; i++) {

      t_delta_change = deltat;
      total_t_delta = 0.;
      total_err = 0.;
      normw = 0.;

      Amat = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat1 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat2 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat3 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat4 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat5 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);
      Amat6 = phgMatAXPBY(1.0, Amat_ini, 0.0, NULL);

      phgMatDisassemble(Amat);
      phgMatDisassemble(Amat1);
      phgMatDisassemble(Amat2);
      phgMatDisassemble(Amat3);
      phgMatDisassemble(Amat4);
      phgMatDisassemble(Amat5);
      phgMatDisassemble(Amat6);
      t1 = t_now + 1.0 * t_delta_change / 5.0;
      t2 = t_now + 3.0 * t_delta_change / 10.0;
      t3 = t_now + 4.0 * t_delta_change / 5.0;
      t4 = t_now + 8.0 * t_delta_change / 9.0;
      t5 = t_now + t_delta_change;
      t6 = t5;
      if (phgRank == 0) {
        for (j = 0; j < Nbasis; j++) {
          //*(I + j) = j;
          for (k = 0; k < Nbasis; k++) {
            *(AddMat + k) = Cos(omega *
                                (1.0 / ((k + 1.0) * (k + 1.0)) -
                                 1.0 / ((j + 1.0) * (j + 1.0))) *
                                t_now);
            *(AddMat1 + k) = Cos(omega *
                                 (1.0 / ((k + 1.0) * (k + 1.0)) -
                                  1.0 / ((j + 1.0) * (j + 1.0))) *
                                 t1);
            *(AddMat2 + k) = Cos(omega *
                                 (1.0 / ((k + 1.0) * (k + 1.0)) -
                                  1.0 / ((j + 1.0) * (j + 1.0))) *
                                 t2);
            *(AddMat3 + k) = Cos(omega *
                                 (1.0 / ((k + 1.0) * (k + 1.0)) -
                                  1.0 / ((j + 1.0) * (j + 1.0))) *
                                 t3);
            *(AddMat4 + k) = Cos(omega *
                                 (1.0 / ((k + 1.0) * (k + 1.0)) -
                                  1.0 / ((j + 1.0) * (j + 1.0))) *
                                 t4);
            *(AddMat5 + k) = Cos(omega *
                                 (1.0 / ((k + 1.0) * (k + 1.0)) -
                                  1.0 / ((j + 1.0) * (j + 1.0))) *
                                 t5);
          }
          // if (i == (tN - 1)) {
          // phgPrintf("MAT= %e,\n", *AddMat);
          // phgPrintf("I=%d \n", *(I + j));
          // }

          phgMatAddGlobalEntries(Amat, 1, I + j, Nbasis, I, AddMat);

          phgMatAddGlobalEntries(Amat1, 1, I + j, Nbasis, I, AddMat1);
          phgMatAddGlobalEntries(Amat2, 1, I + j, Nbasis, I, AddMat2);
          phgMatAddGlobalEntries(Amat3, 1, I + j, Nbasis, I, AddMat3);
          phgMatAddGlobalEntries(Amat4, 1, I + j, Nbasis, I, AddMat4);
          phgMatAddGlobalEntries(Amat5, 1, I + j, Nbasis, I, AddMat5);
          phgMatAddGlobalEntries(Amat6, 1, I + j, Nbasis, I, AddMat5);
        }
      }

      phgMatAssemble(Amat);
      // phgMatDumpMATLAB(Amat, mat_name, output_mat_name);
      phgMatAssemble(Amat1);

      phgMatAssemble(Amat2);

      phgMatAssemble(Amat3);

      phgMatAssemble(Amat4);

      phgMatAssemble(Amat5);

      phgMatAssemble(Amat6);
      /*
          phgVecDestroy(&c1);
          phgVecDestroy(&wc1);
          phgVecDestroy(&c2);
          phgVecDestroy(&wc2);
          phgVecDestroy(&c3);
          phgVecDestroy(&wc3);
          phgVecDestroy(&c4);
          phgVecDestroy(&wc4);
          phgVecDestroy(&c5);
          phgVecDestroy(&wc5);
          phgVecDestroy(&c6);
          phgVecDestroy(&wc6);
          phgVecDestroy(&c7);
          */

      c1 = phgMatVec(0, 1.0, Amat, c_pre, 0.0, NULL);
      // phgMatVec(0, 1.0, Amat, c_pre, 0.0, &c1);
      //      check;
      phgVecCopy(c1, &wc1);
      //    check;
      phgVecAXPBY(1.0, c_pre, 1.0 * t_delta_change / 5.0, &wc1);

      c2 = phgMatVec(0, 1.0, Amat1, wc1, 0.0, NULL);
      phgVecCopy(c1, &wc2);
      phgVecAXPBY(1.0, c_pre, 3.0 * t_delta_change / 40.0, &wc2);

      phgVecAXPBY(9.0 * t_delta_change / 40.0, c2, 1.0, &wc2);

      c3 = phgMatVec(0, 1.0, Amat2, wc2, 0.0, NULL);
      phgVecCopy(c1, &wc3);
      phgVecAXPBY(1.0, c_pre, 44.0 * t_delta_change / 45.0, &wc3);
      phgVecAXPBY(-56.0 * t_delta_change / 15.0, c2, 1.0, &wc3);
      phgVecAXPBY(32.0 * t_delta_change / 9.0, c3, 1.0, &wc3);

      c4 = phgMatVec(0, 1.0, Amat3, wc3, 0.0, NULL);
      phgVecCopy(c1, &wc4);
      phgVecAXPBY(1.0, c_pre, 19372.0 * t_delta_change / 6561.0, &wc4);
      phgVecAXPBY(-25360.0 * t_delta_change / 2187.0, c2, 1.0, &wc4);
      phgVecAXPBY(64448.0 * t_delta_change / 6561.0, c3, 1.0, &wc4);
      phgVecAXPBY(-212.0 * t_delta_change / 729.0, c4, 1.0, &wc4);

      c5 = phgMatVec(0, 1.0, Amat4, wc4, 0.0, NULL);
      phgVecCopy(c1, &wc5);
      phgVecAXPBY(1.0, c_pre, 9017.0 * t_delta_change / 3168.0, &wc5);
      phgVecAXPBY(-355.0 * t_delta_change / 33.0, c2, 1.0, &wc5);
      phgVecAXPBY(46732.0 * t_delta_change / 5247.0, c3, 1.0, &wc5);
      phgVecAXPBY(49.0 * t_delta_change / 176.0, c4, 1.0, &wc5);
      phgVecAXPBY(-5103.0 * t_delta_change / 18656.0, c5, 1.0, &wc5);

      c6 = phgMatVec(0, 1.0, Amat5, wc5, 0.0, NULL);
      phgVecCopy(c1, &wc6);
      phgVecAXPBY(1.0, c_pre, 35.0 * t_delta_change / 384.0, &wc6);
      // wc6 = phgVecAXPBY(1.0, wc5, -355.0 * deltat / 33.0, &c2);
      phgVecAXPBY(500.0 * t_delta_change / 1113.0, c3, 1.0, &wc6);
      phgVecAXPBY(125.0 * t_delta_change / 192.0, c4, 1.0, &wc6);
      phgVecAXPBY(-2187.0 * t_delta_change / 6784.0, c5, 1.0, &wc6);
      phgVecAXPBY(11.0 * t_delta_change / 84.0, c6, 1.0, &wc6);
      c7 = phgMatVec(0, 1.0, Amat6, wc6, 0.0, NULL);

      phgVecCopy(c1, &c_err);
      phgVecAXPBY(-71.0 / 16695.0, c3, 71.0 / 57600.0, &c_err);
      phgVecAXPBY(71.0 / 1920.0, c4, 1.0, &c_err);
      phgVecAXPBY(-17253. / 339200.0, c5, 1.0, &c_err);
      phgVecAXPBY(22.0 / 525.0, c6, 1.0, &c_err);

      phgVecAXPBY(1.0 / 40.0, c7, 1.0, &c_err);
      c_err_norm2 = phgVecNorm2(c_err, 0, NULL);

      c_err_norm2 = c_err_norm2 * t_delta_change;
      normw = phgVecNorm2(wc6, 0, NULL);
      // phgPrintf("%e \n", c_err_norm2);
      err = c_err_norm2 / normw;
      t_now = t_now + t_delta_change;

      if (t_delta_change != tmin) {
        if (err > tol) {
          t_delta_change = t_delta_change * alpha * Pow(tol / err, 1. / 5.);
          // printf("%e\n", t_delta_change);
          if (t_delta_change < tmin) {
            t_delta_change = tmin;
            i = i - 1;
            continue;
          } else {
            i = i - 1;
            continue;
          }
        }
      }
      if ((total_t_delta + t_delta_change) == deltat) {

        phgVecCopy(wc6, &c_pre);
        break;
      }
      if ((total_t_delta + t_delta_change) > deltat) {
        t_delta_change = deltat - total_t_delta;
      }
      total_t_delta += t_delta_change;
      // phgVecCopy(wc6, &c_pre);
      check;
      // phgMatDestroy(&Amat);
    }
  }
  phgVecCopy(c_pre, &c);
  // phgVecDumpMATLAB(wc6, wc_mat_name1, output_wc_mat_name1);
  phgPrintf("%e \n", t_now - deltat);
  phgVecDumpMATLAB(wc6, wc_mat_name1, output_wc_mat_name1);
  phgVecDumpMATLAB(c_pre, vec_name, outputname);
  // phgMatDumpMATLAB(Amat, wc_mat_name1, output_wc_mat_name1);
  /*
    phgMatDumpMATLAB(Amat1, wc_mat_name1, output_wc_mat_name1);
    phgMatDumpMATLAB(Amat2, wc_mat_name2, output_wc_mat_name2);
    phgMatDumpMATLAB(Amat3, wc_mat_name3, output_wc_mat_name3);
    phgMatDumpMATLAB(Amat4, wc_mat_name4, output_wc_mat_name4);
    phgMatDumpMATLAB(Amat5, wc_mat_name5, output_wc_mat_name5);
    // phgVecDumpMATLAB(c6, wc_mat_name6, output_wc_mat_name6);
  */

  // phgVecDumpMATLAB(c1, wc_mat_name1, output_wc_mat_name1);
  /*
    phgVecDumpMATLAB(c2, wc_mat_name2, output_wc_mat_name2);
    phgVecDumpMATLAB(c3, wc_mat_name3, output_wc_mat_name3);
    phgVecDumpMATLAB(c4, wc_mat_name4, output_wc_mat_name4);
    phgVecDumpMATLAB(c5, wc_mat_name5, output_wc_mat_name5);
    phgVecDumpMATLAB(c6, wc_mat_name6, output_wc_mat_name6);
  */
  /*
  phgVecDumpMATLAB(wc1, wc_mat_name1, output_wc_mat_name1);
  phgVecDumpMATLAB(wc2, wc_mat_name2, output_wc_mat_name2);
  phgVecDumpMATLAB(wc3, wc_mat_name3, output_wc_mat_name3);
  phgVecDumpMATLAB(wc4, wc_mat_name4, output_wc_mat_name4);
  phgVecDumpMATLAB(wc5, wc_mat_name5, output_wc_mat_name5);
  phgVecDumpMATLAB(wc6, wc_mat_name6, output_wc_mat_name6);
*/
  // check;
  phgFree(AddMat);
  phgFree(AddMat1);
  phgFree(AddMat2);
  phgFree(AddMat3);
  phgFree(AddMat4);
  phgFree(AddMat5);
  // phgMatDumpMATLAB(Amat, mat_name, output_mat_name);
  phgVecDestroy(&cini);
  // check;
  phgMatDestroy(&Amat);
  phgMatDestroy(&Amat1);
  phgMatDestroy(&Amat2);
  phgMatDestroy(&Amat3);
  phgMatDestroy(&Amat4);
  phgMatDestroy(&Amat5);
  phgMatDestroy(&Amat6);
  // check;
  // phgMapDestroy(&vecmap);
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
  ode45_solver_parallel(vecmap, c, cini, Nbasis, max_step, tol, alpha, tN);
  phgFinalize();

  return 0;
}