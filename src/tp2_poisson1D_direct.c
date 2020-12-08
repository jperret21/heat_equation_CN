/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la; // nb ligne de matrice A
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double *B;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=102;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  B=(double *) malloc(sizeof(double)*la);

  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_dense_RHS_DBC_1D(B,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  write_vec(RHS, &la, "RHS.dat");       //b
  write_vec(EX_SOL, &la, "EX_SOL.dat"); // sol exact
  write_vec(X, &la, "X_grid.dat");      // subdicision

  kv=1; //une ligne de 0
  ku=1; // nd diag uper
  kl=1; // nb diag lower
  lab=kv+kl+ku+1; // taille matrice



  AB = (double *) malloc(sizeof(double)*lab*la);

  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //

  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");

    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);

  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
  }


  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;

  printf("\nThe relative residual error is relres = %e\n",relres);

//----------------------------------------------------------------------------//

  // Test de dgbmv  EXO 4 :
  double *CB;
  CB = (double *) malloc(sizeof(double)*(lab-1)*la);
  for (int i=1; i<la-1; i++)
  {
    CB[i*(lab-1)]=-1.0;
    CB[i*(lab-1) + 1]=2.0;
    CB[i*(lab-1) + 2]=-1.0;
  }
  CB[1]=2.0;
  CB[2]=-1.0;
  CB[(lab-1)*la - 3] = -1.0;
  CB[(lab-1)*la - 2] = 2.0;


  //dgbmv en row major :
  double *RHS_2;
  RHS_2=(double *) malloc(sizeof(double)*la);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

  //dgbmv en Col major :
  cblas_dgbmv(CblasColMajor, CblasNoTrans, la, la, kl, ku, 1.0, CB, lab-1, EX_SOL, 1, 0, RHS_2, 1);
  //dgbmv en Row major :
  //cblas_dgbmv(CblasRowMajor, CblasNoTrans, la, la, kl, ku, 1.0, AB, lab, EX_SOL, 1, 0, RHS, 1);

  write_vec(RHS_2, &la, "my_dgbmv.dat");
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, RHS_2, 1);
  relres = cblas_ddot(la, RHS_2, 1, RHS_2,1);
  relres = sqrt(relres);
  relres = relres / temp;

  printf("\nThe relative residual error is relres (cblas_dgbmv) = %e\n",relres);


  //

/*
  //test de dgbmv //
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,1,1,1.0,AB,kl+ku+1,RHS,1,0,RHS,1 );
  //cblas_dgbmv(CblasColMajor,CblasNoTrans,lab,la,1,1,1.0,AB,3,EX_SOL,1,0,RHS,1 );
  //write_vec(RHS, &la, "dgbmv.dat"); // rhs doit egale au rhs init
  write_xy(RHS, B, &la, "test_result.dat");
*/
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(B);
  free(ipiv);

  printf("\n\n--------- End -----------\n");



}
