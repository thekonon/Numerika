//make ukol4
//mpiexec -n 1 ./ukol4
#include <iostream>
#include <cmath>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
using namespace std;
const int n = 1000;

int main(int argc,char **args)
{
  /* Inicializace */
  CHKERRQ( PetscInitialize( &argc , &args , (char *)0 , 0 ) );

  /* Vytvoreni vektoru a nastaveni jeho celkove velikosti na 10 a hodnot na 0 */
  Vec y;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &y);
  VecSet(y, 2.0);

  // Vytvoreni matice o velkosti 10x10
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(A);
  MatSetUp(A);

  // Nastaveni prvku na 1, -2, 1
  PetscInt istart, iend;
  MatGetOwnershipRange(A, &istart, &iend);
  int m = round(sqrt(n));
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++){
        if(i==j){
          MatSetValue(A, i, j, 4, INSERT_VALUES);
        }
        if(j == i+1 || j == i-1){
          MatSetValue(A, i, j, -1, INSERT_VALUES);
        }
        if(j == i+m || j == i-m){
          MatSetValue(A, i, j, -1, INSERT_VALUES);
        }
    }
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // Reseni soustavy Ax = y
  Vec x;  
  VecDuplicate(y, &x);

  KSP solver;
  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A);
  
  
  KSPSetType(solver, KSPCG);
  PC prec;
  KSPGetPC(solver,&prec);
  PCSetType(prec,PCJACOBI);

  KSPSetFromOptions(solver);
  KSPSetUp(solver);

  KSPSolve(solver, y, x);

  //VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  PetscViewer lab;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"ukol4.txt",&lab);
  VecView(x, lab);
  KSPDestroy(&solver);
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&y);

  cout<<"Řešení bylo uloženo do souboru ukol4.txt"<<endl;

  CHKERRQ( PetscFinalize() );

  return 0;
}
