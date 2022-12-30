//g++ -I . ukol3.cpp -o ukol3 -lumfpack
#include <iostream>
#include <cmath>
#include <fstream>
#include "suitesparse/umfpack.h"
#include<eigen/Eigen/Sparse>
using namespace std;
using Eigen::SparseMatrix;

int main(){
    unsigned int N = 10000;
    int m = (int) sqrt(N);

    void *Symbolic, *Numeric ;
    double *null = (double *) NULL;

    double x[N];
    double b[N];
    // double b [ ] = {8., 45., -3., 3., 19.} ;

    SparseMatrix<double> A(N, N);

    for (size_t i=0; i<N; i++) {
        for (size_t j=0; j<N; j++){
            if(i==j){
                A.coeffRef(i, j) = 4;
            }
            if(j == i+1 || j == i-1){
                A.coeffRef(i, j) = -1;
            }
                if(j == i+m || j == i-m){
                A.coeffRef(i, j) = -1;
            }
        }
        b[i] = 2.0;
    }


    A.makeCompressed();
    umfpack_di_symbolic (N, N, A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr(), &Symbolic, null, null) ;
    umfpack_di_numeric (A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr(), Symbolic, &Numeric, null, null) ;
    umfpack_di_free_symbolic (&Symbolic) ;
    umfpack_di_solve (UMFPACK_A, A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr(), x, b, Numeric, null, null);
    umfpack_di_free_numeric (&Numeric) ;

    //Výpis prvních 10 výsledků
    for (int i=0; i<10; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }
    //Uložení do souborů
  ofstream myfile;
  myfile.open ("Ukol3.txt");
  for(int i = 0; i<N; i++){
    myfile << "Reseni "<< i << ": "; 
    myfile <<x[i];
    myfile << "\n";
  }
  myfile.close();


}
    
