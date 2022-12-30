//g++ ukol12.cpp -o Ukol1 -llapacke -llapack
#include <iostream>
#include <lapacke.h>
#include <fstream>
#include <cmath>
using namespace std;

class Matrix {
public:

  Matrix(size_t rows, size_t cols):
    rows_(rows),
    cols_(cols)
  { data_ = new double[rows*cols]; }

  ~Matrix() { delete[] data_; }

  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  size_t lda() const { return cols_; }

  double* data() { return data_; }

  double* operator[](int row) { return data_ + row*lda(); }

private:
  size_t rows_;
  size_t cols_;
  double* data_;
};

int main() {
  int n = 1e3;

  double* x = new double[n];
  double* y = new double[n];

  //Matice A je pro úkol 2 a AA je úkol 1
  Matrix A(n, n);
  Matrix AA(n,n);

  //Vyplnění matici
  int m = round(sqrt(n));
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++){
        if(j == i+1 || j == i-1){
        A[i][j] = -1;
      }else{
        A[i][j] = 0;
      }
      if(j == i+m || j == i-m){
        A[i][j] = -1;
      }
    }
    A[i][i] = 4;
    y[i] = 2.0;
  }
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++){
        if(j == i+1 || j == i-1){
        AA[i][j] = -1;
      }else{
        AA[i][j] = 0;
      }
      if(j == i+m || j == i-m){
        AA[i][j] = -1;
      }
    }
    AA[i][i] = 4;
  }
  cout<<"Velikost m: "<<m<<endl;

  ///Příprava proměnných
  int* ipiv = new int[n];
  double* trA = new double[n];
  double* my_eigs = new double[n];
  char equed[2] = "N"; 
  double r[n], c[n], rcond, ferr, berr, rpivot;
  int LWORK;

  //Úkol 1 - výpočet
  int info2 = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'L', n, AA.data(), AA.lda(), my_eigs);

  //Výpočet úkol2 x = A \ y
  int info = LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 
          'E',
          'N', 
			    A.rows(), 
          1, 
          A.data(), 
          A.lda(),
			    A.data(), 
          A.lda(),
			    ipiv, 
          equed, 
          r, 
          c,
			    y, 
          1, 
			    x, 
          1,
			    &rcond, 
          &ferr, 
          &berr, 
          &rpivot);

  //Výpis vlastních čísel
  for(int i = 0; i<10; i++)
  {
    cout<<"Eig["<<i<<"]: "<<my_eigs[i]<<endl;
  }
  //Výpis řešení Ax=b
  cout<<"Prvních 10 vektoru výsledků x[i]:"<<endl;
  for (int i = 0; i<10; i++){
    cout<<"x[i] = "<<x[i]<<endl;
  }

  //Uložení do souborů
  ofstream myfile;
  myfile.open ("Ukol1.txt");
  for(int i = 0; i<n; i++){
    myfile << "Vl. cislo "<< i << ": "; 
    myfile << my_eigs[i];
    myfile << "\n";
  }
  myfile.close();
  myfile.open ("Ukol2.txt");
  for(int i = 0; i<n; i++){
    myfile << "Vysledek "<< i << ": "; 
    myfile << x[i];
    myfile << "\n";
  }
  myfile.close();
  cout<<"Výsledky byly uloženy do souborů Ukol1.txt a Ukol2.txt"<<endl;
  //Uvolnění paměti
  delete[] ipiv;
  delete[] y;
  delete[] x;

  return 0;
}
