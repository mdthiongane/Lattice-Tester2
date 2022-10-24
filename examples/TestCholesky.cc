#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>
#include <fstream>
#include <cmath>
#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/WriterRes.h"
#include "Examples.h"

using namespace std;
using namespace LatticeTester;

namespace {

void transposeMatrix(IntMat &mat1, IntMat &mat2)
{  int n=mat1.NumRows();
   for(int i=0; i<n;i++){
    for(int j=0;j<n;j++)
      mat2[j][i]=mat1[i][j];
    //{ std::cout << i<<j<<" "<<std::endl;}
   //   std::cout<<""<<std::endl;
   }

}

void produit(IntMat &mat1, IntMat &mat2, IntMat &mat3)
{   int n=mat1.NumRows();
    for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
      { Int z(0);
        mat3[i][j]=z;
        for(int k = 0; k < n; k++)
        {
          mat3[i][j] =mat3[i][j]+ mat1[i][k] * mat2[k][j];
         }
       }
    }

}

void Cholesky_Decomposition(IntMat &matrix, IntMat &lower)
{   int n=matrix.NumRows();
    // lower[n][n];
   // memset(lower, 0, sizeof(lower));
 
    // Decomposing a matrix into Lower Triangular
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            Int sum(0);
 
            if (j == i) // summation for diagonals
            {
                for (int k = 0; k < j; k++)
                    sum =sum+ lower[j][k]*lower[j][k];
               // lower[j][j] = sqrt(matrix[j][j] - sum);
                 lower[j][j] = SqrRoot(matrix[j][j] - sum);
            } else {
 
                // Evaluating L(i, j) using L(j, j)
                for (int k = 0; k < j; k++)
                    sum =sum+ (lower[i][k] * lower[j][k]);
                lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
            }
        }

    }  

}



void printMat(IntMat &mat){
    int n=mat.NumRows();
   for(int i=0;i<n;i++)
     {for(int j=0;j<n;j++){
       std::cout << mat[i][j]<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}

void copy(IntMat &b1, IntMat &b2){
 
     for(int i=0;i<b1.size1();i++)
     { for(int j=0;j<b1.size2();j++){
          b2(i,j)=b1(i,j);
         }   
     }

}


}

int main() {

   /** int mat1[][5] = { { -199, -146, -19,   1,   69 },
                        { 61,   -129,   -125,   -204,   51  },
                        {30,   203,   244,   44,   626},
                        {-170,   214,   507,   -602,   -50},
                        { 488, -521, 436, 357, 149  } };   
    int mat2[5][5];
    int mat3[5][5];
    
    */
     IntLatticeBase<Int, Real, RealRed>* basis;
      Reducer<Int, Real, RealRed>* red;

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1, matrix2;
      IntMat mat1,mat2,mat3,lower,uper;
      unsigned int ln;
      
        std::string prime = primes[0];
        std::string s1("cholesky");

      //! Reader shenanigans
      //name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      name = "bench/" + prime+ "_2" + "_001" ;
       // name = "bench/" + prime+ "_4" + "_002" ;
      //  name = "bench/" + prime+ "_5" + "_4" ;



      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);
      mat1.SetDims(numlines, numlines);
      mat2.SetDims(numlines, numlines);
      mat3.SetDims(numlines, numlines);
      lower.SetDims(numlines, numlines);
      uper.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      std::cout << " ########### The initial basis ###############"<<std::endl;
      printMat(matrix1);
 
      // BKZ reduction before shortest vector search
      Int m(101);
      CalcDual (matrix1, matrix2,numlines, m);
      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1,matrix2,m, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);
      red->redBKZ();
      basis->updateVecNorm();
     //  std::cout << " ########### Before copy ###############"<<std::endl;
      copy((red->getIntLatticeBase())->getBasis(), mat1);
      
      std::cout << " ########### The reducer matrice basis ###############"<<std::endl;
       printMat(mat1);

       transposeMatrix(mat1, mat2);
       std::cout << " ########### The Transpose matrice of the basis ###############"<<std::endl;
       printMat(mat2);
       produit(mat1, mat2, mat3);
       std::cout << " ########### The product matrice of the basis and it transpose ###############"<<std::endl;
       printMat(mat3);

      
       Cholesky_Decomposition(mat3,lower);
       std::cout << " ########### After  Cholesky the lower triangular matrice ###############"<<std::endl; 
       printMat(lower);

      
       std::cout << " ########### After  Cholesky the uper triangular matrice ###############"<<std::endl;  
        transposeMatrix( lower, uper);
        printMat(uper);

      if (!red->shortestVector(L2NORM,s1)) {
        std::cout << " shortestVector failed"<<std::endl; 
      }
   
  return 0;
}
