#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/WriterRes.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"

#include "Examples.h"

using namespace LatticeTester;

namespace {
  // Returns the average of the length of this vector
  Real average(RealVec vector) {
    Real sum(0);
    for (int i = 0; i<vector.length(); i++) {
      sum += vector[i];
    }
    return sum/Real(vector.length());
  }

void printRes (RealMat mat, int lin, int col){
   std::ofstream out("ResultatUtilTriang.csv");

    for (int j=0;j<lin;j++) {
      for (int k=0;k<col;j++) 
        out << j <<',';
     out << '\n';
    }
   out.close();
}



void printBase(IntMat bas_mat){
     int L=bas_mat.NumRows();
     int C=bas_mat.NumCols();

     for(int i=0;i<L;i++)
     {for(int j=0;j<C;j++){
       std::cout <<  bas_mat(i,j)<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}


void ModuloMatrix(IntMat &bas_mat, Int mod){
     int L=bas_mat.NumRows();
     int C=bas_mat.NumCols();

     for(int i=0;i<L;i++)
     { for(int j=0;j<C;j++){
        Modulo(bas_mat(i,j),mod, bas_mat(i,j));
       }
   
     }

}

}

int main() {
 
      std::string prime = primes[0];

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1,matrix2;
      unsigned int ln;
      BasisConstruction<Int> constr;
      Int m(1021);

      name = "bench/" + prime+ "_5" + "_2" ;
      std::cout << name<<std::endl;
      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);
     
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);


      std::cout << " The base before triangularization\n"; 
      printBase(matrix1); 

      //constr.GCDTriangularBasis(matrix1,m);
      constr.GCDTriangularBasis(matrix1);
      ModuloMatrix(matrix1,m);
    
      std::cout << " The base after triangularization\n"; 
      printBase(matrix1); 

 
  return 0;
}
