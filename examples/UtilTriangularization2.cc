//An example of program to test the speed of triangularization
//with Util::LowerTriangularization (Lecuyer method)


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
    int l=bas_mat.size1();
    int c=bas_mat.size2();
     for(int i=0;i<l;i++)
     {for(int j=0;j<c;j++){
       std::cout <<  bas_mat(i,j)<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}

}

int main() {


      std::string prime = primes[0];

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1,matrix2,matrix3;
      unsigned int ln;
      IntLatticeBase<Int, Real, RealRed>* basis;
      Reducer<Int, Real, RealRed>* red;

      Int m(1021);

      std::cout << name<<std::endl;
      name = "bench/" + prime+ "_5" + "_2" ;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);
      matrix3.SetDims(numlines, numlines);
     
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);
      
      std::cout << " The base before reduction\n"; 
      printBase(matrix1);

      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1,matrix1,m, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);
      //red->redLLL(0.999,1000000,numlines);
      basis->updateVecNorm();


      std::cout << " The base after reduction\n"; 
      
       printBase((red->getIntLatticeBase())->getBasis()); 
       CopyMatr(matrix2,(red->getIntLatticeBase())->getBasis(), numlines, numlines);
       //Triangularization2<IntMat,IntVec, Int> (matrix2, matrix3, m);
       TriangularizationLower<IntMat,IntVec,Int>(matrix1, matrix2 ,m);
       printBase(matrix2) ;

 
  return 0;
}
