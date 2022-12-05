//An example to use BasisConstruction::GCDTriangular
//The triangular method implement by Mark-Antoine 
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


}

int main() {
  clock_t timer = clock();
  int max_dim = 15; //! Actual max dim is 5*max_dim

  clock_t tmp;

  RealMat matTriangularTime;
  matTriangularTime.resize(15, 10);

  std::string prime = primes[2];

  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1,matrix2;
      unsigned int ln;
      BasisConstruction<Int> constr;
      Int m(1021);

    
      //! Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      std::cout << name<<std::endl;
      //name = "bench/" + prime+ "_2" + "_001" ;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);
     
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      tmp = clock(); 

     // constr.GCDTriangularBasis(matrix1,m);
      constr.GCDTriangularBasis(matrix1);

      double tps=(double)(clock() - tmp); //(CLOCKS_PER_SEC);
      matTriangularTime(j,k)=tps;
     // delete constr;
  
    }
  }

  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";
  for(int i=0;i<15;i++)
   {for(int j=0;j<10;j++){
       std::cout << matTriangularTime(i,j)<< "   ";
     }
     std::cout << ""<< std::endl;
   }
 
  return 0;
}
