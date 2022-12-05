/*
**An example programm to compute a shortest vector using 
* a triangular decomposition method. The implemented method
* is discribe in LatticeGuide section 8.2
*/

#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/WriterRes.h"
#include "latticetester/Util.h"

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

}

int main() {
  clock_t timer = clock();
  int max_dim = 1; //! Actual max dim is 5*max_dim
  //! This is basically the C method of timing a program. We time globally, but
  //! also for eache dimension and for each size of integers in the matrix.
  clock_t sho_bkz[max_dim], tmp;
 // clock_t total_times[1];
  for (int i = 0; i < max_dim; i++){
     sho_bkz[i] = 0;
  }
  int bkz_fails=0;
 // Real vec_length[3];
 // vec_length[0] = vec_length[1] = vec_length[2] = 0;

  RealMat matShortVecLeng;
  matShortVecLeng.resize(10, 10);

  std::string prime = primes[0];

  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 1; k++) {
      // We dynamically allocate memory to these two pointers every time we need to
      // create an object of their type. This is because of the OOP approach
      // to lattice reduction.
      IntLatticeBase<Int, Real, RealRed>* basis;
      Reducer<Int, Real, RealRed>* red;

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1;
      unsigned int ln;
      
 
       std::string s4("UtilTriangular2");

      //name = "bench/" + prime+ "_4" + "_002" ;
      // name = "bench/" + prime+ "_5" + "_2" ;
         name = "bench/" + prime+ "_4" + "_001" ;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);


      // BKZ reduction before shortest vector search
      Int m(1021);
      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1,matrix1,m, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);

      std::cout << " The base before reduction\n"; 
      printBase((red->getIntLatticeBase())->getBasis()); 

     // red->redBKZ(0.9, 10, QUADRUPLE, 0);

      std::cout << " The base after reduction\n"; 
      printBase((red->getIntLatticeBase())->getBasis()); 

      basis->updateVecNorm();
    ///  vec_length[2] += average(basis->getVecNorm());
       std::cout << "Short vector in initial base " <<  red->getMinLength() << "\n";

      tmp = clock();
      if (!red->shortestVector(L2NORM,s4)) {
        bkz_fails++;
      }
      sho_bkz[j] += clock() - tmp;
      matShortVecLeng(j,k)=red->getMinLength();
      delete red;
      //std::cout << "BKZ: " << average(basis->getVecNorm()) << "\n";
      delete basis;
    }
  }

  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";
  for(int i=0;i<10;i++)
   {for(int j=0;j<10;j++){
       std::cout <<  matShortVecLeng(i,j)<< "   ";
     }
     std::cout << ""<< std::endl;
   }
  return 0;
}
