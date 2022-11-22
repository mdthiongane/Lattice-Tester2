//An example of program to test the speed of triangularization
//with Util::Triangularization (Couture method)

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
  //! This is basically the C method of timing a program. We time globally, but
  //! also for eache dimension and for each size of integers in the matrix.
  clock_t tmp;
 // clock_t total_times[1];
 // for (int i = 0; i < max_dim; i++){
 //    sho_bkz[i] = 0;
 // }
 // int bkz_fails=0;
 // Real vec_length[3];
 // vec_length[0] = vec_length[1] = vec_length[2] = 0;

  RealMat matTriangularTime;
  matTriangularTime.resize(15, 10);

  std::string prime = primes[2];

  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 10; k++) {
      // We dynamically allocate memory to these two pointers every time we need to
      // create an object of their type. This is because of the OOP approach
      // to lattice reduction.
    //  IntLatticeBase<Int, Real, RealRed>* basis;
     // Reducer<Int, Real, RealRed>* red;

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1,matrix2;
      unsigned int ln;
     // BasisConstruction<Int> constr;
      Int m(1021);
    // std::string s2("GCDTriangular");

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

        Triangularization(matrix1 ,matrix2, numlines,numlines,m);

      double tps=(double)(clock() - tmp); //(CLOCKS_PER_SEC);
      matTriangularTime(j,k)=tps;
     // delete constr;
  
    }
  }

  //! Printing the results in a somewhat formated way.
/*  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  std::cout << "          ";
  int width6 = getWidth(sho_bkz, max_dim, "SV BKZ and Util Trangularization", total_times, 0);
  std::cout << std::endl;

  std::cout << "Total time" << std::setw(width6) << total_times[0]
     << total_times[5] << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 
      << std::setw(width6) << sho_bkz[i] 
      << std::endl;
  }
  std::cout << "Fails     "
    << std::setw(width6) << bkz_fails 
    << std::endl;
**/


  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";
  for(int i=0;i<15;i++)
   {for(int j=0;j<10;j++){
       std::cout << matTriangularTime(i,j)<< "   ";
     }
     std::cout << ""<< std::endl;
   }
 
  return 0;
}
