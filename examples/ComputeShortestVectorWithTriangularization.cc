/*
**An example programm to compute a shortest vector using 
* a triangular decomposition method. The implemented method
* is discribe in LatticeGuide section 8.2. 
* The basis have 5x5 dimension, the modulo 'm' is 1021.
 * The name of dat file is '1021_5_2.dat'. It is the 'examples' folder of LatticeTester.
 * The absolue path is 'examples/bench/1021_5_2'.  
*/

#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>
#include "latticetester/Types.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/Reducer.h"
#include "latticetester/Const.h"
#include "Examples.h"
#include "latticetester/WriterRes.h"



using namespace LatticeTester;
namespace
{
  const std::string prime = primes[0];
}


int main() {
 // clock_t total_times[1];
  int bkz_fails=0;
  std::string prime = primes[0];
         // to lattice reduction.
      IntLatticeBase<Int, Real, RealRed>* basis;
      Reducer<Int, Real, RealRed>* red;

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1,matrix2;
      unsigned int ln;
      std::string s1("triangular");
      name = "bench/" + prime+ "_5" + "_2" ;
      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);
      Int m(1021);
      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1,matrix2,m, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);
      std::cout << " The base before reduction\n"; 
       printBase((red->getIntLatticeBase())->getBasis()); 
       //The BKZ reduction 
       red->redBKZ(0.9, 10, QUADRUPLE, 0);
       std::cout << " The base after reduction\n"; 
       printBase((red->getIntLatticeBase())->getBasis()); 
       basis->updateVecNorm();
       if (!red->shortestVector(L2NORM,s1)) {
         bkz_fails++;
       }
  
       std::cout << "Short vector length " <<  red->getMinLength() << "\n";

  return 0;
}
