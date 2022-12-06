/**An example of program to test the speed of
 *Util::Triangularization @author Couture.
 * We use 150 basis. We begin with 5x5 dimension
 *to 75x75 dimension. 10 different basis for each dimension
 **/

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
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"
#include "latticetester/NTLWrap.h"
#include "Examples.h"

using namespace LatticeTester;

namespace
{
  const std::string prime = primes[0];
}



int main()
{
  clock_t timer = clock();
  int max_dim = 15; //! Actual max dim is 5*max_dim
  clock_t tmp;

  RealMat matTriangularTime;
  matTriangularTime.resize(15, 10);

  std::string prime = primes[2];

  for (int j = 0; j < max_dim; j++)
  {
    for (int k = 0; k < 10; k++)
    {

      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1, matrix2;
      unsigned int ln;
      // BasisConstruction<Int> constr;
      Int m(1021);
      // std::string s2("GCDTriangular");

      //! Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5 * (j + 1)) + "_" + std::to_string(k);
      std::cout << name << std::endl;
      // name = "bench/" + prime+ "_2" + "_001" ;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      matrix2.SetDims(numlines, numlines);

      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      tmp = clock();

      Triangularization(matrix1, matrix2, numlines, numlines, m);

      double tps = (double)(clock() - tmp); //(CLOCKS_PER_SEC);
      matTriangularTime(j, k) = tps;
      // delete constr;
    }
  }

  std::cout << "Total time: " << (double)(clock() - timer) / (CLOCKS_PER_SEC * 60) << " minutes\n";
  for (int i = 0; i < 15; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      std::cout << matTriangularTime(i, j) << "   ";
    }
    std::cout << "" << std::endl;
  }

  return 0;
}
