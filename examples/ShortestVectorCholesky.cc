/*An example programm to compute a basis shortest vector using the
**Cholesky decomposition.
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


/*namespace
{
  // Returns the average of the length of this vector
  Real average(RealVec vector)
  {
    Real sum(0);
    for (int i = 0; i < vector.length(); i++)
    {
      sum += vector[i];
    }
    return sum / Real(vector.length());
  }

  void printBase(IntMat bas_mat)
  {
    int l = bas_mat.size1();
    int c = bas_mat.size2();
    for (int i = 0; i < l; i++)
    {
      for (int j = 0; j < c; j++)
      {
        std::cout << bas_mat(i, j) << "   ";
      }
      std::cout << "" << std::endl;
    }
  }

}*/

int main()
{
  clock_t timer = clock();
  int max_dim = 1; //! Actual max dim is 5*max_dim

  int bkz_fails = 0;

  RealMat matShortVecLeng;
  matShortVecLeng.resize(10, 10);

  std::string prime = primes[0];

  for (int j = 0; j < max_dim; j++)
  {
    for (int k = 0; k < 1; k++)
    {

      IntLatticeBase<Int, Real, RealRed> *basis;
      Reducer<Int, Real, RealRed> *red;

      //! Variables definition
      ParamReader<Int, RealRed> reader, reader2;
      std::string name;
      int numlines;
      IntMat matrix1;
      unsigned int ln;

      std::string s1("cholesky");

      name = "bench/" + prime + "_5" + "_2";
      //  name = "bench/" + prime+ "_2" + "_002" ;

      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);

      // BKZ reduction before shortest vector search
      Int m(1021);
      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, matrix1, m, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);

      std::cout << " The base before reduction\n";
      printBase((red->getIntLatticeBase())->getBasis());

      // red->redBKZ();
      red->redLLL(0.999, 1000000, numlines);
      basis->updateVecNorm();

      std::cout << " The base after reduction\n";
      printBase((red->getIntLatticeBase())->getBasis());

      if (!red->shortestVector(L2NORM, s1))
      {
        bkz_fails++;
      }

      matShortVecLeng(j, k) = red->getMinLength();
      delete red;
      delete basis;
    }
  }

  //! Printing the results in a somewhat formated way.

  std::cout << "Total time: " << (double)(clock() - timer) / (CLOCKS_PER_SEC * 60) << " minutes\n";
  for (int i = 0; i < 10; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      std::cout << matShortVecLeng(i, j) << "   ";
    }
    std::cout << "" << std::endl;
  }

  return 0;
}
