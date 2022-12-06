// An example of program to test the speed of triangularization
// with Util::LowerTriangularization (Lecuyer method)

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
  const std::string prime = primes[1];
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

  void printRes(RealMat mat, int lin, int col)
  {
    std::ofstream out("ResultatUtilTriang.csv");

    for (int j = 0; j < lin; j++)
    {
      for (int k = 0; k < col; j++)
        out << j << ',';
      out << '\n';
    }
    out.close();
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

}
*/

int main()
{

  //  clock_t timer = clock();
  IntLatticeBase<Int, Real, RealRed> *basis;
  BasisConstruction<Int> constr; // The basis constructor we will use
  Reducer<Int, Real, RealRed> *red;
  clock_t tmps;
  std::string prime = primes[0];

  //! Variables definition
  ParamReader<Int, RealRed> reader;
  std::string name;
  int numlines;
  IntMat matrix1, matrix2, matrix3;
  unsigned int ln;


  Int m(1021);

  std::cout << name << std::endl;
  name = "bench/" + prime + "_75" + "_2";

  reader = ParamReader<Int, RealRed>(name + ".dat");
  reader.getLines();
  reader.readInt(numlines, 0, 0);
  matrix1.SetDims(numlines, numlines);
  matrix2.SetDims(numlines, numlines);
  matrix3.SetDims(numlines, numlines);

  ln = 1;
  reader.readBMat(matrix1, ln, 0, numlines);

  // std::cout << " The base before reduction\n";
  // printBase(matrix1);

  basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, matrix1, m, numlines);
  red = new Reducer<Int, Real, RealRed>(*basis);
  // red->redLLL(0.999,1000000,numlines);
  basis->updateVecNorm();

  //  std::cout << " The base after reduction\n";

  // printBase((red->getIntLatticeBase())->getBasis());
  copy(matrix2, (red->getIntLatticeBase())->getBasis());
  double tps = 0;
  for (int i = 0; i < 500; i++)
  {
    tmps = clock();
    TriangularizationLower<IntMat, IntVec, Int>(matrix1, matrix2, m);
    tps = tps + (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
    copy(matrix2, (red->getIntLatticeBase())->getBasis());
  }
  std::cout << " The triangular compute time: " << tps << std::endl;

  std::cout << " #############################################################\n";
  tps = 0;
  for (int i = 0; i < 500; i++)
  {
    tmps = clock();
    // TriangularizationLower<IntMat,IntVec,Int>(matrix1, matrix2 ,m);
    red->redLLL(0.999, 1000000, numlines);
    tps = tps + (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
    // CopyMatr(matrix2,(red->getIntLatticeBase())->getBasis(), numlines, numlines);
    red = new Reducer<Int, Real, RealRed>(*basis);
  }
  std::cout << " The LLL basis compute time: " << tps << std::endl;

  return 0;
}
