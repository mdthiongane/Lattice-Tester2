/**An example of program to test the speed of
 *Util::LowerTriangularization @author Lecuyer.
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

  std::string prime = primes[0];

  //! Variables definition
  ParamReader<Int, RealRed> reader;
  std::string name;
  int numlines;
  IntMat matrix1, matrix2, matrix3;
  unsigned int ln;
  IntLatticeBase<Int, Real, RealRed> *basis;
  Reducer<Int, Real, RealRed> *red;

  Int m(1021);

  std::cout << name << std::endl;
  name = "bench/" + prime + "_5" + "_2";

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

  basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, matrix1, m, numlines);
  red = new Reducer<Int, Real, RealRed>(*basis);
  // red->redLLL(0.999,1000000,numlines);
  basis->updateVecNorm();

  std::cout << " The base after reduction\n";

  printBase((red->getIntLatticeBase())->getBasis());
  copy(matrix2, (red->getIntLatticeBase())->getBasis());
  // Triangularization2<IntMat,IntVec, Int> (matrix2, matrix3, m);
  TriangularizationLower<IntMat, IntVec, Int>(matrix1, matrix2, m);
  printBase(matrix2);

  return 0;
}
