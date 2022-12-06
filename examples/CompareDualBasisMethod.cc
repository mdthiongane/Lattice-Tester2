
/* The code compare the speed of Util::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and Util::CalcDualUpper method which compute an m-dual basis
 * with an upper triangular basis. We triangularize the basis before calling
 * the Util::CalcDualUpper method.
 * The used basis is in the 'examples/bench' folder of LatticeTester. 
 * In 'bench' sub-folder  of 'examples'you can find many basis with different dimension. 
 * Each basis in a '.dat' file.
 * The dimension  of 5x5 to 75x75. Here is the contain of the contain of file '1021_5_1.dat'
 * 
 *  5
 *  30 203 244 44 626 
 *  965 905 182 890 975 
 *  245 882 657 654 335 
 *  794 232 968 807 535 
 *  107 607 513 451 455 
 * 
 * This file contain information of 5x5 basis. There are six line of data. 
 * The first give the dimension of the basis (5 that mean 5x5 basis). Line 2 to line
 * 6 give the basis vector
 * 
 *  
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

using namespace LatticeTester;

namespace
{
  const std::string prime = primes[1];
}

/*void printBase(IntMat bas_mat)
{
  int L = bas_mat.NumRows();
  int C = bas_mat.NumCols();

  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < C; j++)
    {
      std::cout << bas_mat(i, j) << "   ";
    }
    std::cout << "" << std::endl;
  }
}

void copy(IntMat &b1, IntMat &b2)
{

  for (int i = 0; i < b1.size1(); i++)
  {
    for (int j = 0; j < b1.size2(); j++)
    {
      b2(i, j) = b1(i, j);
    }
  }
}*/

int main()
{

  //  clock_t timer = clock();
  clock_t tmps;
  IntLatticeBase<Int, Real, RealRed> *lattice;
  // IntLatticeBase<Int, Real, RealRed> *m_latCopie;
  Reducer<Int, Real, RealRed> *red;
  IntMat bas_mat, bas_mat2, dua_mat;
  IntMat m_v, m_v2, m_v3, m_v4;
  // Int m(101), G;
  Int m(1048573);
  //  Int m(1021);
  Int G;
  IntVec vec, coeff, vl, tmp;

  std::string name = "bench/" + prime + "_30_4";
  ParamReader<Int, RealRed> reader(name + ".dat");
  std::cout << name << std::endl;

  reader.getLines();
  int numlines;
  unsigned int ln;
  reader.readInt(numlines, 0, 0);

  bas_mat.SetDims(numlines, numlines);
  bas_mat2.SetDims(numlines, numlines);
  dua_mat.SetDims(numlines, numlines);
  ln = 1;
  //! Filling the matrix
  reader.readBMat(bas_mat, ln, 0, numlines);

  // Creating a lattice basis
  lattice = new IntLatticeBase<Int, Real, RealRed>(bas_mat, bas_mat, m, numlines);
  // IntLatticeBase<Int, Real, RealRed> lattice(bas_mat,bas_mat,m, numlines);

  red = new Reducer<Int, Real, RealRed>(*lattice);

  // BKZ reduction before shortest vector search
  //   red->redBKZ();

  // We copy the base in m_v and m_v2
  m_v.SetDims(numlines, numlines);
  m_v2.SetDims(numlines, numlines);
  m_v3.SetDims(numlines, numlines);
  m_v4.SetDims(numlines, numlines);

  // copy base to m_v
  copy((red->getIntLatticeBase())->getBasis(), m_v);
  copy((red->getIntLatticeBase())->getBasis(), m_v3);

  double tps = 0;
  Triangularization2<IntMat, IntVec, Int>(m_v, m_v2, m);
  tmps = clock();
  for (int i = 0; i < 200; i++)
  {
    calcDual(m_v2, m_v3, numlines, m);
    //  m_v2(i,i)= m_v2(i,i)+1;
  }
  tps = (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
  std::cout << " Time clock calcul m-dual with upper triangular basis: " << tps << std::endl;

  tmps = clock();
  for (int i = 0; i < 200; i++)
  {
    CalcDual2(m_v3, m_v4, m);
    //   m_v2(i,i)= m_v2(i,i)+1;
  }
  tps = (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
  std::cout << " Time clock calcul m-dual without traingular basis: " << tps << std::endl;

  return 0;
}