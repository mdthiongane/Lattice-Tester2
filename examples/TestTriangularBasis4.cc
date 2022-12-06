/*
 *An example of program to use the
 *Util::TriangularizationLower  @author Lecuyer.
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

  IntLatticeBase<Int, Real, RealRed> *lattice;
  // IntLatticeBase<Int, Real, RealRed> *m_latCopie;
  Reducer<Int, Real, RealRed> *red;
  IntMat bas_mat, bas_mat2, dua_mat;
  IntMat m_v, m_v2;
  Int m(7);

  // std::string name = "bench/"  + prime+"_5_0" ;

  std::string name = "bench/" + prime + "_4" + "_001";
  ParamReader<Int, RealRed> reader(name + ".dat");

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

  std::cout << " print initial Base \n";
  printBase(bas_mat);

  // Creating a lattice basis
  lattice = new IntLatticeBase<Int, Real, RealRed>(bas_mat, bas_mat, m, numlines);
  // IntLatticeBase<Int, Real, RealRed> lattice(bas_mat,bas_mat,m, numlines);
  red = new Reducer<Int, Real, RealRed>(*lattice);
  std::cout << " The initial base\n";
  printBase(bas_mat);

  // BKZ reduction before shortest vector search
  red->redBKZ();

  std::cout << " The base after reduction\n";
  printBase((red->getIntLatticeBase())->getBasis());

  m_v.SetDims(numlines, numlines);
  m_v2.SetDims(numlines, numlines);

  // copy base to m_v
  copy((red->getIntLatticeBase())->getBasis(), m_v);

  std::cout << " Print the copy basis" << std::endl;
  printBase(m_v);

  // Triangularization2<IntMat,IntVec, Int> (m_v, m_v2, m);
  TriangularizationLower<IntMat, IntVec, Int>(m_v, m_v2, m);

  std::cout << " Print Base after triangularization:" << std::endl;

  printBase(m_v2);

  return 0;
}