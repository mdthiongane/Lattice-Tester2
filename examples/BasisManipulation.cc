
/**An example of program to use the
 *Util::Triangularization @author Couture.
 **/

//#define NTL_TYPES_CODE 2

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
  const std::string prime = primes[0];
}

/*void printBase(IntMat bas_mat)
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

  IntLatticeBase<Int, Real, RealRed> *lattice;
  // IntLatticeBase<Int, Real, RealRed> *m_latCopie;
  BasisConstruction<Int> constr; // The basis constructor we will use
  Reducer<Int, Real, RealRed> *red;
  IntMat bas_mat, dua_mat;
  IntMat w_copie, m_v, m_v2, m_dual;
  Int m(1021);

  std::string name = "bench/" + prime + "_5_0";
  ParamReader<Int, RealRed> reader(name + ".dat");

  reader.getLines();
  int numlines;
  unsigned int ln;
  reader.readInt(numlines, 0, 0);

  bas_mat.SetDims(numlines, numlines);
  dua_mat.SetDims(numlines, numlines);
  ln = 1;
  //! Filling the matrix
  reader.readBMat(bas_mat, ln, 0, numlines);

  // Creating a lattice basis
  lattice = new IntLatticeBase<Int, Real, RealRed>(bas_mat, bas_mat, m, numlines);

 // red = new Reducer<Int, Real, RealRed>(*lattice);

  std::cout << " The initial base\n";
  printBase(bas_mat);

  // BKZ reduction before shortest vector search
 // red->redBKZ();

  //std::cout << " The base after reduction\n";
 // printBase((red->getIntLatticeBase())->getBasis());

  //for the upper triangular  basis
  m_v.SetDims(numlines, numlines);
  //for the lower triangular  basis
  m_v2.SetDims(numlines, numlines);
  //for the m-dual  basis
  m_dual.SetDims(numlines, numlines);
  //A copie of lattice basis
  w_copie.SetDims(numlines, numlines);

  // copy the lattice in w_copie 
  // copy((red->getIntLatticeBase())->getBasis(), w_copie);
  copy(bas_mat, w_copie);

  std::cout << "Print the copie which is a copy of the basis \n";
  printBase(w_copie);

  constr.upperTriangular(w_copie, m_v, m);
 // Triangularization(w_copie, m_v, numlines, numlines, m);
  std::cout << " The upper triangular basis \n";
  printBase(m_v);
  //std::cout << " The base W \n";
  //printBase(red->getIntLatticeBase()->getBasis());
 // std::cout << " The copie W_copie after Triangularization \n";

  std::cout << " The construction of the lower triangular basis \n";
  constr.lowerTriangular(w_copie, m_v2, m);
  // Triangularization(m_v, m_v2, numlines, numlines, m);
  std::cout << " Print the lower triangular basis \n";
  printBase(m_v2);

  std::cout << " The construction of m-dual basis \n";
  constr.calcDual(w_copie, m_dual, m);
  // Triangularization(m_v, m_v2, numlines, numlines, m);
  std::cout << " Print the m-dual basis \n";
  printBase(m_v2);

  constr.LLLConstruction(lattice->getBasis());
 


  return 0;
}