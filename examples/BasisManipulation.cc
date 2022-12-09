
/**An example of program to use the
 *Util::Triangularization @author Couture.
 **/
#define NTL_TYPES_CODE 2
#include <iostream>
#include <ctime>
#include <NTL/mat_GF2.h>
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


int main()
{
   IntLatticeBase<Int, Real, RealRed> *lattice;
  BasisConstruction<Int> constr; // The basis constructor we will use
 // Reducer<Int, Real, RealRed> *red;
  IntMat bas_mat, dua_mat;
  IntMat w_copie, m_v, m_v2;
  NTL::Mat<NTL::ZZ> m_dual, w_copie2;
  Int m(1021);
  clock_t tmps;

  std::string name = "bench/" + prime + "_5_1";
  ParamReader<Int, RealRed> reader(name + ".dat");

  reader.getLines();
  int numlines;
  unsigned int ln;
  reader.readInt(numlines, 0, 0);

  bas_mat.SetDims(numlines, numlines);
  dua_mat.SetDims(numlines, numlines);
    //for the upper triangular  basis
  m_v.SetDims(numlines, numlines);
  //for the lower triangular  basis
  m_v2.SetDims(numlines, numlines);
  //for the m-dual  basis
  m_dual.SetDims(numlines, numlines);
  //A copie of lattice basis
  w_copie.SetDims(numlines, numlines);
  //for the m-dual  basis
   w_copie2.SetDims(numlines, numlines);
 
  ln = 1;
 
  //! Filling the matrix
  reader.readBMat(bas_mat, ln, 0, numlines);

  // Creating a lattice basis
  lattice = new IntLatticeBase<Int, Real, RealRed>(bas_mat, bas_mat, m, numlines);
  std::cout << " The initial base\n";
  printBase(lattice->getBasis()); //The primal basi is 'bas_mat'

 // We copy the primal basis  in 'w_copie'
  copy(bas_mat, w_copie);
  //construction of the lower triangular basis 
  constr.upperTriangular(w_copie, m_v, m);
  std::cout << " Print the upper triangular basis \n";
  printBase(m_v);

  
  copy(bas_mat, w_copie);
  // The construction of the lower triangular basis 
  constr.lowerTriangular(w_copie, m_v2, m);
  std::cout << " Print the lower triangular basis \n";
  printBase(m_v2);

  std::cout << " The LLL reduction basis with delta=0.8 \n";
  copy(bas_mat, w_copie);
  constr.LLLConstruction(w_copie, 0.8);
  printBase(w_copie);

  std::cout << " The LLL reduction basis with delta=0.99 \n";
  copy(bas_mat, w_copie);
  constr.LLLConstruction(w_copie, 0.99999);
  printBase(w_copie);
 
 // create basis a with NTL::Mat<NTL::ZZ> matrix 
 // because in BasisConstruction class matrix is define with NTL::matrix<Int>,  
 // It is define with constant IntMat by the instruction :	
 //              typedef NTL::matrix<Int> IntMat;
 NTL::ZZ mm(1021);
   for(int i=0;i<numlines;i++){
    for(int j=0;j<numlines;j++)
        w_copie2[i][j]=NTL::conv<NTL::ZZ>(bas_mat[i][j]);
  }  
  
  copyMatrixToMat(bas_mat, w_copie2);
  //std::cout << " The basis initial \n";
  //printBase(bas_mat);
  std::cout << " The copy w_copie 'NTL::Mat' basis for which we compute the m-dual \n";
  printBase2(w_copie2);

  constr.calcDual(w_copie2, m_dual, mm);
  // Triangularization(m_v, m_v2, numlines, numlines, m);
  std::cout << " The m-dual of the primal basis \n";
  printBase2(m_dual);


  /* To compare the speed of BasisConstruction::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and BasisConstruction::calcDualUpperTriangular method which compute an m-dual basis
 * with an upper triangular basis. We triangularize the basis before calling
 * the Util::CalcDualUpper method.*/
  double tps = 0;
  
  //we work with a copy of bas_mat
  copy(bas_mat, w_copie);
  constr.upperTriangular(w_copie, m_v2, m);
  tmps = clock();
  for (int i = 0; i < 100; i++) {
    constr.calcDualUpperTriangular(m_v2, m_v, numlines, m);
  
  }
  tps = (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
  std::cout << " Time (in second) to compute 100 m-dual from upper triangular basis: " << tps << std::endl;

  tmps = clock();
  for (int i = 0; i < 100; i++){
     constr.calcDual(w_copie2, m_dual, m); 
  }
  tps = (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
  std::cout << " Time (in second) to compute 100 m-dual from non-traingular basis: " << tps << std::endl;

   
   /*
    * To compare the speed of triangular 'BasisConstruction::upperTriangular'
    * and the speed of 'BasisConstruction::LLLConstruction'
   */
  copy(w_copie, bas_mat);
  tps = 0;
  for (int i = 0; i < 100; i++){
    tmps = clock();
    constr.upperTriangular(w_copie, m_v2, m);
    tps = tps + (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
    copy(w_copie, bas_mat);
  }
  std::cout << " The triangular compute time: " << tps << std::endl;

  std::cout << " #############################################################\n";
  tps = 0;
  copy(w_copie, bas_mat);
  for (int i = 0; i < 500; i++) {
    tmps = clock();
    constr.LLLConstruction(w_copie, 0.99999);
    tps = tps + (double)(clock() - tmps) / (CLOCKS_PER_SEC * 60);
    copy(w_copie, bas_mat);
  }
  std::cout << " The LLL basis compute time: " << tps << std::endl;
 

  return 0;
}