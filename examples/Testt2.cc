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

namespace {
  const std::string prime = primes[0];
}


void printBase(IntMat bas_mat, int i, int j){
 
     for(int i=0;i<4;i++)
     {for(int j=0;j<4;j++){
       std::cout <<  bas_mat(i,j)<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}

int main() {
 // clock_t timer = clock();
  // The different clocks we will use for benchmarking
  // We use ctime for implementation simplicity
 // int max_dim = 1; // Actual max dim is 5*max_dim
  

  // Defining constants for the execution of the algorithms
  BasisConstruction<Int> constr; // The basis constructor we will use
  IntLatticeBase<Int, Real, RealRed> *lattice;
 // IntLatticeBase<Int, Real, RealRed> *m_latCopie; 
  Reducer<Int, Real, RealRed>* red;
  IntMat bas_mat, dua_mat;
  IntMat m_v,m_v2;
  Int m(1021); 
  //clock_t tmp;
 
      //! Reader shenanigans
     // std::string name = "bench/" + prime[0] + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      //std::string name = "bench/" + prime+ "_4" + "_001" ;
     // std::string name = "bench/" + prime+ "_4" + "_002" ;
      std::string name = "bench/"  + prime+"_5_0" ;
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
      lattice= new IntLatticeBase<Int, Real, RealRed>(bas_mat,bas_mat,m, numlines);
      //IntLatticeBase<Int, Real, RealRed> lattice(bas_mat,bas_mat,m, numlines);
      
      red = new Reducer<Int, Real, RealRed>(*lattice);
       
     
     std::cout << " The initial base\n"; 
     printBase(bas_mat, 4, 4);

     
       // BKZ reduction before shortest vector search
     red->redBKZ();

     std::cout << " The base after reduction\n"; 
     printBase((red->getIntLatticeBase())->getBasis(), 4, 4); 

      
     
    // Tringular GCD basis 
   
    
     constr.GCDTriangularBasis((red->getIntLatticeBase())->getBasis());
     std::cout << " The base after GCD triangularization\n";  
     printBase((red->getIntLatticeBase())->getBasis(), 4, 4);
     
     // Tringular GCD basis agian
     constr.GCDTriangularBasis((red->getIntLatticeBase())->getBasis()); 
     std::cout << " The base after second GCD triangularization\n";  
     printBase((red->getIntLatticeBase())->getBasis(), 4, 4);
     
    
 

  return 0;
}