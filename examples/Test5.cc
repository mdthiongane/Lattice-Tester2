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

  IntLatticeBase<Int, Real, RealRed> *lattice;
  IntLatticeBase<Int, Real, RealRed> *m_latCopie; 
  Reducer<Int, Real, RealRed>* red;
  IntMat bas_mat, dua_mat;
  IntMat m_v,m_v2;
  Int m(7), G; 
  IntVec vec, coeff,vl,tmp;
  int pc=0,pl=0,K;

 
 
      //std::string name = "bench/" + prime+ "_4" + "_001" ;
      std::string name = "bench/" + prime+ "_4" + "_002" ;
      //  std::string name = "bench/"  + prime+"_5_0" ;
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
	      m_v.SetDims(numlines, numlines);
        m_v2.SetDims(numlines, numlines);
		    m_latCopie = new IntLatticeBase<Int, Real, RealRed>((red->getIntLatticeBase())->getBasis(),(red->getIntLatticeBase())->getBasis(),
        (red->getIntLatticeBase())->getModulo(),(red->getIntLatticeBase())->getDim());
             // Tringular basis from util
        // Triangularization(m_latCopie->getBasis() ,m_v, numlines,numlines,m_latCopie->getModulo());
         Triangularization2<IntMat,IntVec, Int> (m_latCopie->getBasis(), m, vec,  coeff, vl, G, K,tmp,pc,pl);

        std::cout << " The base V after Util triangularization\n";  
        printBase(m_v, 4, 4);
         std::cout << " The base W \n";  
        printBase(red->getIntLatticeBase()->getBasis(), 4, 4);
         std::cout << " The copie W' after Triangularization \n";  
        printBase(m_latCopie->getBasis(), 4, 4);

       // Triangularization2(m_v,m_v2, numlines,numlines,m_latCopie->getModulo());
        Triangularization2<IntMat,IntVec, Int> (m_v, m, vec,  coeff, vl, G, K,tmp,pc,pl);
        std::cout << " The base after second Util triangularization\n";  
         printBase(m_v, 4, 4);
    

 
  //std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";

  return 0;
}