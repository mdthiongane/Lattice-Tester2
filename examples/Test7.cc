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
  const std::string prime = primes[1];
}


void printBase(IntMat bas_mat){
     int L=bas_mat.NumRows();
     int C=bas_mat.NumCols();

     for(int i=0;i<L;i++)
     {for(int j=0;j<C;j++){
       std::cout <<  bas_mat(i,j)<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}


void copy(IntMat &b1, IntMat &b2){
 
     for(int i=0;i<b1.size1();i++)
     { for(int j=0;j<b1.size2();j++){
          b2(i,j)=b1(i,j);
         }   
     }

}



int main() {

  IntLatticeBase<Int, Real, RealRed> *lattice;
  //IntLatticeBase<Int, Real, RealRed> *m_latCopie; 
  Reducer<Int, Real, RealRed>* red;
  IntMat bas_mat, bas_mat2, dua_mat;
  IntMat m_v,m_v2;
  //Int m(101), G; 
   Int m(1021), G; 
  IntVec vec, coeff,vl,tmp;
  //int pc=0,pl=0,K;

 
 
      //std::string name = "bench/" + prime+ "_4" + "_001" ;
     // std::string name = "bench/" + prime+ "_4" + "_002" ;
      std::string name = "bench/"  + prime+"_5_0" ;
     // std::string name = "bench/" + prime+ "_2" + "_001" ;
    //  std::string name = "bench/" + prime+ "_10" + "_1" ;
      ParamReader<Int, RealRed> reader(name + ".dat");
     std::cout <<name<<std::endl; 
                      
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
      lattice= new IntLatticeBase<Int, Real, RealRed>(bas_mat,bas_mat,m, numlines);
      //IntLatticeBase<Int, Real, RealRed> lattice(bas_mat,bas_mat,m, numlines);
      
       red = new Reducer<Int, Real, RealRed>(*lattice);
       
   
  
       std::cout << " The initial base\n"; 
       printBase(bas_mat);

     
       // BKZ reduction before shortest vector search
        red->redBKZ();

        std::cout << " The base after reduction\n"; 
        printBase((red->getIntLatticeBase())->getBasis()); 

       // We copy the base in m_v and m_v2
	        m_v.SetDims(numlines, numlines);
          m_v2.SetDims(numlines, numlines);
		
         //copy base to m_v
         copy((red->getIntLatticeBase())->getBasis(), m_v);
    

       
        std::cout << " The base m_v before triangularization\n";  
        printBase(m_v);
       
       
        //Triangularization of m_v
       
       // Triangularization2<IntMat,IntVec, Int> (m_v, m_v2, m);
         Triangularization(m_v ,m_v2, numlines,numlines,m);

        std::cout << " The base m_v after  triangularization2\n";  
        printBase(m_v);
        std::cout << " The base i \n";  
 

      
     //   Triangularization2<IntMat,IntVec, Int> (m_v, m, vec,  coeff, vl, G, K,tmp,pc,pl);
        std::cout << " The base m_v2 after second  triangularization\n";  
        printBase(m_v2);
    

  return 0;
}