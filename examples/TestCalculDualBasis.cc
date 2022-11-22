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
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_GF2.h>


#include "Examples.h"


using namespace LatticeTester;

namespace {
  const std::string prime = primes[0];
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


  IntMat bas_mat, bas_mat2, dua_mat;
  IntMat m_v,m_v2,m_v3;
  //Int m(101), G; 
  Int m(100), G, d; 
  IntVec vec, coeff,vl,tmp;
  std::string name = "bench/"  + prime+"_5_001" ;

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

  
     
             // We copy the base in m_v and m_v2
	     m_v.SetDims(numlines, numlines);
       m_v2.SetDims(numlines, numlines);
       m_v3.SetDims(numlines, numlines);
   
  
       std::cout << " The initial base\n"; 
       printBase(bas_mat);


       inv(d,m_v,bas_mat);
       std::cout << " The inverse base with determinanr d="<<d<<std::endl;  
       printBase(m_v);

       std::cout << " The inverse base 2 \n";  
      // m_v2=inv(bas_mat);
       printBase(m_v2); 
        

       CalcDual2 (bas_mat, m_v3, m) ;

       std::cout << " The m-dual basis \n"; 
       printBase(m_v3); 


  return 0;
}