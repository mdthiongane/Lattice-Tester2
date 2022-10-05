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

#include "latticetester/Const.h"
#include "latticetester/NTLWrap.h"



#include "Examples.h"

using namespace LatticeTester;

namespace {
  const std::string prime = primes[0];
}


void printBase(IntMat bas_mat, int lin, int col){
 
     for(int i=0;i<lin;i++)
     {for(int j=0;j<col;j++){
       std::cout <<  bas_mat(i,j)<< "   ";
     }
      std::cout << ""<< std::endl;
     }

}

void printVector(IntVec vec){
 
     for(int i=0;i<vec.length();i++)
    
       std::cout <<  vec[i]<< "    ";  
     std::cout << " "<<std::endl;            
}


void getMatColumnVec(IntMat mat,int lin, int col, int numCol, int pos, IntVec &vec){
  int k=0;
  vec.SetLength(lin-pos);
  for(int i=pos;i<lin;i++)
    vec[k++]=mat(i,numCol);    
}

void getMatRowVec(IntMat mat,int lin, int col, int numRo, int pos, IntVec &vec){
  int k=0;
  vec.SetLength(col-pos);
  for(int i=pos;i<col;i++)
    vec[k++]=mat(numRo,i);    
}


int main() {
 
  IntMat bas_mat, dua_mat;
  IntMat m_v,m_v2;
  Int m(1021); 
  //clock_t tmp;
 
      //! Reader shenanigans
     // std::string name = "bench/" + prime[0] + "_" + std::to_string((j+1)*5) + "_" + std::to_string(k);
      //std::string name = "bench/" + prime+ "_4" + "_001" ;
     // std::string name = "bench/" + prime+ "_4" + "_002" ;
     // std::string name = "bench/"  + prime+"_5_0" ;

      std::string name = "bench/" + prime+ "_4" + "_003" ;
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
      
   
     IntVec l, c;
     Int mod(1021);
    // getMatColumnVec<std::vector<std::int64_t,>(bas_mat,numlines, numlines, 0, 0, c);
     getMatRowVec(bas_mat,numlines, numlines, 0,0, l);

     std::cout << " print Base \n"; 
     printBase(bas_mat, 4, 4);


     std::cout << " print vector colonne \n"; 
     printVector(c);     

     Int G;
     IntVec coeff;
     IntVec vl;

     //coeff.SetLength(4);
     GCDvect (c, coeff, G);
     
     //std::vector<std::int64_t> t=c;
    // GCD2vect (c, 0, c.length());
  
     
     //vl.SetLength(4);
    std::cout << " G="<<G<<std::endl; 
     std::cout << " Coeeficient: \n"; 
     printVector(coeff);    

    // for(int i=0;i<=4;i++)
     //  coeff[i]=4-i;
     int x=0;
     computeVl(bas_mat, coeff,  vl,x,mod) ;

     std::cout << " vecteur vl:"<<std::endl; 
     printVector(vl); 

     int K=0; 
     int pl=0;
     int pc=0;   
     replaceVector(bas_mat,vl,pl, pc,K) ;
     std::cout << " La ligne remplacer K="<<K<<std::endl; 

    std::cout << " print Base after set vl in first line \n"; 
     printBase(bas_mat, 4, 4);
     //int pl=0; 
     //int pc=0;
     updateMatrive(bas_mat, vl, pl,pc, G ,mod);
     
     //std::cout << " Print matrice relplace first line:"<<std::endl; 

     // swapVector<IntMat,IntVec>(bas_mat, 0 , 3);
     
    std::cout << " Print matrice after update:"<<std::endl; 

     printBase(bas_mat, 4, 4);



    IntVec tmp;
    int l1=0,l2=3;
    swapVector(bas_mat, l1 , l2, tmp);

   std::cout << " Print matrice after swap l=0 et l=3:"<<std::endl; 

     printBase(bas_mat, 4, 4);

    /***
     coeff.SetLength(c.length());
     GCDvect (c, coeff,  G) ;

     std::cout << "G="<<G<<std::endl; 
     std::cout << " column vec"<<std::endl; 
     printVector(c); 
     std::cout << " coefficient:"<<std::endl; 
     printVector(coeff); 
     Int som(0);
     for(int i=0;i<c.length();i++)
        som=som+ c[i]*coeff[i];
     std::cout << "Som="<<som<<std::endl;   
     **/ 

  return 0;
}