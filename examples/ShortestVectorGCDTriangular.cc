/**
 * This is an example to the usage of the reducer class to perform reduction
 * of lattices. This serves as a comparison between the different reduction
 * methods when used before searching for the shortest vector. The output is
 * formated to include the number of clock ticks spent on each algorithm as well
 * as the number of times the program failed to find the shortest vector for
 * each pre-reduction.
 *
 * This is an example ouput for the program:
 * ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS
 *             Dieter    LLL     BKZ   SV Dieter     SV LLL     SV BKZ
 * Total time 7479736 804050 2414321 11960009667 2665878917 2444318476
 * Dim      5    3889    873     945         634        583        579
 * Dim     10   25622   3104    4025        3012       2603       2582
 * Dim     15   52247   6485   10449        8268       6081       6023
 * Dim     20   94343  12015   23153       20410      12513      12014
 * Dim     25  137098  22926   48952      262480      36137      29878
 * Dim     30  195888  30963   87739      604659     141528     105142
 * Dim     35  350096  44862  149397     9703739    1469959     968583
 * Dim     40  468672  54571  176943    28913921    5581257    3571159
 * Dim     45  619861  80332  229724   607680675   66124526   53250411
 * Dim     50  839252  98698  269042  1260636786  110130393  102112537
 * Dim     55 1119413 125894  364620  2351365585  186791907  204776212
 * Dim     60 1591436 120304  420965  3423588351  416827286  371013095
 * Dim     65 1981919 203023  628367  4277221147 1878754144 1708470261
 * Fails           14      2       2
 * Total time: 284.743 minutes
 * */

// We define the numeric types.
// It is possible to use this example with TYPES 2 and 3. For now 1 calls the
// same function for both execution and we look forward to change that.
#define NTL_TYPES_CODE 2

#include <iostream>
#include <ctime>

#include "latticetester/ParamReader.h"
#include "latticetester/Types.h"
#include "latticetester/Reducer.h"
#include "latticetester/IntLatticeBase.h"
#include "latticetester/WriterRes.h"

#include "Examples.h"

using namespace LatticeTester;

namespace {
  // Returns the average of the length of this vector
  Real average(RealVec vector) {
    Real sum(0);
    for (int i = 0; i<vector.length(); i++) {
      sum += vector[i];
    }
    return sum/Real(vector.length());
  }
}

int main() {
  clock_t timer = clock();
  int max_dim = 3; //! Actual max dim is 5*max_dim
  //! This is basically the C method of timing a program. We time globally, but
  //! also for eache dimension and for each size of integers in the matrix.
  clock_t sho_bkz[max_dim], tmp;
  clock_t total_times[1];
  for (int i = 0; i < max_dim; i++){
     sho_bkz[i] = 0;
  }
  int bkz_fails=0;
  Real vec_length[3];
  vec_length[0] = vec_length[1] = vec_length[2] = 0;

  std::string prime = primes[0];

  for (int j = 0; j < max_dim; j++) {
    for (int k = 0; k < 1; k++) {
      // We dynamically allocate memory to these two pointers every time we need to
      // create an object of their type. This is because of the OOP approach
      // to lattice reduction.
      IntLatticeBase<Int, Real, RealRed>* basis;
      Reducer<Int, Real, RealRed>* red;

      //! Variables definition
      ParamReader<Int, RealRed> reader;
      std::string name;
      int numlines;
      IntMat matrix1;
      unsigned int ln;
 
      
       std::string s2("GCDTriangular");

      //! Reader shenanigans
      name = "bench/" + prime + "_" + std::to_string(5*(j+1)) + "_" + std::to_string(k);
      std::cout << name<<std::endl;
      reader = ParamReader<Int, RealRed>(name + ".dat");
      reader.getLines();
      reader.readInt(numlines, 0, 0);
      matrix1.SetDims(numlines, numlines);
      ln = 1;
      reader.readBMat(matrix1, ln, 0, numlines);


      // BKZ reduction before shortest vector search
      basis = new IntLatticeBase<Int, Real, RealRed>(matrix1, numlines);
      red = new Reducer<Int, Real, RealRed>(*basis);
      red->redBKZ();
      basis->updateVecNorm();
      vec_length[2] += average(basis->getVecNorm());
      tmp = clock();
      if (!red->shortestVector(L2NORM,s2)) {
        bkz_fails++;
      }
      sho_bkz[j] += clock() - tmp;
      delete red;
      //std::cout << "BKZ: " << average(basis->getVecNorm()) << "\n";ss
      delete basis;
    }
  }

  //! Printing the results in a somewhat formated way.
  std::cout << "ALL THE RESULTS ARE NUMBERED IN TERMS OF SYSTEM CLOCK TICKS\n";
  std::cout << "          ";
  int width6 = getWidth(sho_bkz, max_dim, "SV BKZ and Cholesky", total_times, 0);
  std::cout << std::endl;

  std::cout << "Total time" << std::setw(width6) << total_times[0]
     << total_times[5] << std::endl;
  for (int i = 0; i < max_dim; i++) {
    std::cout << "Dim " << std::setw(6) << (i+1)*5 
      << std::setw(width6) << sho_bkz[i] 
      << std::endl;
  }
  std::cout << "Fails     "
    << std::setw(width6) << bkz_fails 
    << std::endl;

  std::cout << std::fixed << std::setprecision(2) << "Averages: " << vec_length[0]/vec_length[1]
    << 1.0 << std::setw(width6) << vec_length[2]/vec_length[1]
    <<std::endl;

  std::cout << "Total time: " << (double)(clock()-timer)/(CLOCKS_PER_SEC*60) << " minutes\n";
  
  return 0;
}
