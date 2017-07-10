//
//  main.cpp
//  Lattice Tester
//
//  Created by Erwan Bourceret on 18/04/2017.
//  Copyright © 2017 DIRO. All rights reserved.
//
/*
 * The purpose of this example is to compare the time
 *  of each reduction method. We compute several
 * random matrix and apply algorithms. The parameters of this
 * analysis are defined at the beginning of this file,
 * after include part. The output is a graphic generated by R.
 * The User need R distribution and the header RInside.h
 */


// Include Header
#include <iostream>
#include <map>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

// Include LatticeTester Header
#include "latticetester/Util.h"
#include "latticetester/Const.h"
#include "latticetester/Types.h"
#include "latticetester/IntFactor.h"
#include "latticetester/IntLatticeBasis.h"
#include "latticetester/Reducer.h"
#include "latticetester/Types.h"
#include "latticetester/ParamReader.h"
#include "latticetester/LatticeTesterConfig.h"

// Include NTL Header
#include <NTL/tools.h>
#include <NTL/ctools.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "NTL/vec_ZZ.h"
#include "NTL/vec_ZZ_p.h"
#include <NTL/vec_vec_ZZ.h>
#include <NTL/vec_vec_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>

// Include Boost Header
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>

// Include Random Generator of MRG Matrix an tools
#include "SimpleMRG.h"
#include "Tools.h"

using namespace std;
using namespace NTL;
using namespace LatticeTester;

int main (int argc, char *argv[])
{

   if (argc < 2) {
      cerr << "\n*** Usage:\n   "
           << argv[0] << " data_file1 data_file2 ...." << endl
           << "or\n   "
           << argv[0] << " dir1 dir2 ...." << endl
           << endl;
      return -1;
   }


   for (int j = 1; j < argc; j++) {
      // Do the test for each data file or directory on the command line

      //if (0 != S_ISDIR(buf.st_mode))         // directory
      //   status |= testall.doTestDir (argv[j]);
      string fname (argv[j]);
      fname += ".dat";
      ParamReader paramRdr (fname.c_str ());
      fname.clear ();
      LatticeTesterConfig config;
      paramRdr.read (config);

      config.write();
   }


   return 0;
}

