/**
ParamReader.cc for ISO C++
version
modified 28/05/06 10:36
authors: Hicham Wahbi
       Frederik Rozon
       Richard Simard
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cassert>
#include <cstdlib>

#include "latticetester/Types.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/LatticeTesterConfig.h"

using namespace std;
using namespace LatticeTester;


namespace LatticeTester
{

//===========================================================================

ParamReader::ParamReader()
{
   m_fileName.reserve(MAX_WORD_SIZE);
}


//===========================================================================

ParamReader::ParamReader(string fileName)
{
   m_fileName.reserve(fileName.length());
   m_fileName = fileName;
}


//===========================================================================

ParamReader::~ParamReader()
{
   for (int i = 0; i < (int) m_lines.size(); i++)
      m_lines[i].clear();
   m_lines.clear();
}


//===========================================================================

void ParamReader::getLines()
{
   ifstream inFile(m_fileName.c_str());
   if (inFile.fail()) {
      cerr << "An error occurred. Unable to read input file:" <<
               m_fileName << endl;
      exit(1);
   }
   m_lines.clear();
   const string blancs = " \t";
   string line;
   string::size_type i;

   while (getline(inFile, line, '\n')) {
      i = line.find_first_not_of(blancs);  // find first non blank char
      if (i != string::npos && line[i] != '#')  // comment line begins with a #
         m_lines.push_back(line);
   }

   inFile.close();
}

//===========================================================================

void ParamReader::getToken(string& field, unsigned int ln, unsigned int pos)
{
   std::vector<string> tokens;
   tokens.reserve(20);
   tokens.clear();

   tokenize(tokens, ln - 1);
   if (pos > tokens.size()) {
      cerr << "Warning: position " << pos << " exceeds number of params at line "
      << ln << endl;
      field.clear();
   } else
      field = tokens[pos];
   return;
}


//===========================================================================

int ParamReader::tokenize(std::vector<string>& tokens, unsigned int ln)
{
   if (ln >= m_lines.size()) {
      cerr << "Warning: line " << ln << " doesn't exist in param file\n";
      exit(0);
   }
   string str = m_lines[ln];
   string word;
   int wnum = 0;
   bool firstFound = false;
   char c = 0;
   for (unsigned int i = 0; i < str.length(); i++) {
      c = str[i];
      if (!isDelim(c)) {
         string tmp(1, c);
         word.append(tmp);
         //tokens[wnum] += c;
         firstFound = true;

      } else {
         if (firstFound) {
            tokens.push_back(word);
            word.clear();
            wnum++;
         }
         firstFound = false;
      }
   }

   if (firstFound) {
      tokens.push_back(word);
      word.clear();
      wnum++;
   }
   if (isDelim(c))
      wnum--;
   return wnum;
}


//===========================================================================

bool ParamReader::isDelim(char c)
{
   const string delim = " ?!,\t\n";
   bool bRetVal = 0;
   for (int i = 0; delim[i] != 0; ++i) {
      if (delim[i] == c) {
         bRetVal = 1;
         break;
      }
   }

   return bRetVal;
}


//===========================================================================

void ParamReader::readBool(bool & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (strcasecmp(val.c_str(), "true") == 0)
      field = true;
   else if (strcasecmp(val.c_str(), "false") == 0)
      field = false;
   else
      MyExit(1, "readBool:   NO SUCH CASE");
}


//===========================================================================

void ParamReader::readString(string & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = string(val);
}


//===========================================================================

void ParamReader::readChar(char & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = val[0];
}


//===========================================================================

void ParamReader::readNumber3 (MScal & m, long & m1, long & m2, long & m3,
                               unsigned int ln, unsigned int pos)
{
   m2 = m3 = 0;
   string val;
   getToken(val, ln, pos++);
   istringstream in (val);
   in >> m;
   getToken(val, ln, pos++);
   if (val.empty())
      return;
   istringstream in2 (val);
   conv (m1, m);
   in2 >> m2;
   if (! in2)
      return;
   assert (m2 > 0);
   getToken(val, ln, pos++);
   if (val.empty())
      return;
   istringstream in3 (val);
   in3 >> m3;
   if (! in3)
      return;
   int sign;
   if (m < 0) {
      sign = -1;
      m = -m;
   } else
      sign = 1;
   m = sign*(power (m, m2) + m3);
}


//===========================================================================

void ParamReader::readInt(int& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = (int) strtol(val.c_str(), (char **)NULL, 10);
//   field = atoi(val.c_str());
}

//===========================================================================

void ParamReader::readLong(long& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = strtol(val.c_str(), (char **)NULL, 10);
//   field = atol(val.c_str());
}

//===========================================================================

void ParamReader::readZZ(NTL::ZZ & field, unsigned int ln, int pos)
{
   string val;
   getToken(val, ln, pos);
   field = to_ZZ(val.c_str());
}

//===========================================================================

void ParamReader::readDouble(double& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   field = strtod(val.c_str(), (char **)NULL);
//   field = atof(val.c_str());
}

//===========================================================================

void ParamReader::readMScal(MScal & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   conv(field, val.c_str());
}

//===========================================================================

void ParamReader::readBScal(BScal& field, unsigned int ln, int pos)
{
   string val;
   getToken(val, ln, pos);
   conv(field, val.c_str());
}

//===========================================================================

void ParamReader::readMVect(MVect & fields, unsigned int & ln, unsigned int pos,
   unsigned int numPos, int j)
{
   for (unsigned int i = pos; i < numPos; i++) {
      readMScal(fields[j], ln, i);
      j++;
   }
}

//===========================================================================

void ParamReader::readBMat(BMat & fields, unsigned int & ln, unsigned int pos,
   unsigned int numPos)
{
   for (unsigned int j = pos; j <= numPos; j++){
      for (unsigned int i = pos; i <= numPos; i++) {
         readBScal(fields[j][i], ln, i);
      }
      ++ln;
   }
}

//===========================================================================

void ParamReader::readIntVect (int* fields, unsigned int ln, unsigned int pos,
                               unsigned int num, int j)
{
   for (unsigned int i = pos; i < num; i++) {
      readInt(fields[j], ln, i);
      j++;
   }
}

//===========================================================================

void ParamReader::readDoubleVect(double* fields, unsigned int ln,
   unsigned int pos, unsigned int numPos, int j)
{
   for (unsigned int i = pos; i <= numPos; i++) {
      readDouble(fields[j], ln, i);
      j++;
   }
}


//===========================================================================

void ParamReader::readInterval (MVect & B, MVect & C, unsigned int & ln, int k)
{
   long m1, m2, m3;
   for (int i = 1; i <= k; i++) {
      readNumber3 (B[i], m1, m2, m3, ++ln, 1);
      readNumber3 (C[i], m1, m2, m3, ++ln, 1);
      assert (C[i] >= B[i]);
   }
}

//===========================================================================

void ParamReader::readNormType (NormType & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);
   cout << "valeur : " << val.c_str() << endl;
   if (0 == strcasecmp(val.c_str(), "SUPNORM"))
      field = SUPNORM;
   else if (0 == strcasecmp(val.c_str(), "L1NORM"))
      field = L1NORM;
   else if (0 == strcasecmp(val.c_str(), "L2NORM"))
      field = L2NORM;
   else if (0 == strcasecmp(val.c_str(), "ZAREMBANORM"))
      field = ZAREMBANORM;
    else
      MyExit(1, "readNormType:   NO SUCH CASE");
}


//===========================================================================


void ParamReader::readNormaType(NormaType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "BESTLAT"))
      field = BESTLAT;
   else if (0 == strcasecmp(val.c_str(), "LAMINATED"))
      field = LAMINATED;
   else if (0 == strcasecmp(val.c_str(), "ROGERS"))
      field = ROGERS;
   else if (0 == strcasecmp(val.c_str(), "MINKL1"))
      field = MINKL1;
   else if (0 == strcasecmp(val.c_str(), "MINKOWSKI"))
      field = MINKOWSKI;
   else if (0 == strcasecmp(val.c_str(), "NORMA_GENERIC"))
      field = NORMA_GENERIC;
   else
      MyExit(1, "readNormaType:   NO SUCH CASE");
}


//===========================================================================


void ParamReader::readOutputType(OutputType & field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "TERMINAL"))
      field = TERMINAL;
   else if (0 == strcasecmp(val.c_str(), "RES"))
      field = RES;
   else if (0 == strcasecmp(val.c_str(), "GEN")) {
      field = GEN;
      MyExit(1, "readOutputType:   GEN case not ready");
   } else if (0 == strcasecmp(val.c_str(), "TEX")) {
      field = TEX;
      MyExit(1, "readOutputType:   TEX case not ready");
   } else
      MyExit(1, "readOutputType:   NO SUCH CASE");
}


//===========================================================================

void ParamReader::readPreRed(PreReductionType& field, unsigned int ln, unsigned int pos)
{
   string val;
   getToken(val, ln, pos);

   if (0 == strcasecmp(val.c_str(), "BKZ"))
      field = BKZ;
   else if (0 == strcasecmp(val.c_str(), "PreRedDieter"))
      field = PreRedDieter;
   else if (0 == strcasecmp(val.c_str(), "LenstraLL"))
      field = LenstraLL;
   else
      MyExit(1, "readPreRed:   NO SUCH CASE");
}

//===========================================================================


void ParamReader::read (LatticeTesterConfig & config)
{
   getLines ();
   unsigned int ln = 1;

   readNormType (config.norm, ln, 1);
   readNormaType (config.normalizer, ++ln, 1);
   readPreRed (config.prereduction, ++ln, 1);

   if(config.prereduction == BKZ){
      readDouble (config.fact, ++ln, 1);
      readInt (config.blocksize, ++ln, 1);
   }
   else if(config.prereduction == PreRedDieter){
   }

   else if(config.prereduction == LenstraLL){
      readDouble (config.fact, ++ln, 1);
   }

   readInt (config.dimension, ++ln, 1);

   readBMat(config.basis, ++ln, 1, config.dimension);

   readLong (config.maxNodesBB, ++ln, 1);

   readOutputType(config.outputType, ++ln, 1);
}


//===========================================================================

} // End namespace LatMRG
