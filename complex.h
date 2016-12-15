#ifndef COMPLEX_H
#define COMPLEX_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <string>
#include <map>

#include <cstdio>
#include <cstdlib>

using namespace std;

#include "data_types_xyzt.h"
#include "instructions.h"
extern std::string ARCH_NAME;

void movCVec(InstVector& ivector, FVec *r, FVec *s1, const string &mask);
void addCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, const string &mask);
void subCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, const string &mask);
void addiCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, const string &mask);
void subiCVec(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, const string &mask);
void mulCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask="");
void fmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask="");
void fmsmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask);
void fnmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask="");
void fmsnmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask);
void mulConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask="");
void fmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask="");
void fnmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask="");
void fmsmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask);
void fmsnmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask);

void movCVec(InstVector& ivector, string r, string s1, const string &mask);
void addCVec(InstVector& ivector, string r, string s1, string s2, const string &mask);
void subCVec(InstVector& ivector, string r, string s1, string s2, const string &mask);
void addiCVec(InstVector& ivector, string r, string s1, string s2, const string &mask);
void subiCVec(InstVector& ivector, string r, string s1, string s2, const string &mask);
void mulCVec(InstVector& ivector, string out, string in1, string in2, const string &mask="");
void fmaddCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask="");
void fmsmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void fnmaddCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask="");
void fmsnmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void fmcsmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void fmcsnmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void mulConjCVec(InstVector& ivector, string out, string in1, string in2, const string &mask="");
void fmaddConjCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask="");
void fnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask="");
void fmsmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask="");
void fmsnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask="");
void fmcsmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void fmcsnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask);
void Conj_CrossProd(InstVector& ivector, string r, string s1, string s2, string s3, string s4, const string &mask="");

#endif
