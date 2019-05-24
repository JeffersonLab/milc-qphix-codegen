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

#include "data_utils.h"
#include "instructions.h"
//#include "data_types.h"


extern std::string ARCH_NAME;


void declare_HalfSpinor(InstVector& ivector, FVec h_spinor[2][3][2]);
void declare_Gauge(InstVector& ivector, FVec u_gauge[3][3][2]);
void declare_WilsonSpinor(InstVector& ivector, FVec spinor[4][3][2]);
void declare_KSSpinor(InstVector& ivector, FVec ks_spinor[3][2]);
void declare_Clover(InstVector& ivector, FVec diag[6], FVec offdiag[15][2]);
// r[RE] = s1[RE]-beta_vec*s2[RE] = fnmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]-beta_vec*s2[IM] = fnamdd(beta_vec,s2[IM],s1[IM])
void addCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] + beta_vec*s2[RE] = fmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[IM] = fmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] + beta_vec * s2[IM] = fmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] - beta_vec * s2[RE] = fnmadd(beta_vec, s2[RE], s1[IM])

void addiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] - beta_vec*s2[IM] = fnmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[RE] = fmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE]+beta_vec*s2[RE] = fmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]+beta_vec*s2[IM] = fmadd(beta_vec,s2[IM],s1[IM])
void addCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] - beta_vec*s2[RE] = fnmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[IM] = fnmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] - beta_vec * s2[IM] = fnmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] + beta_vec * s2[RE] = fmadd(beta_vec, s2[RE], s1[IM])
void addiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r[RE] = s1[RE] + beta_vec*s2[IM] = fmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[RE] = fnmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask);
// r = (s1*s2-s3*s4)'
//r[RE] = (s1[RE]*s2[RE])-(s1[IM]*s2[IM])-(s3[RE]*s4[RE])+(s3[IM]*s4[IM])
//r[IM] = (s3[RE]*s4[IM])+(s3[IM]*s4[RE])-(s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void Conj_CrossProd(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec *s3, FVec *s4, string &mask);
// Merge L2 prefetches with another instruction stream
void mergeIvectorWithL2Prefetches(InstVector& ivector, InstVector& l2prefs);
// Dump an instruction stream into a file
void dumpIVector(InstVector& ivector, string filename);
