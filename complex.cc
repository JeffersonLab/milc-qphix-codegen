/************ complex.cc ***************/
/*
 * Complex linear algrithms
 */
/***************************************/   
#include "complex.h"
using namespace std;

void movCVec(InstVector& ivector, vComplex out, vComplex in1, const string &mask)
{
    movFVec(ivector, out[RE], in1[RE], mask);
    movFVec(ivector, out[IM], in1[IM], mask);
}

void addCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    addFVec(ivector, out[RE], in1[RE], in2[RE], mask);
    addFVec(ivector, out[IM], in1[IM], in2[IM], mask);
}

void subCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    subFVec(ivector, out[RE], in1[RE], in2[RE], mask);
    subFVec(ivector, out[IM], in1[IM], in2[IM], mask);
}

void addiCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    subFVec(ivector, out[RE], in1[RE], in2[IM], mask);
    addFVec(ivector, out[IM], in1[IM], in2[RE], mask);
}

void subiCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    addFVec(ivector, out[RE], in1[RE], in2[IM], mask);
    subFVec(ivector, out[IM], in1[IM], in2[RE], mask);
}

/* out = in1 * in2 */
void mulCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    mulFVec(ivector, out[RE], in1[RE], in2[RE], mask);
    fnmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
    mulFVec(ivector, out[IM], in1[RE], in2[IM], mask);
    fmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask); 
}

/* out = in1 * in2 + in3 */
void fmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask)
{
    fmaddFVec(ivector, out[RE], in1[RE], in2[RE], in3[RE], mask);
    fnmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
    fmaddFVec(ivector, out[IM], in1[RE], in2[IM], in3[IM], mask);
    fmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask);
}

/* out = s(Real scalar) * in1 * in2 + in3 */
void fmsmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask)
{
    fmsmaddFVec(ivector, out[RE], in1[RE], in2[RE], s, in3[RE], mask);
    fmsnmaddFVec(ivector, out[RE], in1[IM], in2[IM], s, out[RE], mask);
    fmsmaddFVec(ivector, out[IM], in1[RE], in2[IM], s, in3[IM], mask);
    fmsmaddFVec(ivector, out[IM], in1[IM], in2[RE], s, out[IM], mask);
}

/* out = -in1 * in2 + in3 */
void fnmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask)
{
    fnmaddFVec(ivector, out[RE], in1[RE], in2[RE], in3[RE], mask);
    fmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
    fnmaddFVec(ivector, out[IM], in1[RE], in2[IM], in3[IM], mask);
    fnmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask);
}

/* out = -s(Real scalar) * in1 * in2 + in3 */
void fmsnmaddCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask)
{
    fmsnmaddFVec(ivector, out[RE], in1[RE], in2[RE], s, in3[RE], mask);
    fmsmaddFVec(ivector, out[RE], in1[IM], in2[IM], s, out[RE], mask);
    fmsnmaddFVec(ivector, out[IM], in1[RE], in2[IM], s, in3[IM], mask);
    fmsnmaddFVec(ivector, out[IM], in1[IM], in2[RE], s, out[IM], mask);
}

/* out = in1' * in2 */
void mulConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &mask)
{
    mulFVec(ivector, out[RE], in1[RE], in2[RE], mask);
    fmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
    mulFVec(ivector, out[IM], in1[RE], in2[IM], mask);
    fnmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask);
}

/* out = in1' * in2 + in3 */
void fmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask)
{
    fmaddFVec(ivector, out[RE], in1[RE], in2[RE], in3[RE], mask);
    fmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
    fmaddFVec(ivector, out[IM], in1[RE], in2[IM], in3[IM], mask);
    fnmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask);
}

/* out = -in1' * in2 + in3 */
void fnmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, vComplex in3, const string &mask)
{
  fnmaddFVec(ivector, out[RE], in1[RE], in2[RE], in3[RE], mask);
  fnmaddFVec(ivector, out[RE], in1[IM], in2[IM], out[RE], mask);
  fnmaddFVec(ivector, out[IM], in1[RE], in2[IM], in3[IM], mask);
  fmaddFVec(ivector, out[IM], in1[IM], in2[RE], out[IM], mask);
}

/* out = s * in1' * in2 + in3 */
void fmsmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask)
{
    fmsmaddFVec(ivector, out[RE], in1[RE], in2[RE], s, in3[RE], mask);
    fmsmaddFVec(ivector, out[RE], in1[IM], in2[IM], s, out[RE], mask);
    fmsmaddFVec(ivector, out[IM], in1[RE], in2[IM], s, in3[IM], mask);
    fmsnmaddFVec(ivector, out[IM], in1[IM], in2[RE], s, out[IM], mask);
}

/* out = -s * in1' * in2 + in3 */
void fmsnmaddConjCVec(InstVector& ivector, vComplex out, vComplex in1, vComplex in2, const string &s, vComplex in3, const string &mask)
{
    fmsnmaddFVec(ivector, out[RE], in1[RE], in2[RE], s, in3[RE], mask);
    fmsnmaddFVec(ivector, out[RE], in1[IM], in2[IM], s, out[RE], mask);
    fmsnmaddFVec(ivector, out[IM], in1[RE], in2[IM], s, in3[IM], mask);
    fmsmaddFVec(ivector, out[IM], in1[IM], in2[RE], s, out[IM], mask);
}

