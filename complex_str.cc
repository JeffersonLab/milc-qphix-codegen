/************ complex.cc ***************/
/*
 * Complex linear algrithms
 */
/***************************************/   
#include "complex.h"
using namespace std;

void movCVec(InstVector& ivector, string out, string in1, const string &mask)
{
    movFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), mask);
    movFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), mask);
}

void addCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    addFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), mask);
    addFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), mask);
}

void subCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    subFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), mask);
    subFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), mask);
}

void addiCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    subFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), mask);
    addFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), mask);
}

void subiCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    addFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), mask);
    subFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), mask);
}

/* out = in1 * in2 */
void mulCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    mulFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), mask);
    fnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
    mulFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), mask);
    fmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask); 
}

/* out = in1 * in2 + in3 */
void fmaddCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask)
{
    fmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), FVec(in3+"[RE]"), mask);
    fnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
    fmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), FVec(in3+"[IM]"), mask);
    fmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask);
}

/* out = s(Real scalar) * in1 * in2 + in3 */
void fmsmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[IM]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[IM]"), mask);
}

/* out = -in1 * in2 + in3 */
void fnmaddCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask)
{
    fnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), FVec(in3+"[RE]"), mask);
    fmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
    fnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), FVec(in3+"[IM]"), mask);
    fnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask);
}

/* out = -s(Real scalar) * in1 * in2 + in3 */
void fmsnmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[IM]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[IM]"), mask);
}

/* out = Ic * s * in1 * in2 + in3 */
void fmcsmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[IM]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[IM]"), mask);
}

/* out = -Ic * s * in1 * in2 + in3 */
void fmcsnmaddCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[IM]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[IM]"), mask);
}

/* out = in1' * in2 */
void mulConjCVec(InstVector& ivector, string out, string in1, string in2, const string &mask)
{
    mulFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), mask);
    fmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
    mulFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), mask);
    fnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask);
}

/* out = in1' * in2 + in3 */
void fmaddConjCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask)
{
    fmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), FVec(in3+"[RE]"), mask);
    fmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
    fmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), FVec(in3+"[IM]"), mask);
    fnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask);
}

/* out = -in1' * in2 + in3 */
void fnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, string in3, const string &mask)
{
  fnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), FVec(in3+"[RE]"), mask);
  fnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), FVec(out+"[RE]"), mask);
  fnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), FVec(in3+"[IM]"), mask);
  fmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), FVec(out+"[IM]"), mask);
}

/* out = s * in1' * in2 + in3 */
void fmsmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[IM]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[IM]"), mask);
}

/* out = -s * in1' * in2 + in3 */
void fmsnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[IM]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[IM]"), mask);
}

/* out = Ic * s * in1' * in2 + in3 */
void fmcsmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[RE]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[IM]"), mask);
    fmsmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[IM]"), mask);
}

/* out = -Ic * s * in1' * in2 + in3 */
void fmcsnmaddConjCVec(InstVector& ivector, string out, string in1, string in2, const string &s, string in3, const string &mask)
{
    fmsmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[RE]"), FVec(in2+"[IM]"), s, FVec(in3+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[RE]"), FVec(in1+"[IM]"), FVec(in2+"[RE]"), s, FVec(out+"[RE]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[RE]"), FVec(in2+"[RE]"), s, FVec(in3+"[IM]"), mask);
    fmsnmaddFVec(ivector, FVec(out+"[IM]"), FVec(in1+"[IM]"), FVec(in2+"[IM]"), s, FVec(out+"[IM]"), mask);    
}

void Conj_CrossProd(InstVector& ivector, string r, string s1, string s2, string s3, string s4, const string &mask)
{
    mulFVec(ivector, FVec(r+"[RE]"), FVec(s1+"[RE]"), FVec(s2+"[RE]"), mask);
    fnmaddFVec(ivector, FVec(r+"[RE]"), FVec(s1+"[IM]"), FVec(s2+"[IM]"), FVec(r+"[RE]"), mask);
    fnmaddFVec(ivector, FVec(r+"[RE]"), FVec(s3+"[RE]"), FVec(s4+"[RE]"), FVec(r+"[RE]"), mask);
    fmaddFVec(ivector, FVec(r+"[RE]"), FVec(s3+"[IM]"), FVec(s4+"[IM]"), FVec(r+"[RE]"), mask);

    mulFVec(ivector, FVec(r+"[IM]"), FVec(s3+"[RE]"), FVec(s4+"[IM]"), mask);
    fmaddFVec(ivector, FVec(r+"[IM]"), FVec(s3+"[IM]"), FVec(s4+"[RE]"), FVec(r+"[IM]"), mask);
    fnmaddFVec(ivector, FVec(r+"[IM]"), FVec(s1+"[RE]"), FVec(s2+"[IM]"), FVec(r+"[IM]"), mask);
    fnmaddFVec(ivector, FVec(r+"[IM]"), FVec(s1+"[IM]"), FVec(s2+"[RE]"), FVec(r+"[IM]"), mask);
}

