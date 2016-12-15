/*************** mat_utils.cc ***************/
/*
 * linear algebra of SU3 and Hermitian matrix 
 */
/********************************************/ 
#include "mat_utils.h"

extern string tmpMtr3str;

/* Reserved indices: i, j, k, c, r, c1, ir */

#define Q(M) #M
#define QUOTE(M) Q(M)

#define MATRLOOP(i, j, k) \
  forLoopInc(ivector, QUOTE(i), QUOTE(0), QUOTE(3));\
  forLoopInc(ivector, QUOTE(j), QUOTE(0), QUOTE(3));\
  forLoopInc(ivector, QUOTE(k), QUOTE(0), QUOTE(2));\

#define ENDMATRLOOP \
  endScope(ivector);\
  endScope(ivector);\
  endScope(ivector);\


/* SU3 */

void decompressGauge(InstVector& ivector, string u_gauge, bool compress12, string mask)
{
    if( compress12 ) 
    {
	forLoopInc(ivector, "c", "0", "3");
	{
            Conj_CrossProd(ivector, u_gauge+"[2][c]", u_gauge+"[0][(c+1)%3]", u_gauge+"[1][(c+2)%3]", u_gauge+"[0][(c+2)%3]", u_gauge+"[1][(c+1)%3]", mask);
	}
	endScope(ivector);
    }
}

void decompressGauge(InstVector& ivector, string u_gauge, FVec norm, bool compress12, string mask)
{
    if( compress12 ) {
        mulFVec(ivector, norm, FVec(u_gauge+"[0][0][0]"), FVec(u_gauge+"[0][0][0]"), mask);
        fmaddFVec(ivector, norm, FVec(u_gauge+"[0][0][1]"), FVec(u_gauge+"[0][0][1]"), norm, mask);
        forLoopInc(ivector, "c", "1", "3");
 	forLoopInc(ivector, "ir", "0", "2");
	{
            fmaddFVec(ivector, norm, FVec(u_gauge+"[0][c][ir]"), FVec(u_gauge+"[0][c][ir]"), norm, mask);
        }
	endScope(ivector);
	endScope(ivector);
        sqrtFVec(ivector, norm, norm, mask);
	forLoopInc(ivector, "c", "0", "3"); 
	{
            Conj_CrossProd(ivector, u_gauge+"[2][c]", u_gauge+"[0][(c+1)%3]", u_gauge+"[1][(c+2)%3]", u_gauge+"[0][(c+2)%3]", u_gauge+"[1][(c+1)%3]", mask);
            divFVec(ivector, FVec(u_gauge+"[2][c][0]"), FVec(u_gauge+"[2][c][0]"), norm, mask);
            divFVec(ivector, FVec(u_gauge+"[2][c][1]"), FVec(u_gauge+"[2][c][1]"), norm, mask);
        }
	endScope(ivector);
    }
}

void movMat(InstVector& ivector, string out, string in, string mask)
{
  MATRLOOP(i, j, k)
  {
    movFVec(ivector, FVec(out+"[i][j][k]"), FVec(in+"[i][j][k]"), mask);
  }
  ENDMATRLOOP
}

/* Hermitian */

/* Calculate the hermition part of a 3x3 matrix */
/* scale factor 2.0 */
void matHerm(InstVector& ivector, string herm, string matr, string mask)
{
  forLoopInc(ivector, "r", "0", "3");
  forLoopStatement(ivector, "int i=0, c=r+1", "c<3", "++c, ++i");
  {
    addFVec(ivector, FVec(herm+"[2*i]"), FVec(matr+"[r][c][0]"), FVec(matr+"[c][r][0]"), mask);
    subFVec(ivector, FVec(herm+"[2*i+1]"), FVec(matr+"[r][c][1]"), FVec(matr+"[c][r][1]"), mask); 
  }
  endScope(ivector);
  endScope(ivector);
  forLoopInc(ivector, "r", "6", "9");
  {
    addFVec(ivector, FVec(herm+"[r]"), FVec(matr+"[r-6][r-6][0]"), FVec(matr+"[r-6][r-6][0]"), mask);
  }
  endScope(ivector);
}

/* Calculate the anti-hermition part of a 3x3 matrix */
/* scale factor 2.0 */
void matAntiHerm(InstVector& ivector, string aherm, string matr, string mask)
{
  forLoopStatement(ivector, "int i=0, r=0", "r<3", "++r");
  forLoopStatement(ivector, "int c=r+1", "c<3", "++c, ++i");
  {
    subFVec(ivector, FVec(aherm+"[2*i]"), FVec(matr+"[r][c][0]"), FVec(matr+"[c][r][0]"), mask);
    addFVec(ivector, FVec(aherm+"[2*i+1]"), FVec(matr+"[r][c][1]"), FVec(matr+"[c][r][1]"), mask);
  }
  endScope(ivector);
  endScope(ivector);
  forLoopInc(ivector, "r", "6", "9");
  {
    addFVec(ivector, FVec(aherm+"[r]"), FVec(matr+"[r-6][r-6][1]"), FVec(matr+"[r-6][r-6][1]"), mask);
  }
  endScope(ivector);
}

/* Calculate the traceless part of a 3x3 matrix */
void matTrls(InstVector& ivector, string trls, string matr, string mask)
{
  movFVec(ivector, FVec(trls+"[0][0][0]"), FVec(matr+"[0][0][0]"), mask);
  movFVec(ivector, FVec(trls+"[0][0][1]"), FVec(matr+"[0][0][1]"), mask);
  forLoopInc(ivector, "i", "1", "3");
  {
    addFVec(ivector, FVec(trls+"[0][0][0]"), FVec(trls+"[0][0][0]"), FVec(matr+"[i][i][0]"), mask);
    addFVec(ivector, FVec(trls+"[0][0][1]"), FVec(trls+"[0][0][1]"), FVec(matr+"[i][i][1]"), mask);
  }
  endScope(ivector);
  smulFVec(ivector, FVec(trls+"[0][0][0]"), FVec(trls+"[0][0][0]"), "1.0/3.0", mask);
  smulFVec(ivector, FVec(trls+"[0][0][1]"), FVec(trls+"[0][0][1]"), "1.0/3.0", mask);
  MATRLOOP(r, c, k)
  {
    ifStatement(ivector, "r+c!=0");
    {
      movFVec(ivector, FVec(trls+"[r][c][k]"), FVec(matr+"[r][c][k]"), mask);
    }
    endScope(ivector);
  }
  ENDMATRLOOP
  forLoopInc(ivector, "i", "1", "3");
  {
    subFVec(ivector, FVec(trls+"[i][i][0]"), FVec(matr+"[i][i][0]"), FVec(trls+"[0][0][0]"), mask);
    subFVec(ivector, FVec(trls+"[i][i][1]"), FVec(matr+"[i][i][1]"), FVec(trls+"[0][0][1]"), mask);
  }
  endScope(ivector);
  subFVec(ivector, FVec(trls+"[0][0][0]"), FVec(matr+"[0][0][0]"), FVec(trls+"[0][0][0]"), mask);
  subFVec(ivector, FVec(trls+"[0][0][1]"), FVec(matr+"[0][0][1]"), FVec(trls+"[0][0][1]"), mask);
}

/* Calculate the traceless part of a 3x3 (anti-)hermition matrix */
void hermTrls(InstVector& ivector, string trls, string tmp, string mask)
{
  addFVec(ivector, FVec(tmp+"[0]"), FVec(trls+"[6]"), FVec(trls+"[7]"), mask);
  addFVec(ivector, FVec(tmp+"[0]"), FVec(tmp+"[0]"), FVec(trls+"[8]"), mask);
  smulFVec(ivector, FVec(tmp+"[0]"), FVec(tmp+"[0]"), "1.0/3.0", mask);
  forLoopInc(ivector, "i", "0", "3");
  {
    subFVec(ivector, FVec(trls+"[6+i]"), FVec(trls+"[6+i]"), FVec(tmp+"[0]"), mask);
  }
  endScope(ivector);
}

/* Calculate the anti-hermition 
 * and traceless part of a 3x3 matrix */
void matAntiHermTrls(InstVector& ivector, string trls, string matr, string mask)
{
  matAntiHerm(ivector, trls, matr, mask);
  string tmp = "(&"+tmpMtr3str+"[0][0][0])";
  hermTrls(ivector, trls, tmp, mask);
}

/**************/
/* Algorithms */
/**************/
/* out = in1 + in2 */
void addMat(InstVector& ivector, string out, string in1, string in2, string mask)
{
/*
  forLoopInc(ivector, "i", "0", "3");
  forLoopInc(ivector, "j", "0", "3");
  forLoopInc(ivector, "k", "0", "2");
*/
  MATRLOOP(i, j, k)
  {
    addFVec(ivector, FVec(out+"[i][j][k]"), FVec(in1+"[i][j][k]"), FVec(in2+"[i][j][k]"), mask);
  }
  ENDMATRLOOP
/*
  endScope(ivector);
  endScope(ivector);
  endScope(ivector);
*/
}

/* out = s(Real scalar) * in */
void sMulMat(InstVector& ivector, string out, string in, string s, string mask)
{
  MATRLOOP(i, j, k)
  {
    smulFVec(ivector, FVec(out+"[i][j][k]"), FVec(in+"[i][j][k]"), s, mask);
  }
  ENDMATRLOOP
}

/* out = in1 + s(Real scalar) * in2 */
void addSMulMat(InstVector& ivector, string out, string in1, string in2, string s, string mask)
{
  MATRLOOP(i, j, k)
  {
    fsmaddFVec(ivector, FVec(out+"[i][j][k]"), FVec(in2+"[i][j][k]"), s, FVec(in1+"[i][j][k]"), mask);
  }
  ENDMATRLOOP
}

/* r[3][1] = u_mat[3][3](') * s_vec[3][1] */
void matMultVec(InstVector& ivector, string r, string u_mat, string s_vec, bool adjMul, string mask)
{
  forLoopInc(ivector, "c1", "0", "3");
  {
    if(!adjMul) {
	mulCVec(ivector, r+"[c1]", u_mat+"[c1][0]", s_vec+"[0]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[c1][1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[c1][2]", s_vec+"[2]", r+"[c1]", mask);
    }
    else {
    	mulConjCVec(ivector, r+"[c1]", u_mat+"[0][c1]", s_vec+"[0]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[1][c1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[2][c1]", s_vec+"[2]", r+"[c1]", mask);
    }
  }
  endScope(ivector);
}

/* r[1][3] = s_vec[1][3] * u_mat[3][3](')
 * Transposed matrix vector multiplication
 */
void matMultVecT(InstVector& ivector, string r, string s_vec, string u_mat, bool adjMul, string mask)
{
  forLoopInc(ivector, "c1", "0", "3");
  {
    if(!adjMul) {
	mulCVec(ivector, r+"[c1]", u_mat+"[0][c1]", s_vec+"[0]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[1][c1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[2][c1]", s_vec+"[2]", r+"[c1]", mask);
    }
    else {
	mulConjCVec(ivector, r+"[c1]", u_mat+"[c1][0]", s_vec+"[0]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][2]", s_vec+"[2]", r+"[c1]", mask);
    }
  }
  endScope(ivector);
}

/* r[1][3] = add[1][3] + s_vec[1][3] * u_mat[3][3](')
 *
 */
void addMatMultVecT(InstVector& ivector, string r, string add, string s_vec, string u_mat, bool adjMul, string mask)
{
  forLoopInc(ivector, "c1", "0", "3");
  {
    if(!adjMul) {
	fmaddCVec(ivector, r+"[c1]", u_mat+"[0][c1]", s_vec+"[0]", add+"[c1]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[1][c1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddCVec(ivector, r+"[c1]", u_mat+"[2][c1]", s_vec+"[2]", r+"[c1]", mask);
    }
    else {
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][0]", s_vec+"[0]", add+"[c1]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][1]", s_vec+"[1]", r+"[c1]", mask);
	fmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][2]", s_vec+"[2]", r+"[c1]", mask);
    }
  }
  endScope(ivector);
}

/* r[1][3] = add[1][3] + s * s_vec[1][3] * u_mat[3][3](') */
void addSMulMatMultVecT(InstVector& ivector, string r, string add, string s_vec, string u_mat, string s, bool adjMul, string mask)
{
  forLoopInc(ivector, "c1", "0", "3");
  {
    if(!adjMul) {
 	fmsmaddCVec(ivector, r+"[c1]", u_mat+"[0][c1]", s_vec+"[0]", s, add+"[c1]", mask);
	fmsmaddCVec(ivector, r+"[c1]", u_mat+"[1][c1]", s_vec+"[1]", s, r+"[c1]", mask);
	fmsmaddCVec(ivector, r+"[c1]", u_mat+"[2][c1]", s_vec+"[2]", s, r+"[c1]", mask);
    }
    else {
	fmsmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][0]", s_vec+"[0]", s, add+"[c1]", mask);
	fmsmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][1]", s_vec+"[1]", s, r+"[c1]", mask);
	fmsmaddConjCVec(ivector, r+"[c1]", u_mat+"[c1][2]", s_vec+"[2]", s, r+"[c1]", mask);
    }
  }
  endScope(ivector);
}

/* u_ret = u_mat1 * u_mat2(') */
void matMultMat(InstVector& ivector, string u_ret, string u_mat1, string u_mat2, bool adjMul, string mask)
{
    forLoopInc(ivector, "i", "0", "3");
    {
        matMultVecT(ivector, u_ret+"[i]", u_mat1+"[i]", u_mat2, adjMul, mask);
    }
    endScope(ivector);
}

/* u_ret = u_mat1(') * u_mat2 */
void adjMatMultMat(InstVector& ivector, string u_ret, string u_mat1, string u_mat2, bool adjMul, string mask)
{
  if(adjMul) {
    forLoopInc(ivector, "i", "0", "3");
    {
        forLoopInc(ivector, "j", "0", "3");
	{
            mulConjCVec(ivector, u_ret+"[i][j]", u_mat1+"[0][i]", u_mat2+"[0][j]", mask);
            fmaddConjCVec(ivector, u_ret+"[i][j]", u_mat1+"[1][i]", u_mat2+"[1][j]", u_ret+"[i][j]", mask);
            fmaddConjCVec(ivector, u_ret+"[i][j]", u_mat1+"[2][i]", u_mat2+"[2][j]", u_ret+"[i][j]", mask);
        }
	endScope(ivector);
    }
    endScope(ivector);
  }
  else
    matMultMat(ivector, u_ret, u_mat1, u_mat2, false, mask);
}

/* matr = su3 * matr1^* */
void matMultMatSconj(InstVector& ivector, string matr, string su3, string matr1, string mask)
{
  forLoopInc(ivector, "r", "0", "3");
  forLoopInc(ivector, "c", "0", "3");
  {
    mulConjCVec(ivector, matr+"[r][c]", matr1+"[0][c]", su3+"[r][0]", mask);
    fmaddConjCVec(ivector, matr+"[r][c]", matr1+"[1][c]", su3+"[r][1]", matr+"[r][c]", mask);
    fmaddConjCVec(ivector, matr+"[r][c]", matr1+"[2][c]", su3+"[r][2]", matr+"[r][c]", mask);
  }
  endScope(ivector);
  endScope(ivector); 
}

/* out = add + in1 * in2(') */
void addMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, bool conj, string mask)
{
    forLoopInc(ivector, "i", "0", "3");
    {
      addMatMultVecT(ivector, out+"[i]", add+"[i]", in1+"[i]", in2, conj, mask);
    }
    endScope(ivector);
}

/* out = add + s * in1 * in2(') */
void addSMulMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, string s, bool conj, string mask)
{
    forLoopInc(ivector, "i", "0", "3");
    {
      addSMulMatMultVecT(ivector, out+"[i]", add+"[i]", in1+"[i]", in2, s, conj, mask);
    }
    endScope(ivector);
}

/* out = add + in1(') * in2 */
void addAdjMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, bool conj, string mask)
{
  if(conj) {
    forLoopInc(ivector, "i", "0", "3");
    {
	forLoopInc(ivector, "j", "0", "3");
	{
            fmaddConjCVec(ivector, out+"[i][j]", in1+"[0][i]", in2+"[0][j]", add+"[i][j]", mask);
	    fmaddConjCVec(ivector, out+"[i][j]", in1+"[1][i]", in2+"[1][j]", out+"[i][j]", mask);
	    fmaddConjCVec(ivector, out+"[i][j]", in1+"[2][i]", in2+"[2][j]", out+"[i][j]", mask);
        }
	endScope(ivector);
    }
    endScope(ivector);
  }
  else
    addMatMultMat(ivector, out, add, in1, in2, false, mask);  
}

/* out = add + s * in1(') * in2 */
void addSMulAdjMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, string s, bool conj, string mask)
{
  if(conj) {
    forLoopInc(ivector, "i", "0", "3");
    {
        forLoopInc(ivector, "j", "0", "3");
        {
            fmsmaddConjCVec(ivector, out+"[i][j]", in1+"[0][i]", in2+"[0][j]", s, add+"[i][j]", mask);
            fmsmaddConjCVec(ivector, out+"[i][j]", in1+"[1][i]", in2+"[1][j]", s, out+"[i][j]", mask);
            fmsmaddConjCVec(ivector, out+"[i][j]", in1+"[2][i]", in2+"[2][j]", s, out+"[i][j]", mask);
	}
	endScope(ivector);
    }
    endScope(ivector);
  }
  else
    addSMulMatMultMat(ivector, out, add, in1, in2, s, false, mask);
}

/* out = in1 + in2 */
void addHerm(InstVector& ivector, string out, string in1, string in2, string mask)
{
  forLoopInc(ivector, "i", "0", "9");
  {
    addFVec(ivector, FVec(out+"[i]"), FVec(in1+"[i]"), FVec(in2+"[i]"), mask);
  }
  endScope(ivector);
}

void addHermTl(InstVector& ivector, string out, string in1, string in2, string mask, bool add)
{
  forLoopInc(ivector, "i", "0", "8");
  {
    addFVec(ivector, FVec(out+"[i]"), FVec(in1+"[i]"), FVec(in2+"[i]"), mask);
  }
  endScope(ivector);
  if(add) naddFVec(ivector, FVec(out+"[8]"), FVec(out+"[6]"), FVec(out+"[7]"));
}

/* out = in1 + Ic * scalar(Real) * in2 */
void hermAddSMultAntiHermIc(InstVector& ivector, string out, string in1, string scalar, string in2, string mask)
{
  forLoopInc(ivector, "i", "0", "3");
  {
    /* off-diagonal part */
    fnmaddFVec(ivector, FVec(out+"[2*i]"), scalar, FVec(in2+"[2*i+1]"), FVec(in1+"[2*i]"), mask);
    fmaddFVec(ivector, FVec(out+"[2*i+1]"), scalar, FVec(in2+"[2*i]"), FVec(in1+"[2*i+1]"), mask);
    /* diagonal part */
    fnmaddFVec(ivector, FVec(out+"[6+i]"), scalar, FVec(in2+"[6+i]"), FVec(in1+"[6+i]"), mask);
  }  
  endScope(ivector);
}
 
/* out = Ic * scalar(Real) * in1 */
void sMultAntiHermIc(InstVector& ivector, string out, string scalar, string in1, string mask)
{
  forLoopInc(ivector, "i", "0", "3");
  {
    smulFVec(ivector, FVec(out+"[2*i+1]"), FVec(in1+"[2*i]"), scalar, mask);
    nsmulFVec(ivector, FVec(out+"[2*i]"), FVec(in1+"[2*i+1]"), scalar, mask);
    nsmulFVec(ivector, FVec(out+"[6+i]"), FVec(in1+"[6+i]"), scalar, mask);
  }
  endScope(ivector);
}

/* out = in1 + Ic * scalar(Real) * in2 */
void addSMultAntiHermIc(InstVector& ivector, string out, string in1, string scalar, string in2, string mask)
{
  forLoopInc(ivector, "i", "0", "3");
  {
    fsnmaddFVec(ivector, FVec(out+"[2*i]"), FVec(in2+"[2*i+1]"), scalar, FVec(in1+"[2*i]"), mask);
    fsmaddFVec(ivector, FVec(out+"[2*i+1]"), FVec(in2+"[2*i]"), scalar, FVec(in1+"[2*i+1]"), mask);
    fsnmaddFVec(ivector, FVec(out+"[6+i]"), FVec(in2+"[6+i]"), scalar, FVec(in1+"[6+i]"), mask);
  }
  endScope(ivector);
}

/* out = in2 + Ic * scalar(Real) * herm1 * in2 */
void addSMultHermMultMatIc(InstVector& ivector, string out, string scalar, string herm1, string in2, string mask)
{
  forLoopInc(ivector, "r", "0", "3");
  {
    forLoopInc(ivector, "c", "0", "3");
    {
      fmsnmaddFVec(ivector, FVec(out+"[r][c][0]"), FVec(herm1+"[6+r]"), FVec(in2+"[r][c][1]"), scalar, FVec(in2+"[r][c][0]"), mask);
      fmsmaddFVec(ivector, FVec(out+"[r][c][1]"), FVec(herm1+"[6+r]"), FVec(in2+"[r][c][0]"), scalar, FVec(in2+"[r][c][1]"), mask);
      forLoopInc(ivector, "i", "0", "3");
      ifStatement(ivector, "r+i!=2");
      {
	ifStatement(ivector, "r<i+1-r");
	{
	  fmcsmaddCVec(ivector, out+"[r][c]", "(&"+herm1+"[2*i])", in2+"[i+1-r][c]", scalar, out+"[r][c]", mask);
	  //fmsnmaddConjCVec(ivector, out+"[r][c]", "(&"+herm1+"[2*i])", in2+"[i+1-r][c]", scalar, out+"[r][c]", mask);
	}
	elseStatement(ivector);
	{
	  fmcsmaddConjCVec(ivector, out+"[r][c]", "(&"+herm1+"[2*i])", in2+"[i+1-r][c]", scalar, out+"[r][c]", mask);
	  //fmsnmaddCVec(ivector, out+"[r][c]", "(&"+herm1+"[2*i]", in2+"[i+1-r][c]", scalar, out+"[r][c]", mask);
	}
	endScope(ivector);
      }
      endScope(ivector);
      endScope(ivector);
    }
    endScope(ivector);
  }
  endScope(ivector);
}


