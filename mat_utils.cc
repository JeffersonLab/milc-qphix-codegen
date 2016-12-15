/*************** mat_utils.cc ***************/
/*
 * linear algebra of SU3 and Hermitian matrix 
 */
/********************************************/ 
#include "mat_utils.h"

extern vMatrix tmpMtr3;

static FVec tmpMtr[3][3][2] = { { {FVec("_mat_utils_tmpMtr_00R"), FVec("_mat_utils_tmpMtr_00I")}, {FVec("_mat_utils_tmpMtr_01R"), FVec("_mat_utils_tmpMtr_01I")}, {FVec("_mat_utils_tmpMtr_02R"), FVec("_mat_utils_tmpMtr_02I")} },
				{ {FVec("_mat_utils_tmpMtr_10R"), FVec("_mat_utils_tmpMtr_10I")}, {FVec("_mat_utils_tmpMtr_11R"), FVec("_mat_utils_tmpMtr_11I")}, {FVec("_mat_utils_tmpMtr_12R"), FVec("_mat_utils_tmpMtr_12I")} }, 
				{ {FVec("_mat_utils_tmpMtr_20R"), FVec("_mat_utils_tmpMtr_20I")}, {FVec("_mat_utils_tmpMtr_21R"), FVec("_mat_utils_tmpMtr_21I")}, {FVec("_mat_utils_tmpMtr_22R"), FVec("_mat_utils_tmpMtr_22I")} } };

static int Init = 0;

/* Initiate temporary variables */
static inline 
void init_mat_utils(InstVector& ivector)
{
  if(Init) return;
  declare_Gauge(ivector, tmpMtr);
  Init=1;
}

/* SU3 */

void decompressGauge(InstVector& ivector, vMatrix u_gauge, bool compress12, string mask)
{
    if( compress12 ) 
	for(int c = 0; c < 3; c++) {
            Conj_CrossProd(ivector, u_gauge[2][c], u_gauge[0][(c+1)%3], u_gauge[1][(c+2)%3], u_gauge[0][(c+2)%3], u_gauge[1][(c+1)%3], mask);
	}
}

void decompressGauge(InstVector& ivector, vMatrix u_gauge, FVec norm, bool compress12, string mask)
{
    if( compress12 ) {
        mulFVec(ivector, norm, u_gauge[0][0][0], u_gauge[0][0][0], mask);
        fmaddFVec(ivector, norm, u_gauge[0][0][1], u_gauge[0][0][1], norm, mask);
        for(int c=1; c<3; c++) for(int ir=0; ir<2; ++ir) {
            fmaddFVec(ivector, norm, u_gauge[0][c][ir], u_gauge[0][c][ir], norm, mask);
        }
        sqrtFVec(ivector, norm, norm, mask);
        for(int c = 0; c < 3; c++) {
            Conj_CrossProd(ivector, u_gauge[2][c], u_gauge[0][(c+1)%3], u_gauge[1][(c+2)%3], u_gauge[0][(c+2)%3], u_gauge[1][(c+1)%3], mask);
            divFVec(ivector, u_gauge[2][c][0], u_gauge[2][c][0], norm, mask);
            divFVec(ivector, u_gauge[2][c][1], u_gauge[2][c][1], norm, mask);
        }
    }
}

void movMat(InstVector& ivector, vMatrix out, vMatrix in, string mask)
{
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) for(int k=0; k<2; ++k)
    movFVec(ivector, out[i][j][k], in[i][j][k], mask);
}

/* Hermitian */

/* Calculate the hermition part of a 3x3 matrix */
/* scale factor 2.0 */
void matHerm(InstVector& ivector, FVec *herm, vMatrix matr, string mask)
{
  int i=0;
  for(int r=0; r<3; ++r) for(int c=r+1; c<3; ++c) 
  {
    addFVec(ivector, herm[2*i], matr[r][c][0], matr[c][r][0], mask);
    subFVec(ivector, herm[2*i+1], matr[r][c][1], matr[c][r][1], mask); 
    i++;
  }
  for(i=6; i<9; ++i)
    addFVec(ivector, herm[i], matr[i-6][i-6][0], matr[i-6][i-6][0], mask);
}

/* Calculate the anti-hermition part of a 3x3 matrix */
/* scale factor 2.0 */
void matAntiHerm(InstVector& ivector, FVec *aherm, vMatrix matr, string mask)
{
  int i=0;
  for(int r=0; r<3; ++r) for(int c=r+1; c<3; ++c)                       
  {
    subFVec(ivector, aherm[2*i], matr[r][c][0], matr[c][r][0], mask);
    addFVec(ivector, aherm[2*i+1], matr[r][c][1], matr[c][r][1], mask);
    i++;
  }
  for(i=6; i<9; ++i)
    addFVec(ivector, aherm[i], matr[i-6][i-6][1], matr[i-6][i-6][1], mask);
}

/* Calculate the traceless part of a 3x3 matrix */
void matTrls(InstVector& ivector, vMatrix trls, vMatrix matr, string mask)
{
  movFVec(ivector, trls[0][0][0], matr[0][0][0], mask);
  movFVec(ivector, trls[0][0][1], matr[0][0][1], mask);
  for(int i=1; i<3; ++i)
  {
    addFVec(ivector, trls[0][0][0], trls[0][0][0], matr[i][i][0], mask);
    addFVec(ivector, trls[0][0][1], trls[0][0][1], matr[i][i][1], mask);
  }
  smulFVec(ivector, trls[0][0][0], trls[0][0][0], "1.0/3.0", mask);
  smulFVec(ivector, trls[0][0][1], trls[0][0][1], "1.0/3.0", mask);
  for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
  if(r+c!=0) {
    movFVec(ivector, trls[r][c][ir], matr[r][c][ir], mask);
  }
  for(int i=1; i<3; ++i) for(int ir=0; ir<2; ++ir)
    subFVec(ivector, trls[i][i][ir], matr[i][i][ir], trls[0][0][ir], mask);
  subFVec(ivector, trls[0][0][0], matr[0][0][0], trls[0][0][0], mask);
  subFVec(ivector, trls[0][0][1], matr[0][0][1], trls[0][0][1], mask);
}

/* Calculate the traceless part of a 3x3 (anti-)hermition matrix */
void hermTrls(InstVector& ivector, FVec *trls, FVec *tmp, string mask)
{
  addFVec(ivector, tmp[0], trls[6], trls[7], mask);
  addFVec(ivector, tmp[0], tmp[0], trls[8], mask);
  smulFVec(ivector, tmp[0], tmp[0], "1.0/3.0", mask);
/*
  for(int i=0; i<3; ++i) for(int ir=0; ir<2; ++ir)
    movFVec(ivector, trls[2*i+ir], herm[2*i+ir], mask);
*/
  for(int i=0; i<3; ++i) 
    subFVec(ivector, trls[6+i], trls[6+i], tmp[0], mask);
}

/* Calculate the anti-hermition 
 * and traceless part of a 3x3 matrix */
void matAntiHermTrls(InstVector& ivector, FVec *trls, vMatrix matr, string mask)
{
  matAntiHerm(ivector, trls, matr, mask);
  FVec *tmp = &tmpMtr3[0][0][0];
  hermTrls(ivector, trls, tmp, mask);
}

/**************/
/* Algorithms */
/**************/
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

/* out = in1 + in2 */
void addMat(InstVector& ivector, vMatrix out, vMatrix in1, vMatrix in2, string mask)
{
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) for(int k=0; k<2; ++k)
    addFVec(ivector, out[i][j][k], in1[i][j][k], in2[i][j][k], mask);
}

/* out = s(Real scalar) * in */
void sMulMat(InstVector& ivector, vMatrix out, vMatrix in, string s, string mask)
{
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) for(int k=0; k<2; ++k)
    smulFVec(ivector, out[i][j][k], in[i][j][k], s, mask);
}

/* out = in1 + s(Real scalar) * in2 */
void addSMulMat(InstVector& ivector, vMatrix out, vMatrix in1, vMatrix in2, string s, string mask)
{
  for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) for(int k=0; k<2; ++k)
    fsmaddFVec(ivector, out[i][j][k], in2[i][j][k], s, in1[i][j][k], mask);
}

/* r[3][1] = u_mat[3][3](') * s_vec[3][1] */
void matMultVec(InstVector& ivector, FVec r[3][2], FVec u_mat[3][3][2], FVec s_vec[3][2], bool adjMul, string mask)
{
    for(int c1 = 0; c1 < 3; c1++) {
        if(!adjMul) {
            mulCVec(ivector, r[c1], u_mat[c1][0], s_vec[0], mask);
            fmaddCVec(ivector, r[c1], u_mat[c1][1], s_vec[1], r[c1], mask);
            fmaddCVec(ivector, r[c1], u_mat[c1][2], s_vec[2], r[c1], mask);
        }
        else {
            mulConjCVec(ivector, r[c1], u_mat[0][c1], s_vec[0], mask);
            fmaddConjCVec(ivector, r[c1], u_mat[1][c1], s_vec[1], r[c1], mask);
            fmaddConjCVec(ivector, r[c1], u_mat[2][c1], s_vec[2], r[c1], mask);
        }
    }
}

/* r[1][3] = s_vec[1][3] * u_mat[3][3](')
 * Transposed matrix vector multiplication
 */
void matMultVecT(InstVector& ivector, FVec r[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], bool adjMul, string mask)
{
    for(int c1 = 0; c1 < 3; c1++) {
        if(!adjMul) {
            mulCVec(ivector, r[c1], u_mat[0][c1], s_vec[0], mask);
            fmaddCVec(ivector, r[c1], u_mat[1][c1], s_vec[1], r[c1], mask);
            fmaddCVec(ivector, r[c1], u_mat[2][c1], s_vec[2], r[c1], mask);
        }
        else {
            mulConjCVec(ivector, r[c1], u_mat[c1][0], s_vec[0], mask);
            fmaddConjCVec(ivector, r[c1], u_mat[c1][1], s_vec[1], r[c1], mask);
            fmaddConjCVec(ivector, r[c1], u_mat[c1][2], s_vec[2], r[c1], mask);
        }
    }
}

/* r[1][3] = add[1][3] + s_vec[1][3] * u_mat[3][3](')
 *
 */
void addMatMultVecT(InstVector& ivector, FVec r[3][2], FVec add[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], bool adjMul, string mask)
{
  for(int c1 = 0; c1 < 3; c1++) {
        if(!adjMul) {
	  fmaddCVec(ivector, r[c1], u_mat[0][c1], s_vec[0], add[c1], mask);
	  fmaddCVec(ivector, r[c1], u_mat[1][c1], s_vec[1], r[c1], mask);
	  fmaddCVec(ivector, r[c1], u_mat[2][c1], s_vec[2], r[c1], mask);
	}
	else {
	  fmaddConjCVec(ivector, r[c1], u_mat[c1][0], s_vec[0], add[c1], mask);
	  fmaddConjCVec(ivector, r[c1], u_mat[c1][1], s_vec[1], r[c1], mask);
	  fmaddConjCVec(ivector, r[c1], u_mat[c1][2], s_vec[2], r[c1], mask);
	}
  }
}

/* r[1][3] = add[1][3] + s * s_vec[1][3] * u_mat[3][3](') */
void addSMulMatMultVecT(InstVector& ivector, FVec r[3][2], FVec add[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], string s, bool adjMul, string mask)
{
  for(int c1 = 0; c1 < 3; c1++) {
        if(!adjMul) {
	  fmsmaddCVec(ivector, r[c1], u_mat[0][c1], s_vec[0], s, add[c1], mask);
	  fmsmaddCVec(ivector, r[c1], u_mat[1][c1], s_vec[1], s, r[c1], mask);
	  fmsmaddCVec(ivector, r[c1], u_mat[2][c1], s_vec[2], s, r[c1], mask);
	}
	else {
	  fmsmaddConjCVec(ivector, r[c1], u_mat[c1][0], s_vec[0], s, add[c1], mask);
	  fmsmaddConjCVec(ivector, r[c1], u_mat[c1][1], s_vec[1], s, r[c1], mask);
	  fmsmaddConjCVec(ivector, r[c1], u_mat[c1][2], s_vec[2], s, r[c1], mask);
	}
  }
}

/* u_ret = u_mat1 * u_mat2(') */
void matMultMat(InstVector& ivector, FVec u_ret[3][3][2], FVec u_mat1[3][3][2], FVec u_mat2[3][3][2], bool adjMul, string mask)
{
    for(int i = 0; i < 3; i++) {
        matMultVecT(ivector, u_ret[i], u_mat1[i], u_mat2, adjMul, mask);
    }
}

/* u_ret = u_mat1(') * u_mat2 */
void adjMatMultMat(InstVector& ivector, FVec u_ret[3][3][2], FVec u_mat1[3][3][2], FVec u_mat2[3][3][2], bool adjMul, string mask)
{
  if(adjMul)
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            mulConjCVec(ivector, u_ret[i][j], u_mat1[0][i], u_mat2[0][j], mask);
            fmaddConjCVec(ivector, u_ret[i][j], u_mat1[1][i], u_mat2[1][j], u_ret[i][j], mask);
            fmaddConjCVec(ivector, u_ret[i][j], u_mat1[2][i], u_mat2[2][j], u_ret[i][j], mask);
        }
    }
  else
    matMultMat(ivector, u_ret, u_mat1, u_mat2, false, mask);
}

/* matr = su3 * matr1^* */
void matMultMatSconj(InstVector& ivector, vMatrix matr, vMatrix su3, vMatrix matr1, string mask)
{
  for(int r=0; r<3; ++r) for(int c=0; c<3; ++c)
    for(int k=0; k<3; ++k)
      if(k==0)
	mulConjCVec(ivector, matr[r][c], matr1[k][c], su3[r][k], mask);
      else
	fmaddConjCVec(ivector, matr[r][c], matr1[k][c], su3[r][k], matr[r][c], mask);
}

/* out = add + in1 * in2(') */
void addMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, bool conj, string mask)
{
    for(int i=0; i<3; ++i) 
      addMatMultVecT(ivector, out[i], add[i], in1[i], in2, conj, mask);
}

/* out = add + s * in1 * in2(') */
void addSMulMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, string s, bool conj, string mask)
{
    for(int i=0; i<3; ++i)
      addSMulMatMultVecT(ivector, out[i], add[i], in1[i], in2, s, conj, mask);
}

/* out = add + in1(') * in2 */
void addAdjMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, bool conj, string mask)
{
  if(conj)
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            fmaddConjCVec(ivector, out[i][j], in1[0][i], in2[0][j], add[i][j], mask);
	    fmaddConjCVec(ivector, out[i][j], in1[1][i], in2[1][j], out[i][j], mask);
	    fmaddConjCVec(ivector, out[i][j], in1[2][i], in2[2][j], out[i][j], mask);
        }
    }
  else
    addMatMultMat(ivector, out, add, in1, in2, false, mask);  
}

/* out = add + s * in1(') * in2 */
void addSMulAdjMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, string s, bool conj, string mask)
{
  if(conj)
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            fmsmaddConjCVec(ivector, out[i][j], in1[0][i], in2[0][j], s, add[i][j], mask);
	    fmsmaddConjCVec(ivector, out[i][j], in1[1][i], in2[1][j], s, out[i][j], mask);
	    fmsmaddConjCVec(ivector, out[i][j], in1[2][i], in2[2][j], s, out[i][j], mask);
	}
    }
  else
    addSMulMatMultMat(ivector, out, add, in1, in2, s, false, mask);
}

/* out = in1 + in2 */
void addHerm(InstVector& ivector, FVec *out, FVec *in1, FVec *in2, string mask)
{
  for(int i=0; i<9; ++i)
  {
    addFVec(ivector, out[i], in1[i], in2[i], mask);
  }
}

void addHermTl(InstVector& ivector, FVec *out, FVec *in1, FVec *in2, string mask)
{
  for(int i=0; i<8; ++i)
  {
    addFVec(ivector, out[i], in1[i], in2[i], mask);
  }
  naddFVec(ivector, out[8], out[6], out[7]);
}

/* out = in1 + Ic * scalar(Real) * in2 */
void hermAddSMultAntiHermIc(InstVector& ivector, FVec *out, FVec *in1, FVec scalar, FVec *in2, string mask)
{
  for(int i=0; i<3; ++i) {
    /* off-diagonal part */
    fnmaddFVec(ivector, out[2*i], scalar, in2[2*i+1], in1[2*i], mask);
    fmaddFVec(ivector, out[2*i+1], scalar, in2[2*i], in1[2*i+1], mask);
    /* diagonal part */
    fnmaddFVec(ivector, out[6+i], scalar, in2[6+i], in1[6+i], mask);
  }  
}
 
/* out = Ic * scalar(Real) * in1 */
void sMultAntiHermIc(InstVector& ivector, FVec *out, string scalar, FVec *in1, string mask)
{
  for(int i=0; i<3; ++i) {
    smulFVec(ivector, out[2*i+1], in1[2*i], scalar, mask);
    nsmulFVec(ivector, out[2*i], in1[2*i+1], scalar, mask);
    nsmulFVec(ivector, out[6+i], in1[6+i], scalar, mask);
  }
}

/* out = in1 + Ic * scalar(Real) * in2 */
void addSMultAntiHermIc(InstVector& ivector, FVec *out, FVec *in1, string scalar, FVec *in2, string mask)
{
  for(int i=0; i<3; ++i) {
    fsnmaddFVec(ivector, out[2*i], in2[2*i+1], scalar, in1[2*i], mask);
    fsmaddFVec(ivector, out[2*i+1], in2[2*i], scalar, in1[2*i+1], mask);
    fsnmaddFVec(ivector, out[6+i], in2[6+i], scalar, in1[6+i], mask);
  }
}

#if 0
/* out = in2 + Ic * scalar(Real) * herm1 * in2 */
void addSMultHermMultMatIc(InstVector& ivector, vMatrix out, FVec scalar, FVec *herm1, vMatrix in2, string mask)
{
  FVec tmp2 = tmpMtr3[0][0][0];
  FVec *tmpc = &tmpMtr3[0][1][0];
  for(int r=0; r<3; ++r) {
    mulFVec(ivector, tmp2, scalar, herm1[6+r], mask);
    for(int c=0; c<3; ++c)
    {
      fnmaddFVec(ivector, out[r][c][0], tmp2, in2[r][c][1], in2[r][c][0], mask); 
      fmaddFVec(ivector, out[r][c][1], tmp2, in2[r][c][0], in2[r][c][1], mask);
      for(int i=0; i<3; ++i)
      if(r+i!=2) {
	mulFVec(ivector, tmpc[0], scalar, herm1[2*i+1], mask);
	mulFVec(ivector, tmpc[1], scalar, herm1[2*i], mask);
	if(r<i+1-r)  
	  fnmaddConjCVec(ivector, out[r][c], tmpc, in2[i+1-r][c], out[r][c], mask);
	else
	  fnmaddCVec(ivector, out[r][c], tmpc, in2[i+1-r][c], out[r][c], mask);
      }
    }
  }
}
#endif

/* out = in2 + Ic * scalar(Real) * herm1 * in2 */
void addSMultHermMultMatIc(InstVector& ivector, vMatrix out, string scalar, FVec *herm1, vMatrix in2, /*FVec *tmp,*/ string mask)
{
  //FVec *tmpc = &tmp[2];
  for(int r=0; r<3; ++r) {
    for(int c=0; c<3; ++c)
    {
      fmsnmaddFVec(ivector, out[r][c][0], herm1[6+r], in2[r][c][1], scalar, in2[r][c][0], mask);
      fmsmaddFVec(ivector, out[r][c][1], herm1[6+r], in2[r][c][0], scalar, in2[r][c][1], mask);
      for(int i=0; i<3; ++i)
      if(r+i!=2) {
	//smulFVec(ivector, tmpc[0], herm1[2*i+1], scalar, mask);
	//smulFVec(ivector, tmpc[1], herm1[2*i], scalar, mask);
	if(r<i+1-r)
	  fmsnmaddConjCVec(ivector, out[r][c], herm1+2*i, in2[i+1-r][c], scalar, out[r][c], mask);
	else
	  fmsnmaddCVec(ivector, out[r][c], herm1+2*i, in2[i+1-r][c], scalar, out[r][c], mask);
      }
    }
  }
}

/* out = in2 + Ic * scalar(Real) * in2 * herm1 */
void addSMultMatMultHermIc(InstVector& ivector, vMatrix out, FVec scalar, vMatrix in2, FVec *herm1, string mask)
{
  FVec tmp2 = tmpMtr3[0][0][1];
  FVec *tmpc = &tmpMtr3[0][1][0];
  for(int r=0; r<3; ++r) {
    mulFVec(ivector, tmp2, scalar, herm1[6+r], mask);
    for(int c=0; c<3; ++c)
    {
      fnmaddFVec(ivector, out[c][r][0], in2[c][r][1], tmp2, in2[c][r][0], mask);
      fmaddFVec(ivector, out[c][r][1], in2[c][r][0], tmp2, in2[c][r][1], mask);
      for(int i=0; i<3; ++i)
      if(r+i!=2) {
	mulFVec(ivector, tmpc[0], scalar, herm1[2*i+1], mask);
	mulFVec(ivector, tmpc[1], scalar, herm1[2*i], mask);
	if(r<i+1-r)
	  fnmaddConjCVec(ivector, out[c][r], tmpc, in2[c][i+1-r], out[c][r], mask);
	else
	  fnmaddCVec(ivector, out[c][r], tmpc, in2[c][i+1-r], out[c][r], mask);	
      }
    }
  }
}


