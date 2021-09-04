/*************** staples.cc ****************************/
/*

Calculate three-link staples on each site

*/
/********************************************************/

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <string>

using namespace std;

#include "data_types_xyzt.h"
//#define _STAPLES_C
#include "staples.h"

extern vMatrix tmpMtr3;
extern FVec mOne;
extern string mask;

/*
 * Staple = U1(dag) * U2 * U3(!dag)
 */
//static 
void staple( InstVector& ivector, vMatrix Staple, vMatrix U1, vMatrix U2, vMatrix U3, int updown, string mask )
{
  bool dag = (updown==DN);
  adjMatMultMat(ivector, tmpMtr3, U1, U2, dag, mask);
  matMultMat(ivector, Staple, tmpMtr3, U3, !dag, mask);
}

/*
 *  ____      __        __
 * |    | =  |  |   *  |  | .
 *
 * Rect0 += Staple1 * Staple2 (add==1)
 * Rect0  = Staple1 * Staple2 (add==0)
 *
 */ 
void rectangle0( InstVector& ivector, vMatrix Rect0, vMatrix Staple1, vMatrix Staple2, int add, string mask )
{
  bentchair0(ivector, Rect0, Staple1, Staple2, false, false, add, mask);
}

/*
 *  ____       ____	    
 * |  __|  =  |    |  * __  ,
 *  ____              ____
 * |__  |  =   __ *  |    | .
 *
 * Rect1 (+)= s * Rect0 * U1dag ,
 * Rect2 (+)= s * U2dag * Rect0 .
 *
 * updown = 1 (down), 2 (up), 3 (both).
 *
 */
//static 
void rectangle1( InstVector& ivector, vMatrix Rect1, vMatrix Rect2, vMatrix Rect0, vMatrix U1, vMatrix U2, /* bool order, */ /*int updown,*/ string kappa, int add, string mask )
{
  if(add==1)
  {
    addSMulMatMultMat(ivector, Rect1, Rect1, Rect0, U1, kappa, true, mask);
    addSMulAdjMatMultMat(ivector, Rect2, Rect2, U2, Rect0, kappa, true, mask);
  }
  else
  {
    matMultMat(ivector, Rect1, Rect0, U1, true, mask);
    adjMatMultMat(ivector, Rect2, U2, Rect0, true, mask);
    sMulMat(ivector, Rect1, Rect1, kappa, mask);
    sMulMat(ivector, Rect2, Rect2, kappa, mask);
  }
}

/*
 * ____            __    __
 * ____| =  __  *  __| *    ,
 *  ____           __    __
 * |____ =  __  * |__  *    .    
 *
 * Rect += s * U1 * Staple * U2dag (UP)  or  s * U1dag * Staple * U2 (DN)
 *
 * updown = UP or DN.
 *
 */
//static 
void rectangle2( InstVector& ivector, vMatrix Rect, vMatrix U1, vMatrix Staple, vMatrix U2, vMatrix tmpMtr, string s, int updown, /*int add=1,*/ string mask ) 
{
  if(updown==UP) for(int r=0; r<3; ++r)
  {
      matMultVecT(ivector, tmpMtr[0], U1[r], Staple, false, mask);
      addSMulMatMultVecT(ivector, Rect[r], Rect[r], tmpMtr[0], U2, s, true, mask);
  }
  else {
      adjMatMultMat(ivector, tmpMtr, U1, Staple, true, mask);
      addSMulMatMultMat(ivector, Rect, Rect, tmpMtr, U2, s, false, mask);
  }
}

/*  __     __
 * |  \   |  |   |\
 *  __| =      *  |  *  __ .
 *
 * Bent (+)= U(dag) * Staple1 * Staple2 (backward & UP), 
 * 	     or Staple1(dag) * Staple2 * U (forward & UP),
 * 	     or U * Staple1 * Staple2(dag) (backward & DN),
 * 	     or Staple1 * Staple2 * U(dag) (forward & DN).
 *
 * backward = 0, forward = 1.
 */
//static 
void bentchair( InstVector& ivector, vMatrix Bent, vMatrix Staple1, vMatrix Staple2, vMatrix U, bool backforw, bool updown, int add )
{

}

/*  __     __
 * |  \   |  |   |\
 *    | =      *  | .
 *
 * Bent = Staple1(dag1) * Staple2(dag2) .
 *
 */ 
void bentchair0( InstVector& ivector, vMatrix Bent, vMatrix Staple1, vMatrix Staple2, bool dag1, bool dag2, int add, string mask )
{
  if(dag1&&dag2) {
    if(add==1) 
      for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) for(int k=0; k<3; ++k) {
        fmaddFVec(ivector, Bent[i][j][0], Staple2[j][k][0], Staple1[k][i][0], Bent[i][j][0], mask);
 	fnmaddFVec(ivector, Bent[i][j][0], Staple2[j][k][1], Staple1[k][i][1], Bent[i][j][0], mask);
	fnmaddFVec(ivector, Bent[i][j][1], Staple2[j][k][0], Staple1[k][i][1], Bent[i][j][1], mask);
	fnmaddFVec(ivector, Bent[i][j][1], Staple2[j][k][1], Staple1[k][i][0], Bent[i][j][1], mask);
      }
    else
      for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) { 
	for(int k=0; k<3; ++k) {
	  if(k==0) mulFVec(ivector, Bent[i][j][0], Staple2[j][0][0], Staple1[0][i][0], mask);
	  else fmaddFVec(ivector, Bent[i][j][0], Staple2[j][k][0], Staple1[k][i][0], Bent[i][j][0], mask);
	  fnmaddFVec(ivector, Bent[i][j][0], Staple2[j][k][1], Staple1[k][i][1], Bent[i][j][0], mask);
	  if(k==0) mulFVec(ivector, Bent[i][j][1], Staple2[j][0][0], Staple1[0][i][1], mask);
	  else fmaddFVec(ivector, Bent[i][j][1], Staple2[j][k][0], Staple1[k][i][1], Bent[i][j][1], mask);
	  fmaddFVec(ivector, Bent[i][j][1], Staple2[j][k][1], Staple1[k][i][0], Bent[i][j][1], mask);
        }
	smulFVec(ivector, Bent[i][j][1], Bent[i][j][1], "-1.0", mask);
      }
  }
  else if(dag1) {
    if(add==1) addAdjMatMultMat(ivector, Bent, Bent, Staple1, Staple2, dag1, mask);
    else adjMatMultMat(ivector, Bent, Staple1, Staple2, dag1, mask);
  }
  else {
    if(add==1) addMatMultMat(ivector, Bent, Bent, Staple1, Staple2, dag2, mask);
    else matMultMat(ivector, Bent, Staple1, Staple2, dag2, mask);
  }
}


