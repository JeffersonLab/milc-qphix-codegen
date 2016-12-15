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

extern string tmpMtr3str;
extern string mask;

/*
 * Staple = U1(dag) * U2 * U3(!dag)
 */
//static 
void staple( InstVector& ivector, string Staple, string U1, string U2, string U3, int updown, string mask )
{
  bool dag = (updown==DN);
  adjMatMultMat(ivector, tmpMtr3str, U1, U2, dag, mask);
  matMultMat(ivector, Staple, tmpMtr3str, U3, !dag, mask);
}

/*
 *  ____      __        __
 * |    | =  |  |   *  |  | .
 *
 * Rect0 += Staple1 * Staple2 (add==1)
 * Rect0  = Staple1 * Staple2 (add==0)
 *
 */ 
void rectangle0( InstVector& ivector, string Rect0, string Staple1, string Staple2, int add, string mask )
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
void rectangle1( InstVector& ivector, string Rect1, string Rect2, string Rect0, string U1, string U2, /* bool order, */ /*int updown,*/ string kappa, int add, string mask )
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
void rectangle2( InstVector& ivector, string Rect, string U1, string Staple, string U2, string tmpMtr, string s, int updown, /*int add=1,*/ string mask ) 
{
  if(updown==UP) 
  {
    forLoopInc(ivector, "sr", "0", "3");
    {
      matMultVecT(ivector, tmpMtr+"[0]", U1+"[sr]", Staple, false, mask);
      addSMulMatMultVecT(ivector, Rect+"[sr]", Rect+"[sr]", tmpMtr+"[0]", U2, s, true, mask);
    }
    endScope(ivector);
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
void bentchair( InstVector& ivector, string Bent, string Staple1, string Staple2, string U, bool backforw, bool updown, int add )
{

}

/*  __     __
 * |  \   |  |   |\
 *    | =      *  | .
 *
 * Bent = Staple1(dag1) * Staple2(dag2) .
 *
 */ 
void bentchair0( InstVector& ivector, string Bent, string Staple1, string Staple2, bool dag1, bool dag2, int add, string mask )
{
  if(dag1&&dag2) {
    if(add==1) 
    {
      forLoopInc(ivector,"si", "0", "3");
      forLoopInc(ivector,"sj", "0", "3");
      forLoopInc(ivector,"sk", "0", "3");
      {
        fmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][sk][0]"), FVec(Staple1+"[sk][si][0]"), FVec(Bent+"[si][sj][0]"), mask);
 	fnmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][sk][1]"), FVec(Staple1+"[sk][si][1]"), FVec(Bent+"[si][sj][0]"), mask);
	fnmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][sk][0]"), FVec(Staple1+"[sk][si][1]"), FVec(Bent+"[si][sj][1]"), mask);
	fnmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][sk][1]"), FVec(Staple1+"[sk][si][0]"), FVec(Bent+"[si][sj][1]"), mask);
      }
      endScope(ivector);
      endScope(ivector);
      endScope(ivector);
    }
    else
    {
      forLoopInc(ivector,"si", "0", "3");
      forLoopInc(ivector,"sj", "0", "3");
      {
      	{
	  mulFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][0][0]"), FVec(Staple1+"[0][si][0]"), mask);
	  fnmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][0][1]"), FVec(Staple1+"[0][si][1]"), FVec(Bent+"[si][sj][0]"), mask);
	  mulFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][0][0]"), FVec(Staple1+"[0][si][1]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][0][1]"), FVec(Staple1+"[0][si][0]"), FVec(Bent+"[si][sj][1]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][1][0]"), FVec(Staple1+"[1][si][0]"), FVec(Bent+"[si][sj][0]"), mask);
	  fnmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][1][1]"), FVec(Staple1+"[1][si][1]"), FVec(Bent+"[si][sj][0]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][1][0]"), FVec(Staple1+"[1][si][1]"), FVec(Bent+"[si][sj][1]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][1][1]"), FVec(Staple1+"[1][si][0]"), FVec(Bent+"[si][sj][1]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][2][0]"), FVec(Staple1+"[2][si][0]"), FVec(Bent+"[si][sj][0]"), mask);
	  fnmaddFVec(ivector, FVec(Bent+"[si][sj][0]"), FVec(Staple2+"[sj][2][1]"), FVec(Staple1+"[2][si][1]"), FVec(Bent+"[si][sj][0]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][2][0]"), FVec(Staple1+"[2][si][1]"), FVec(Bent+"[si][sj][1]"), mask);
	  fmaddFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Staple2+"[sj][2][1]"), FVec(Staple1+"[2][si][0]"), FVec(Bent+"[si][sj][1]"), mask);
        }
	smulFVec(ivector, FVec(Bent+"[si][sj][1]"), FVec(Bent+"[si][sj][1]"), "-1.0", mask);
      }
      endScope(ivector);
      endScope(ivector);
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


