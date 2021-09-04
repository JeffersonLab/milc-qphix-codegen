#ifndef STAPLES_H
#define STAPLES_H

/*
#undef EXTERN
#ifdef _STAPLES_C
#define EXTERN 
#else
#define EXTERN extern
#endif
*/

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

#include "mat_utils.h"
#include "complex.h"

//		#2Dface
//		    |
//		    ^
//EXTERN FVec staple1[6][16][3][3][2]; 

//EXTERN FVec staple2[6][16][3][3][2];
void staple( InstVector& ivector, vMatrix Staple, vMatrix U1, vMatrix U2, vMatrix U3, int updown, string mask="" );
void rectangle0( InstVector& ivector, vMatrix Rect0, vMatrix Staple1, vMatrix Staple2, int add=1, string mask="" );
void rectangle1( InstVector& ivector, vMatrix Rect1, vMatrix Rect2, vMatrix Rect0, vMatrix U1, vMatrix U2, string s, int add=1, string mask="" );
void rectangle2( InstVector& ivector, vMatrix Rect, vMatrix U1, vMatrix Staple, vMatrix U2, vMatrix tmpMtr, string s, int updown, string mask="" );
//void bentchair( InstVector& ivector, vMatrix Bent, vMatrix Staple1, vMatrix Staple2, vMatrix U, bool backforw, bool updown, int add=1 );
void bentchair0( InstVector& ivector, vMatrix Bent, vMatrix Staple1, vMatrix Staple2, bool dag1, bool dag2, int add=1, string mask="" );

void staple( InstVector& ivector, string Staple, string U1, string U2, string U3, int updown, string mask="" );
void rectangle0( InstVector& ivector, string Rect0, string Staple1, string Staple2, int add=1, string mask="" );
void rectangle1( InstVector& ivector, string Rect1, string Rect2, string Rect0, string U1, string U2, string s, int add=1, string mask="" );
void rectangle2( InstVector& ivector, string Rect, string U1, string Staple, string U2, string tmpMtr, string s, int updown, string mask="" );
void bentchair0( InstVector& ivector, string Bent, string Staple1, string Staple2, bool dag1, bool dag2, int add=1, string mask="" );
#endif
