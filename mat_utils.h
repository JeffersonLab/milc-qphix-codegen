#ifndef MAT_UTILS_H
#define MAT_UTILS_H

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

#include "complex.h"
#include "dslash_common.h"

void decompressGauge(InstVector& ivector, vMatrix u_gauge, bool compress12, string mask="");
void decompressGauge(InstVector& ivector, vMatrix u_gauge, FVec norm, bool compress12, string mask="");
void movMat(InstVector& ivector, vMatrix out, vMatrix in, string mask="");
void addMat(InstVector& ivector, vMatrix out, vMatrix in1, vMatrix in2, string mask="");
void sMulMat(InstVector& ivector, vMatrix out, vMatrix in, string s, string mask="");
void addSMulMat(InstVector& ivector, vMatrix out, vMatrix in1, vMatrix in2, string s, string mask="");
void matHerm(InstVector& ivector, FVec *herm, vMatrix matr, string mask="");
void matAntiHerm(InstVector& ivector, FVec *aherm, vMatrix matr, string mask="");
void matTrls(InstVector& ivector, vMatrix trls, vMatrix matr, string mask="");
void hermTrls(InstVector& ivector, FVec *trls, FVec *tmp, string mask="");
void matAntiHermTrls(InstVector& ivector, FVec *trls, vMatrix matr, string mask="");
void matMultVec(InstVector& ivector, FVec r[3][2], FVec u_mat[3][3][2], FVec s_vec[3][2], bool adjMul, string mask="");
void matMultVecT(InstVector& ivector, FVec r[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], bool adjMul, string mask="");
void addMatMultVecT(InstVector& ivector, FVec r[3][2], FVec add[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], bool adjMul, string mask="");
void addSMulMatMultVecT(InstVector& ivector, FVec r[3][2], FVec add[3][2], FVec s_vec[3][2], FVec u_mat[3][3][2], string s, bool adjMul, string mask="");
void matMultMat(InstVector& ivector, FVec u_ret[3][3][2], FVec u_mat1[3][3][2], FVec u_mat2[3][3][2], bool adjMul, string mask="");
void adjMatMultMat(InstVector& ivector, FVec u_ret[3][3][2], FVec u_mat1[3][3][2], FVec u_mat2[3][3][2], bool adjMul, string mask="");
void matMultMatSconj(InstVector& ivector, vMatrix matr, vMatrix su3, vMatrix matr1, string mask="");
void addMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, bool conj, string mask="");
void addSMulMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, string s, bool conj, string mask="");
void addAdjMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, bool conj, string mask="");
void addSMulAdjMatMultMat(InstVector& ivector, vMatrix out, vMatrix add, vMatrix in1, vMatrix in2, string s, bool conj, string mask="");
void addHerm(InstVector& ivector, FVec *out, FVec *in1, FVec *in2, string mask="");
void addHermTl(InstVector& ivector, FVec *out, FVec *in1, FVec *in2, string mask="");
void addSMultAntiHermIc(InstVector& ivector, FVec *out, FVec *in1, string scalar, FVec *in2, string mask="");
void hermAddSMultAntiHermIc(InstVector& ivector, FVec *out, FVec *in1, string scalar, FVec *in2, string mask="");
//void addSMultHermMultMatIc(InstVector& ivector, vMatrix out, FVec scalar, FVec *herm1, vMatrix in2, string mask="");
void addSMultHermMultMatIc(InstVector& ivector, vMatrix out, string scalar, FVec *herm1, vMatrix in2, string mask="");
void addSMultMatMultHermIc(InstVector& ivector, vMatrix out, FVec scalar, vMatrix in2, FVec *herm1, string mask="");
void sMultAntiHermIc(InstVector& ivector, FVec *out, string scalar, FVec *in1, string mask="");

void decompressGauge(InstVector& ivector, string u_gauge, bool compress12, string mask="");
void decompressGauge(InstVector& ivector, string u_gauge, FVec norm, bool compress12, string mask="");
void movMat(InstVector& ivector, string out, string in, string mask="");
void addMat(InstVector& ivector, string out, string in1, string in2, string mask="");
void sMulMat(InstVector& ivector, string out, string in, string s, string mask="");
void addSMulMat(InstVector& ivector, string out, string in1, string in2, string s, string mask="");
void matHerm(InstVector& ivector, string herm, string matr, string mask="");
void matAntiHerm(InstVector& ivector, string aherm, string matr, string mask="");
void matTrls(InstVector& ivector, string trls, string matr, string mask="");
void hermTrls(InstVector& ivector, string trls, string tmp, string mask="");
void matAntiHermTrls(InstVector& ivector, string trls, string matr, string mask="");
void matMultVec(InstVector& ivector, string r, string u_mat, string s_vec, bool adjMul, string mask="");
void matMultVecT(InstVector& ivector, string r, string s_vec, string u_mat, bool adjMul, string mask="");
void addMatMultVecT(InstVector& ivector, string r, string add, string s_vec, string u_mat, bool adjMul, string mask="");
void addSMulMatMultVecT(InstVector& ivector, string r, string add, string s_vec, string u_mat, string s, bool adjMul, string mask="");
void matMultMat(InstVector& ivector, string u_ret, string u_mat1, string u_mat2, bool adjMul, string mask="");
void adjMatMultMat(InstVector& ivector, string u_ret, string u_mat1, string u_mat2, bool adjMul, string mask="");
void matMultMatSconj(InstVector& ivector, string matr, string su3, string matr1, string mask="");
void addMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, bool conj, string mask="");
void addSMulMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, string s, bool conj, string mask="");
void addAdjMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, bool conj, string mask="");
void addSMulAdjMatMultMat(InstVector& ivector, string out, string add, string in1, string in2, string s, bool conj, string mask="");
void addHerm(InstVector& ivector, string out, string in1, string in2, string mask="");
void addHermTl(InstVector& ivector, string out, string in1, string in2, string mask="", bool add=true);
void addSMultAntiHermIc(InstVector& ivector, string out, string in1, string scalar, string in2, string mask="");
void hermAddSMultAntiHermIc(InstVector& ivector, string out, string in1, string scalar, string in2, string mask="");
void addSMultHermMultMatIc(InstVector& ivector, string out, string scalar, string herm1, string in2, /*string tmp, */string mask="");
void sMultAntiHermIc(InstVector& ivector, string out, string scalar, string in1, string mask="");

#endif
