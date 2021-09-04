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

#include "mat_utils.h"
#include "complex.h"

extern std::string ARCH_NAME;

string getTypeName(size_t s);

void declare_Hermit(InstVector &ivector, vHerm herm);
void declare_Hermit7W(InstVector &ivector, vHerm *herm);
void declare_Gauge12(InstVector &ivector, FVec matr[2][3][2]);
void declare_Gauge7Way(InstVector &ivector, FVec matr[7][2][3][2], bool compress12);
void loadGaugeDir(InstVector &ivector, vMatrix out, string inGBase, int dir, bool compress12, string mask);
void loadGauge12Dir(InstVector &ivector, FVec out[2][3][2], string &inGBase, int dir, string mask);
void loadGaugeDir(InstVector &ivector, vMatrix out, FVec norm, string &inGBase, int dir, bool compress12, string mask);
void loadHermitDir(InstVector &ivector, vHerm out, string &inHBase, int dir, string mask);
void loadHermitHelperDir(InstVector &ivector, vHerm out, string &inHBase, int dir1, int dir2, string mask);
void loadHermitYZTDir(InstVector &ivector, vHerm out, string &inHBase, int dir1, int dir2, string mask);
void storeGaugeDir(InstVector &ivector, string &outGBase, vMatrix in, int dir, string mask, bool compress12, int isStreaming);
void storeHermitDir(InstVector &ivector, string &outHBase, vHerm in, int dir, string mask, int isStreaming);
void storeHermitDir(InstVector &ivector, string &outHBase, vHerm in, int dir, string mask);
void storeHermitHelperDir(InstVector &ivector, string &outHBase, vHerm in, int dir1, int dir2, string mask, int isStreaming);
void unpackGaugeDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, bool compress12, string mask);
void unpackGauge7WayDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, bool compress12, string mask);
void unpackHermitDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, string mask);
void unpackHermitDir(InstVector &ivector, string ret, string rbuf, string dir, string mask);

void declare_Gauge(InstVector &ivector, FVec gauge);
void declare_Hermit(InstVector &ivector, FVec herm);
void declare_HermTl(InstVector &ivector, FVec herm);
void declare_Hermit7Way(InstVector &ivector, FVec herm);
void declare_Gauge12(InstVector &ivector, FVec matr);
void declare_Gauge7Way(InstVector &ivector, FVec matr, bool compress12);
void loadGaugeDir(InstVector &ivector, string out, string inGBase, string dir, bool compress12, string mask);
void loadGauge12Dir(InstVector &ivector, string out, string &inGBase, string dir, string mask);
void loadGaugeDir(InstVector &ivector, string out, FVec norm, string inGBase, string dir, bool compress12, string mask);
void loadHermitDir(InstVector &ivector, string out, string inHBase, string dir, string mask);
void loadHermitYZTDir(InstVector &ivector, string out, string inHBase, string dir1, string dir2, string mask);
void loadHermitHelperDir(InstVector &ivector, string out, string inHBase, string dir1, string dir2, string mask);
void storeGaugeDir(InstVector &ivector, string outGBase, string in, string dir, string mask, bool compress12, int isStreaming);
void storeHermitDir(InstVector &ivector, string outHBase, string in, string dir, string mask, int isStreaming);
void storeHermitDir(InstVector &ivector, string outHBase, string in, string dir, string mask);
void storeHermitHelperDir(InstVector &ivector, string outHBase, string in, string dir1, string dir2, string mask, int isStreaming);

