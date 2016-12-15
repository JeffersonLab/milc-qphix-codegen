#include "data_utils.h"

string getTypeName(size_t s)
{
    if(s == 2) return "half";
    else if(s == 4) return "float";
    else if(s == 8) return "double";
    else return "Unknown";
}

void declare_Hermit(InstVector &ivector, vHerm herm)
{
  for(int k=0; k<9; ++k)
    declareFVecFromFVec(ivector, herm[k]);
}

void declare_Hermit7Way(InstVector &ivector, vHerm *herm)
{
  for(int i=0; i<7; ++i)
    declare_Hermit(ivector, herm[i]);
}

void declare_Gauge12(InstVector &ivector, FVec matr[2][3][2])
{
  for(int r=0; r<2; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
    declareFVecFromFVec(ivector, matr[r][c][ir]);
}

void declare_Gauge7Way(InstVector &ivector, FVec matr[7][2][3][2], bool compress12)
{
  for(int r=0; r<7; ++r) for(int i=0; i<3-compress12; ++i) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir)
    declareFVecFromFVec(ivector, matr[r][i][c][ir]);
}

void loadGaugeDir(InstVector &ivector, vMatrix out, string inGBase, int dir, bool compress12, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, compress12);
  decompressGauge(ivector, out, compress12, mask);
}

void loadGauge12Dir(InstVector &ivector, FVec out[2][3][2], string &inGBase, int dir, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, true);
}

void loadGaugeDir(InstVector &ivector, vMatrix out, FVec norm, string &inGBase, int dir, bool compress12, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, compress12);
  decompressGauge(ivector, out, norm, compress12, mask);
}

void loadHermitDir(InstVector &ivector, vHerm out, string &inHBase, int dir, string mask="")
{
  LoadFullHermitDir(ivector, out, inHBase, dir, mask);
}

void loadHermitYZTDir(InstVector &ivector, vHerm out, string &inHBase, int dir1, int dir2, string mask="")
{
  LoadFullHermitHelperDir(ivector, out, inHBase, dir1, dir2-((dir1+2)/2>dir2/2 ? 0 : ((dir1+2)/2<dir2/2 ? 1 : dir2%2)));
}

void loadHermitHelperDir(InstVector &ivector, vHerm out, string &inHBase, int dir1, int dir2, string mask="")
{
  LoadFullHermitHelperDir(ivector, out, inHBase, dir1, dir2);
}

void storeGaugeDir(InstVector &ivector, string &outGBase, vMatrix in, int dir, string mask, bool compress12, int isStreaming)
{
  StoreFullGaugeDir(ivector, in, outGBase, dir, compress12, isStreaming);
}

void storeHermitDir(InstVector &ivector, string &outHBase, vHerm in, int dir, string mask, int isStreaming)
{
  StoreFullHermitDir(ivector, in, outHBase, dir, isStreaming);
}

void storeHermitDir(InstVector &ivector, string &outHBase, vHerm in, int dir, string mask)
{
  StoreFullHermitDir(ivector, in, outHBase, dir, mask);
}

/*
void storeHermitXYZTDir(InstVector &ivector, string &outHBase, vHerm in, int dir1, int dir2, string mask, int isStreaming)
{
  StoreHermitTDir(ivector, in, outHBase, dir1, dir2, isStreaming);
}
*/

void storeHermitHelperDir(InstVector &ivector, string &outHBase, vHerm in, int dir1, int dir2, string mask, int isStreaming)
{
  StoreFullHermitHelperDir(ivector, in, outHBase, dir1, dir2, isStreaming);
}

void unpackGaugeDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, bool compress12, string mask)
{
  UnpackGaugeDir(ivector, ret, lbuf, rbuf, dir, compress12);
  decompressGauge(ivector, ret, compress12, mask);
}

void unpackGauge7WayDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, bool compress12, string mask)
{
  UnpackGauge7WayDir(ivector, ret, lbuf, rbuf, dir, compress12);
  if(compress12) 
  {
    forLoopInc(ivector, "d", "0", "7");
    {
	decompressGauge(ivector, ret+"[d]", compress12, mask);
    }
    endScope(ivector);
  }
}

void unpackHermitDir(InstVector &ivector, string ret, string rbuf, string lbuf, string dir, string mask)
{
  UnpackHermitDir(ivector, ret, lbuf, rbuf, dir);
}

void unpackHermitDir(InstVector &ivector, string ret, string rbuf, string dir, string mask)
{
  UnpackHermitDir(ivector, ret, rbuf, dir);
}

void declare_Gauge(InstVector &ivector, FVec gauge)
{
  int len[3]={3, 3, 2};
  declareFVecArrayFromFVec(ivector, gauge, 3, len);
}

void declare_Hermit(InstVector &ivector, FVec herm)
{
  int len[1] = {9};
  declareFVecArrayFromFVec(ivector, herm, 1, len);
}

void declare_HermTl(InstVector &ivector, FVec herm)
{
  int len[1] = {8};
  declareFVecArrayFromFVec(ivector, herm, 1, len);
}

void declare_Hermit7Way(InstVector &ivector, FVec herm)
{
  int len[2] = {7, 9};
  declareFVecArrayFromFVec(ivector, herm, 2, len);
}

void declare_Gauge12(InstVector &ivector, FVec matr)
{
  int len[3] = {2, 3, 2};
  declareFVecArrayFromFVec(ivector, matr, 3, len);
}

void declare_Gauge7Way(InstVector &ivector, FVec matr, bool compress12)
{
  int len[4] = {7, 3, 3, 2};
  if(compress12) len[1]--;
  declareFVecArrayFromFVec(ivector, matr, 4, len);
}

void loadGaugeDir(InstVector &ivector, string out, string inGBase, string dir, bool compress12, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, compress12);
  decompressGauge(ivector, out, compress12, mask);
}

void loadGauge12Dir(InstVector &ivector, string out, string &inGBase, string dir, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, true);
}

void loadGaugeDir(InstVector &ivector, string out, FVec norm, string &inGBase, string dir, bool compress12, string mask="")
{
  LoadFullGaugeDir(ivector, out, inGBase, dir, compress12);
  decompressGauge(ivector, out, norm, compress12, mask);
}

void loadHermitDir(InstVector &ivector, string out, string inHBase, string dir, string mask="")
{
  LoadFullHermitDir(ivector, out, inHBase, dir, mask);
}

void loadHermitYZTDir(InstVector &ivector, string out, string inHBase, string dir1, string dir2, string mask="")
{
  LoadFullHermitHelperDir(ivector, out, inHBase, dir1, dir2+"-((("+dir1+"+2)>>1)>("+dir2+">>1) ? 0 : ((("+dir1+"+2)>>1)<("+dir2+">>1) ? 1 : ("+dir2+"&1)))");
}

void loadHermitHelperDir(InstVector &ivector, string out, string inHBase, string dir1, string dir2, string mask="")
{
  LoadFullHermitHelperDir(ivector, out, inHBase, dir1, dir2);
}

void storeGaugeDir(InstVector &ivector, string outGBase, string in, string dir, string mask, bool compress12, int isStreaming)
{
  StoreFullGaugeDir(ivector, in, outGBase, dir, compress12, isStreaming);
}

void storeHermitDir(InstVector &ivector, string outHBase, string in, string dir, string mask, int isStreaming)
{
  StoreFullHermitDir(ivector, in, outHBase, dir, isStreaming);
}

void storeHermitDir(InstVector &ivector, string outHBase, string in, string dir, string mask)
{
  StoreFullHermitDir(ivector, in, outHBase, dir, mask);
}

void storeHermitHelperDir(InstVector &ivector, string outHBase, string in, string dir1, string dir2, string mask, int isStreaming)
{
  StoreFullHermitHelperDir(ivector, in, outHBase, dir1, dir2, isStreaming);
}


