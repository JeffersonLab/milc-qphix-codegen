#ifndef _DATA_TYPES_H_
#define _DATA_TYPES_H_

#include "macros.h"
#include "instructions.h"

#if PRECISION == 1
#define FPTYPE float
#define LPTYPE short
#elif PRECISION == 2
#define FPTYPE double
#define LPTYPE float
#else
#error "Invalid FP PRECISION"
#endif

#ifndef USE_LP_SPINOR
typedef FPTYPE SpinorBaseType;
#define SpinorType	0
#else
typedef LPTYPE SpinorBaseType;
#define SpinorType	1
#endif

#ifndef USE_LP_GAUGE
typedef FPTYPE GaugeBaseType;
#define GaugeType	0
#else
typedef LPTYPE GaugeBaseType;
#define GaugeType	1
#endif

#ifndef USE_LP_HERMIT
typedef FPTYPE HermitBaseType;
typedef FPTYPE HermitHelperBaseType;
#define HermitType	0
#define HermitHelperType 0
#else
typedef LPTYPE HermitBaseType;
typedef LPTYPE HermitHelperBaseType;
#define HermitType	1
#define HermitHelperType 1
#endif

#ifndef USE_LP_CLOVER
typedef FPTYPE CloverBaseType;
#define CloverType	0
#else
typedef LPTYPE CloverBaseType;
#define CloverType	1
#endif

template<typename FT, int VL, bool compress12>
class data_types
{
public:
    typedef FT Spinor[4][3][2][VL];
    typedef FT Gauge[8][(compress12? 2 : 3)][3][2][VL];
    typedef FT Hermit[8][HERM_ELEM][VL];
    typedef struct {
        FT diag1[6][VL];
        FT off_diag1[15][2][VL];
        FT diag2[6][VL];
        FT off_diag2[15][2][VL];
    } Clover;
};

typedef FVec vComplex[2];
typedef FVec vMatrix[3][3][2];
typedef FVec vHerm[9];

void LoadSpinorElement(InstVector& ivector, const FVec& ret, const string& base, int spin, int col, int reim, int dir = -1);
void LoadFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base);
void StoreFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base, int isStreaming = 0);
void StreamFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base);
void PackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], const string& lBase, const string& rBase, int dir);
void UnpackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], const string& lBase, const string& rBase, int dir);
void LoadFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12);
void LoadFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12);
void StoreFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12, int isStreaming);
void StoreFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12, int isStreaming);
void StreamFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12);
void StreamFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12);
void LoadFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, const string& mask);
void LoadFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, const string& mask);
void StoreFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, int isStreaming);
void StoreFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, string mask);
void StoreFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, int isStreaming);
void StoreFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, string mask);
void StreamFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir);
void StreamFullHermitDir(InstVector& ivector, string ret, const string& base, string dir);
void LoadFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2);
void LoadFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2);
void StoreFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2, int isStreaming);
void StoreFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2, int isStreaming);
void StreamFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2);
void StreamFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2);
void LoadFullCloverBlock(InstVector& ivector, const FVec diag[6], const FVec off_diag[15][2], const string& base, int block);
void LoadKSElement(InstVector& ivector, const FVec& ret, const string& base, int col, int reim, int dir);
void LoadFullKS(InstVector& ivector, const FVec ret[3][2], const string& base);
void StoreFullKS(InstVector& ivector, const FVec ret[3][2], const string& base, int isStreaming=0);
void StreamFullKS(InstVector& ivector, const FVec ret[3][2], const string& base);
void PackKSSpinor(InstVector& ivector, const FVec ret[3][2], const string& lBase, const string& rBase, int dir);
void UnpackKSSpinor(InstVector& ivector, const FVec ret[3][2], const string& lBase, const string& rBase, int dir);
void PackGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& rBase, int dir, bool compress12);
void PackGaugeDir(InstVector& ivector, string ret, const string& rBase, string dir, bool compress12);
void PackGauge7WayDir(InstVector& ivector, const FVec ret[7][2][3][2], /*const string& lBase,*/ const string& rBase, int dir);
void PackGauge7WayDir(InstVector& ivector, const FVec ret[7][2][3][2], const string& lBase, const string& rBase, int dir);
void PackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir);
//void PackGauge7WayDir(InstVector& ivector, string ret, const string& rBase, string dir, bool compress12);
//void PackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir, bool compress12);
void UnpackGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& lBase, const string& rBase, int dir, bool compress12);
void UnpackGaugeDir(InstVector& ivector, const string& ret, const string& lBase, const string& rBase, const string& dir, bool compress12);
void UnpackGaugeDir(InstVector& ivector, const string& ret, const string& rBase, const string& dir, bool compress12);
void UnpackGauge7WayDir(InstVector& ivector, const FVec *ret[3][2], /*const string& lBase,*/ const string& rBase, int dir, bool compress12);
void UnpackGauge7WayDir(InstVector& ivector, const FVec *ret[3][2], const string& lBase, const string& rBase, int dir, bool compress12);
void UnpackGauge7WayDir(InstVector& ivector, const string ret, const string& lBase, const string& rBase, string dir, bool compress12);
//void UnpackGauge7WayDir(InstVector& ivector, string ret, /*const string& lBase,*/ const string& rBase, string dir, bool compress12);
//void UnpackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir, bool compress12);
void PackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& lBase, const string& rBase, int dir);
void PackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& rBase, int dir);
//void PackHermit7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir);
//void PackHermit7WayDir(InstVector& ivector, string ret, const string& rBase, string dir);
void PackHermitDir(InstVector& ivector, const FVec ret[9], const string& lBase, const string& rBase, int dir);
void PackHermitDir(InstVector& ivector, const FVec ret[9], const string& rBase, int dir);
void PackHermitDir(InstVector& ivector, const string& ret, const string& lBase, const string& rBase, const string& dir);
void PackHermitDir(InstVector& ivector, const string& ret, const string& rBase, const string& dir);
void UnpackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& lBase, const string& rBase, int dir);
void UnpackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& rBase, int dir);
void UnpackHermitDir(InstVector& ivector, const FVec ret[9], const string& rBase, int dir);
void UnpackHermitDir(InstVector& ivector, const FVec ret[9], const string& lBase, const string& rBase, int dir);
void UnpackHermitDir(InstVector& ivector, string ret, const string& rBase, string dir);
void UnpackHermitDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir);
//void UnpackHermit7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir);
//void UnpackHermit7WayDir(InstVector& ivector, string ret, const string& rBase, string dir);

void PrefetchL1FullSpinorDirIn(InstVector& ivector, const string& base, int dir = -1, int type = 0);
void PrefetchL1FullSpinorOut(InstVector& ivector, const string& base);
void PrefetchL2FullSpinorDirIn(InstVector& ivector, const string& base, const string& pref_dist, int dir=-1);
void PrefetchL2FullSpinorOut(InstVector& ivector, const string& base, const string& pref_dist);
void PrefetchL1FullKSDirIn(InstVector& ivector, const string& base, int dir=-1, int type=0);
void PrefetchL1FullKSOut(InstVector& ivector, const string& base);
void PrefetchL2FullKSDirIn(InstVector& ivector, const string& base, const string& pref_dist, int dir=-1);
void PrefetchL2FullKSOut(InstVector& ivector, const string& base, const string& pref_dist);
void PrefetchL1FullGaugeDirIn(InstVector& ivector, const string& base, int dir, bool compress12, int type = 0);
void PrefetchL2FullGaugeDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, bool compress12, int type = 0);
void PrefetchL2FullGaugeIn(InstVector& ivector, const string& base, const string& pref_dist, bool compress12, int type = 0);
void PrefetchL1FullHermitDirIn(InstVector& ivector, const string& base, int dir, int type = 0);
void PrefetchL2FullHermitDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, int type = 0);
void PrefetchL2FullHermitIn(InstVector& ivector, const string& base, const string& pref_dist, int type = 0);
void PrefetchL1FullHermitHelperDirIn(InstVector& ivector, const string& base, int dir1, int dir2, int type = 0);
void PrefetchL2FullHermitHelperDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, int type = 0);
void PrefetchL2FullHermitHelperIn(InstVector& ivector, const string& base, const string& pref_dist, int type = 0);
void PrefetchL1FullCloverBlockIn(InstVector& ivector, const string& base, int block);
void PrefetchL2FullCloverIn(InstVector& ivector, const string& base, const string& pref_dist);
void PrefetchL1HalfSpinorDir(InstVector& ivector, const string& lBase, const string& rBase, int dir, bool isPrefforWrite, int type);
void PrefetchL2HalfSpinorDir(InstVector& ivector, const string& lBase, const string& rBase, const string& pref_dist, int dir, bool isPrefforWrite, int type);

#endif // _DATA_TYPES_H_
