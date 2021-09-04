
#include "data_types_xyzt.h"
#include "instructions.h"

#if 0
typedef SpinorBaseType Spinor[4][3][2][VECLEN];
typedef GaugeBaseType Gauge[8][3][3][2][VECLEN];
typedef struct {
    CloverBaseType diag1[6][VECLEN];
    CloverBaseType off_diag1[15][2][VECLEN];
    CloverBaseType diag2[6][VECLEN];
    CloverBaseType off_diag2[15][2][VECLEN];
} Clover;
#endif


void readFVecSpinor(InstVector& ivector, const FVec& ret, const string& base, int spin, int color, int reim)
{
    loadFVec(ivector, ret, new SpinorAddress(base,spin,color,reim,SpinorType), string(""));
}

void writeFVecSpinor(InstVector& ivector, const FVec& ret, const string& base, int spin, int color, int reim, int isStreaming)
{
    storeFVec(ivector, ret, new SpinorAddress(base,spin,color,reim,SpinorType), isStreaming);
}

void readFVecKS(InstVector& ivector, const FVec& ret, const string& base, int color, int reim)
{
    loadFVec(ivector, ret, new KSAddress(base,color,reim,SpinorType), string(""));
}

void writeFVecKS(InstVector& ivector, const FVec& ret, const string& base, int color, int reim, int isStreaming)
{
    storeFVec(ivector, ret, new KSAddress(base,color,reim,SpinorType), isStreaming);
}

void readFVecGauge(InstVector& ivector, const FVec& ret, const string& base, int dir, int c1, int c2, int reim)
{
    loadFVec(ivector, ret, new GaugeAddress(base,dir,c1,c2,reim,GaugeType), string(""));
}

void readFVecGauge(InstVector& ivector, const FVec& ret, const string& base, const string& dir, const string& c1, const string& c2, int reim)
{
    loadFVec(ivector, ret, new GaugeAddressStr(base,dir,c1,c2,reim,GaugeType), string(""));
}

void writeFVecGauge(InstVector& ivector, const FVec& ret, const string& base, int dir, int c1, int c2, int reim, int isStreaming)
{
    storeFVec(ivector, ret, new GaugeAddress(base,dir,c1,c2,reim,GaugeType), isStreaming);
}

void writeFVecGauge(InstVector& ivector, const FVec& ret, const string& base, const string& dir, const string& c1, const string& c2, int reim, int isStreaming)
{
    storeFVec(ivector, ret, new GaugeAddressStr(base,dir,c1,c2,reim,GaugeType), isStreaming);
}

void readFVecHermit(InstVector& ivector, const FVec& ret, const string& base, int dir, int k, const string& mask="")
{
    loadFVec(ivector, ret, new HermitAddress(base,dir,k,HermitType), mask);
}

void readFVecHermit(InstVector& ivector, const FVec& ret, const string& base, const string& dir, const string& k, const string& mask="")
{
    loadFVec(ivector, ret, new HermitAddressStr(base,dir,k,HermitType), mask);
}

void readFVecHermitHelper(InstVector& ivector, const FVec& ret, const string& base, int dir1, int dir2, int k)
{
    loadFVec(ivector, ret, new HermitHelperAddress(base,dir1,dir2,k,HermitHelperType), string(""));
}

void readFVecHermitHelper(InstVector& ivector, const FVec& ret, const string& base, const string& dir1, const string& dir2, const string& k)
{
    loadFVec(ivector, ret, new HermitHelperAddressStr(base,dir1,dir2,k,HermitHelperType), string(""));
}

void writeFVecHermit(InstVector& ivector, const FVec& ret, const string& base, int dir, int k, int isStreaming)
{
    storeFVec(ivector, ret, new HermitAddress(base,dir,k,HermitType), isStreaming);
}

void writeFVecHermit(InstVector& ivector, const FVec& ret, const string& base, const string& dir, const string& k, int isStreaming)
{
    storeFVec(ivector, ret, new HermitAddressStr(base,dir,k,HermitType), isStreaming);
}

void writeFVecHermit(InstVector& ivector, const FVec& ret, const string& base, int dir, int k, string mask)
{
    storeFVec(ivector, ret, new HermitAddress(base,dir,k,HermitType), mask);
}

void writeFVecHermit(InstVector& ivector, const FVec& ret, const string& base, const string& dir, const string& k, string mask)
{
    storeFVec(ivector, ret, new HermitAddressStr(base,dir,k,HermitType), mask);
}

void writeFVecHermitHelper(InstVector& ivector, const FVec& ret, const string& base, int dir1, int dir2, int k, int isStreaming)
{
    storeFVec(ivector, ret, new HermitHelperAddress(base,dir1,dir2,k,HermitHelperType), isStreaming);
}

void writeFVecHermitHelper(InstVector& ivector, const FVec& ret, const string& base, const string& dir1, const string& dir2, const string& k, int isStreaming)
{
    storeFVec(ivector, ret, new HermitHelperAddressStr(base,dir1,dir2,k,HermitHelperType), isStreaming);
}

void readFVecClovDiag(InstVector& ivector, const FVec& ret, const string& base, int block, int c)
{
    loadFVec(ivector, ret, new ClovDiagAddress(base,block,c,CloverType), string(""));
}

void readFVecClovOffDiag(InstVector& ivector, const FVec& ret, const string& base, int block, int c, int reim)
{
    loadFVec(ivector, ret, new ClovOffDiagAddress(base,block,c,reim,CloverType), string(""));
}

void LoadSpinorElement(InstVector& ivector, const FVec& ret, const string& base, int spin, int col, int reim, int dir)
{
    readFVecSpinor(ivector, ret, base, spin, col, reim);
}

void LoadFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base)
{
    for(int col=0; col < 3; col++) {
        for(int spin=0; spin < 4; spin++) {
            LoadSpinorElement(ivector, ret[spin][col][RE], base, spin, col, RE, -1);
            LoadSpinorElement(ivector, ret[spin][col][IM], base, spin, col, IM, -1);
        }
    }
}

void StoreFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base, int isStreaming)
{

#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif

    for(int spin=0; spin < 4; spin++) {
        for(int col=0; col < 3; col++) {
            writeFVecSpinor(ivector, ret[spin][col][RE], base, spin, col, RE, isStreaming);
            writeFVecSpinor(ivector, ret[spin][col][IM], base, spin, col, IM, isStreaming);
        }
    }
}

void StreamFullSpinor(InstVector& ivector, const FVec ret[4][3][2], const string& base)
{
    StoreFullSpinor(ivector, ret, base, 1);
}

void PackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], const string& lBase, const string& rBase, int dir)
{
	int lpackoffset = 0, rpackoffset = 0;
    for(int spin=0; spin < 2; spin++) {
        for(int col=0; col < 3; col++) {
			int rsz = packXYZTFVec(ivector, ret[spin][col], new AddressImm(new GenericAddress(lBase, SpinorType), lpackoffset),
				new AddressImm(new GenericAddress(rBase, SpinorType), rpackoffset), dir);
			rpackoffset += rsz;
			if(rsz == VECLEN) lpackoffset += rsz;
		}
	}
}

void UnpackHalfSpinor(InstVector& ivector, const FVec ret[2][3][2], const string& lBase, const string& rBase, int dir)
{
	int lunpackoffset = 0, runpackoffset = 0;
    for(int spin=0; spin < 2; spin++) {
        for(int col=0; col < 3; col++) {
			int rsz = unpackXYZTFVec(ivector, ret[spin][col], new AddressImm(new GenericAddress(lBase, SpinorType), lunpackoffset),
				new AddressImm(new GenericAddress(rBase, SpinorType), runpackoffset), dir);
			runpackoffset += rsz;
			if(rsz == VECLEN) lunpackoffset += rsz;
		}
	}
}

void LoadFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12)
{
    int nrows = 3;

    if(compress12 == true) {
        nrows = 2;
    }

    for(int c1=0; c1 < nrows; c1++) {
        for(int c2=0; c2 < 3; c2++) {
            readFVecGauge(ivector, ret[c1][c2][RE], base, dir, c1, c2, RE);
            readFVecGauge(ivector, ret[c1][c2][IM], base, dir, c1, c2, IM);
        }
    }
}

void LoadFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12)
{
    if(compress12 == true) forLoopInc(ivector, "c1", "0", "2");
    else forLoopInc(ivector, "c1", "0", "3");
    {
        forLoopInc(ivector, "c2", "0", "3");
	{
	    readFVecGauge(ivector, FVec(ret+"[c1][c2][0]"), base, dir, "c1", "c2", RE);
	    readFVecGauge(ivector, FVec(ret+"[c1][c2][1]"), base, dir, "c1", "c2", IM);
	}
	endScope(ivector);
    }
    endScope(ivector);
}

void StoreFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    int nrows = 3;
    if(compress12 == true) {
        nrows = 2;
    }
    for(int c1=0; c1 < nrows; c1++) 
        for(int c2=0; c2 < 3; c2++) {
        writeFVecGauge(ivector, ret[c1][c2][RE], base, dir, c1, c2, RE, isStreaming);
        writeFVecGauge(ivector, ret[c1][c2][IM], base, dir, c1, c2, IM, isStreaming);
    }
}

void StoreFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    if(compress12 == true) forLoopInc(ivector, "c1", "0", "2");
    else forLoopInc(ivector, "c1", "0", "3");
    {
        forLoopInc(ivector, "c2", "0", "3");
        {
	    writeFVecGauge(ivector, FVec(ret+"[c1][c2][0]"), base, dir, "c1", "c2", RE, isStreaming);
	    writeFVecGauge(ivector, FVec(ret+"[c1][c2][1]"), base, dir, "c1", "c2", IM, isStreaming);
	}
	endScope(ivector);
    }
    endScope(ivector);
}

void StreamFullGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& base, int dir, bool compress12)
{
    StoreFullGaugeDir(ivector, ret, base, dir, compress12, 1);
}

void StreamFullGaugeDir(InstVector& ivector, string ret, const string& base, string dir, bool compress12)
{
    StoreFullGaugeDir(ivector, ret, base, dir, compress12, 1);
}

void LoadFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, const string& mask="")
{
    for(int k=0; k<HERM_ELEM; ++k)
	readFVecHermit(ivector, ret[k], base, dir, k, mask);
}

void LoadFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, const string& mask="")
{
    forLoopInc(ivector, "k", "0", "8");
    {
	readFVecHermit(ivector, FVec(ret+"[k]"), base, dir, "k", mask);
    }
    endScope(ivector);
}

void StoreFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    for(int col=0; col < HERM_ELEM; col++) {
        writeFVecHermit(ivector, ret[col], base, dir, col, isStreaming);
    }
}

void StoreFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir, string mask)
{
    for(int col=0; col < HERM_ELEM; col++) {
        writeFVecHermit(ivector, ret[col], base, dir, col, mask);
    }
}

void StoreFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    forLoopInc(ivector, "k", "0", "8");
    {
	writeFVecHermit(ivector, FVec(ret+"[k]"), base, dir, "k", isStreaming);
    }
    endScope(ivector);
}

void StoreFullHermitDir(InstVector& ivector, string ret, const string& base, string dir, string mask)
{
    forLoopInc(ivector, "k", "0", "8");
    {
        writeFVecHermit(ivector, FVec(ret+"[k]"), base, dir, "k", mask);
    }
    endScope(ivector);
}

void StreamFullHermitDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir)
{
    StoreFullHermitDir(ivector, ret, base, dir, 1);
}

void StreamFullHermitDir(InstVector& ivector, string ret, const string& base, string dir)
{
    StoreFullHermitDir(ivector, ret, base, dir, 1);
}

void LoadFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2)
{
    for(int k=0; k<HERM_ELEM; ++k)
        readFVecHermitHelper(ivector, ret[k], base, dir1, dir2, k);
}

void LoadFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2)
{
    forLoopInc(ivector, "k", "0", "8");
    {
	readFVecHermitHelper(ivector, FVec(ret+"[k]"), base, dir1, dir2, "k");
    }
    endScope(ivector);
}

void StoreFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    for(int col=0; col < HERM_ELEM; col++) {
	writeFVecHermitHelper(ivector, ret[col], base, dir1, dir2, col, isStreaming);
    }
}

void StoreFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2, int isStreaming)
{
#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif
    forLoopInc(ivector, "k", "0", "8");
    {
	writeFVecHermitHelper(ivector, FVec(ret+"[k]"), base, dir1, dir2, "k", isStreaming);
    }
    endScope(ivector);
}

void StreamFullHermitHelperDir(InstVector& ivector, const FVec ret[HERM_ELEM], const string& base, int dir1, int dir2)
{
    StoreFullHermitHelperDir(ivector, ret, base, dir1, dir2, 1);
}

void StreamFullHermitHelperDir(InstVector& ivector, string ret, const string& base, string dir1, string dir2)
{
    StoreFullHermitHelperDir(ivector, ret, base, dir1, dir2, 1);
}

void LoadFullCloverBlock(InstVector& ivector, const FVec diag[6], const FVec off_diag[15][2], const string& base, int block)
{
    for(int c=0; c < 6; c++) {
        readFVecClovDiag(ivector, diag[c], base, block, c);
    }

    for(int c=0; c < 15; c++) {
        readFVecClovOffDiag(ivector, off_diag[c][RE], base, block, c, RE);
        readFVecClovOffDiag(ivector, off_diag[c][IM], base, block, c, IM);
    }
}


void LoadKSElement(InstVector& ivector, const FVec& ret, const string& base, int col, int reim, int dir)
{
    readFVecKS(ivector, ret, base, col, reim);
}
void LoadFullKS(InstVector& ivector, const FVec ret[3][2], const string& base)
{
    for(int col=0; col < 3; col++) {
        LoadKSElement(ivector, ret[col][RE], base, col, RE, -1);
        LoadKSElement(ivector, ret[col][IM], base, col, IM, -1);
    }
}


void StoreFullKS(InstVector& ivector, const FVec ret[3][2], const string& base, int isStreaming)
{

#ifndef ENABLE_STREAMING_STORES
    isStreaming = 0;
#endif

    for(int col=0; col < 3; col++) {
        writeFVecKS(ivector, ret[col][RE], base, col, RE, isStreaming);
        writeFVecKS(ivector, ret[col][IM], base, col, IM, isStreaming);
    }
}

void StreamFullKS(InstVector& ivector, const FVec ret[3][2], const string& base)
{
    StoreFullKS(ivector, ret, base, 1);
}

void PackKSSpinor(InstVector& ivector, const FVec ret[3][2], const string& lBase, const string& rBase, int dir)
{
	int lpackoffset = 0, rpackoffset = 0;
    for(int col=0; col < 3; col++) {
		int rsz = packXYZTFVec(ivector, ret[col], new AddressImm(new GenericAddress(lBase, SpinorType), lpackoffset),
			new AddressImm(new GenericAddress(rBase, SpinorType), rpackoffset), dir);
		rpackoffset += rsz;
		if(rsz == VECLEN) lpackoffset += rsz;
	}
}

void UnpackKSSpinor(InstVector& ivector, const FVec ret[3][2], const string& lBase, const string& rBase, int dir)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int col=0; col < 3; col++) {
		int rsz = unpackXYZTFVec(ivector, ret[col], new AddressImm(new GenericAddress(lBase, SpinorType), lunpackoffset),
			new AddressImm(new GenericAddress(rBase, SpinorType), runpackoffset), dir);
		runpackoffset += rsz;
		if(rsz == VECLEN) lunpackoffset += rsz;
    }
}

void PackGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& rBase, int dir, bool compress12)
{
    int rpackoffset = 0;
    for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
	int rsz = packXYZTFVec(ivector, ret[r][c], 
		new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir);
	rpackoffset += rsz;
    }
}

void PackGaugeDir(InstVector& ivector, string ret, const string& rBase, string dir, bool compress12)
{
    if(compress12)
        forLoopStatement(ivector, "int r=0, rsz=0, rpackoffset = 0", "r<2", "++r");
    else
	forLoopStatement(ivector, "int r=0, rsz=0, rpackoffset = 0", "r<3", "++r");
    {
            forLoopStatement(ivector, "int c=0", "c<3", "++c, rpackoffset += rsz");
            {
                packXYZTFVec(ivector, ret+"[r][c]", new AddressImmStr(new GenericAddress(rBase, GaugeType), "rpackoffset"), dir, "rsz");
            }
            endScope(ivector);
    }
    endScope(ivector); 
}

void PackGauge7WayDir(InstVector& ivector, const FVec ret[7][2][3][2], /*const string& lBase,*/ const string& rBase, int dir)
{
    int rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<2; r++) for(int c=0; c<3; ++c) 
    {
	int rsz = packXYZTFVec(ivector, ret[d][r][c], /*new AddressImm(new GenericAddress(lBase, GaugeType), lpackoffset),*/
		new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir);
	rpackoffset += rsz;
    }
}

void PackGauge7WayDir(InstVector& ivector, const FVec ret[7][2][3][2], const string& lBase, const string& rBase, int dir)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<2; r++) for(int c=0; c<3; ++c)
    {        int rsz = packXYZTFVec(ivector, ret[d][r][c], new AddressImm(new GenericAddress(lBase, GaugeType), lpackoffset),
		new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir);
        rpackoffset += rsz;
        if(rsz == VECLEN) lpackoffset += rsz;
    }
}

void PackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir)
{
    forLoopStatement(ivector, "int d=0, rsz=0, lpackoffset = 0, rpackoffset = 0", "d<7", "++d");
    {
	forLoopInc(ivector, "r", "0", "2");
	{
	    forLoopStatement(ivector, "int c=0", "c<3", "++c, rpackoffset += rsz");
	    {
		packXYZTFVec(ivector, ret+"[d][r][c]", new AddressImmStr(new GenericAddress(lBase, GaugeType), "lpackoffset"), new AddressImmStr(new GenericAddress(rBase, GaugeType), "rpackoffset"), dir, "rsz");
		ifStatement(ivector, "rsz == VECLEN");
		{
		    initInt(ivector, "lpackoffset", "lpackoffset+VECLEN");
		}
		endScope(ivector);
	    }
	    endScope(ivector);
	}
	endScope(ivector);
    }
    endScope(ivector);
}

void UnpackGaugeDir(InstVector& ivector, const FVec ret[3][3][2], const string& lBase, const string& rBase, int dir, bool compress12)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int r=0; r<3-compress12; ++r) for(int c=0; c<3; ++c)
    {
	int rsz = unpackXYZTFVec(ivector, ret[r][c], new AddressImm(new GenericAddress(lBase, GaugeType), lpackoffset),
		new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir);
	rpackoffset += rsz;
	if(rsz == VECLEN) lpackoffset += rsz;
    }
}

void UnpackGaugeDir(InstVector& ivector, const string& ret, const string& lBase, const string& rBase, const string& dir, bool compress12)
{
    int lpackoffset = 0, rpackoffset = 0;
    if(compress12) forLoopStatement(ivector, "int r=0, rsz=0, lpackoffset = 0, rpackoffset = 0", "r<2", "++r");
    else forLoopStatement(ivector, "int r=0, lpackoffset = 0, rpackoffset = 0", "r<3", "++r");
    forLoopStatement(ivector, "int c=0", "c<3", "++c, rpackoffset += rsz");
    {
	unpackXYZTFVec(ivector, ret+"[r][c]", new AddressImmStr(new GenericAddress(lBase, GaugeType), "lpackoffset"),
		new AddressImmStr(new GenericAddress(rBase, GaugeType), "rpackoffset"), dir, "rsz");
	//if(rsz==VECLEN) lpackoffset+=VECLEN;
	ifStatement(ivector, "rsz==VECLEN");
	{
		initInt(ivector, "lpackoffset", "lpackoffset+VECLEN");
	}
	endScope(ivector);
    }
    endScope(ivector);
    endScope(ivector);
}

void UnpackGaugeDir(InstVector& ivector, const string& ret, const string& rBase, const string& dir, bool compress12)
{
    if(compress12) forLoopStatement(ivector, "int r=0, rsz=0, rpackoffset = 0", "r<2", "++r");
    else forLoopStatement(ivector, "int r=0, rpackoffset = 0", "r<3", "++r");
    forLoopStatement(ivector, "int c=0", "c<3", "++c, rpackoffset += rsz");
    {
	unpackXYZTFVec(ivector, ret+"[r][c]", new AddressImmStr(new GenericAddress(rBase, GaugeType), "rpackoffset"), dir, "rsz");
    }
    endScope(ivector);
    endScope(ivector);
}

#if 0
void PackHermit7WayDir(InstVector& ivector, const FVec *ret[3][2], const string& lBase, const string& rBase, int dirl, int dirr, bool compress12)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
        int rsz = unpackXYZTFVec(ivector, ret[d*(3-compress12)+r][c], new AddressImm(new GenericAddress(lBase, GaugeType), lunpackoffset),
                new AddressImm(new GenericAddress(rBase, GaugeType), runpackoffset), dir);
        runpackoffset += rsz;
        if(rsz == VECLEN) lunpackoffset += rsz;
    }
}

void PackGauge7WayDir(InstVector& ivector, string ret, const string& rBase, string dir, bool compress12)
{
    int rpackoffset = 0;
    for(int d=0; d<7; ++d)
    for(int c1=0; c1<3-compress12; ++c1) for(int c2=0; c2<3; ++c2)
    {
	ostringstream tmp;
	tmp << "[" << d << "][" << c1 << "][" << c2 << "]";
        int rsz = packXYZTFVec(ivector, ret+tmp.str(), /*new AddressImm(new GenericAddress(lBase, GaugeType), lpackoffset),*/
                new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir); 
	rpackoffset += rsz;
    }
}

void PackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir, bool compress12)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {        
	ostringstream tmp;
	tmp << "[" << d << "][" << r << "][" << c << "]";
	int rsz = packXYZTFVec(ivector, ret+tmp.str(), new AddressImm(new GenericAddress(lBase, GaugeType), lpackoffset),
                new AddressImm(new GenericAddress(rBase, GaugeType), rpackoffset), dir);
        rpackoffset += rsz;
        if(rsz == VECLEN) lpackoffset += rsz;
    }
}
#endif

void UnpackGauge7WayDir(InstVector& ivector, const FVec *ret[3][2], /*const string& lBase,*/ const string& rBase, int dir, bool compress12)
{
    int runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
	int rsz = unpackXYZTFVec(ivector, ret[d*(3-compress12)+r][c], /*new AddressImm(new GenericAddress(lBase, GaugeType), lunpackoffset),*/
                new AddressImm(new GenericAddress(rBase, GaugeType), runpackoffset), dir);
        runpackoffset += rsz;
    }
}

void UnpackGauge7WayDir(InstVector& ivector, const FVec *ret[3][2], const string& lBase, const string& rBase, int dir, bool compress12)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
        int rsz = unpackXYZTFVec(ivector, ret[d*(3-compress12)+r][c], new AddressImm(new GenericAddress(lBase, GaugeType), lunpackoffset),
		new AddressImm(new GenericAddress(rBase, GaugeType), runpackoffset), dir);
        runpackoffset += rsz;
        if(rsz == VECLEN) lunpackoffset += rsz;
    }
}

void UnpackGauge7WayDir(InstVector& ivector, const string ret, const string& lBase, const string& rBase, string dir, bool compress12)
{
    forLoopStatement(ivector, "int d=0, rsz=0, lpackoffset = 0, rpackoffset = 0", "d<7", "++d");
    {
        if(compress12) forLoopInc(ivector, "r", "0", "2");
	else forLoopInc(ivector, "r", "0", "3");
        {
            forLoopStatement(ivector, "int c=0", "c<3", "++c, rpackoffset += rsz");
            {
                unpackXYZTFVec(ivector, ret+"[d][r][c]", new AddressImmStr(new GenericAddress(lBase, GaugeType), "lpackoffset"), new AddressImmStr(new GenericAddress(rBase, GaugeType), "rpackoffset"), dir, "rsz");
                ifStatement(ivector, "rsz == VECLEN");
                {
                    initInt(ivector, "lpackoffset", "lpackoffset+VECLEN");
                }
                endScope(ivector);
            }
            endScope(ivector);
        }
        endScope(ivector);
    }
    endScope(ivector);
}

#if 0
void UnpackGauge7WayDir(InstVector& ivector, string ret, /*const string& lBase,*/ const string& rBase, string dir, bool compress12)
{
    int runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
	ostringstream tmp;
        tmp << "[" << d << "][" << r << "][" << c << "]";
        int rsz = unpackXYZTFVec(ivector, ret+tmp.str(), /*new AddressImm(new GenericAddress(lBase, GaugeType), lunpackoffset),*/
                new AddressImm(new GenericAddress(rBase, GaugeType), runpackoffset), dir);
        runpackoffset += rsz;
    }
}

void UnpackGauge7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir, bool compress12)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<3-compress12; r++) for(int c=0; c<3; ++c)
    {
        ostringstream tmp;
	tmp << "[" << d << "][" << r << "][" << c << "]";
	int rsz = unpackXYZTFVec(ivector, ret+tmp.str(), new AddressImm(new GenericAddress(lBase, GaugeType), lunpackoffset),
                new AddressImm(new GenericAddress(rBase, GaugeType), runpackoffset), dir);
        runpackoffset += rsz;
        if(rsz == VECLEN) lunpackoffset += rsz;
    }
}
#endif

void PackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& lBase, const string& rBase, int dir)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<4; r++)
    {
        int rsz = packXYZTFVec(ivector, ret[d]+2*r, new AddressImm(new GenericAddress(lBase, HermitType), lpackoffset),
		new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
        rpackoffset += rsz;
	if(rsz == VECLEN) lpackoffset += rsz;
    }
}

void PackHermitDir(InstVector& ivector, const FVec ret[9], const string& lBase, const string& rBase, int dir)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int r=0; r<4; r++)
    {
	int rsz = packXYZTFVec(ivector, ret+2*r, new AddressImm(new GenericAddress(lBase, HermitType), lpackoffset),
		new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
	rpackoffset += rsz;
	if(rsz == VECLEN) lpackoffset += rsz;
    }
}

void PackHermitDir(InstVector& ivector, const FVec ret[9], const string& rBase, int dir)
{
    int rpackoffset = 0;
    for(int r=0; r<4; r++)
    {
        int rsz = packXYZTFVec(ivector, ret+2*r, new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
	rpackoffset += rsz;
    }
}

void PackHermitDir(InstVector& ivector, const string& ret, const string& lBase, const string& rBase, const string& dir)
{
/*
    forLoopStatement(ivector, "int r=0, rsz=0, lpackoffset = 0, rpackoffset = 0", "r<4", "++r, rpackoffset += rsz");
    {
	int rsz = packXYZTFVec(ivector, ret+"+2*r", new AddressImm(new GenericAddress(lBase, HermitType), lpackoffset),
		new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
	if(rsz==VECLEN) initInt(ivector, "rsz", "VECLEN");
	else initInt(ivector, "rsz", "2*VECLEN");
	if(rsz == VECLEN) initInt(ivector, "lpackoffset", "lpackoffset+VECLEN");
    }
    endScope(ivector);
*/
    //int lpackoffset = 0, rpackoffset = 0;
    //for(int r=0; r<4; r++)
    forLoopStatement(ivector, "int r=0, rsz=0, lpackoffset = 0, rpackoffset = 0", "r<4", "++r, rpackoffset += rsz");
    {
	packXYZTFVec(ivector, "("+ret+"+2*r)", new AddressImmStr(new GenericAddress(lBase, HermitType), "lpackoffset"),
		new AddressImmStr(new GenericAddress(rBase, HermitType), "rpackoffset"), dir, "rsz");
	//rpackoffset += rsz;
	ifStatement(ivector, "rsz == VECLEN");
	{
		initInt(ivector, "lpackoffset", "lpackoffset+rsz");
	}
	endScope(ivector);
    }
    endScope(ivector);
}

void PackHermitDir(InstVector& ivector, const string& ret, const string& rBase, const string& dir)
{
    forLoopStatement(ivector, "int r=0, rsz=0, rpackoffset = 0", "r<4", "++r, rpackoffset += rsz");
    {
	packXYZTFVec(ivector, "("+ret+"+2*r)", new AddressImmStr(new GenericAddress(rBase, HermitType), "rpackoffset"), dir, "rsz");
    }
    endScope(ivector);
}

void PackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& rBase, int dir)
{
    int rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<4; r++)
    {
        int rsz = packXYZTFVec(ivector, ret[d]+2*r, new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
        rpackoffset += rsz;
    }
}

void UnpackHermitDir(InstVector& ivector, const FVec ret[9], const string& rBase, int dir)
{
    int rpackoffset = 0;
    for(int r=0; r<4; r++)
    {
	int rsz = unpackXYZTFVec(ivector, ret+2*r, new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
	rpackoffset += rsz;
    }
}

void UnpackHermitDir(InstVector& ivector, const FVec ret[9], const string& lBase, const string& rBase, int dir)
{
    int rpackoffset = 0, lpackoffset = 0;
    for(int r=0; r<4; r++)
    {
	int rsz = unpackXYZTFVec(ivector, ret+2*r, new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset),
		new AddressImm(new GenericAddress(lBase, HermitType), lpackoffset), dir);
	rpackoffset += rsz;
	if(rsz==VECLEN) lpackoffset += rsz;
    }
}

void UnpackHermitDir(InstVector& ivector, string ret, const string& rBase, string dir)
{
    int rpackoffset = 0;
    //for(int r=0; r<4; r++)
    forLoopStatement(ivector, "int r=0, rpackoffset=0, rsz=0", "r<4", "++r, rpackoffset += rsz");
    {
	/*std::stringstream Str;
	Str << (2*r);
	string r2 = Str.str();
	*/
	unpackXYZTFVec(ivector, "("+ret+"+2*r)", new AddressImmStr(new GenericAddress(rBase, HermitType), "rpackoffset"), dir, "rsz");
	//rpackoffset += rsz;
    }
    endScope(ivector);
}

void UnpackHermitDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir)
{
    int rpackoffset = 0, lpackoffset = 0;
    //for(int r=0; r<4; r++)
    forLoopStatement(ivector, "int r=0, rpackoffset=0, lpackoffset=0, rsz=0", "r<4", "++r, rpackoffset += rsz");
    {
	/*std::stringstream Str;
        Str << (2*r);
        string r2 = Str.str();
	*/
	unpackXYZTFVec(ivector, "("+ret+"+2*r)", new AddressImmStr(new GenericAddress(rBase, HermitType), "rpackoffset"),
                new AddressImmStr(new GenericAddress(lBase, HermitType), "lpackoffset"), dir, "rsz");
        //rpackoffset += rsz;
        //if(rsz==VECLEN) lpackoffset += rsz;
        ifStatement(ivector, "rsz==VECLEN");
	{
		initInt(ivector, "lpackoffset", "lpackoffset+VECLEN");
	}
	endScope(ivector);
    }
    endScope(ivector);
}

#if 0
void PackHermit7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir)
{
    int lpackoffset = 0, rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<4; r++)
    {
	ostringstream tmp;
        tmp << "[" << d << "][" << 2*r << "]";
        int rsz = packXYZTFVec(ivector, "((vector *)&"+ret+tmp.str()+")", new AddressImm(new GenericAddress(lBase, HermitType), lpackoffset),
                new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
        rpackoffset += rsz;
        if(rsz == VECLEN) lpackoffset += rsz;
    }
}

void PackHermit7WayDir(InstVector& ivector, string ret, const string& rBase, string dir)
{
    int rpackoffset = 0;
    for(int d=0; d<7; ++d) for(int r=0; r<4; r++)
    {
	ostringstream tmp;
	tmp << "[" << d << "][" << 2*r << "]";
        int rsz = packXYZTFVec(ivector, "((vector *)&"+ret+tmp.str()+")", new AddressImm(new GenericAddress(rBase, HermitType), rpackoffset), dir);
        rpackoffset += rsz;
    }
}
#endif

void UnpackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& lBase, const string& rBase, int dir)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int k=0; k<4; ++k)
    {
	int rsz = unpackXYZTFVec(ivector, ret[d]+2*k, new AddressImm(new GenericAddress(lBase, HermitType), lunpackoffset),
		new AddressImm(new GenericAddress(rBase, HermitType), runpackoffset), dir);
	runpackoffset += rsz;
	if(rsz == VECLEN) lunpackoffset += rsz;
    }
}

void UnpackHermit7WayDir(InstVector& ivector, const FVec *ret[9], const string& rBase, int dir)
{
    int runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int k=0; k<4; ++k)
    {
	int rsz = unpackXYZTFVec(ivector, ret[d]+2*k, new AddressImm(new GenericAddress(rBase, HermitType), runpackoffset), dir);
	runpackoffset += rsz;
    }
}

#if 0
void UnpackHermit7WayDir(InstVector& ivector, string ret, const string& lBase, const string& rBase, string dir)
{
    int lunpackoffset = 0, runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int k=0; k<4; ++k)
    {
	ostringstream tmp;
        tmp << "[" << d << "][" << 2*k << "]";
        int rsz = unpackXYZTFVec(ivector, "((vector *)&"+ret+tmp.str()+")", new AddressImm(new GenericAddress(lBase, HermitType), lunpackoffset),
                new AddressImm(new GenericAddress(rBase, HermitType), runpackoffset), dir);
        runpackoffset += rsz;
        if(rsz == VECLEN) lunpackoffset += rsz;
    }
}

void UnpackHermit7WayDir(InstVector& ivector, string ret, const string& rBase, string dir)
{
    int runpackoffset = 0;
    for(int d=0; d<7; ++d) for(int k=0; k<4; ++k)
    {
	ostringstream tmp;
        tmp << "[" << d << "][" << 2*k << "]";
        int rsz = unpackXYZTFVec(ivector, "((vector *)&"+ret+tmp.str()+")", new AddressImm(new GenericAddress(rBase, HermitType), runpackoffset), dir);
        runpackoffset += rsz;
    }
}
#endif

// Prefetches

void prefetchL1SpinorIn(InstVector& ivector, string base, int imm, int dir, int type)
{
#ifdef PREF_L1_SPINOR_IN
    PrefetchL1(ivector, new AddressImm(new SpinorAddress(base,0,0,0,SpinorType), imm), type, dir);
#endif
}

void prefetchL2SpinorIn(InstVector& ivector, string base, const string& pref_dist, int imm, int dir)
{
#ifdef PREF_L2_SPINOR_IN
    PrefetchL2(ivector, new AddressImm(new AddressOffset(new SpinorAddress(base,0,0,0,SpinorType), pref_dist), imm));
#endif
}

void prefetchL1SpinorOut(InstVector& ivector, string base, int imm)
{
#ifdef PREF_L1_SPINOR_OUT
    PrefetchL1(ivector, new AddressImm(new SpinorAddress(base,0,0,0,SpinorType), imm), 3 /*NT & EX */);
#endif
}

void prefetchL2SpinorOut(InstVector& ivector, string base, const string& pref_dist, int imm)
{
#ifdef PREF_L2_SPINOR_OUT
    PrefetchL2(ivector, new AddressImm(new AddressOffset(new SpinorAddress(base,0,0,0,SpinorType), pref_dist), imm), 3 /*NT & EX */);
#endif
}

void prefetchL1KSIn(InstVector& ivector, string base, int imm, int dir, int type)
{
#ifdef PREF_L1_SPINOR_IN
    PrefetchL1(ivector, new AddressImm(new KSAddress(base,0,0,SpinorType), imm), type, dir);
#endif
}

void prefetchL2KSIn(InstVector& ivector, string base, const string& pref_dist, int imm, int dir)
{
#ifdef PREF_L2_SPINOR_IN
    PrefetchL2(ivector, new AddressImm(new AddressOffset(new KSAddress(base,0,0,SpinorType), pref_dist), imm));
#endif
}

void prefetchL1KSOut(InstVector& ivector, string base, int imm)
{
#ifdef PREF_L1_SPINOR_OUT
    PrefetchL1(ivector, new AddressImm(new KSAddress(base,0,0,SpinorType), imm), 3 /*NT & EX */);
#endif
}

void prefetchL2KSOut(InstVector& ivector, string base, const string& pref_dist, int imm)
{
#ifdef PREF_L2_SPINOR_OUT
    PrefetchL2(ivector, new AddressImm(new AddressOffset(new KSAddress(base,0,0,SpinorType), pref_dist), imm), 3 /*NT & EX */);
#endif
}

void prefetchL1GuageDirIn(InstVector& ivector, string base, int dir, int imm, int type)
{
#ifdef PREF_L1_GAUGE
    prefetchL1(ivector, new AddressImm(new GaugeAddress(base,dir,0,0,0,GaugeType), imm), type);
#endif
}

void prefetchL2GuageDirIn(InstVector& ivector, string base, int dir, const string& pref_dist, int imm, int type)
{
#ifdef PREF_L2_GAUGE
    prefetchL2(ivector, new AddressImm(new AddressOffset(new GaugeAddress(base,dir,0,0,0,GaugeType), pref_dist), imm), type);
#endif
}

void prefetchL1HermitDirIn(InstVector& ivector, string base, int dir, int imm, int type)
{
#ifdef PREF_L1_HERMIT
    prefetchL1(ivector, new AddressImm(new HermitAddress(base,dir,0,HermitType), imm), type);
#endif
}

void prefetchL2HermitDirIn(InstVector& ivector, string base, int dir, const string& pref_dist, int imm, int type)
{
#ifdef PREF_L2_HERMIT
    prefetchL2(ivector, new AddressImm(new AddressOffset(new HermitAddress(base,dir,0,HermitType), pref_dist), imm), type);
#endif
}

void prefetchL1HermitHelperDirIn(InstVector& ivector, string base, int dir1, int dir2, int imm, int type)
{
#ifdef PREF_L1_HERMIT
    prefetchL1(ivector, new AddressImm(new HermitHelperAddress(base,dir1, dir2, 0,HermitHelperType), imm), type);
#endif
}

void prefetchL2HermitHelperDirIn(InstVector& ivector, string base, int dir, const string& pref_dist, int imm, int type)
{
#ifdef PREF_L2_HERMIT
    prefetchL2(ivector, new AddressImm(new AddressOffset(new HermitHelperAddress(base,dir,0,0,HermitHelperType), pref_dist), imm), type);
#endif
}

void prefetchL1CloverBlockIn(InstVector& ivector, string base, int block, int imm)
{
#ifdef PREF_L1_CLOVER
    prefetchL1(ivector, new AddressImm(new ClovDiagAddress(base,block,0,CloverType), imm), 0);
#endif
}

void prefetchL2CloverBlockIn(InstVector& ivector, string base, int block, const string& pref_dist, int imm)
{
#ifdef PREF_L2_CLOVER
    prefetchL2(ivector, new AddressImm(new AddressOffset(new ClovDiagAddress(base,block,0,CloverType), pref_dist), imm), 0);
#endif
}

void PrefetchL1FullSpinorDirIn(InstVector& ivector, const string& base, int dir, int type)
{
    // for now we ignore direction but it can be used for specialization
    for(int i = 0; i < (24*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1SpinorIn(ivector, base, i*(64/sizeof(SpinorBaseType)), dir, type);
    }
}

void PrefetchL1FullSpinorOut(InstVector& ivector, const string& base, const string& off)
{
    for(int i = 0; i < (24*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1SpinorOut(ivector, base, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL2FullSpinorDirIn(InstVector& ivector, const string& base, const string& pref_dist, int dir)
{
    // for now we ignore direction but itcan be used for specialization
    for(int i = 0; i < (24*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2SpinorIn(ivector, base, pref_dist, i*(64/sizeof(SpinorBaseType)), dir);
    }
}

void PrefetchL2FullSpinorOut(InstVector& ivector, const string& base, const string& pref_dist)
{
    for(int i = 0; i < (24*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2SpinorOut(ivector, base, pref_dist, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL1FullKSDirIn(InstVector& ivector, const string& base, int dir, int type)
{
    // for now we ignore direction but it can be used for specialization
    for(int i = 0; i < (6*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1KSIn(ivector, base, i*(64/sizeof(SpinorBaseType)), dir, type);
    }
}

void PrefetchL1FullKSOut(InstVector& ivector, const string& base, const string& off)
{
    for(int i = 0; i < (6*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1KSOut(ivector, base, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL2FullKSDirIn(InstVector& ivector, const string& base, const string& pref_dist, int dir)
{
    // for now we ignore direction but itcan be used for specialization
    for(int i = 0; i < (6*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2KSIn(ivector, base, pref_dist, i*(64/sizeof(SpinorBaseType)), dir);
    }
}

void PrefetchL2FullKSOut(InstVector& ivector, const string& base, const string& pref_dist)
{
    for(int i = 0; i < (6*VECLEN*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2KSOut(ivector, base, pref_dist, i*(64/sizeof(SpinorBaseType)));
    }
}

void PrefetchL1FullGaugeDirIn(InstVector& ivector, const string& base, int dir, bool compress12, int type)
{
    int nSites = VECLEN;
    int g_size=0;

    if( compress12 ) {
        g_size=2*3*2;
    }
    else {
        g_size=3*3*2;
    }

    for(int i = 0; i < ((g_size*nSites*sizeof(GaugeBaseType)+63)/64); i++) {
        prefetchL1GuageDirIn(ivector, base, dir, i*(64/sizeof(GaugeBaseType)), type);
    }
}

void PrefetchL2FullGaugeDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, bool compress12, int type)
{
    int nSites = VECLEN;
    int g_size=0;

    if( compress12 ) {
        g_size=2*3*2;
    }
    else {
        g_size=3*3*2;
    }

    for(int i = 0; i < ((g_size*nSites*sizeof(GaugeBaseType)+63)/64); i++) {
        prefetchL2GuageDirIn(ivector, base, dir, pref_dist, i*(64/sizeof(GaugeBaseType)), type);
    }
}

void PrefetchL2FullGaugeIn(InstVector& ivector, const string& base, const string& pref_dist, bool compress12, int type)
{
    for(int dir = 0; dir < 8; dir++) {
        PrefetchL2FullGaugeDirIn(ivector, base, dir, pref_dist, compress12, type);
    }
}

void PrefetchL1FullHermitDirIn(InstVector& ivector, const string& base, int dir, int type)
{
    int nSites = VECLEN;
    int h_size=HERM_ELEM;

    for(int i = 0; i < ((h_size*nSites*sizeof(HermitBaseType)+63)/64); i++) {
        prefetchL1HermitDirIn(ivector, base, dir, i*(64/sizeof(HermitBaseType)), type);
    }
}

void PrefetchL2FullHermitDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, int type)
{
    int nSites = VECLEN;
    int h_size=HERM_ELEM;
    for(int i = 0; i < ((h_size*nSites*sizeof(HermitBaseType)+63)/64); i++) {
        prefetchL2HermitDirIn(ivector, base, dir, pref_dist, i*(64/sizeof(HermitBaseType)), type);
    }
}

void PrefetchL2FullHermitIn(InstVector& ivector, const string& base, const string& pref_dist, int type)
{
    for(int dir = 0; dir < 8; dir++) 
	PrefetchL2FullHermitDirIn(ivector, base, dir, pref_dist, type);
}

void PrefetchL1FullHermitHelperDirIn(InstVector& ivector, const string& base, int dir1, int dir2, int type)
{
    int nSites = VECLEN;
    int h_size=HERM_ELEM;

    for(int i = 0; i < ((h_size*nSites*sizeof(HermitBaseType)+63)/64); i++) {
        prefetchL1HermitHelperDirIn(ivector, base, dir1, dir2, i*(64/sizeof(HermitBaseType)), type);
    }
}

void PrefetchL2FullHermitHelperDirIn(InstVector& ivector, const string& base, int dir, const string& pref_dist, int type)
{
    int nSites = VECLEN;
    int h_size=HERM_ELEM*7;
    for(int i = 0; i < ((h_size*nSites*sizeof(HermitBaseType)+63)/64); i++) {
        prefetchL2HermitHelperDirIn(ivector, base, dir, pref_dist, i*(64/sizeof(HermitHelperBaseType)), type);
    }
}

void PrefetchL2FullHermitHelperIn(InstVector& ivector, const string& base, const string& pref_dist, int type)
{
    for(int dir = 0; dir < 8; dir++)
	PrefetchL2FullHermitHelperDirIn(ivector, base, dir, pref_dist, type);
}

void PrefetchL1FullCloverBlockIn(InstVector& ivector, const string& base, int block)
{
    int nSites = VECLEN;

    for(int i = 0; i < ((36*nSites*sizeof(CloverBaseType)+63)/64); i++) {
        prefetchL1CloverBlockIn(ivector, base, block, i*(64/sizeof(CloverBaseType)));
    }
}

void PrefetchL2FullCloverIn(InstVector& ivector, const string& base, const string& pref_dist)
{
    int nSites = VECLEN;

    for(int i = 0; i < ((2*36*nSites*sizeof(CloverBaseType)+63)/64); i++) {
        prefetchL2CloverBlockIn(ivector, base, 0, pref_dist, i*(64/sizeof(CloverBaseType)));
    }
}

void PrefetchL1HalfSpinorDir(InstVector& ivector, const string& lBase, const string& rBase, int dir, bool isPrefforWrite, int type)
{
	int split = 0;
	if((1 << dir/2) < VECLEN) split = 1;
	int size = (split == 1 ? VECLEN / 2 : VECLEN);

#ifdef ENABLE_STREAMING_STORES

    if(isPrefforWrite) {
        return;
    }

#endif

    for(int i = 0; i < (12*size*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL1(ivector, new AddressImm(new GenericAddress(rBase, SpinorType), i*(64/sizeof(SpinorBaseType))), type);
		if(split == 1) 
	        prefetchL1(ivector, new AddressImm(new GenericAddress(lBase, SpinorType), i*(64/sizeof(SpinorBaseType))), type);
    }
}

void PrefetchL2HalfSpinorDir(InstVector& ivector, const string& lBase, const string& rBase, const string& pref_dist, int dir, bool isPrefforWrite, int type)
{
	int split = 0;
	if((1 << dir/2) < VECLEN) split = 1;
	int size = (split == 1 ? VECLEN / 2 : VECLEN);

#ifdef ENABLE_STREAMING_STORES

    if(isPrefforWrite) {
        return;
    }

#endif

    for(int i = 0; i < (12*size*sizeof(SpinorBaseType)+63)/64; i++) {
        prefetchL2(ivector, new AddressImm(new AddressOffset(new GenericAddress(rBase, SpinorType), pref_dist), i*(64/sizeof(SpinorBaseType))), type);
		if(split == 1) 
	        prefetchL2(ivector, new AddressImm(new AddressOffset(new GenericAddress(lBase, SpinorType), pref_dist), i*(64/sizeof(SpinorBaseType))), type);
    }
}

