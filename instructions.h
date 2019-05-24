#ifndef __INSTRUCTIONS_H__
#define __INSTRUCTIONS_H__

#include <string>
#include <sstream>
#include <vector>

#include <fstream>
#include <iostream>
#include <typeinfo>
#include <map>

#include <cstdio>
#include <cstdlib>

using namespace std;

#ifndef PRECISION
#define PRECISION 1
#endif

#include "address_types.h"

using namespace std;

#define RE 0
#define IM 1
#define MAXARRAYDIM 10

class FVec
{
public:
    FVec(const string& name_);
    FVec(const FVec& v_) : name(v_.getName()), type(v_.getType()) {}
    const string& getName() const
    {
        return name;
    }
    const string& getType() const
    {
        return type;
    }

private:
    const string name;
    const string type;
};

class Instruction
{
public:
    // string class return empty string
    virtual string serialize() const
    {
        return string("");
    }
    virtual bool hasAddress() const
    {
        return false;
    }
    virtual int numArithmeticInst() const
    {
        return 0;
    }
    virtual int numDeclarations() const
    {
        return 0;
    }
    virtual int numScopes() const
    {
        return 0;
    }
    virtual int numIfs() const
    {
        return 0;
    }
    virtual int numLoops() const
    {
	return 0;
    }
};

typedef vector<Instruction *> InstVector;

class BeginScope : public Instruction
{
public:
    string serialize() const
    {
        return "{";
    }
    int numScopes() const
    {
        return 1;
    }
};

class EndScope : public Instruction
{
public:
    string serialize() const
    {
        return "}";
    }
    int numScopes() const
    {
        return 1;
    }
};

class ElseStatement : public Instruction
{
public:
    string serialize() const
    {
        return "} else {";
    }
};

class ElseIfStatement : public Instruction {
public:
    ElseIfStatement(const string& condition_) : condition(condition_) {}
    string serialize() const {
        return "} else if ( "+condition+" ) { ";
    }
    int numIfs() const {
        return 1;
    }
private:
    const string condition;
};

class IfStringCond : public Instruction
{
public:
    IfStringCond(const string& condition_) : condition(condition_) {}
    string serialize() const
    {
        return " if ( "+condition+" ) { ";
    }
    int numIfs() const
    {
        return 1;
    }
    int numScopes() const
    {
	return 1;
    }
private:
    const string condition;
};

class ForLoopInc : public Instruction
{
public:
    ForLoopInc(const string& index_, const string& begin_, const string& end_, const string& inc_) : index(index_), begin(begin_), end(end_), inc(inc_) {}
    string serialize() const
    {
	return " for ( int "+index+"="+begin+"; "+index+"<"+end+"; "+index+"+="+inc+" ) { ";
    }
    int numLoops() const
    {
	return 1;
    }
private:
    const string index;
    const string begin;
    const string end;
    const string inc;
};

class ForLoopDec : public Instruction
{
public:
    ForLoopDec(const string& index_, const string& begin_, const string& end_, const string& dec_) : index(index_), begin(begin_), end(end_), dec(dec_) {}
    string serialize() const
    {
        return " for ( int "+index+"="+begin+"; "+index+">"+end+"; "+index+"-="+dec+" ) { ";
    }
    int numLoops() const
    {
        return 1;
    }
private:
    const string index;
    const string begin;
    const string end;
    const string dec;
};

class ForLoopStatement : public Instruction
{
public:
    ForLoopStatement(const string& begin_, const string& end_, const string& inc_) : begin(begin_), end(end_), inc(inc_) {}
    string serialize() const
    {
	return " for ( "+begin+"; "+end+"; "+inc+ ") { ";
    }
    int numLoops() const
    {
        return 1;
    }
private:
    const string begin;
    const string end;
    const string inc;
};

class InlineCode : public Instruction
{
public:
    InlineCode(const string& code_) : code(code_) {}
    string serialize() const
    {
        return code;
    }
private:
    const string code;
};

enum MemRefType { LOAD_ALIGNED_VEC, LOAD_UNALIGNED_VEC, LOAD_MASKED_VEC, STORE_VEC, STREAM_VEC, STORE_MASKED_VEC, LOAD_NONVEC, L1_PREFETCH, NTA_PREFETCH, L2_PREFETCH, L1_EVICT, L2_EVICT, GATHER_VEC, SCATTER_VEC, GATHER_PREFETCH_L1, GATHER_PREFETCH_L2, GATHER_PREFETCH_NTA };

class MemRefInstruction : public Instruction
{
public:
    // Override virtual
    virtual bool hasAddress() const
    {
        return true;
    }
    virtual const Address* getAddress() const = 0;
    virtual MemRefType getType() const = 0;
    virtual bool hasGSAddress() const
    {
        return false;
    }
};

class DeclareFVec : public Instruction
{
public:
    DeclareFVec( const FVec &v_ ) : v(v_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }

private:
    const FVec v;
};

class InitFVec : public Instruction
{
public:
    InitFVec( const FVec &v_ ) : v(v_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return 0;
    }

private:
    const FVec v;
};

class DeclareFVecArray : public Instruction
{
public:
    DeclareFVecArray( const FVec &v_, const int dim_, const int *len_ ) : v(v_), dim(dim_) { 
	for(int i=0; i<dim; ++i) len[i] = len_[i]; 
    }
    string serialize() const
    {
    	//string tmp = v.getType()+" "+v.getName();
    	ostringstream tmp;
	tmp << v.getType() << " " << v.getName();
	for(int i=0; i<dim; ++i) {
	    //std::ostringstream itos;
	    //itos << len[i];
	    //tmp += "["+itos.str()+"]";
	    tmp << "[" << len[i] << "]";
	}
	tmp << ";";
    	return tmp.str();
    };
    int numDeclarations() const
    {
        return 1;
    }

private:
    const FVec v;
    const int dim;
    int len[MAXARRAYDIM];
};

/*
class Declare_MM_PERM : public Instruction
{
public:
    Declare_MM_PERM() {}
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }
};
*/

class DeclarePACKMASK : public Instruction
{
public:
    DeclarePACKMASK() {}
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }
};

/*
class Declare_FVec : public Instruction
{
public:
    Declare_FVec(string &vect_) : vect(vect_) {}
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }
private:
    const string vect;
};
*/
#if 0
class DeclareFVecArray2D : public Instruction
{
public:
    DeclareFVecArray( const FVec &v_, const int len1_, const int len2_ ) : v(v_), len1(len1_), len2(len2_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }

private:
    const FVec v;
    const int len1;
    const int len2;
};

class DeclareFVecArray3D : public Instruction
{
public:
    DeclareFVecArray( const FVec &v_, const int len1_, const int len2_, const int len3_ ) : v(v_), len1(len1_), len2(len2_), len3(len3_) { }   
    string serialize() const;
    int numDeclarations() const
    {
        return 1;
    }

private:
    const FVec v;
    const int len1;
    const int len2;
    const int len3;
};
#endif

class DeclareMask : public Instruction
{
public:
    DeclareMask(string name_, string value_="") : name(name_), value(value_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return (1);
    }
private:
    const string name;
    const string value;
};

class DeclareInt : public Instruction
{
public:
    DeclareInt(string name_, string value_="") : name(name_), value(value_) { }
    string serialize() const
    {
    	ostringstream outbuf;

    	if(value.empty()) {
            outbuf << "int " << name << ";" << endl;
    	}
    	else {
            outbuf << "int " << name << " = " << value << ";" << endl;
    	}

    	return outbuf.str();
    }
    int numDeclarations() const
    {
        return (1);
    }
private:
    const string name;
    const string value;
};

class IntToMask : public Instruction
{
public:
    IntToMask(const string maskname, const string valname) : mask(maskname), value(valname) {}
    string serialize() const;
private:
    const string mask, value;
};

class DeclareOffsets : public Instruction
{
public:
    DeclareOffsets(string pname_, string vname_) : pname(pname_), vname(vname_) { }
    string serialize() const;
    int numDeclarations() const
    {
        return (1);
    }
private:
    const string vname;
    const string pname;
};

class IfAllOneCond : public Instruction
{
public:
    IfAllOneCond(const string& condition_) : condition(condition_) {}
    string serialize() const;
    int numIfs() const
    {
        return 1;
    }
private:
    const string condition;
};

class LoadFVec : public MemRefInstruction
{
public:
    LoadFVec( const FVec& v_, const Address* a_, const string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return LOAD_ALIGNED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;

};

class StoreFVec : public MemRefInstruction
{
public:
    StoreFVec( const FVec& v_, const Address* a_, int isStreaming_) : v(v_), a(a_), isStreaming(isStreaming_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_VEC;
    }

private:
    const FVec v;
    const Address* a;
    int isStreaming;
};

class StoreFVecMask : public MemRefInstruction
{
public:
    StoreFVecMask( const FVec& v_, const Address* a_, string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_VEC;
    }

private:
    const FVec v;
    const Address* a;
    const string mask;
};

class GatherFVec : public MemRefInstruction
{
public:
    GatherFVec( const FVec& v_, const GatherAddress* a_, const string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_VEC;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const FVec v;
    const GatherAddress* a;
    const string mask;
};

class ScatterFVec : public MemRefInstruction
{
public:
    ScatterFVec( const FVec& v_, const GatherAddress* a_) : v(v_), a(a_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return SCATTER_VEC;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const FVec v;
    const GatherAddress* a;
};

class LoadBroadcast : public MemRefInstruction
{
public:
    LoadBroadcast( const FVec& v_, const Address* a_) : v(v_), a(a_) {}
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return LOAD_NONVEC;
    }
private:
    const FVec v;
    const Address* a;
};


class PrefetchL1 : public MemRefInstruction
{
public:
    PrefetchL1( const Address* a_, int type = 0);

    string serialize() const
    {
        ostringstream stream;
        stream << " _mm_prefetch((const char *)( " << a->serialize() << " ), " << hint << ");" << endl;
        return stream.str();
    }

    MemRefType getType() const
    {
        return L1_PREFETCH;
    }
    const Address* getAddress() const
    {
        return a;
    }

private:
    const Address* a;
    string hint;
};

class PrefetchL2 : public MemRefInstruction
{
public:
    PrefetchL2( const Address* a_, int type = 0);
    string serialize() const
    {
        ostringstream stream;
        stream << " _mm_prefetch((const char *)( " << a->serialize() << " ), " << hint << ");" << endl;
        return stream.str();
    }

    MemRefType getType() const
    {
        return L2_PREFETCH;
    }
    const Address* getAddress() const
    {
        return a;
    }

private:
    const Address* a;
    string hint;
};

class GatherPrefetchL1 : public MemRefInstruction
{
public:
    GatherPrefetchL1( const GatherAddress* a_, int type = 0);
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_PREFETCH_L1;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const GatherAddress* a;
    string hint;
};

class GatherPrefetchL2 : public MemRefInstruction
{
public:
    GatherPrefetchL2( const GatherAddress* a_, int type = 0);
    string serialize() const;
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return GATHER_PREFETCH_L2;
    }
    virtual bool hasGSAddress() const
    {
        return true;
    }
private:
    const GatherAddress* a;
    string hint;
};


// Arithmetic Instructions
class SetZero : public Instruction
{
public:
    SetZero( const FVec& ret_) : ret(ret_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
};

class Set1Const : public Instruction {
public:
    Set1Const( const FVec& ret_, const double val_) : ret(ret_), val(val_) {}
    string serialize() const;
    int numArithmeticInst() const {
        return 1;
    }
private:
    const FVec ret;
    const double val;
};

class Mul : public Instruction
{
public:
    Mul( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {};
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class SMul : public Instruction
{
public:
    SMul( const FVec& ret_, const FVec& a_, const string &scalar_, const string& mask_) : ret(ret_), a(a_), scalar(scalar_), mask(mask_) {};
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const string scalar;
    const string mask;
};

class NSMul : public Instruction
{
public:
    NSMul( const FVec& ret_, const FVec& a_, const string &scalar_, const string& mask_) : ret(ret_), a(a_), scalar(scalar_), mask(mask_) {};
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const string scalar;
    const string mask;
};

class FnMAdd : public Instruction
{
public:
    FnMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const FVec c;
    const string mask;
};

class FSnMAdd : public Instruction
{
public:
    FSnMAdd( const FVec& ret_, const FVec& a_, const string& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const string b;
    const FVec c;
    const string mask;
};

class FMSnMAdd : public Instruction
{
public:
    FMSnMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const string& c_, const FVec& d_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), d(d_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string c;
    const FVec d;
    const string mask;
};
class FMAdd : public Instruction
{
public:
    FMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const FVec c;
    const string mask;
};

class FSMAdd : public Instruction
{
public:
    FSMAdd( const FVec& ret_, const FVec& a_, const string& b_, const FVec& c_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const string b;
    const FVec c;
    const string mask;
};

class FMSMAdd : public Instruction
{
public:
    FMSMAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const string& c_, const FVec& d_, const string& mask_) : ret(ret_), a(a_), b(b_), c(c_), d(d_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string c;
    const FVec d;
    const string mask;
};

class Add : public Instruction
{
public:
    Add( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class NAdd : public Instruction
{
public:
    NAdd( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class Sub : public Instruction
{
public:
    Sub( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class Div : public Instruction
{
public:
    Div( const FVec& ret_, const FVec& a_, const FVec& b_, const string& mask_) : ret(ret_), a(a_), b(b_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    const string mask;
};

class Sqrt : public Instruction
{
public:
    Sqrt( const FVec& ret_, const FVec& a_, const string& mask_) : ret(ret_), a(a_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 1;
    }
private:
    const FVec ret;
    const FVec a;
    const string mask;
};

class MovFVec : public Instruction
{
public:
    MovFVec( const FVec& ret_, const FVec& a_, const string mask_) : ret(ret_), a(a_), mask(mask_) {}
    string serialize() const;
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    const string mask;
};


void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen);
void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen);

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward);
void storeSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward);
void loadSplit3SOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, const Address *a3, int soanum, int soalen, int forward);
void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask=0xFFFF);
void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask=0xFFFF);
inline void packFVec(InstVector& ivector, const FVec& ret, Address *a, int mask, int possibleMask)
{
    stringstream msk_str;
    msk_str << "0x" << hex << mask;
    packFVec(ivector, ret, a, msk_str.str(), possibleMask);
}


void gatherFVec(InstVector& ivector, const FVec& ret, GatherAddress *a, string mask);
void scatterFVec(InstVector& ivector, const FVec& ret, GatherAddress *a);
void gatherPrefetchL1(InstVector& ivector, GatherAddress *a, int type = 0);
void gatherPrefetchL2(InstVector& ivector, GatherAddress *a, int type = 0);

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen);

void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir);
void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir);
void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir);
void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir);

int packXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir);
void packXYZTFVec(InstVector& ivector, const string& r, const Address*lAddr, const Address*rAddr, string dir, string srz);
int packXYZTFVec(InstVector& ivector, const FVec r[2], const Address*rAddr, int dir);
void packXYZTFVec(InstVector& ivector, const string& r, const Address*rAddr, string dir, string srz);
int unpackXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir);
void unpackXYZTFVec(InstVector& ivector, const string &r, const Address*lAddr, const Address*rAddr, string dir, string srz);
int unpackXYZTFVec(InstVector& ivector, const FVec r[2], const Address*rAddr, int dir);
void unpackXYZTFVec(InstVector& ivector, const string &r, const Address*rAddr, string dir, string srz);

inline void movFVec(InstVector& ivector, const FVec& ret, const FVec& a, string mask)
{
    ivector.push_back(new MovFVec(ret, a, mask));
}

inline void beginScope(InstVector& ivector)
{
    ivector.push_back(new BeginScope());
}

inline void endScope(InstVector& ivector)
{
    ivector.push_back(new EndScope());
}

inline void ifStatement(InstVector& ivector, string condition)
{
    ivector.push_back( new IfStringCond(condition));
}

inline void elseStatement(InstVector& ivector)
{
    ivector.push_back( new ElseStatement());
}

inline void elseIfStatement(InstVector& ivector, string condition)
{
    ivector.push_back( new ElseIfStatement(condition));
}

inline void forLoopInc(InstVector& ivector, string index, string begin, string end, string inc="1")
{
    ivector.push_back( new ForLoopInc(index, begin, end, inc));
}

inline void forLoopDec(InstVector& ivector, string index, string begin, string end, string dec="1")
{
    ivector.push_back( new ForLoopDec(index, begin, end, dec));
}

inline void forLoopStatement(InstVector& ivector, string begin, string end, string inc)
{
    ivector.push_back( new ForLoopStatement(begin, end, inc));
}

inline void ifAllOneStatement(InstVector& ivector, string condition)
{
    ivector.push_back( new IfAllOneCond(condition));
}

inline void inlineCode(InstVector& ivector, string code)
{
    ivector.push_back( new InlineCode(code));
}

inline void declareFVecFromFVec(InstVector& ivector, const FVec& v)
{
    ivector.push_back(new DeclareFVec(v));
}

inline FVec declareFVec(InstVector& ivector, const std::string name)
{
    FVec tmp(name);
    ivector.push_back(new DeclareFVec( tmp ));
    return tmp;
}

inline void initFVec(InstVector& ivector, const FVec& ret)
{
    ivector.push_back(new InitFVec(ret));
}

inline void declareInt(InstVector& ivector, const std::string name, const std::string value="")
{
    ivector.push_back(new DeclareInt(name, value));
}

inline void initInt(InstVector& ivector, const string name, const string value)
{
    ivector.push_back(new IntToMask(name, value));
}

inline void declareFVecArray(InstVector& ivector, const std::string name, int dim, int *len)
{
    FVec tmp(name);
    ivector.push_back(new DeclareFVecArray(tmp, dim, len));
}

inline void declareFVecArrayFromFVec(InstVector& ivector, const FVec& v, int dim, int *len)
{
    ivector.push_back(new DeclareFVecArray(v, dim, len));
}

/*
inline void declare_MM_PERM(InstVector& ivector)
{
    ivector.push_back(new Declare_MM_PERM());
}
*/

inline void declarePACKMASK(InstVector& ivector)
{
    ivector.push_back(new DeclarePACKMASK());
}

/*
inline void declare_FVec(InstVector& ivector, string vect)
{
    ivector.push_back(new Declare_FVec(vect));
}
*/

inline void setZero(InstVector& ivector, const FVec& ret)
{
    ivector.push_back(new SetZero(ret));
}

inline void set1Const(InstVector& ivector, const FVec& ret, const double val)
{
    ivector.push_back(new Set1Const(ret, val));
}

inline void mulFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Mul(ret, a, b, mask));
}

inline void smulFVec(InstVector& ivector, const FVec& ret, const FVec& a, string s, string mask = "")
{
    ivector.push_back(new SMul(ret, a, s, mask));
}

inline void nsmulFVec(InstVector& ivector, const FVec& ret, const FVec& a, string s, string mask = "")
{
    ivector.push_back(new NSMul(ret, a, s, mask));
}

inline void addFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Add(ret, a, b, mask));
}

inline void naddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new NAdd(ret, a, b, mask));
}

inline void subFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Sub(ret, a, b, mask));
}

inline void fmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const FVec& c, string mask = "")
{
    ivector.push_back(new FMAdd(ret, a, b, c, mask));
}

inline void fsmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, string s, const FVec& c, string mask = "")
{
    ivector.push_back(new FSMAdd(ret, a, s, c, mask));
}

inline void fmsmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string s, const FVec& c, string mask = "")
{
    ivector.push_back(new FMSMAdd(ret, a, b, s, c, mask));
}

inline void fnmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const FVec& c, string mask = "")
{
    ivector.push_back(new FnMAdd(ret, a, b, c, mask));
}

inline void fsnmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, string s, const FVec& c, string mask = "")
{
    ivector.push_back(new FSnMAdd(ret, a, s, c, mask));
}

inline void fmsnmaddFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string s, const FVec& c, string mask = "")
{
    ivector.push_back(new FMSnMAdd(ret, a, b, s, c, mask));
}

inline void divFVec(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, string mask = "")
{
    ivector.push_back(new Div(ret, a, b, mask));
}

inline void sqrtFVec(InstVector& ivector, const FVec& ret, const FVec& a, string mask = "")
{
    ivector.push_back(new Sqrt(ret, a, mask));
}

inline void loadFVec(InstVector& ivector, const FVec& ret, const Address *a, string mask)
{
    ivector.push_back( new LoadFVec(ret, a, mask));
}

inline void storeFVec(InstVector& ivector, const FVec& ret, const Address *a, string mask)
{
    ivector.push_back( new StoreFVecMask(ret, a, mask));
}

inline void storeFVec(InstVector& ivector, const FVec& ret, const Address *a, int isStreaming)
{
    ivector.push_back( new StoreFVec(ret, a, isStreaming));
}

inline void loadBroadcastScalar(InstVector& ivector, const FVec& ret, string scalar_name, int type = 0)
{
    ivector.push_back( new LoadBroadcast(ret, new AddressOfScalar(scalar_name, type)));
}

inline void prefetchL1(InstVector& ivector, const Address *a, int type = 0)
{
    ivector.push_back( new PrefetchL1(a, type));
}

inline void prefetchL2(InstVector& ivector, const Address *a, int type = 0)
{
    ivector.push_back( new PrefetchL2(a, type));
}

inline void declareMask(InstVector& ivector, const string name, const string value="")
{
    ivector.push_back(new DeclareMask(name, value));
}

inline void intToMask(InstVector& ivector, const string maskname, const string intname)
{
    ivector.push_back(new IntToMask(maskname, intname));
}

inline void intToMask(InstVector& ivector, const string maskname, const int msk) {
    stringstream msk_str;
    msk_str << "0x" << hex << msk;
    ivector.push_back(new IntToMask(maskname, msk_str.str()));
}

#endif
