#include <stdlib.h>
#include <stdio.h>
#include "instructions.h"

#if PRECISION == 2 && VECLEN == 8
#pragma message "Using double Precision"

#define FVECTYPE "__m512d"

const string fullMask("0xFF");

FVec::FVec(const string& name_) : name(name_), type(FVECTYPE) {}

string DeclareFVec::serialize() const
{
#if 0
    return v.getType()+" "+v.getName()+ " = _mm512_undefined_pd(); ";
#else
    return v.getType()+" "+v.getName()+ " = _mm512_setzero_pd(); ";
#endif
}

string InitFVec::serialize() const
{
#if 1
    return v.getName()+ " = _mm512_undefined_pd(); ";
#else
    return v.getName()+ " = _mm512_setzero_pd(); ";
#endif
}

/*
string Declare_MM_PERM::serialize() const
{
#ifndef AVX512
    return "_MM_PERM_ENUM MM_PERM[8] = { _MM_SWIZ_REG_CDAB, _MM_SWIZ_REG_CDAB, _MM_SWIZ_REG_BADC, _MM_SWIZ_REG_BADC, _MM_PERM_BADC, _MM_PERM_BADC, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00 };";
#else
    return "_MM_PERM_ENUM MM_PERM[8] = { (_MM_PERM_ENUM)0x55, (_MM_PERM_ENUM)0x55, (_MM_PERM_ENUM)0xB1, (_MM_PERM_ENUM)0xB1, (_MM_PERM_ENUM)0x4E, (_MM_PERM_ENUM)0x4E, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00 };";
#endif
}
*/

string DeclarePACKMASK::serialize() const
{
    return " __mmask8 PACKMASK[2][4] = {{0x55, 0x33, 0x0F, 0xFF}, {0xAA, 0xCC, 0xF0, 0xFF}};";
}

string DeclareMask::serialize() const
{
    ostringstream outbuf;

    if(value.empty()) {
        outbuf << "__mmask8 " << name << ";" << endl;
    }
    else {
        outbuf << "__mmask8 " << name << " = " << value << ";" << endl;
    }

    return outbuf.str();
}

string IntToMask::serialize() const
{
    ostringstream outbuf;
    outbuf << mask << " = " << value << ";" << endl;
    return outbuf.str();
}

string DeclareOffsets::serialize() const
{
    ostringstream outbuf;

    outbuf << "__m512i " << vname
           << " = _mm512_mask_load_epi32(_mm512_setzero_epi32(), "
           << fullMask << ", " << pname << ");" << endl;
    return outbuf.str();
}

string IfAllOneCond::serialize() const
{
    return " if ((" + condition + " & " + fullMask + ") == " + fullMask + ") { ";
}

string LoadFVec::serialize() const
{
    std::ostringstream buf;
    string lmask = mask;

    if(mask.empty()) {
        lmask = fullMask;
    }


    if(!a->isHalfType()) {
        buf << v.getName() << " = _mm512_mask_load_pd(" << v.getName() << ", " << lmask << ", "  << a->serialize() << ");" <<endl;
    }
    else {
        buf << v.getName() << " = _mm512_mask_cvtpslo_pd(" << v.getName() << ", " << lmask << ", "  <<
            "_mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), " << fullMask << ", " << a->serialize() << "));" <<endl;
    }

    return buf.str();
}

string StoreFVec::serialize() const
{
    ostringstream buf;
    int streaming = isStreaming;

    if(a->isHalfType()) {
        streaming = 0;
    }

    if(streaming) {
#ifdef AVX512
        buf << "_mm512_stream_pd((void*)" << a->serialize() << "," << v.getName() <<  ");" <<endl;
#else
        buf << "_mm512_storenrngo_pd((void*)" << a->serialize() << "," << v.getName() <<  ");" <<endl;
#endif
    }
    else {
        if(!a->isHalfType()) {
            buf << "_mm512_store_pd(" << a->serialize() << "," << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm512_mask_packstorelo_ps(" << a->serialize() << ", " << fullMask <<  ", _mm512_cvtpd_pslo(" << v.getName() <<  "));" <<endl;
        }
    }

    return buf.str();
}

string StoreFVecMask::serialize() const
{
    ostringstream buf;
        if(!a->isHalfType()) {
            buf << "_mm512_mask_store_pd(" << a->serialize() << ", " << mask << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm512_mask_packstorelo_ps(" << a->serialize() << ", " << mask <<  ", _mm512_cvtpd_pslo(" << v.getName() <<  "));" <<endl;
        }
    return buf.str();
}

string GatherFVec::serialize() const
{
    std::ostringstream buf;
    string lmask = mask;

    if(mask.empty()) {
        lmask = fullMask;
    }


    if(!a->isHalfType()) {
        buf << v.getName() << " = _mm512_mask_i32logather_pd(" << v.getName()
            << ", " << lmask << ", " << a->getOffsets() << ", (void*)"
            << a->getBase()  << ", _MM_SCALE_8);" <<endl;
    }
    else {
        buf << v.getName() << " = _mm512_mask_cvtpslo_pd(" << v.getName()
            << ", " << lmask << ", "
            << "_mm512_mask_i32gather_ps(_mm512_undefined_ps(), " << ", "
            << lmask << ", " <<  a->getOffsets() << ", (void*)"
            << a->getBase()  << ", _MM_SCALE_4));" <<endl;
    }

    return buf.str();

}

string ScatterFVec::serialize() const
{
    std::ostringstream buf;

    string downConv = "_MM_DOWNCONV_PS_NONE";
    string scale = "_MM_SCALE_4";

    if(a->isHalfType()) {
        downConv = "_MM_DOWNCONV_PS_FLOAT16";
        scale = "_MM_SCALE_2";
    }


    if(!a->isHalfType()) {
        buf << "_mm512_i32loscatter_pd((void*)" << a->getBase() << ", "
            << a->getOffsets() << ", " << v.getName() << ", _MM_SCALE_4);" <<endl;
    }
    else {
        buf << "_mm512_mask_i32scatter_ps((void*)" << a->getBase() << ", "
            << fullMask << ", " << a->getOffsets() << ", _mm512_cvtpd_pslo("
            << v.getName() << "), _MM_SCALE_4);" <<endl;
    }

    return buf.str();
}

string LoadBroadcast::serialize() const
{
    std::ostringstream buf;


    if(!a->isHalfType()) {
        buf << v.getName() << " = _mm512_extload_pd(" << a->serialize()
            << ", _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);" << endl;
    }
    else {
        buf << v.getName() << " = _mm512_cvtpslo_pd(_mm512_extload_ps("
            << a->serialize() << ", _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE));"
            << endl;
    }

    return buf.str();
}

class LoadUnpackFVec : public MemRefInstruction
{
public:
    LoadUnpackFVec( const FVec& v_, const Address* a_, const string mask_) : v(v_), a(a_), mask(mask_) {}

    string serialize() const
    {
        std::ostringstream buf;
        string lmask = mask;

        if(mask.empty()) {
            lmask = fullMask;
        }

#ifndef AVX512

        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm512_mask_loadunpacklo_pd(" << v.getName()
                << ", " << lmask << ", " << a->serialize() << ");" << endl;
        }
        else {
            buf << v.getName() << " = _mm512_mask_cvtpslo_pd(" << v.getName()
                << ", " << lmask << ", " <<"_mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), "
                << lmask << ", " << a->serialize() << "));" << endl;
        }

#else

        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm512_mask_expandloadu_pd(" << v.getName()
                << ", " << lmask << ", " << a->serialize() << ");" << endl;
        }
        else {
            buf << v.getName() << " = _mm512_mask_cvtpslo_pd(" << v.getName()
                << ", " << lmask << ", " <<"_mm512_mask_expandloadu_ps(_mm512_undefined_ps(), "
                << lmask << ", " << a->serialize() << "));" << endl;
        }

#endif
        return buf.str();
    }
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return LOAD_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
};

class PackStoreFVec : public MemRefInstruction
{
public:
    PackStoreFVec( const FVec& v_, const Address* a_, string mask_) : v(v_), a(a_), mask(mask_) {}

    string serialize() const
    {
        std::ostringstream buf;
        string lmask= mask;

        if(mask.empty()) {
            lmask = fullMask;
        }

#ifndef AVX512

        if(!a->isHalfType()) {
            buf << "_mm512_mask_packstorelo_pd((void*)" << a->serialize() << ", "
                << lmask << ", " << v.getName() << ");" << endl;
        }
        else {
            buf << "_mm512_mask_packstorelo_ps((void*)" << a->serialize() << ", "
                << lmask << ", _mm512_cvtpd_pslo(" << v.getName() << "));" << endl;
        }

#else

        if(!a->isHalfType()) {
            buf << "_mm512_mask_compressstoreu_pd(" << a->serialize() << ", "
                << lmask << ", " << v.getName() << ");" << endl;
        }
        else {
            buf << "_mm512_mask_compressstoreu_ps(" << a->serialize() << ", "
                << lmask << ", _mm512_cvtpd_pslo(" << v.getName() << "));" << endl;
        }

#endif
        return buf.str();

    }
    const Address* getAddress() const
    {
        return a;
    }
    MemRefType getType() const
    {
        return STORE_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
};

PrefetchL1::PrefetchL1( const Address* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch (type) {
    case 0:
        hint = "_MM_HINT_T0";
        break;

    case 1:
        hint = "_MM_HINT_NTA";
        break;

    case 2:
        hint = "_MM_HINT_ET0";
        break;

    case 3:
        hint = "_MM_HINT_ENTA";
        break;
    }
}

PrefetchL2::PrefetchL2( const Address* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch(type) {
    case 0:
        hint = "_MM_HINT_T1";
        break;

    case 1:
        hint = "_MM_HINT_T2";
        break;

    case 2:
        hint = "_MM_HINT_ET1";
        break;

    case 3:
        hint = "_MM_HINT_ET2";
        break;
    }
}

GatherPrefetchL1::GatherPrefetchL1( const GatherAddress* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch (type) {
    case 0:
        hint = "_MM_HINT_T0";
        break;

    case 1:
        hint = "_MM_HINT_NTA";
        break;

    case 2:
        hint = "_MM_HINT_ET0";
        break;

    case 3:
        hint = "_MM_HINT_ENTA";
        break;
    }
}

string GatherPrefetchL1::serialize() const
{
    std::ostringstream buf;
    string scale = "_MM_SCALE_8";

    if(a->isHalfType()) {
        scale = "_MM_SCALE_4";
    }

    buf << "_mm512_mask_prefetch_i32gather_ps(" << a->getOffsets() << ", "
        << fullMask << ", (void*)" << a->getBase() << ", " << scale << ", "
        << hint << ");" <<endl;

    return buf.str();
}


GatherPrefetchL2::GatherPrefetchL2( const GatherAddress* a_, int type) : a(a_)
{
    // Type: 0 - none, 1 - NT, 2 - Ex, 3 - NT+Ex
    switch(type) {
    case 0:
        hint = "_MM_HINT_T1";
        break;

    case 1:
        hint = "_MM_HINT_T2";
        break;

    case 2:
        hint = "_MM_HINT_ET1";
        break;

    case 3:
        hint = "_MM_HINT_ET2";
        break;
    }
}

string GatherPrefetchL2::serialize() const
{
    std::ostringstream buf;
    string scale = "_MM_SCALE_4";

    if(a->isHalfType()) {
        scale = "_MM_SCALE_2";
    }

    buf << "_mm512_mask_prefetch_i32gather_ps(" << a->getOffsets() << ", "
        << fullMask << ", (void*)" << a->getBase() << ", " << scale << ", "
        << hint << ");" <<endl;

    return buf.str();

}

string SetZero::serialize() const
{
    return  ret.getName()+" = _mm512_setzero_pd(); ";
}

string Set1Const::serialize() const {
    std::ostringstream buf;
    buf << ret.getName() << " = _mm512_set1_pd(" << val << "); " << endl;
    return  buf.str();
}

string Mul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_mul_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+ " = _mm512_mask_mul_pd( " + ret.getName() + "," + mask + "," + a.getName() + " , " + b.getName()+" );" ;
    }
}

string SMul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_mul_pd( "+a.getName()+" , _mm512_set1_pd( "+ scalar +" ));" ;
    }
    else {
	return  ret.getName()+ " = _mm512_mask_mul_pd( " + ret.getName() + "," + mask + "," + a.getName() + " , _mm512_set1_pd( " + scalar + " ));" ;
    }
}

string NSMul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_mul_pd( "+a.getName()+" , _mm512_set1_pd( -("+ scalar +") ));" ;
    }
    else {
        return  ret.getName()+ " = _mm512_mask_mul_pd( " + ret.getName() + "," + mask + "," + a.getName() + " , _mm512_set1_pd( -(" + scalar + ") ));" ;
    }
}

string FnMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fnmadd_pd("+a.getName()+", "+b.getName()+", "+c.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fnmadd_pd(" +a.getName()+", "+b.getName()+", "+c.getName()+"));" ;
    }
}

string FSnMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fnmadd_pd("+a.getName()+", _mm512_set1_pd( "+b+"), "+c.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fnmadd_pd(" +a.getName()+", _mm512_set1_pd( "+b+"), "+c.getName()+"));" ; 
    }
}

string FMSnMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fnmadd_pd("+a.getName()+", _mm512_mul_pd( "+b.getName()+", _mm512_set1_pd( "+c+")), "+d.getName()+");" ;
    }
    else {
	return  ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fnmadd_pd("+a.getName()+", _mm512_mul_pd( "+b.getName()+", _mm512_set1_pd("+c+")), "+d.getName()+"));" ;
    }
}

string FMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fmadd_pd("+a.getName()+", "+b.getName()+", "+c.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fmadd_pd(" +a.getName()+", "+b.getName()+", "+c.getName()+"));" ;
    }
}

string FSMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fmadd_pd("+a.getName()+", _mm512_set1_pd( "+b+"), "+c.getName()+");" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fmadd_pd(" +a.getName()+", _mm512_set1_pd( "+b+"), "+c.getName()+"));" ;
    }
}

string FMSMAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_fmadd_pd("+a.getName()+", _mm512_mul_pd("+b.getName()+", _mm512_set1_pd( "+c+" )), "+d.getName()+");" ;
    }
    else {
	return ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", _mm512_fmadd_pd("+a.getName()+", _mm512_mul_pd("+b.getName()+", _mm512_set1_pd( "+c+" )), "+d.getName()+"));" ;
    }
}

string Add::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_add_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_add_pd( " + ret.getName() + "," + mask + "," +a.getName()+" , "+b.getName()+" );" ;
    }
}

string NAdd::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_sub_pd( _mm512_set1_pd(0.0), _mm512_add_pd( "+a.getName()+" , "+b.getName()+" ));" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_sub_pd( " + ret.getName() + "," + mask + ", _mm512_set1_pd(0.0), _mm512_add_pd( " +a.getName()+" , "+b.getName()+" ));" ;
    }
}

string Sub::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_sub_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_sub_pd( " + ret.getName() + "," + mask + "," +a.getName()+" , "+b.getName()+" );" ;
    }
}

string Div::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_div_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_div_pd( " + ret.getName() + "," + mask + "," +a.getName()+" , "+b.getName()+" );" ;
    }
}

string Sqrt::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm512_sqrt_pd( "+a.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm512_mask_sqrt_pd( " + ret.getName() + "," + mask + "," +a.getName()+" );" ;
    }
}

string MovFVec::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = " + a.getName()+";" ;
    }
    else {
        return ret.getName()+" = _mm512_mask_mov_pd(" + ret.getName() + ", " + mask + ", " + a.getName() + ");";
    }
}

class Perm64x2 : public Instruction
{
public:
    Perm64x2(const FVec& ret_, const string mk_, const FVec& a_, const int imm_) : ret(ret_), mk(mk_), a(a_), imm(imm_) {}

    string serialize() const
    {
        ostringstream stream;
        stream << ret.getName()
               << " = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps("
               << ret.getName() << "), " << mk
               << ", _mm512_castpd_ps(" << a.getName()
               << "), (_MM_PERM_ENUM)" << imm << "));" ;

        return stream.str();
    }

    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    const string mk;
    const int imm;
};

class PermuteXYZTFVec : public Instruction
{
public:
    PermuteXYZTFVec(const FVec& ret_, const FVec& a_, int dir_) : ret(ret_), a(a_), dir(dir_/2) {}
    string serialize() const
    {
        ostringstream stream;
#ifndef AVX512
		if(dir < 2) {
			string imm = (dir == 0 ? "_MM_SWIZ_REG_CDAB" : "_MM_SWIZ_REG_BADC");
	        stream << ret.getName() << " = _mm512_swizzle_pd(" << a.getName() << ", "  << imm << ");";
		}
		else if(dir == 2) {
			string imm = "_MM_PERM_BADC";
			stream << ret.getName() << " = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(" << a.getName() << "), " << imm << "));";
		}
		else {
			// nothing needs t be done for dir 3
		}
#else // AVX512 defined
		if(dir == 0) {
			string imm = "0x55";
	        stream << ret.getName() << " = _mm512_permute_pd(" << a.getName() << ", "  << imm << ");";
		} 
		else if(dir < 3) {
			string imm = (dir == 1 ? "0xB1" : "0x4E");
			stream << ret.getName() << " = _mm512_shuffle_f64x2(" << a.getName() << ", "  << a.getName() << ", " << imm << ");" ;
		}
		else {
			// nothing needs t be done for dir 3
		}
#endif // AVX512
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    int dir;
};

class PermuteXYZTFVecDir : public Instruction
{
public:
    PermuteXYZTFVecDir(const FVec& ret_, const FVec& a_, const string& dir_) : ret(ret_), a(a_), dir(dir_) {}
    string serialize() const
    {
        ostringstream stream;
#ifndef AVX512
                stream << "if(((" << dir << ")>>1) == 0) " << ret.getName() << " = _mm512_swizzle_pd(" << a.getName() << ", "  << " _MM_SWIZ_REG_CDAB);" << endl
			<< "else if(((" << dir << ")>>1) == 1) " << ret.getName() << " = _mm512_swizzle_pd(" << a.getName() << ", "  << " _MM_SWIZ_REG_BADC);" << endl
                	<< "else if(((" << dir << ")>>1)== 2) " << ret.getName() << " = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(" << a.getName() << "), _MM_PERM_BADC));";
#else
                stream << "if(((" << dir << ")>>1) == 0) " << ret.getName() << " = _mm512_permute_pd(" << a.getName() << ", (_MM_PERM_ENUM)0x55);" << endl
			<< "else if(((" << dir << ")>>1) == 1) " << ret.getName() << " = _mm512_shuffle_f64x2(" << a.getName() << ", "  << a.getName() << ", (_MM_PERM_ENUM)0xB1);" << endl
                        << "else if(((" << dir << ")>>1) == 2) " << ret.getName() << " = _mm512_shuffle_f64x2(" << a.getName() << ", "  << a.getName() << ", (_MM_PERM_ENUM)0x4E);" ;
#endif
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    const string dir;
};

class APermuteXYZTFVec : public Instruction
{
public:
    APermuteXYZTFVec(const FVec& ret_, const FVec& a_, int dir_) : ret(ret_), a(a_), dir(dir_/2) {}
    string serialize() const
    {
        ostringstream stream;
#ifndef AVX512
                if(dir < 2) {
                        string imm = (dir == 0 ? "_MM_SWIZ_REG_CDAB" : "_MM_SWIZ_REG_BADC");
                stream << ret.getName() << " = _mm512_swizzle_pd(" << a.getName() << ", "  << imm << ");";
                }
                else if(dir == 2) {
                        string imm = "_MM_PERM_BADC";
                        stream << ret.getName() << " = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(" << a.getName() << "), " << imm << "));";
                }
                else {
			// nothing needs t be done for dir 3
                }
#else // AVX512 defined
                if(dir == 0) {
                        string imm = "0x55";
                stream << ret.getName() << " = _mm512_permute_pd(" << a.getName() << ", "  << imm << ");";
                }
                else if(dir < 3) {
                        string imm = (dir == 1 ? "0xB1" : "0x4E");
                        stream << ret.getName() << " = _mm512_shuffle_f64x2(" << a.getName() << ", "  << a.getName() << ", " << imm << ");" ;
                }
                else {
                        // nothing needs t be done for dir 3
                }
#endif 
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    int dir;
};

void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    int mskbits = (((1 << soalen)-1) << (soanum*soalen));
    stringstream mk;
    mk << "0x" << hex << mskbits;
    string localmask = mk.str();
    ivector.push_back( new LoadUnpackFVec(ret, a, localmask));
}

void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    int mskbits = (((1 << soalen)-1) << (soanum*soalen));
    stringstream mk;
    mk << "0x" << hex << mskbits;
    string localmask = mk.str();
    ivector.push_back( new PackStoreFVec(ret, a, localmask));
}

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward)
{
    int mskbits1, mskbits2;

    if(forward) {
        mskbits1 = (((1 << (soalen-1))-1) << (soanum*soalen));
        mskbits2 = (1 << ((soanum+1)*soalen-1));
    }
    else {
        mskbits1 = (1 << (soanum*soalen));
        mskbits2 = (((1 << (soalen-1))-1) << (soanum*soalen+1));
    }

    stringstream mk1;
    mk1 << "0x" << hex << mskbits1;
    string localmask1 = mk1.str();
    stringstream mk2;
    mk2 << "0x" << hex << mskbits2;
    string localmask2 = mk2.str();
    ivector.push_back( new LoadUnpackFVec(ret, a1, localmask1));
    ivector.push_back( new LoadUnpackFVec(ret, a2, localmask2));
}

void storeSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward)
{
    int mskbits1, mskbits2;
    if(forward) {
        mskbits1 = (((1 << (soalen-1))-1) << (soanum*soalen));
        mskbits2 = (1 << ((soanum+1)*soalen-1));
    }
    else {
        mskbits1 = (1 << (soanum*soalen));
        mskbits2 = (((1 << (soalen-1))-1) << (soanum*soalen+1));
    }

    stringstream mk1;
    mk1 << "0x" << hex << mskbits1;
    string localmask1 = mk1.str();
    stringstream mk2;
    mk2 << "0x" << hex << mskbits2;
    string localmask2 = mk2.str();
    ivector.push_back( new PackStoreFVec(ret, a1, localmask1));
    ivector.push_back( new PackStoreFVec(ret, a2, localmask2));
}

void loadSplit3SOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, const Address *a3, int soanum, int soalen, int forward)
{
    int mskbits1, mskbits2, mskbits3;
    if(forward) {
        mskbits1 = (((1 << (soalen-2))-1) << (soanum*soalen));
        mskbits2 = (1 << ((soanum+1)*soalen-2));
        mskbits3 = (1 << ((soanum+1)*soalen-1));
    }
    else {
        mskbits1 = (1 << (soanum*soalen));
        mskbits2 = (1 << (soanum*soalen+1));
        mskbits3 = (((1 << (soalen-2))-1) << (soanum*soalen+2));
    }

    stringstream mk1;
    mk1 << "0x" << hex << mskbits1;
    string localmask1 = mk1.str();
    stringstream mk2;
    mk2 << "0x" << hex << mskbits2;
    string localmask2 = mk2.str();
    stringstream mk3;
    mk3 << "0x" << hex << mskbits3;
    string localmask3 = mk3.str();
    ivector.push_back( new LoadUnpackFVec(ret, a1, localmask1));
    ivector.push_back( new LoadUnpackFVec(ret, a2, localmask2));
    ivector.push_back( new LoadUnpackFVec(ret, a3, localmask3));
}

void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    ivector.push_back( new LoadUnpackFVec(ret, a, mask));
}

void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    ivector.push_back( new PackStoreFVec(ret, a, mask));
}

void gatherFVec(InstVector& ivector, const FVec& ret, GatherAddress *a, string mask)
{
    ivector.push_back( new GatherFVec(ret, a, mask));
}

void scatterFVec(InstVector& ivector, const FVec& ret, GatherAddress *a)
{
    ivector.push_back( new ScatterFVec(ret, a));
}

void gatherPrefetchL1(InstVector& ivector, GatherAddress *a, int type)
{
    ivector.push_back( new GatherPrefetchL1(a, type));
}

void gatherPrefetchL2(InstVector& ivector, GatherAddress *a, int type)
{
    ivector.push_back( new GatherPrefetchL2(a, type));
}

void fperm64x2(InstVector& ivector, const FVec& ret, const string mk, const FVec& a, const int imm)
{
    ivector.push_back(new Perm64x2(ret, mk, a, imm));
}

void transpose4x4(InstVector& ivector, const FVec r[4], const FVec f[4])
{
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            int mskbits = (((1 << 4)-1) << (j*4));
            stringstream mk;
            mk << "0x" << hex << mskbits;
            string localmask = mk.str();
            fperm64x2(ivector, r[i], localmask, f[j], 0x55*i);
        }
    }
}

void transpose2x2(InstVector& ivector, const FVec r[2], const FVec f[2])
{
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            int mskbits = (((1 << 8)-1) << (j*8));
            stringstream mk;
            mk << "0x" << hex << mskbits;
            string localmask = mk.str();
            fperm64x2(ivector, r[i], localmask, f[j], 0x44+0xAA*i);
        }
    }
}

void transpose1x1(InstVector& ivector, const FVec r[1], const FVec f[1])
{
    movFVec(ivector, r[0], f[0], string(""));
}

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen)
{
    switch (soalen) {
    case 2:
        transpose4x4(ivector, r, f);
        break;

    case 4:
        transpose2x2(ivector, r, f);
        break;

    case 8:
        transpose1x1(ivector, r, f);
        break;

    default:
        printf("SOALEN = %d Not Supported (only SOALEN = 2, 4 & 8 are supported)\n", soalen);
    }
}

void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir)
{
	ivector.push_back(new PermuteXYZTFVec(r, f, dir));
}

void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir)
{
	ivector.push_back(new APermuteXYZTFVec(r, f, dir));
}

void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir)
{
	ivector.push_back(new PermuteXYZTFVecDir(r, f, dir));
}

void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir)
{
	ivector.push_back(new PermuteXYZTFVecDir(r, f, dir));
}

int packXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir) 
{
	int dim = dir/2;
	int fb = dir % 2;
	if(dim < 3) {
		string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
		ivector.push_back( new PackStoreFVec(r[0], rAddr, masks[fb][dim]));
		ivector.push_back( new PackStoreFVec(r[0], lAddr, masks[1-fb][dim]));

		ivector.push_back( new PackStoreFVec(r[1], new AddressImm(rAddr, VECLEN/2), masks[fb][dim]));
		ivector.push_back( new PackStoreFVec(r[1], new AddressImm(lAddr, VECLEN/2), masks[1-fb][dim]));
		return VECLEN;
	}
	else {
		ivector.push_back( new StoreFVec( r[0], rAddr, 1));
		ivector.push_back( new StoreFVec( r[1], new AddressImm(rAddr, VECLEN), 1));
		return 2*VECLEN;
	}
}

void packXYZTFVec(InstVector& ivector, const string &r, const Address*lAddr, const Address*rAddr, string dir, string srz)
{
/*
        int dim = dir/2;
        int fb = dir % 2;
*/
        ifStatement(ivector, "(("+dir+")>>1) < 3");
	{
//                string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
                ivector.push_back( new PackStoreFVec(FVec(r+"[0]"), rAddr, "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new PackStoreFVec(FVec(r+"[0]"), lAddr, "PACKMASK[1-(("+dir+")&1)][("+dir+")>>1]"));

                ivector.push_back( new PackStoreFVec(FVec(r+"[1]"), new AddressImm(rAddr, VECLEN/2), "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new PackStoreFVec(FVec(r+"[1]"), new AddressImm(lAddr, VECLEN/2), "PACKMASK[1-(("+dir+")&1)][("+dir+")>>1]"));
                //return VECLEN;
                ivector.push_back( new IntToMask(srz, "VECLEN") );
        }
        elseStatement(ivector);
	{
                ivector.push_back( new StoreFVec( FVec(r+"[0]"), rAddr, 1));
                ivector.push_back( new StoreFVec( FVec(r+"[1]"), new AddressImm(rAddr, VECLEN), 1));
                //return 2*VECLEN;
                ivector.push_back( new IntToMask(srz, "2*VECLEN") );
        }
	endScope(ivector);
}

int packXYZTFVec(InstVector& ivector, const FVec r[2], const Address*rAddr, int dir)
{
        int dim = dir/2;
        int fb = dir % 2;
        if(dim < 3) {
                string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
                ivector.push_back( new PackStoreFVec(r[0], rAddr, masks[fb][dim]));
		ivector.push_back( new PackStoreFVec(r[1], new AddressImm(rAddr, VECLEN/2), masks[fb][dim]));
		return VECLEN;
	}
	else {
		ivector.push_back( new StoreFVec( r[0], rAddr, 1));
		ivector.push_back( new StoreFVec( r[1], new AddressImm(rAddr, VECLEN), 1));
		return 2*VECLEN;
	}
}

void packXYZTFVec(InstVector& ivector, const string &r, const Address*rAddr, string dir, string srz)
{
        ifStatement(ivector, "(("+dir+")>>1) < 3");
	{
		ivector.push_back( new PackStoreFVec(FVec(r+"[0]"), rAddr, "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new PackStoreFVec(FVec(r+"[1]"), new AddressImm(rAddr, VECLEN/2), "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
		ivector.push_back( new IntToMask(srz, "VECLEN") );
               // return VECLEN;
        }
        elseStatement(ivector);
 	{
                ivector.push_back( new StoreFVec( FVec(r+"[0]"), rAddr, 1));
                ivector.push_back( new StoreFVec( FVec(r+"[1]"), new AddressImm(rAddr, VECLEN), 1));
                //return 2*VECLEN;
                ivector.push_back( new IntToMask(srz, "2*VECLEN") );
        }
	endScope(ivector);
}

int unpackXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir) {
	int dim = dir/2;
	int fb = dir % 2;
	if(dim < 3) {
		string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
		ivector.push_back( new LoadUnpackFVec(r[0], rAddr, masks[fb][dim]));
		ivector.push_back( new LoadUnpackFVec(r[0], lAddr, masks[1-fb][dim]));

		ivector.push_back( new LoadUnpackFVec(r[1], new AddressImm(rAddr, VECLEN/2), masks[fb][dim]));
		ivector.push_back( new LoadUnpackFVec(r[1], new AddressImm(lAddr, VECLEN/2), masks[1-fb][dim]));

		return VECLEN;
	}
	else {
		ivector.push_back( new LoadFVec( r[0], rAddr, string("")));
		ivector.push_back( new LoadFVec( r[1], new AddressImm(rAddr, VECLEN), string("")));
		return 2*VECLEN;
	}
}

void unpackXYZTFVec(InstVector& ivector, const string &r, const Address*lAddr, const Address*rAddr, string dir, string srz) 
{
/*
        int dim = dir/2;
        int fb = dir % 2;
*/
        ifStatement(ivector, "(("+dir+")>>1) < 3"); 
	{
//                string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
                ivector.push_back( new LoadUnpackFVec(FVec(r+"[0]"), rAddr, "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new LoadUnpackFVec(FVec(r+"[0]"), lAddr, "PACKMASK[1-(("+dir+")&1)][("+dir+")>>1]"));

                ivector.push_back( new LoadUnpackFVec(FVec(r+"[1]"), new AddressImm(rAddr, VECLEN/2), "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new LoadUnpackFVec(FVec(r+"[1]"), new AddressImm(lAddr, VECLEN/2), "PACKMASK[1-(("+dir+")&1)][("+dir+")>>1]"));

                //return VECLEN;
                ivector.push_back( new IntToMask(srz, "VECLEN") );
        }
        elseStatement(ivector); 
	{
                ivector.push_back( new LoadFVec(FVec(r+"[0]"), rAddr, string("")));
                ivector.push_back( new LoadFVec(FVec(r+"[1]"), new AddressImm(rAddr, VECLEN), string("")));
                //return 2*VECLEN;
                ivector.push_back( new IntToMask(srz, "2*VECLEN") );
        }
	endScope(ivector);
}

int unpackXYZTFVec(InstVector& ivector, const FVec r[2], const Address*rAddr, int dir) {
        int dim = dir/2;
        int fb = dir % 2;
        if(dim < 3) {
                string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
                ivector.push_back( new LoadUnpackFVec(r[0], rAddr, masks[fb][dim]));
		ivector.push_back( new LoadUnpackFVec(r[1], new AddressImm(rAddr, VECLEN/2), masks[fb][dim]));
		return VECLEN;
	}
	else {
		ivector.push_back( new LoadFVec( r[0], rAddr, string("")));
		ivector.push_back( new LoadFVec( r[1], new AddressImm(rAddr, VECLEN), string("")));
		return 2*VECLEN;
	}
}

void unpackXYZTFVec(InstVector& ivector, const string& r, const Address*rAddr, string dir, string srz) {
/*
        int dim = dir/2;
        int fb = dir % 2;
*/
        ifStatement(ivector, "(("+dir+")>>1) < 3"); 
	{
                string masks[2][3] = {{"0x55", "0x33", "0x0F"}, {"0xAA", "0xCC", "0xF0"}};
                ivector.push_back( new LoadUnpackFVec(FVec(r+"[0]"), rAddr, "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                ivector.push_back( new LoadUnpackFVec(FVec(r+"[1]"), new AddressImm(rAddr, VECLEN/2), "PACKMASK[("+dir+")&1][("+dir+")>>1]"));
                //return VECLEN;
                ivector.push_back( new IntToMask(srz, "VECLEN") );
        }
        elseStatement(ivector); 
	{
                ivector.push_back( new LoadFVec( FVec(r+"[0]"), rAddr, string("")));
                ivector.push_back( new LoadFVec( FVec(r+"[1]"), new AddressImm(rAddr, VECLEN), string("")));
                //return 2*VECLEN;
                ivector.push_back( new IntToMask(srz, "2*VECLEN") );
        }
	endScope(ivector);
}

#endif // PRECISION == 2
