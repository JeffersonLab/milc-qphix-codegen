#include <stdlib.h>
#include <stdio.h>
#include "instructions.h"

#if PRECISION == 1 && VECLEN == 8
#pragma message "Using AVX Single Precision"

#define FVECTYPE "__m256"

const string fullIntMask("0xFF");
const string fullMask("0xFF");

FVec::FVec(const string& name_) : name(name_), type(FVECTYPE) {}

string DeclareFVec::serialize() const
{
#if 0
    return v.getType()+" "+v.getName()+ " = _mm256_undefined_ps();";
#else
    return v.getType()+" "+v.getName()+ " = _mm256_setzero_ps();";
#endif
}

string InitFVec::serialize() const
{
#if 0
    return v.getName()+ " = _mm256_undefined_ps();";
#else
    return v.getName()+ " = _mm256_setzero_ps();";
#endif
}

string Declare_MM_PERM::serialize() const
{
  return "_MM_PERM_ENUM MM_PERM[8] = { (_MM_PERM_ENUM)0xB1, (_MM_PERM_ENUM)0xB1, (_MM_PERM_ENUM)0x4E, (_MM_PERM_ENUM)0x4E, (_MM_PERM_ENUM)0x01, (_MM_PERM_ENUM)0x01, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00 };";
}

string DeclareMask::serialize() const
{
    ostringstream outbuf;

    if(value.empty()) {
        outbuf << "__m256 " << name << ";" << endl;
    }
    else {
        outbuf << "__m256 " << name << " = " << value << ";" << endl;
    }

    return outbuf.str();
}

string IntToMask::serialize() const
{
    ostringstream outbuf;
    outbuf << mask << " = _mm256_int2mask_ps(" << value << ");" << endl;
    return outbuf.str();
}

string DeclareOffsets::serialize() const
{
    ostringstream outbuf;
    outbuf << "__m256i " << vname << " = _mm256_load_si256((__m256i const *)" << pname << ");" << endl;
    return outbuf.str();
}

string IfAllOneCond::serialize() const
{
    return " if ((" + condition + " & " + fullIntMask + ") == " + fullIntMask + ") { ";
}

string LoadFVec::serialize() const
{
    std::ostringstream buf;

    if(mask.empty()) {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_load_ps(" << a->serialize() << ");" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtph_ps(_mm_load_si128((__m128i*)" << a->serialize() << "));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_maskload_ps(" << a->serialize() << ", _mm256_castps_si256(" << mask << "));" <<endl;
        }
        else {
            printf("Warning: AVX Masked load ignoring mask for half type\n");
            buf << v.getName() << " = _mm256_cvtph_ps(_mm_load_si128((__m128i*)" << a->serialize() << "));" <<endl;
        }
    }

    return buf.str();

}

string StoreFVec::serialize() const
{
    ostringstream buf;
    int streaming = isStreaming;

    if(streaming) {
        if(!a->isHalfType()) {
            buf << "_mm256_stream_ps(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm_stream_si128((__m128i *)" << a->serialize() << ", _mm256_cvtps_ph(" << v.getName() <<  ", _MM_FROUND_CUR_DIRECTION));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << "_mm256_store_ps(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm_store_si128((__m128i *)" << a->serialize() << ", _mm256_cvtps_ph(" << v.getName() <<  ", _MM_FROUND_CUR_DIRECTION));" <<endl;
        }
    }

    return buf.str();
}

string GatherFVec::serialize() const
{
    std::ostringstream buf;

    if(mask.empty()) {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_i32gather_ps(" << a->getBase()  << ", " << a->getOffsets() << ", _MM_SCALE_4);" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtph_ps(_mm_i32gather_ps(" << a->getBase()  << ", " << a->getOffsets() << ", _MM_SCALE_2));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_mask_i32gather_ps(" << v.getName() << ", " << a->getBase()  << ", " << a->getOffsets() << ", " << mask << ", _MM_SCALE_4);" <<endl;
        }
        else {
            printf("Error: AVX2 half type gather not supported\n");
            exit(1);
        }
    }

    return buf.str();

}

string ScatterFVec::serialize() const
{
    std::ostringstream buf;
    buf << "#error \"scatter is not supported in AVX/AVX2\"" <<endl;
    printf("scatter is not supported in AVX/AVX2\n");
    exit(1);
    return buf.str();
}

string LoadBroadcast::serialize() const
{
    std::ostringstream buf;

    if(!a->isHalfType()) {
        buf << v.getName() << " = _mm256_broadcast_ss(" << a->serialize() << ");" << endl;
    }
    else {
        buf << v.getName() << " = _mm256_cvtph_ps(_mm_set1_epi16(*" << a->serialize() << "));" << endl;
    }

    return buf.str();
}

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
        hint = "_MM_HINT_T0";
        break;

    case 3:
        hint = "_MM_HINT_NTA";
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
        hint = "_MM_HINT_T1";
        break;

    case 3:
        hint = "_MM_HINT_T2";
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
//		case 2: hint = "_MM_HINT_ET0"; break;
//		case 3: hint = "_MM_HINT_ENTA"; break;
    }
}
string GatherPrefetchL1::serialize() const
{
    std::ostringstream buf;
    buf << "#error \"Gather Prefetch is not supported in AVX/AVX2\"" <<endl;
    printf("Gather Prefetch is not supported in AVX/AVX2\n");
    exit(1);
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
//		case 2: hint = "_MM_HINT_ET1"; break;
//		case 3: hint = "_MM_HINT_ET2"; break;
    }
}
string GatherPrefetchL2::serialize() const
{
    std::ostringstream buf;
    buf << "#error \"Gather Prefetch is not supported in AVX/AVX2\"" <<endl;
    printf("Gather Prefetch is not supported in AVX/AVX2\n");
    exit(1);
    return buf.str();

}

string SetZero::serialize() const
{
    return  ret.getName()+" = _mm256_setzero_ps(); ";
}

string Set1Const::serialize() const {
    std::ostringstream buf;
    buf << ret.getName() << " = _mm256_set1_ps(" << val << "); " << endl;
    return  buf.str();
}

string Mul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_mul_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+ " = _mm256_blendv_ps(" + ret.getName() + ", _mm256_mul_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string FnMAdd::serialize() const
{
#ifndef AVX2

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sub_ps("+c.getName()+", _mm256_mul_ps("+a.getName()+" , "+b.getName()+"));" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", _mm256_sub_ps("+c.getName()+", _mm256_mul_ps("+a.getName()+" , "+b.getName()+")), " + mask + ");" ;
    }

#else

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_fnmadd_ps( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", _mm256_fnmadd_ps( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" ), " + mask + ");" ;
    }

#endif
}

string FMAdd::serialize() const
{
#ifndef AVX2

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_add_ps(_mm256_mul_ps("+a.getName()+", "+b.getName()+"), "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", _mm256_add_ps(_mm256_mul_ps("+a.getName()+", "+b.getName()+"), "+c.getName()+"), " + mask + ");" ;
    }

#else

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_fmadd_ps( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" =_mm256_blendv_ps(" + ret.getName() + ", _mm256_fmadd_ps( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" ), " + mask + ");" ;
    }

#endif
}

string Add::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_add_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", _mm256_add_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Sub::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sub_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", _mm256_sub_ps( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Div::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_div_ps( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_mask_div_ps( " + ret.getName() + "," + mask + "," +a.getName()+" , "+b.getName()+" );" ;
    }
}

string Sqrt::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sqrt_ps( "+a.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_mask_sqrt_ps( " + ret.getName() + "," + mask + "," +a.getName()+" );" ;
    }
}

string MovFVec::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = " + a.getName()+";" ;
    }
    else {
        return ret.getName()+" = _mm256_blendv_ps(" + ret.getName() + ", "+a.getName()+", " + mask + ");" ;
    }
}

class Perm32x4 : public Instruction
{
public:
    Perm32x4(const FVec& ret_, const FVec& a_, const FVec& b_, int imm_) : ret(ret_), a(a_), b(b_), imm(imm_) {}
    string serialize() const
    {
        ostringstream stream;
        stream << ret.getName() << " = _mm256_permute2f128_ps(" << a.getName() << ", "  << b.getName() << ", " << imm << ");" ;
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec ret;
    const FVec a;
    const FVec b;
    int imm;
};

class PermuteXYZTFVec : public Instruction
{
public:
    PermuteXYZTFVec(const FVec& ret_, const FVec& a_, int dir_) : ret(ret_), a(a_), dir(dir_/2) {}
    string serialize() const
    {
        ostringstream stream;
		if(dir < 2) {
			string imm = (dir == 0 ? "0xB1" : "0x4E");
	        stream << ret.getName() << " = _mm256_permute_ps(" << a.getName() << ", "  << imm << ");";
		}
		else if(dir == 2) {
			stream << ret.getName() << " = _mm256_permute2f128_ps(" << a.getName() << ", "  << a.getName() << ", 0x01);" ;
		}
		else {
			// nothing needs to be done for dir 3
		}
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
                stream << "if(((" << dir << ")>>1) == 0) " << ret.getName() << " = _mm256_permute_ps(" << a.getName() << ", (_MM_PERM_ENUM)0xB1);" << endl
			<< "else if(((" << dir << ")>>1) == 1) " << ret.getName() << " = _mm256_permute_ps(" << a.getName() << ", (_MM_PERM_ENUM)0x4E);" << endl
                        << "else if(((" << dir << ")>>1) == 2) " << ret.getName() << " = _mm256_permute2f128_ps(" << a.getName() << ", "  << a.getName() << ", 0x01);" ;
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

class PackXYZTFVec : public Instruction
{
public:
    PackXYZTFVec(const FVec& a_, const FVec& b_, const Address *lAddr_, const Address *rAddr_, int dir_) : a(a_), b(b_), lAddr(lAddr_), rAddr(rAddr_), dir(dir_) {}
    string serialize() const
    {
        ostringstream stream;
		int dim = dir / 2;
		int fb = dir % 2;

		if(dim < 2) {
			string imm[2] = {(dim == 0 ? "0x88" : "0x44"), (dim == 0 ? "0xDD" : "0xEE")};
	        stream << "_mm256_stream_ps(" << rAddr->serialize() << ", _mm256_shuffle_ps(" << a.getName() << ", " << b.getName() << ", "  << imm[fb] << "));" << endl;
	        stream << "_mm256_stream_ps(" << lAddr->serialize() << ", _mm256_shuffle_ps(" << a.getName() << ", " << b.getName() << ", "  << imm[1-fb] << "));";
		}
		else if(dim == 2) {
			string imm[2] = {"0x20", "0x31"};
	        stream << "_mm256_stream_ps(" << rAddr->serialize() << ", _mm256_permute2f128_ps(" << a.getName() << ", " << b.getName() << ", "  << imm[fb] << "));" << endl;
	        stream << "_mm256_stream_ps(" << lAddr->serialize() << ", _mm256_permute2f128_ps(" << a.getName() << ", " << b.getName() << ", "  << imm[1-fb] << "));";
		}
		else {
			// this is taken care outside
		}
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec a;
    const FVec b;
	const Address *lAddr, *rAddr;
    int dir;
};

class UnpackXYZTFVec : public Instruction
{
public:
    UnpackXYZTFVec(const FVec& a_, const FVec& b_, const Address *lAddr_, const Address *rAddr_, int dir_) : a(a_), b(b_), lAddr(lAddr_), rAddr(rAddr_), dir(dir_) {}
    string serialize() const
    {
        ostringstream stream;
		int dim = dir / 2;
		int fb = dir % 2;

		if(dim < 2) {
			const Address *adr[2] = {rAddr, lAddr};
	        stream << a.getName() << " = _mm256_shuffle_ps(_mm256_load_ps(" << adr[fb]->serialize() << "), _mm256_load_ps(" << adr[1-fb]->serialize() << "), 0x44);" << endl;
	        stream << b.getName() << " = _mm256_shuffle_ps(_mm256_load_ps(" << adr[fb]->serialize() << "), _mm256_load_ps(" << adr[1-fb]->serialize() << "), 0xEE);" << endl;
			if(dim == 0) {
				stream << a.getName() << " = _mm256_permute_ps(" << a.getName() << ", 0xD8);" << endl;
				stream << b.getName() << " = _mm256_permute_ps(" << b.getName() << ", 0xD8);" << endl;
			}
		}
		else if(dim == 2) {
			const Address *adr[2] = {rAddr, lAddr};
	        stream << a.getName() << " = _mm256_permute2f128_ps(_mm256_load_ps(" << adr[fb]->serialize() << "), _mm256_load_ps(" << adr[1-fb]->serialize() << "), 0x20);" << endl;
	        stream << b.getName() << " = _mm256_permute2f128_ps(_mm256_load_ps(" << adr[fb]->serialize() << "), _mm256_load_ps(" << adr[1-fb]->serialize() << "), 0x31);" << endl;
		}
		else {
			// this is taken care outside
		}
        return stream.str();
    }
    int numArithmeticInst() const
    {
        return 0;
    }
private:
    const FVec a;
    const FVec b;
	const Address *lAddr, *rAddr;
    int dir;
};

class LoadHalfFVec : public MemRefInstruction
{
public:
    LoadHalfFVec( const FVec& v_, const Address* a_, const int num_) : v(v_), a(a_), num(num_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!a->isHalfType()) {
            buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() << ", _mm_load_ps(" << a->serialize() << "), " << num << ");" << endl;
        }
        else {
            buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() << ", _mm_cvtph_ps(_mm_castpd_si128(_mm_load_sd((double*)" << a->serialize() << "))), " << num << ");" << endl;
        }

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
    const int num;
};

class StoreHalfFVec : public MemRefInstruction
{
public:
    StoreHalfFVec( const FVec& v_, const Address* a_, const int num_) : v(v_), a(a_), num(num_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!a->isHalfType()) {
            buf << "_mm_store_ps(" << a->serialize() <<  ", _mm256_extractf128_ps(" << v.getName() << ", " << num << "));" << endl;
        }
        else {
            buf << "_mm_store_sd((double*)" << a->serialize() <<  ", _mm_castsi128_pd(_mm_cvtps_ph(_mm256_extractf128_ps(" << v.getName() << ", " << num << "), _MM_FROUND_CUR_DIRECTION)));" << endl;
        }

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
    const int num;
};

class LoadSplitSOAFVec : public MemRefInstruction
{
public:
    LoadSplitSOAFVec( const FVec& v_, const Address* a1_, const Address* a2_, const int soanum_, const int soalen_, int forward_) : v(v_), a1(a1_), a2(a2_), soanum(soanum_), soalen(soalen_), forward(forward_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!a1->isHalfType()) {
            if(forward) {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_blend_ps(_mm256_loadu_ps(" << a1->serialize() << "), _mm256_broadcast_ss(" <<
                        a2->serialize() << "), " << (1 << (soalen-1)) << ");" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_blend_ps(_mm_loadu_ps(" << a1->serialize() << "), _mm_broadcast_ss(" << a2->serialize() << "), " << (1 << (soalen-1)) << "), " << soanum << ");" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_blend_ps(_mm256_loadu_ps((" << a2->serialize() << ")-1), _mm256_broadcast_ss(" << a1->serialize() << "), 1);" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_blend_ps(_mm_loadu_ps((" << a2->serialize() << ")-1), _mm_broadcast_ss(" << a1->serialize() << "), 1), " << soanum << ");" << endl;
                }
            }
        }
        else {
            if(forward) {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_cvtph_ps(_mm_blend_epi16(_mm_loadu_si128((__m128i*)" << a1->serialize() << "), _mm_set1_epi16(*" <<
                        a2->serialize() << "), " << (1 << (soalen-1)) << "));" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_cvtph_ps(_mm_blend_epi16(_mm_loadu_si128((__m128i*)" << a1->serialize() << "), _mm_set1_epi16(*" << a2->serialize() << "), " << (1 << (soalen-1)) << ")), " << soanum << ");" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_cvtph_ps(_mm_blend_epi16(_mm_loadu_si128((__m128i*)(" << a2->serialize() << ")-1), _mm_set1_epi16(*" << a1->serialize() << "), 1));" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_cvtph_ps(_mm_blend_epi16(_mm_loadu_si128((__m128i*)(" << a2->serialize() << ")-1), _mm_set1_epi16(*" << a1->serialize() << "), 1)), " << soanum << ");" << endl;
                }
            }
        }

        return buf.str();
    }
    const Address* getAddress() const
    {
        return a1;
    }
    MemRefType getType() const
    {
        return LOAD_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a1;
    const Address* a2;
    const int soalen, soanum;
    const int forward;
};

class StoreSplitSOAFVec : public MemRefInstruction
{
public:
    StoreSplitSOAFVec( const FVec& v_, const Address* a1_, const Address* a2_, const int soanum_, const int soalen_, int forward_) : v(v_), a1(a1_), a2(a2_), soanum(soanum_), soalen(soalen_), forward(forward_) {}
    string serialize() const {
        std::ostringstream buf;
        if(!a1->isHalfType()) {
            if(forward) {
                if(soalen == 8) {
                    buf << "_mm_storeu_ps(" << a1->serialize() << ", _mm256_castps256_ps128(" << v.getName() << "));" << endl;
                    buf << "_mm_store_sd((double*)(" << a1->serialize() << " + 4), _mm_castps_pd(_mm256_extractf128_ps(" << v.getName() << ", 1)));" << endl;
                    buf << "((int*)(" << a1->serialize() << " + 6))[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", 1), 2);" << endl;
                    buf << "((int*)" << a2->serialize() << ")[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", 1), 3);" << endl;
                }
                else {
                    buf << "_mm_store_sd((double*)" << a1->serialize() << ", _mm_castps_pd(_mm256_extractf128_ps(" << v.getName() << ", " << soanum << ")));" << endl;
                    buf << "((int*)(" << a1->serialize() << " + 2))[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", " << soanum << "), 2);" << endl;
                    buf << "((int*)" << a2->serialize() << ")[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", " << soanum << "), 3);" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << "_mm_store_ss(" << a1->serialize() << ", _mm256_castps256_ps128(" << v.getName() << "));" << endl;
                    buf << "((int*)" << a2->serialize() << ")[0] = _mm_extract_ps(_mm256_castps256_ps128(" << v.getName() << "), 1);" << endl;
                    buf << "_mm_storeh_pd((double*)(" << a2->serialize() << " + 1), _mm_castps_pd(_mm256_castps256_ps128(" << v.getName() << ")));" << endl;
                    buf << "_mm_storeu_ps((" << a2->serialize() << " + 3), _mm256_extractf128_ps(" << v.getName() << ", 1));" << endl;
                } else {
                    buf << "_mm_store_ss(" << a1->serialize() << ", _mm256_extractf128_ps(" << v.getName() << ", " << soanum << "));" << endl;
                    buf << "((int*)" << a2->serialize() << ")[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", " << soanum << "), 1);" << endl;
                    buf << "_mm_storeh_pd((double*)(" << a2->serialize() << " + 1), _mm_castps_pd(_mm256_extractf128_ps(" << v.getName() << ", " << soanum << ")));" << endl;
                }
            }
        }
        else {
            buf << "#pragma error \"Half Prec Not Implemented yet!\"" << endl;
        }
        return buf.str();
    }
    const Address* getAddress() const {
        return a1;
    }
    MemRefType getType() const {
        return LOAD_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a1;
    const Address* a2;
    const int soalen, soanum;
    const int forward;
};

class LoadSplit3SOAFVec : public MemRefInstruction
{
public:
    LoadSplit3SOAFVec( const FVec& v_, const Address* a1_, const Address* a2_, const Address* a3_, const int soanum_, const int soalen_, int forward_) : v(v_), a1(a1_), a2(a2_), a3(a3_), soanum(soanum_), soalen(soalen_), forward(forward_) {}
    string serialize() const {
        std::ostringstream buf;
        if(!a1->isHalfType()) {
            if(forward) {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_blend_ps(_mm256_blend_ps(_mm256_loadu_ps(" << a1->serialize() << "), _mm256_broadcast_ss(" <<
                        a2->serialize() << "), " << (1 << (soalen-2)) << "), _mm256_broadcast_ss(" << a3->serialize() << "), " << (1 << (soalen-1)) << ");" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_blend_ps(_mm_blend_ps(_mm_loadu_ps(" << a1->serialize() << "), _mm_broadcast_ss(" << a2->serialize() << "), " << (1 << (soalen-2)) <<
                        "), _mm_broadcast_ss(" << a3->serialize() << "), " << (1 << (soalen-1)) << "), " << soanum << ");" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_blend_ps(_mm256_blend_ps(_mm256_loadu_ps((" << a3->serialize() << ")-2), _mm256_broadcast_ss(" << a2->serialize() <<
                        "), 2), _mm256_broadcast_ss(" << a1->serialize() << "), 1);" << endl;
                } else {
                    buf << v.getName() << " =  _mm256_insertf128_ps(" << v.getName() <<
                        ", _mm_blend_ps(_mm_blend_ps(_mm_loadu_ps((" << a3->serialize() << ")-2), _mm_broadcast_ss(" << a2->serialize() <<
                        "), 2), _mm_broadcast_ss(" << a1->serialize() << "), 1), " << soanum << ");" << endl;
                }
            }
        }
        else {
            buf << "#pragma error \"Half Prec Not Implemented yet!\"" << endl;
        }
        return buf.str();
    }
    const Address* getAddress() const {
        return a1;
    }
    MemRefType getType() const {
        return LOAD_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a1;
    const Address* a2;
    const Address* a3;
    const int soalen, soanum;
    const int forward;
};

class CondInsertFVecElement : public MemRefInstruction
{
public:
    CondInsertFVecElement( const FVec& v_, const Address* a_, const string mask_, int pos_, bool skipCond_) : v(v_), a(a_), mask(mask_), pos(pos_), skipCond(skipCond_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!skipCond) {
            buf << "if(" << mask << " & " << (1 << pos) << ") ";
        }

        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_blend_ps(" << v.getName() << ", _mm256_broadcast_ss(" << a->serialize() << "), " << (1<<pos) << ");" << endl;
        }
        else {
            buf << v.getName() << " = _mm256_blend_ps(" << v.getName() << ", _mm256_cvtph_ps(_mm_set1_epi16(*" << a->serialize() << ")), " << (1<<pos) << ");" << endl;
        }

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
    const int pos;
    const bool skipCond;
};

class CondExtractFVecElement : public MemRefInstruction
{
public:
    CondExtractFVecElement( const FVec& v_, const Address* a_, const string mask_, int pos_, bool skipCond_) : v(v_), a(a_), mask(mask_), pos(pos_), skipCond(skipCond_) {}
    string serialize() const
    {
        std::ostringstream buf;

        if(!skipCond) {
            buf << "if(" << mask << " & " << (1 << pos) << ") ";
        }

        if(!a->isHalfType()) {
            buf << "((int*)" << a->serialize() << ")[0] = _mm_extract_ps(_mm256_extractf128_ps(" << v.getName() << ", " << pos / 4 << "), " << (pos % 4) << ");" << endl;
        }
        else {
            buf << "(" << a->serialize() << ")[0] = _mm_extract_epi16(_mm256_cvtps_ph(" << v.getName() << ", _MM_FROUND_CUR_DIRECTION), " << pos << ");" << endl;
        }

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
    const int pos;
    const bool skipCond;
};

class UnpackFVec : public MemRefInstruction
{
public:
    UnpackFVec( const FVec& v_, const Address* a_, string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const {
        std::ostringstream buf;
        buf << "_mm256_expand_ps(" << a->serialize() << ", " << mask << ", " << v.getName() << ");" << endl;
        return buf.str();
    }
    const Address* getAddress() const {
        return a;
    }
    MemRefType getType() const {
        return STORE_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
};

class PackFVec : public MemRefInstruction
{
public:
    PackFVec( const FVec& v_, const Address* a_, string mask_) : v(v_), a(a_), mask(mask_) {}
    string serialize() const {
        std::ostringstream buf;
        buf << "_mm256_compress_ps(" << a->serialize() << ", " << mask << ", " << v.getName() << ");" << endl;
        return buf.str();
    }
    const Address* getAddress() const {
        return a;
    }
    MemRefType getType() const {
        return STORE_MASKED_VEC;
    }
private:
    const FVec v;
    const Address* a;
    const string mask;
};

void loadSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    if(soalen == 8) {
        ivector.push_back( new LoadFVec(ret, a, string("")));
    }
    else if(soalen == 4) {
        ivector.push_back( new LoadHalfFVec(ret, a, soanum));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    if(soalen == 8) {
        ivector.push_back( new StoreFVec(ret, a, 0));
    }
    else if(soalen == 4) {
        ivector.push_back( new StoreHalfFVec(ret, a, soanum));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void loadSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward)
{
    ivector.push_back( new LoadSplitSOAFVec(ret, a1, a2, soanum, soalen, forward));
}

void storeSplitSOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, int soanum, int soalen, int forward)
{
    ivector.push_back( new StoreSplitSOAFVec(ret, a1, a2, soanum, soalen, forward));
}

void loadSplit3SOAFVec(InstVector& ivector, const FVec& ret, const Address *a1, const Address *a2, const Address *a3, int soanum, int soalen, int forward)
{
    ivector.push_back( new LoadSplit3SOAFVec(ret, a1, a2, a3, soanum, soalen, forward));
}

void unpackFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    int pos = 0, nBits = 0;

    for(int i = 0; i < 8; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    if(possibleMask!=0xFFFF) {
        for(int i = 0; i < 8; i++) if(possibleMask & (1 << i)) {
                ivector.push_back( new CondInsertFVecElement(ret, new AddressImm(a, pos), mask, i, nBits==1));
            }
    }
    else {
        ivector.push_back( new UnpackFVec(ret, a, mask));
    }
}

void packFVec(InstVector& ivector, const FVec& ret, Address *a, string mask, int possibleMask)
{
    int pos = 0, nBits = 0;

    for(int i = 0; i < 8; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    if(possibleMask != 0xFFFF) {
        for(int i = 0; i < 8; i++) if(possibleMask & (1 << i)) {
                ivector.push_back( new CondExtractFVecElement(ret, new AddressImm(a, pos), mask, i, nBits==1));
            }
    }
    else {
        ivector.push_back( new PackFVec(ret, a, mask));
    }
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

void fperm32x4(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const int imm)
{
    ivector.push_back(new Perm32x4(ret, a, b, imm));
}

void transpose2x2(InstVector& ivector, const FVec r[2], const FVec f[2])
{
    for(int i = 0; i < 2; i++) {
        fperm32x4(ivector, r[i], f[0], f[1], 0x20+i*0x11);
    }
}

void transpose1x1(InstVector& ivector, const FVec r[1], const FVec f[1])
{
    movFVec(ivector, r[0], f[0], string(""));
}

void transpose(InstVector& ivector, const FVec r[], const FVec f[], int soalen)
{
    switch (soalen) {
    case 4:
        transpose2x2(ivector, r, f);
        break;

    case 8:
        transpose1x1(ivector, r, f);
        break;

    default:
        printf("SOALEN = %d Not Supported (only SOALEN = 4, 8 supported)\n", soalen);
    }
}

void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir)
{
	ivector.push_back(new PermuteXYZTFVec(r, f, dir));
}

void permuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir)
{
        ivector.push_back(new PermuteXYZTFVecDir(r, f, dir));
}

void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, int dir)
{
        ivector.push_back(new PermuteXYZTFVec(r, f, dir));
}

void aPermuteXYZTFVec(InstVector& ivector, const FVec r, const FVec f, string dir)
{
        ivector.push_back(new PermuteXYZTFVecDir(r, f, dir));
}

int packXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir) 
{
	//string masks[2][3] = {{"0x55", "0xCC", "0x0F"}, {"0xAA", "0x33", "0xF0"}};
	if(dir < 6) {
		ivector.push_back( new PackXYZTFVec(r[0], r[1], lAddr, rAddr, dir));
		return VECLEN;
	}
	else {
		ivector.push_back( new StoreFVec( r[0], rAddr, 1));
		ivector.push_back( new StoreFVec( r[1], new AddressImm(rAddr, VECLEN), 1));
		return 2*VECLEN;
	}
}

int unpackXYZTFVec(InstVector& ivector, const FVec r[2], const Address*lAddr, const Address*rAddr, int dir) {
	if(dir < 6) {
		ivector.push_back( new UnpackXYZTFVec(r[0], r[1], lAddr, rAddr, dir));
		return VECLEN;
	}
	else {
		ivector.push_back( new LoadFVec( r[0], rAddr, string("")));
		ivector.push_back( new LoadFVec( r[1], new AddressImm(rAddr, VECLEN), string("")));
		return 2*VECLEN;
	}
}

#endif // PRECISION == 1
