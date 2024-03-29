#include <stdlib.h>
#include <stdio.h>
#include "instructions.h"

#if PRECISION == 2 && VECLEN == 4
#pragma message "Using AVX double Precision"

#define FVECTYPE "__m256d"

const string fullIntMask("0xF");
const string fullMask("0xF");

FVec::FVec(const string& name_) : name(name_), type(FVECTYPE) {}

string DeclareFVec::serialize() const
{
#if 0
    return v.getType()+" "+v.getName()+ " = _mm256_undefined_pd();";
#else
    return v.getType()+" "+v.getName()+ " = _mm256_setzero_pd();";
#endif
}

string InitFVec::serialize() const
{
#if 0
    return v.getName()+ " = _mm256_undefined_pd();";
#else
    return v.getName()+ " = _mm256_setzero_pd();";
#endif
}

/*
string Declare_MM_PERM::serialize() const
{
    return "_MM_PERM_ENUM MM_PERM[8] = { (_MM_PERM_ENUM)0x05, (_MM_PERM_ENUM)0x05, (_MM_PERM_ENUM)0x01, (_MM_PERM_ENUM)0x01, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00, (_MM_PERM_ENUM)0x00 };";
}
*/

string DeclareMask::serialize() const
{
    ostringstream outbuf;

    if(value.empty()) {
        outbuf << "__m256d " << name << ";" << endl;
    }
    else {
        outbuf << "__m256d " << name << " = " << value << ";" << endl;
    }

    return outbuf.str();
}

string IntToMask::serialize() const
{
    ostringstream outbuf;
    outbuf << mask << " = _mm256_int2mask_pd(" << value << ");" << endl;
    return outbuf.str();
}

string DeclareOffsets::serialize() const
{
    ostringstream outbuf;
    outbuf << "__m128i " << vname << " = _mm_load_si128((__m128i const *)" << pname << ");" << endl;
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
            buf << v.getName() << " = _mm256_load_pd(" << a->serialize() << ");" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtps_pd(_mm_load_ps(" << a->serialize() << "));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_maskload_pd(" << a->serialize() << ", _mm256_castpd_si256(" << mask << "));" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtps_pd(_mm_maskload_ps(" << a->serialize() << ", _mm_castps_si128(_mm256_cvtpd_ps(" << mask << "))));" <<endl;
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
            buf << "_mm256_stream_pd(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm_stream_ps(" << a->serialize() << ", _mm256_cvtpd_ps(" << v.getName() <<  "));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << "_mm256_store_pd(" << a->serialize() << ", " << v.getName() <<  ");" <<endl;
        }
        else {
            buf << "_mm_store_ps(" << a->serialize() << ", _mm256_cvtpd_ps(" << v.getName() <<  "));" <<endl;
        }
    }

    return buf.str();
}

string GatherFVec::serialize() const
{
    std::ostringstream buf;

    if(mask.empty()) {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_i32gather_pd(" << a->getBase()  << ", " << a->getOffsets() << ", _MM_SCALE_8);" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtps_pd(_mm_i32gather_ps(" << a->getBase()  << ", " << a->getOffsets() << ", _MM_SCALE_4));" <<endl;
        }
    }
    else {
        if(!a->isHalfType()) {
            buf << v.getName() << " = _mm256_mask_i32gather_pd(" << v.getName() << ", " << a->getBase()  << ", " << a->getOffsets() << ", " << mask << ", _MM_SCALE_8);" <<endl;
        }
        else {
            buf << v.getName() << " = _mm256_cvtps_pd(_mm__mask_i32gather_ps(" << v.getName() << ", " << a->getBase()  << ", " << a->getOffsets() << ", " << mask << ", _MM_SCALE_4));" <<endl;
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
        buf << v.getName() << " = _mm256_broadcast_sd(" << a->serialize() << ");" << endl;
    }
    else {
        buf << v.getName() << " = _mm256_cvtps_pd(_mm_broadcast_ss(" << a->serialize() << "));" << endl;
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
    return  ret.getName()+" = _mm256_setzero_pd(); ";
}

string Set1Const::serialize() const {
    std::ostringstream buf;
    buf << ret.getName() << " = _mm256_set1_pd(" << val << "); " << endl;
    return  buf.str();
}

string Mul::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_mul_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+ " = _mm256_blendv_pd(" + ret.getName() + ", _mm256_mul_pd( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string FnMAdd::serialize() const
{
#ifndef AVX2

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sub_pd("+c.getName()+", _mm256_mul_pd("+a.getName()+" , "+b.getName()+"));" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", _mm256_sub_pd("+c.getName()+", _mm256_mul_pd("+a.getName()+" , "+b.getName()+")), " + mask + ");" ;
    }

#else

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_fnmadd_pd( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", _mm256_fnmadd_pd( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" ), " + mask + ");" ;
    }

#endif
}

string FMAdd::serialize() const
{
#ifndef AVX2

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_add_pd(_mm256_mul_pd("+a.getName()+", "+b.getName()+"), "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", _mm256_add_pd(_mm256_mul_pd("+a.getName()+", "+b.getName()+"), "+c.getName()+"), " + mask + ");" ;
    }

#else

    if(mask.empty()) {
        return  ret.getName()+" = _mm256_fmadd_pd( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" );" ;
    }
    else {
        return  ret.getName()+" =_mm256_blendv_pd(" + ret.getName() + ", _mm256_fmadd_pd( "+a.getName()+" , "+b.getName()+" , "+c.getName()+" ), " + mask + ");" ;
    }

#endif
}

string Add::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_add_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", _mm256_add_pd( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Sub::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sub_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", _mm256_sub_pd( "+a.getName()+" , "+b.getName()+"), " + mask + ");" ;
    }
}

string Div::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_div_pd( "+a.getName()+" , "+b.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_mask_div_pd( " + ret.getName() + "," + mask + "," +a.getName()+" , "+b.getName()+" );" ;
    }
}

string Sqrt::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = _mm256_sqrt_pd( "+a.getName()+" );" ;
    }
    else {
        return  ret.getName()+" = _mm256_mask_sqrt_pd( " + ret.getName() + "," + mask + "," +a.getName()+" );" ;
    }
}

string MovFVec::serialize() const
{
    if(mask.empty()) {
        return  ret.getName()+" = " + a.getName()+";" ;
    }
    else {
        return ret.getName()+" = _mm256_blendv_pd(" + ret.getName() + ", "+a.getName()+", " + mask + ");" ;
    }
}

class Perm64x2 : public Instruction
{
public:
    Perm64x2(const FVec& ret_, const FVec& a_, const FVec& b_, int imm_) : ret(ret_), a(a_), b(b_), imm(imm_) {}
    string serialize() const
    {
        ostringstream stream;
        stream << ret.getName() << " = _mm256_permute2f128_pd(" << a.getName() << ", "  << b.getName() << ", " << imm << ");" ;
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
			if(dir == 0) {
				stream << ret.getName() << " = _mm256_permute_pd(" << a.getName() << ", 0x05);";
			}
			else if(dir == 1) {
				stream << ret.getName() << " = _mm256_permute2f128_pd(" << a.getName() << ", "  << a.getName() << ", 0x01);" ;
			}
			else {
				// nothing needs t be done for dir 2 & 3
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
                        stream << "if(((" << dir << ")>>1) == 0) " << ret.getName() << " = _mm256_permute_pd(" << a.getName() << ", 0x05);" << endl
                                << "else if(((" << dir << ")>>1) == 1) " << ret.getName() << " = _mm256_permute2f128_pd(" << a.getName() << ", "  << a.getName() << ", 0x01);" ;
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

			if(dim == 0) {
				string imm[2] = {"0x00", "0x0F"};
				stream << "_mm256_stream_pd(" << rAddr->serialize() << ", _mm256_shuffle_pd(" << a.getName() << ", " << b.getName() << ", "  << imm[fb] << "));" << endl;
				stream << "_mm256_stream_pd(" << lAddr->serialize() << ", _mm256_shuffle_pd(" << a.getName() << ", " << b.getName() << ", "  << imm[1-fb] << "));";
			}
			else if(dim == 1) {
				string imm[2] = {"0x20", "0x31"};
				stream << "_mm256_stream_pd(" << rAddr->serialize() << ", _mm256_permute2f128_pd(" << a.getName() << ", " << b.getName() << ", "  << imm[fb] << "));" << endl;
				stream << "_mm256_stream_pd(" << lAddr->serialize() << ", _mm256_permute2f128_pd(" << a.getName() << ", " << b.getName() << ", "  << imm[1-fb] << "));";
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

		if(dim == 0) {
			const Address *adr[2] = {rAddr, lAddr};
	        stream << a.getName() << " = _mm256_shuffle_pd(_mm256_load_pd(" << adr[fb]->serialize() << "), _mm256_load_pd(" << adr[1-fb]->serialize() << "), 0x00);" << endl;
	        stream << b.getName() << " = _mm256_shuffle_pd(_mm256_load_pd(" << adr[fb]->serialize() << "), _mm256_load_pd(" << adr[1-fb]->serialize() << "), 0x0F);" << endl;
		}
		else if(dim == 1) {
			const Address *adr[2] = {rAddr, lAddr};
	        stream << a.getName() << " = _mm256_permute2f128_pd(_mm256_load_pd(" << adr[fb]->serialize() << "), _mm256_load_pd(" << adr[1-fb]->serialize() << "), 0x20);" << endl;
	        stream << b.getName() << " = _mm256_permute2f128_pd(_mm256_load_pd(" << adr[fb]->serialize() << "), _mm256_load_pd(" << adr[1-fb]->serialize() << "), 0x31);" << endl;
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
            buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() << ", _mm_load_pd(" << a->serialize() << "), " << num << ");" << endl;
        }
        else {
            buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() << ", _mm_cvtps_pd(_mm_castpd_ps(_mm_load_sd((double*)" << a->serialize() << "))), " << num << ");" << endl;
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
            buf << "_mm_store_pd(" << a->serialize() <<  ", _mm256_extractf128_pd(" << v.getName() << ", " << num << "));" << endl;
        }
        else {
            buf << "_mm_store_sd((double*)" << a->serialize() <<  ", _mm_castps_pd(_mm_cvtpd_ps(_mm256_extractf128_pd(" << v.getName() << ", " << num << "))));" << endl;
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
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm256_blend_pd(_mm256_loadu_pd(" << a1->serialize() << "), _mm256_broadcast_sd(" <<
                        a2->serialize() << "), " << (1 << (soalen-1)) << ");" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() <<
                        ", _mm_blend_pd(_mm_loadu_pd(" << a1->serialize() << "), _mm_loaddup_pd(" << a2->serialize() << "), " << (1 << (soalen-1)) << "), " << soanum << ");" << endl;
                }
            }
            else {
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm256_blend_pd(_mm256_loadu_pd((" << a2->serialize() << ")-1), _mm256_broadcast_sd(" << a1->serialize() << "), 1);" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() <<
                        ", _mm_blend_pd(_mm_loadu_pd((" << a2->serialize() << ")-1), _mm_loaddup_pd(" << a1->serialize() << "), 1), " << soanum << ");" << endl;
                }
            }
        }
        else {
            if(forward) {
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm256_cvtps_pd(_mm_blend_ps(_mm_loadu_ps(" << a1->serialize() << "), _mm_broadcast_ss(" <<
                        a2->serialize() << "), " << (1 << (soalen-1)) << "));" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() <<
                        ", _mm_cvtps_pd(_mm_blend_ps(_mm_loadu_ps(" << a1->serialize() << "), _mm_broadcast_ss(" << a2->serialize() << "), " << (1 << (soalen-1)) << ")), " << soanum << ");" << endl;
                }
            }
            else {
                if(soalen == 4) {
                    buf << v.getName() << " =  _mm256_cvtps_pd(_mm_blend_ps(_mm_loadu_ps((" << a2->serialize() << ")-1), _mm_broadcast_ss(" << a1->serialize() << "), 1));" << endl;
                }
                else {
                    buf << v.getName() << " =  _mm256_insertf128_pd(" << v.getName() <<
                        ", _mm_cvtps_pd(_mm_blend_ps(_mm_loadu_ps((" << a2->serialize() << ")-1), _mm_broadcast_ss(" << a1->serialize() << "), 1)), " << soanum << ");" << endl;
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
                if(soalen == 4) {
                    buf << "_mm_storeu_ps(" << a1->serialize() << ", _mm_castpd_ps(_mm256_castpd256_pd128(" << v.getName() << ")));" << endl;
                    buf << "_mm_store_sd((" << a1->serialize() << " + 2), _mm256_extractf128_pd(" << v.getName() << ", 1));" << endl;
                    buf << "_mm_storeh_pd(" << a2->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", 1));" << endl;
                }
                else {
                    buf << "_mm_store_sd(" << a1->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << soanum << "));" << endl;
                    buf << "_mm_storeh_pd(" << a2->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << soanum << "));" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << "_mm_store_sd(" << a1->serialize() << ", _mm256_castpd256_pd128(" << v.getName() << "));" << endl;
                    buf << "_mm_storeh_pd(" << a2->serialize() << ", _mm256_castpd256_pd128(" << v.getName() << "));" << endl;
                    buf << "_mm_storeu_pd((" << a2->serialize() << " + 1), _mm256_extractf128_pd(" << v.getName() << ", 1));" << endl;
                } else {
                    buf << "_mm_store_sd(" << a1->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << soanum << "));" << endl;
                    buf << "_mm_storeh_pd(" << a2->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << soanum << "));" << endl;
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
                    buf << v.getName() << " =  _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadu_pd(" << a1->serialize() << "), _mm_loadh_pd(_mm_load_sd(" <<
                        a2->serialize() << "), " <<  a3->serialize() << "), " << "1);" << endl;
                }
                else {
                    buf << "#pragma error \"Should not get here!\"" << endl;
                }
            }
            else {
                if(soalen == 8) {
                    buf << v.getName() << " =  _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_loadh_pd(_mm_load_sd(" << a1->serialize() << "), " <<
                        a2->serialize() << "), _mm_loadu_pd(" <<  a3->serialize() << "), " << "1);" << endl;
                } else {
                    buf << "#pragma error \"Should not get here!\"" << endl;
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
            buf << v.getName() << " = _mm256_blend_pd(" << v.getName() << ", _mm256_broadcast_sd(" << a->serialize() << "), " << (1<<pos) << ");" << endl;
        }
        else {
            buf << v.getName() << " = _mm256_blend_pd(" << v.getName() << ", _mm256_cvtps_pd(_mm_broadcast_ss(" << a->serialize() << ")), " << (1<<pos) << ");" << endl;
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
            if(pos % 2 == 0) {
                buf << "_mm_store_sd(" << a->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << pos / 2 << "));" << endl;
            }
            else {
                buf << "_mm_storeh_pd(" << a->serialize() << ", _mm256_extractf128_pd(" << v.getName() << ", " << pos / 2 << "));" << endl;
            }
        }
        else {
            buf << "((int*)" << a->serialize() << ")[0] = _mm_extract_ps(_mm256_cvtpd_ps(" << v.getName() << "), " << pos << ");" << endl;
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
        buf << "_mm256_expand_pd(" << a->serialize() << ", " << mask << ", " << v.getName() << ");" << endl;
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
        buf << "_mm256_compress_pd(" << a->serialize() << ", " << mask << ", " << v.getName() << ");" << endl;
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
    if(soalen == 4) {
        ivector.push_back( new LoadFVec(ret, a, string("")));
    }
    else if(soalen == 2) {
        ivector.push_back( new LoadHalfFVec(ret, a, soanum));
    }
    else {
        printf("SOALEN = %d not supported\n", soalen);
        exit(1);
    }
}

void storeSOAFVec(InstVector& ivector, const FVec& ret, const Address *a, int soanum, int soalen)
{
    if(soalen == 4) {
        ivector.push_back( new StoreFVec(ret, a, 0));
    }
    else if(soalen == 2) {
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

    for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    if(possibleMask!=0xFFFF) {
        for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
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

    for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
            nBits++;
        }

    if(possibleMask != 0xFFFF) {
        for(int i = 0; i < 4; i++) if(possibleMask & (1 << i)) {
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

void fperm64x2(InstVector& ivector, const FVec& ret, const FVec& a, const FVec& b, const int imm)
{
    ivector.push_back(new Perm64x2(ret, a, b, imm));
}

void transpose2x2(InstVector& ivector, const FVec r[2], const FVec f[2])
{
    for(int i = 0; i < 2; i++) {
        fperm64x2(ivector, r[i], f[0], f[1], 0x20+i*0x11);
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
        transpose2x2(ivector, r, f);
        break;

    case 4:
        transpose1x1(ivector, r, f);
        break;

    default:
        printf("SOALEN = %d Not Supported (only SOALEN = 2, 4 supported)\n", soalen);
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
	if(dir < 4) {
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
	if(dir < 4) {
		ivector.push_back( new UnpackXYZTFVec(r[0], r[1], lAddr, rAddr, dir));
		return VECLEN;
	}
	else {
		ivector.push_back( new LoadFVec( r[0], rAddr, string("")));
		ivector.push_back( new LoadFVec( r[1], new AddressImm(rAddr, VECLEN), string("")));
		return 2*VECLEN;
	}
}

#endif // PRECISION == 2
