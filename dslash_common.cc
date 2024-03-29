#include "dslash_common.h"
using namespace std;

#if PRECISION == 1
#if VECLEN == 16
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 8
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 4
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#elif PRECISION == 2
#if VECLEN == 8
#ifdef AVX512
std::string ARCH_NAME="avx512";
#else
std::string ARCH_NAME="mic";
#endif
#elif VECLEN == 4
#ifdef AVX2
std::string ARCH_NAME="avx2";
#else
std::string ARCH_NAME="avx";
#endif
#elif VECLEN == 2
std::string ARCH_NAME="sse";
#elif VECLEN == 1
std::string ARCH_NAME="scalar";
#endif
#endif //PRECISION

void declare_HalfSpinor(InstVector& ivector, FVec h_spinor[2][3][2])
{
    for(int s=0; s < 2; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, h_spinor[s][c][RE]);
            declareFVecFromFVec(ivector, h_spinor[s][c][IM]);
        }
    }
}

void declare_Gauge(InstVector& ivector, FVec u_gauge[3][3][2]) {
    for(int c1=0; c1 < 3; c1++) {
        for(int c2 = 0; c2 < 3; c2++) {
            declareFVecFromFVec(ivector, u_gauge[c1][c2][RE]);
            declareFVecFromFVec(ivector, u_gauge[c1][c2][IM]);
        }
    }
}

void declare_WilsonSpinor(InstVector& ivector, FVec spinor[4][3][2]) {
    for(int s=0; s < 4; s++) {
        for(int c = 0; c < 3; c++) {
            declareFVecFromFVec(ivector, spinor[s][c][RE]);
            declareFVecFromFVec(ivector, spinor[s][c][IM]);
        }
    }
}

void declare_KSSpinor(InstVector& ivector, FVec ks_spinor[3][2]) {
    for(int c = 0; c < 3; c++) {
        declareFVecFromFVec(ivector, ks_spinor[c][RE]);
        declareFVecFromFVec(ivector, ks_spinor[c][IM]);
    }
}

void declare_Clover(InstVector& ivector, FVec diag[6], FVec offdiag[15][2]) {
    for(int s=0; s < 6; s++) {
        declareFVecFromFVec(ivector, diag[s]);
    }
    for(int s=0; s < 15; s++) {
        declareFVecFromFVec(ivector, offdiag[s][RE]);
        declareFVecFromFVec(ivector, offdiag[s][IM]);
    }
}


// r[RE] = s1[RE]-beta_vec*s2[RE] = fnmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]-beta_vec*s2[IM] = fnamdd(beta_vec,s2[IM],s1[IM])
void addCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] + beta_vec*s2[RE] = fmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[IM] = fmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec * s2[IM] = fmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] - beta_vec * s2[RE] = fnmadd(beta_vec, s2[RE], s1[IM])

void addiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec*s2[IM] = fnmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] + beta_vec*s2[RE] = fmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_mbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE]+beta_vec*s2[RE] = fmadd(beta_vec,s2[RE],s1[RE])
// r[IM] = s1[IM]+beta_vec*s2[IM] = fmadd(beta_vec,s2[IM],s1[IM])
void addCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}
// r[RE] = s1[RE] - beta_vec*s2[RE] = fnmadd(beta_vec, s2[RE], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[IM] = fnmadd(beta_vec, s2[IM], s1[IM]);
void subCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[RE], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[IM], s1[IM], mask);
}

// r[RE] = s1[RE] - beta_vec * s2[IM] = fnmadd(beta_vec,s2[IM], s1[RE])
// r[IM] = s1[IM] + beta_vec * s2[RE] = fmadd(beta_vec, s2[RE], s1[IM])
void addiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fnmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}

// r[RE] = s1[RE] + beta_vec*s2[IM] = fmadd( beta_vec, s2[IM], s1[RE]);
// r[IM] = s1[IM] - beta_vec*s2[RE] = fnmadd ( beta_vec, s2[RE], s1[IM]);
void subiCVec_pbeta(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec &beta_vec, string &mask)
{
    fmaddFVec(ivector, r[RE], beta_vec, s2[IM], s1[RE], mask);
    fnmaddFVec(ivector, r[IM], beta_vec, s2[RE], s1[IM], mask);
}


// r = (s1*s2-s3*s4)'
//r[RE] = (s1[RE]*s2[RE])-(s1[IM]*s2[IM])-(s3[RE]*s4[RE])+(s3[IM]*s4[IM])
//r[IM] = (s3[RE]*s4[IM])+(s3[IM]*s4[RE])-(s1[RE]*s2[IM])-(s1[IM]*s2[RE])
void Conj_CrossProd(InstVector& ivector, FVec *r, FVec *s1, FVec *s2, FVec *s3, FVec *s4, string &mask)
{
    mulFVec(ivector, r[RE], s1[RE], s2[RE], mask);
    fnmaddFVec(ivector, r[RE], s1[IM], s2[IM], r[RE], mask);
    fnmaddFVec(ivector, r[RE], s3[RE], s4[RE], r[RE], mask);
    fmaddFVec(ivector, r[RE], s3[IM], s4[IM], r[RE], mask);

    mulFVec(ivector, r[IM], s3[RE], s4[IM], mask);
    fmaddFVec(ivector, r[IM], s3[IM], s4[RE], r[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[RE], s2[IM], r[IM], mask);
    fnmaddFVec(ivector, r[IM], s1[IM], s2[RE], r[IM], mask);
}

// Merge L2 prefetches with another instruction stream
void mergeIvectorWithL2Prefetches(InstVector& ivector, InstVector& l2prefs)
{
    if( l2prefs.size() == 0 ) {
        cout << "No L2 Prefetches. Returning ivector unchanged" << endl;
    }
    else {

        if( ivector.size() == 0 ) {
            cout << "No actual instructions to merge. Returning " << endl;
            return;
        }


        int ivector_size = ivector.size();
        vector<Instruction*>::iterator it;

        // Skip past first declarations
        it=ivector.begin();
        while ( (*it)->numDeclarations() > 0 ) {
            it++;
            ivector_size--;
        }

        int n_prefs = l2prefs.size();

        //cout << "After declarations Ivector size is " << ivector_size << endl;
        //cout << "PrefetchL2 size is " << n_prefs << endl;

        int spacing_factor = ivector_size / n_prefs;
        //cout << "Spacing Factor is  " << spacing_factor << endl;
        int pref=0;
        for(int i=0; i < n_prefs; i++) {
            it = ivector.insert(it, l2prefs[pref]);
            pref++;
            it++;
            int j=spacing_factor;
            while (j > 0) {
                it++;
                j--;
            }
        }
    }
}

// Dump an instruction stream into a file
void dumpIVector(InstVector& ivector, string filename)
{
    // Check if we have any gather instructions and if we need to convert
    // indices in an offset array to an index register before dumping the instructions
    std::map<string,string> offslist;
    for(int i=0; i < ivector.size(); i++) {
        Instruction *inst = ivector[i];
        if ( inst->hasAddress() ) {
            MemRefInstruction* mr = dynamic_cast< MemRefInstruction* >(inst);
            if(mr->hasGSAddress()) {
                const GatherAddress* ga = dynamic_cast<const GatherAddress *>(mr->getAddress());
                string offs = ga->getOffsets(false);
                string voffs = ga->getOffsets ();
                if(offslist.find(offs) == offslist.end())
                    offslist[offs] = voffs;
            }
        }
    }
    for ( map<string,string>::iterator  it = offslist.begin(); it != offslist.end(); ++it ) {
        ivector.insert(ivector.begin(), new DeclareOffsets(it->first, it->second ));
    }
    ofstream outfile(filename.c_str());
    for(int i=0; i < ivector.size(); i++) {
        outfile << ivector[i]->serialize() << endl;
    }
    outfile.close();
}
