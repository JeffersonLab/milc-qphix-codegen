
#ifndef _ADDRESS_TYPES_H_
#define _ADDRESS_TYPES_H_

using namespace std;

enum AddressType { GAUGE, HERMIT, HERMITHELPER, SPINOR, KS, CLOVER_DIAG, CLOVER_OFFDIAG, ADDRESS_OF_SCALAR, GENERIC_ADDRESS};

class Address
{
public:
    Address(int isHalfType_) : halfType(isHalfType_) { }
    virtual string serialize(void) const=0;
    virtual AddressType getType(void) const = 0;
    int isHalfType(void)
    {
        return halfType;
    }
    int isHalfType(void) const
    {
        return halfType;
    }
private:
    int halfType;
};

class GenericAddress : public Address
{
public:
    GenericAddress(string base_, int isHalfType) : Address(isHalfType), base(base_) {}
    AddressType getType(void) const
    {
        return GENERIC_ADDRESS;
    }
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(" << base << ")";
        return outbuf.str();
    }

private:
    const string base;
};

class AddressOfScalar : public Address
{
public:
    AddressOfScalar(string scalar_name_, int isHalfType) : Address(isHalfType), scalar_name(scalar_name_) {}
    AddressType getType(void) const
    {
        return ADDRESS_OF_SCALAR;
    }
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(&" << scalar_name << ")";
        return outbuf.str();
    }

private:
    const string scalar_name;
};

class GaugeAddress : public Address
{
public:
    GaugeAddress(const string& base_, int dir_, int c1_, int c2_, int reim_, int isHalfType) : Address(isHalfType), base(base_), dir(dir_), c1(c1_), c2(c2_), reim(reim_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir <<"]["<< c1 << "][" << c2 << "][" << reim << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return GAUGE;
    }
protected:
    const string base;
    const int dir;
    const int c1;
    const int c2;
    const int reim;

};

class GaugeAddressStr : public Address
{
public:
    GaugeAddressStr(const string& base_, const string& dir_, const string& c1_, const string& c2_, int reim_, int isHalfType) : Address(isHalfType), base(base_), dir(dir_), c1(c1_), c2(c2_), reim(reim_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir << "][" << c1 << "][" << c2 << "][" << reim << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return GAUGE;
    }
protected:
    const string base;
    const string dir;
    const string c1;
    const string c2;
    const int reim;
};

class HermitAddress : public Address
{
public:
    HermitAddress(const string& base_, int dir_, int k_, int isHalfType) : Address(isHalfType), base(base_), dir(dir_), k(k_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir <<"]["<< k << "]";

	return outbuf.str();
    }
    AddressType getType(void) const
    {
        return HERMIT;
    }
protected:
    const string base;
    const int dir;
    const int k;

};

class HermitAddressStr : public Address
{
public:
    HermitAddressStr(const string& base_, const string& dir_, const string& k_, int isHalfType) : Address(isHalfType), base(base_), dir(dir_), k(k_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir <<"]["<< k << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return HERMIT;
    }
protected:
    const string base;
    const string dir;
    const string k;
};

class HermitHelperAddress : public Address
{
public:
  HermitHelperAddress(const string& base_, int dir1_, int dir2_, int k_, int isHalfType) : Address(isHalfType), base(base_), dir1(dir1_), dir2(dir2_), k(k_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir1 << "][" << dir2 << "]["<< k << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return HERMITHELPER;
    }
protected:
    const string base;
    const int dir1;
    const int dir2;
    const int k;

};

class HermitHelperAddressStr : public Address
{
public:
  HermitHelperAddressStr(const string& base_, const string& dir1_, const string& dir2_, const string& k_, int isHalfType) : Address(isHalfType), base(base_), dir1(dir1_), dir2(dir2_), k(k_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")[" << dir1 << "][" << dir2 << "]["<< k << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return HERMITHELPER;
    }
protected:
    const string base;
    const string dir1;
    const string dir2;
    const string k;

};

class ClovDiagAddress : public Address
{
public:
    ClovDiagAddress(const string& base_, int block_, int component_, int isHalfType) : Address(isHalfType), base(base_), block(block_), component(component_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ").diag"<<block+1<<"[ " << component <<" ] ";
        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return CLOVER_DIAG;
    }
protected:
    const string base;
    const int block;
    const int component;
};

class ClovOffDiagAddress : public Address
{
public:
    ClovOffDiagAddress(const string& base_, int block_, int component_, int reim_, int isHalfType) : Address(isHalfType), base(base_), block(block_), component(component_), reim(reim_) {};

    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ").off_diag"<<block+1<<"[ " << component <<" ][ " << reim << "]";
        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return CLOVER_OFFDIAG;
    }
protected:
    const string base;
    const int block;
    const int component;
    const int reim;
};

class SpinorAddress : public Address
{
public:
    SpinorAddress(const string& base_, int spin_, int color_, int reim_, int isHalfType) : Address(isHalfType), base(base_), spin(spin_), color(color_), reim(reim_) {};
    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")["<< color<< "][" << spin << "][" << reim << "]";

        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return SPINOR;
    }
    //int getOffset(void) const {return (((spin*3+color)*2+reim)*SOALEN);}
protected:
    const string base;
    const int spin;
    const int color;
    const int reim;
};



class KSAddress : public Address {
public:
    KSAddress(const string& base_, int color_, int reim_, int isHalfType) : Address(isHalfType), base(base_), color(color_), reim(reim_) {};
    string serialize(void) const {
        ostringstream outbuf;
        outbuf<< "(*" << base << ")["<< color<< "][" << reim << "]";

        return outbuf.str();
    }
    AddressType getType(void) const {
        return KS;
    }
protected:
    const string base;
    const int color;
    const int reim;
};

class AddressOffset : public Address
{
public:
    AddressOffset(const Address* a_, string offset_var_) : Address(a_->isHalfType()), a(a_), offset_var(offset_var_) {}
    AddressType getType(void) const
    {
        return a->getType();
    }
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(" << a->serialize() << " + " << offset_var <<")";
        return outbuf.str();
    }

private:
    const Address* a;
    const string offset_var;
};

class AddressScaledOffset : public Address {
public:
    AddressScaledOffset(const Address* a_, string offset_var_, int scale_) : Address(a_->isHalfType()), a(a_), offset_var(offset_var_), scale(scale_) {}
    AddressType getType(void) const {
        return a->getType();
    }
    string serialize() const {
        ostringstream outbuf;
        outbuf << "(" << a->serialize() << " + (" << offset_var <<" * "<< scale <<"))";
        return outbuf.str();
    }

private:
    const Address* a;
    const string offset_var;
    const int scale;
};

class AddressImm : public Address
{
public:
    AddressImm(const Address* base_, const int imm_) : Address(base_->isHalfType()), base(base_), imm(imm_) {}
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(" << base->serialize() << "+" << imm << ")";
        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return base->getType();
    }
protected:
    const Address* base;
    const int imm;
};

class AddressImmStr : public Address
{
public:
    AddressImmStr(const Address* base_, string imm_) : Address(base_->isHalfType()), base(base_), imm(imm_) {}
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(" << base->serialize() << "+" << imm << ")";
        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return base->getType();
    }
protected:
    const Address* base;
    const string imm;
};

class IndirectAddress : public Address
{
public:
    IndirectAddress(const Address* a_, string offset_var_, int imm_) : Address(a_->isHalfType()), a(a_), offset_var(offset_var_), imm(imm_) {}
    AddressType getType(void) const
    {
        return a->getType();
    }
    string serialize() const
    {
        ostringstream outbuf;
        outbuf << "(" << a->serialize() << " + " << offset_var <<"[" << imm << "])";
        return outbuf.str();
    }

private:
    const Address* a;
    const string offset_var;
    const int imm;
};

class GatherAddress : public Address
{
public:
    GatherAddress(const Address *base_, const string& offsets_) : Address(base_->isHalfType()), base(base_), offsets(offsets_) {};
    string serialize(void) const
    {
        ostringstream outbuf;
        outbuf<< "(" << base->serialize() << ")["<< offsets<< "[0:VECLEN]]";
        return outbuf.str();
    }
    AddressType getType(void) const
    {
        return base->getType();
    }
    string getBase(void) const
    {
        return base->serialize();
    }
    string getOffsets(bool prefixed = true) const
    {
        if(prefixed == true) {
            return "vec_" + offsets;
        }
        else {
            return offsets;
        }
    }
    const Address* getAddr(int ind) const
    {
        return new IndirectAddress(base, offsets, ind);
    }
protected:
    const Address* base;
    const string offsets;
};

#endif //_ADDRESS_TYPES_H_
