/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Douglas Slotta
 *
 */

/*
 * @file gt8.cpp
 *
 * @author slottad@ncbi.nlm.nih.gov
 *
 * @brief SciDB shared library for encoding genotypes into 8bits
 *
 */

#include <vector>
#include <algorithm>
#include <boost/assign.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "query/Operator.h"
#include "query/FunctionLibrary.h"
#include "query/FunctionDescription.h"
#include "query/TypeSystem.h"
#include "system/ErrorsLibrary.h"

using namespace std;
using namespace scidb;
using namespace boost::assign;

typedef uint8_t gt8_t;

enum {
  GT8_E_CANT_CONVERT_TO_GT8 = SCIDB_USER_ERROR_CODE_START
};

EXPORTED_FUNCTION void GetPluginVersion(uint32_t& major, uint32_t& minor, 
                                        uint32_t& patch, uint32_t& build)
{
    // Provide correct values here. SciDB checks it and does not allow
    // loading plugins that are too new.
    major = scidb::SCIDB_VERSION_MAJOR();
    minor = scidb::SCIDB_VERSION_MINOR();
    patch = scidb::SCIDB_VERSION_PATCH();
    build = scidb::SCIDB_VERSION_BUILD();
}

//string extract_value(const scidb::Value** args)
void extract_value(const scidb::Value** args, scidb::Value* res, void*)
{
    string const& key = args[0]->getString();
    string const& fmtStr = args[1]->getString();
    string const& valStr = args[2]->getString();

    vector<string> keys;
    split(keys, fmtStr, is_any_of(":"));

    vector<string>::iterator it = find(keys.begin(),
                                       keys.end(),
                                       key);
    if (it == keys.end()) return;
    
    size_t pos = it - keys.begin();
    vector<string> values;
    split(values, valStr, is_any_of(":"));

    if (values.size() <= pos) return;
    
    res->setString(values[pos].c_str());
}

void num_csv(const scidb::Value** args, scidb::Value* res, void*)
{
    string const& cell = args[0]->getString();

    uint32_t count = 0;

    if (!cell.empty()) {
        count = 1 + std::count(cell.begin(), cell.end(), ',');
    }

    res->setUint32(count);
}

void gt8_alleleCount(const scidb::Value** args, scidb::Value* res, void* misc)
{
    gt8_t* g = static_cast<gt8_t*>( args[0]->data() );
    int64_t ai = args[1]->getInt64();

    uint32_t ac = 0;

    uint8_t ploidy = *g & 0x80;
    
    if (ploidy == 0) 
    {
        if (*g != 0) 
            if ((*g-1) == ai)
                ac++;
    } 
    else 
    {
        uint8_t a = *g & 0x38;
        a >>= 3;

        uint8_t b = *g & 0x07;

        if (a > 0) 
            if((int64_t)(a-1) == ai)
                ac++;
        if (b > 0) 
            if((int64_t)(b-1) == ai)
                ac++;
    }
    
    res->setUint32(ac);
}

void gt8_ploidy(const scidb::Value** args, scidb::Value* res, void* misc)
{
    gt8_t* g = static_cast<gt8_t*>( args[0]->data() );
    uint8_t pl = *g & 0x80;
    if (pl == 0)
        pl = 1;
    else
        pl = 2;
    res->setUint8(pl);
}

void gt8_phase(const scidb::Value** args, scidb::Value* res, void* misc)
{
    gt8_t* g = static_cast<gt8_t*>( args[0]->data() );
    uint8_t p = *g & 0x40;
    res->setBool(p!=0);
}

void gt8_lessEqualThan(const Value** args, Value* res, void*)
{
    gt8_t& lhs = *(gt8_t*)args[0]->data();
    gt8_t& rhs = *(gt8_t*)args[1]->data();
    bool result = false;

    unsigned loadLhs = lhs;
    unsigned loadRhs = rhs;
    if(lhs & 0x80)
    {
        uint8_t a1 = lhs & 0x38;
        a1 >>= 3;
        uint8_t b1 = lhs & 0x07;
        loadLhs = a1 + b1;
    }
    if (rhs & 0x80)
    {
        uint8_t a2 = rhs & 0x38;
        a2 >>= 3;
        uint8_t b2 = rhs & 0x07;
        loadRhs = a2 + b2;
    }
    // REF < ALT
    if (loadLhs <= loadRhs)
        result = true;
    else if (loadLhs == loadRhs) 
    {
        if ((lhs & 0x40) != (rhs & 0x40))
            // unphased < phased
            result = ((lhs & 0x40) <= (rhs & 0x40));
        else // finaly haploid < diploid
            result = ((lhs & 0x80) <= (rhs & 0x80));
    }
    res->setBool(result);
}

void gt8_lessThan(const Value** args, Value* res, void*)
{
    gt8_t& lhs = *(gt8_t*)args[0]->data();
    gt8_t& rhs = *(gt8_t*)args[1]->data();
    bool result = false;

    unsigned loadLhs = lhs;
    unsigned loadRhs = rhs;
    if(lhs & 0x80)
    {
        uint8_t a1 = lhs & 0x38;
        a1 >>= 3;
        uint8_t b1 = lhs & 0x07;
        loadLhs = a1 + b1;
    }
    if (rhs & 0x80)
    {
        uint8_t a2 = rhs & 0x38;
        a2 >>= 3;
        uint8_t b2 = rhs & 0x07;
        loadRhs = a2 + b2;
    }
    // REF < ALT
    if (loadLhs < loadRhs)
        result = true;
    else if (loadLhs == loadRhs) 
    {
        if ((lhs & 0x40) != (rhs & 0x40))
            // unphased < phased
            result = ((lhs & 0x40) < (rhs & 0x40));
        else // finaly haploid < diploid
            result = ((lhs & 0x80) < (rhs & 0x80));
    }
    res->setBool(result);
}

void gt8_equal(const Value** args, Value* res, void*)
{
    gt8_t& lhs = *(gt8_t*)args[0]->data();
    gt8_t& rhs = *(gt8_t*)args[1]->data();
    bool result = (lhs == rhs);
    res->setBool(result);
}

void gt8_fromString(const scidb::Value** args, scidb::Value* res, void*)
{
    string const& gstr = args[0]->getString();
    size_t phase_loc = gstr.find_first_of("/|");
    bool isDiploid;
    unsigned short int a, b;
    char phase('/');
    if (phase_loc != string::npos) {
        isDiploid = true;
        phase = gstr[phase_loc];
        string strA = gstr.substr(0, phase_loc);
        string strB = gstr.substr(phase_loc+1);
        if (strA[0] == '.') {
            a=0;
        } else {
            if (sscanf(strA.c_str(), "%hu", &a) < 1)
                throw PLUGIN_USER_EXCEPTION("libgt8", scidb::SCIDB_SE_UDO, 
                                             GT8_E_CANT_CONVERT_TO_GT8) << gstr;
            a++;
        }
        if (strB[0] == '.') {
            b=0;
        } else {
            if (sscanf(strB.c_str(), "%hu", &b) < 1)
                throw PLUGIN_USER_EXCEPTION("libgt8", scidb::SCIDB_SE_UDO, 
                                             GT8_E_CANT_CONVERT_TO_GT8) << gstr;
            b++;
        }
    } else {
        isDiploid = false;
        if (gstr[0] == '.') {
            a=0;
        } else {
            if (sscanf(gstr.c_str(), "%hu", &a) < 1)
                throw PLUGIN_USER_EXCEPTION("libgt8", scidb::SCIDB_SE_UDO, 
                                             GT8_E_CANT_CONVERT_TO_GT8) << gstr;
            a++;
        }
    }
    
    gt8_t* g = static_cast<gt8_t*>( res->data() );
    *g = a;
    if (isDiploid) {
        *g <<= 3;
        *g |= b;
        if (phase == '|') {
            *g |= 0xC0;  // 0x80 & 0x40
        } else {
            *g |= 0x80;
        }
    }
}

void gt8_toString(const scidb::Value** args, scidb::Value* res, void*)
{
    gt8_t* g = static_cast<gt8_t*>( args[0]->data() );

    stringstream ss;

    uint8_t ploidy = *g & 0x80;
    
    if (ploidy == 0) {
        if (*g == 0) 
            ss << ".";
        else
            ss << (*g-1);
    } else {
        uint8_t a = *g & 0x38;
        a >>= 3;

        uint8_t b = *g & 0x07;

        uint8_t p = *g & 0x40;
        char phase = (p==0) ? '/' : '|';

        if (a == 0) ss << ".";
        else ss << (static_cast<int>(a)-1);

        ss << phase ;
    
        if (b == 0) ss << ".";
        else ss << (static_cast<int>(b)-1);
    }
    
    res->setString(ss.str().c_str());
}

void construct_gt8(const scidb::Value** args, scidb::Value* res, void*)
{
    *(gt8_t*)res->data() = 0;
}

void gt8_normalize(const scidb::Value** args, scidb::Value* res, void*)
{
    gt8_t* gtIn = static_cast<gt8_t*>( args[0]->data() );
    gt8_t* gtOut = static_cast<gt8_t*>( res->data() );
    
    uint8_t ploidy = *gtIn & 0x80;
    
    if (ploidy == 0) {
        *gtOut = *gtIn;
        return;
    } else {
        uint8_t a = *gtIn & 0x38;
        a >>= 3;
        uint8_t b = *gtIn & 0x07;

        uint8_t p = *gtIn & 0x40;

        if ((a <= b) && (p == 0)) {
            *gtOut = *gtIn;
            return;
        }
        
        *gtOut = b;
        *gtOut <<= 3;
        *gtOut |= a;
        *gtOut |= 0x80;
    }
}

void bitwise_not8(const scidb::Value** args, scidb::Value* res, void*)
{
    uint8_t* left = static_cast<uint8_t*>( args[0]->data() );
    uint8_t* result = static_cast<uint8_t*>( res->data() );

    *result = ~(*left);
}

void bitwise_not16(const scidb::Value** args, scidb::Value* res, void*)
{
    uint16_t* left = static_cast<uint16_t*>( args[0]->data() );
    uint16_t* result = static_cast<uint16_t*>( res->data() );

    *result = ~(*left);
}

void bitwise_not32(const scidb::Value** args, scidb::Value* res, void*)
{
    uint32_t* left = static_cast<uint32_t*>( args[0]->data() );
    uint32_t* result = static_cast<uint32_t*>( res->data() );

    *result = ~(*left);
}

void bitwise_not64(const scidb::Value** args, scidb::Value* res, void*)
{
    uint64_t* left = static_cast<uint64_t*>( args[0]->data() );
    uint64_t* result = static_cast<uint64_t*>( res->data() );

    *result = ~(*left);
}

void bitwise_and8(const scidb::Value** args, scidb::Value* res, void*)
{
    uint8_t* left = static_cast<uint8_t*>( args[0]->data() );
    uint8_t* right = static_cast<uint8_t*>( args[1]->data() );
    uint8_t* result = static_cast<uint8_t*>( res->data() );

    *result = *left & *right;
}

void bitwise_and16(const scidb::Value** args, scidb::Value* res, void*)
{
    uint16_t* left = static_cast<uint16_t*>( args[0]->data() );
    uint16_t* right = static_cast<uint16_t*>( args[1]->data() );
    uint16_t* result = static_cast<uint16_t*>( res->data() );

    *result = *left & *right;
}

void bitwise_and32(const scidb::Value** args, scidb::Value* res, void*)
{
    uint32_t* left = static_cast<uint32_t*>( args[0]->data() );
    uint32_t* right = static_cast<uint32_t*>( args[1]->data() );
    uint32_t* result = static_cast<uint32_t*>( res->data() );

    *result = *left & *right;
}

void bitwise_and64(const scidb::Value** args, scidb::Value* res, void*)
{
    uint64_t* left = static_cast<uint64_t*>( args[0]->data() );
    uint64_t* right = static_cast<uint64_t*>( args[1]->data() );
    uint64_t* result = static_cast<uint64_t*>( res->data() );

    *result = *left & *right;
}

void bitwise_or8(const scidb::Value** args, scidb::Value* res, void*)
{
    uint8_t* left = static_cast<uint8_t*>( args[0]->data() );
    uint8_t* right = static_cast<uint8_t*>( args[1]->data() );
    uint8_t* result = static_cast<uint8_t*>( res->data() );

    *result = *left | *right;
}

void bitwise_or16(const scidb::Value** args, scidb::Value* res, void*)
{
    uint16_t* left = static_cast<uint16_t*>( args[0]->data() );
    uint16_t* right = static_cast<uint16_t*>( args[1]->data() );
    uint16_t* result = static_cast<uint16_t*>( res->data() );

    *result = *left | *right;
}

void bitwise_or32(const scidb::Value** args, scidb::Value* res, void*)
{
    uint32_t* left = static_cast<uint32_t*>( args[0]->data() );
    uint32_t* right = static_cast<uint32_t*>( args[1]->data() );
    uint32_t* result = static_cast<uint32_t*>( res->data() );

    *result = *left | *right;
}

void bitwise_or64(const scidb::Value** args, scidb::Value* res, void*)
{
    uint64_t* left = static_cast<uint64_t*>( args[0]->data() );
    uint64_t* right = static_cast<uint64_t*>( args[1]->data() );
    uint64_t* result = static_cast<uint64_t*>( res->data() );

    *result = *left | *right;
}

void bitwise_xor8(const scidb::Value** args, scidb::Value* res, void*)
{
    uint8_t* left = static_cast<uint8_t*>( args[0]->data() );
    uint8_t* right = static_cast<uint8_t*>( args[1]->data() );
    uint8_t* result = static_cast<uint8_t*>( res->data() );

    *result = *left ^ *right;
}

void bitwise_xor16(const scidb::Value** args, scidb::Value* res, void*)
{
    uint16_t* left = static_cast<uint16_t*>( args[0]->data() );
    uint16_t* right = static_cast<uint16_t*>( args[1]->data() );
    uint16_t* result = static_cast<uint16_t*>( res->data() );

    *result = *left ^ *right;
}

void bitwise_xor32(const scidb::Value** args, scidb::Value* res, void*)
{
    uint32_t* left = static_cast<uint32_t*>( args[0]->data() );
    uint32_t* right = static_cast<uint32_t*>( args[1]->data() );
    uint32_t* result = static_cast<uint32_t*>( res->data() );

    *result = *left ^ *right;
}

void bitwise_xor64(const scidb::Value** args, scidb::Value* res, void*)
{
    uint64_t* left = static_cast<uint64_t*>( args[0]->data() );
    uint64_t* right = static_cast<uint64_t*>( args[1]->data() );
    uint64_t* result = static_cast<uint64_t*>( res->data() );

    *result = *left ^ *right;
}

REGISTER_TYPE(gt8, sizeof(gt8_t));

REGISTER_FUNCTION(extract_value, list_of(TID_STRING)(TID_STRING)(TID_STRING), TID_STRING, extract_value);
REGISTER_FUNCTION(num_csv, list_of(TID_STRING), "uint32", num_csv);
REGISTER_FUNCTION(allele_count, list_of("gt8")("int64"), "uint32", gt8_alleleCount);
REGISTER_FUNCTION(ploidy, list_of("gt8"), "uint8", gt8_ploidy);
REGISTER_FUNCTION(phase, list_of("gt8"), "bool", gt8_phase);
REGISTER_FUNCTION(<=, list_of("gt8")("gt8"), "bool", gt8_lessEqualThan);
REGISTER_FUNCTION(<, list_of("gt8")("gt8"), "bool", gt8_lessThan);
REGISTER_FUNCTION(=, list_of("gt8")("gt8"), "bool", gt8_equal);
REGISTER_FUNCTION(norm, list_of("gt8"), "gt8", gt8_normalize);
REGISTER_FUNCTION(gt8, ArgTypes(), "gt8", construct_gt8);

// Bitwise operators
REGISTER_FUNCTION(bitnot, list_of("uint8"), "uint8", bitwise_not8);
REGISTER_FUNCTION(bitnot, list_of("uint16"), "uint16", bitwise_not16);
REGISTER_FUNCTION(bitnot, list_of("uint32"), "uint32", bitwise_not32);
REGISTER_FUNCTION(bitnot, list_of("uint64"), "uint64", bitwise_not64);

REGISTER_FUNCTION(bitand, list_of("uint8")("uint8"), "uint8", bitwise_and8);
REGISTER_FUNCTION(bitand, list_of("uint16")("uint16"), "uint16", bitwise_and16);
REGISTER_FUNCTION(bitand, list_of("uint32")("uint32"), "uint32", bitwise_and32);
REGISTER_FUNCTION(bitand, list_of("uint64")("uint64"), "uint64", bitwise_and64);

REGISTER_FUNCTION(bitor, list_of("uint8")("uint8"), "uint8", bitwise_or8);
REGISTER_FUNCTION(bitor, list_of("uint16")("uint16"), "uint16", bitwise_or16);
REGISTER_FUNCTION(bitor, list_of("uint32")("uint32"), "uint32", bitwise_or32);
REGISTER_FUNCTION(bitor, list_of("uint64")("uint64"), "uint64", bitwise_or64);

REGISTER_FUNCTION(bitxor, list_of("uint8")("uint8"), "uint8", bitwise_xor8);
REGISTER_FUNCTION(bitxor, list_of("uint16")("uint16"), "uint16", bitwise_xor16);
REGISTER_FUNCTION(bitxor, list_of("uint32")("uint32"), "uint32", bitwise_xor32);
REGISTER_FUNCTION(bitxor, list_of("uint64")("uint64"), "uint64", bitwise_xor64);

REGISTER_CONVERTER(gt8, string, EXPLICIT_CONVERSION_COST, gt8_toString);
REGISTER_CONVERTER(string, gt8, EXPLICIT_CONVERSION_COST, gt8_fromString);

/*
 * Class for registering/unregistering user defined objects
 */
static class Gt8Library
{
public:
    // Registering objects
    Gt8Library()
    {
        Type("gt8", sizeof(gt8_t) * 8);

        _errors[GT8_E_CANT_CONVERT_TO_GT8] = "Cannot convert '%1%' to gt8";
        scidb::ErrorsLibrary::getInstance()->registerErrors("gt8", &_errors);
    }

    ~Gt8Library()
    {
        // This segfaults because it tries to call a mutex that has
        // already been destroyed.  Do we really need it?
        scidb::ErrorsLibrary::getInstance()->unregisterErrors("libgt8");
    }

private:
    scidb::ErrorsLibrary::ErrorsMessages _errors;
} _instance;
