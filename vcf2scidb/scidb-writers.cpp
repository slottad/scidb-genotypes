/*  $Id: scidb-writers.cpp 394139 2013-03-30 03:26:57Z slottad $
 * ===========================================================================
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
 * Authors:  Douglas Slotta
 *
 * File Description:
 *   Implementation of various writer classes for SciDB
 *
 */
#include "scidb-writers.hpp"

#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <sstream>

#include <stdlib.h>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

scidb_writer::scidb_writer(string const& filename)
{
    /* Assumes the file exists and is probably a pipe, pipes required
     * to be read by other processes do not like to have the mode set
     * to be usable by the current user only */
    _out = open(filename.c_str(), O_WRONLY | O_TRUNC);

    if (_out == -1) {
        /* Probably doesn't exist, let try create it and assume the
         * previous failure was not catastrophic */
        _out = open(filename.c_str(), O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR);
    }
}

scidb_writer::~scidb_writer()
{
    close(_out);
}

scidb_text_writer::scidb_text_writer(string const& filename, size_t chunksize)
    : scidb_writer(filename), _chunksize(chunksize), _rowcount(0), _newchunk(false)
{
    write(_out, "[\n", 2);
}

scidb_text_writer::~scidb_text_writer()
{
    if (! _newchunk)
        write(_out, "]\n", 2);
    else
        write(_out, "\n", 1);

}

void scidb_text_writer::set_var(int64_t var)
{
    _prefix =  "(\"" + _chrom + "\"," + _pos + ",";
    _prefix.append(lexical_cast<string>(var));
    _prefix.push_back(',');
}

void scidb_text_writer::put_prefix()
{
    if (_newchunk) {
        write(_out, ";\n[\n", 4);
        _newchunk=false;
    }
    write(_out, _prefix.c_str(), _prefix.size());
}

void scidb_text_writer::put_separator()
{
    write(_out, ",", 1);
}

void scidb_text_writer::put_endrow()
{
    ++_rowcount;
    write(_out, ")\n", 2);
    if ((_rowcount % _chunksize) == 0) {
        write(_out,"]",1);
        _newchunk = true;
    }
}

void scidb_text_writer::put_data(std::string const& data, ENullData nullstatus, EDataType type)
{

    if (data.empty() && (nullstatus == eNullable)) {
        write(_out, "?", 1);
    } else {
        bool isString = ((type == eString) || (type == eGt8));
        if (isString) write(_out, "\"", 1);
        write(_out, data.c_str(), data.size());
        if (isString) write(_out, "\"", 1);
    }
}

void scidb_text_writer::put_uint32(uint32_t data)
{
    string strData = lexical_cast<string>(data);
    write(_out, strData.c_str(), strData.size());
}

void scidb_text_writer::put_int64(int64_t data)
{
    string strData = lexical_cast<string>(data);
    write(_out, strData.c_str(), strData.size());
}

scidb_binary_writer::scidb_binary_writer(string const& filename)
    : scidb_writer(filename)
{}

scidb_binary_writer::~scidb_binary_writer()
{
}

void scidb_binary_writer::set_var(int64_t var)
{
    ostringstream oss;

    uint32_t sz = _chrom.size()+1;
    oss.write(reinterpret_cast<char *>(&sz), sizeof(sz));
    oss.write(_chrom.c_str(), sz);
    int64_t pos = atol(_pos.c_str());
    oss.write(reinterpret_cast<char *>(&pos), sizeof(pos));
    oss.write(reinterpret_cast<char *>(&var), sizeof(var));
    _prefix = oss.str();
}

void scidb_binary_writer::put_prefix()
{
    write(_out, _prefix.c_str(), _prefix.size());
}

void scidb_binary_writer::put_separator()
{
}

void scidb_binary_writer::put_endrow()
{
}

uint8_t str2gt8(string const& gstr)
{
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
                cerr << "Can't convert " << gstr << " to gt8\n";
            a++;
        }
        if (strB[0] == '.') {
            b=0;
        } else {
            if (sscanf(strB.c_str(), "%hu", &b) < 1)
                cerr << "Can't convert " << gstr << " to gt8\n";
            b++;
        }
    } else {
        isDiploid = false;
        if (gstr[0] == '.') {
            a=0;
        } else {
            if (sscanf(gstr.c_str(), "%hu", &a) < 1)
                cerr << "Can't convert " << gstr << " to gt8\n";
            a++;
        }
    }
    
    uint8_t g = a;    
    if (isDiploid) {
        g <<= 3;
        g |= b;
        if (phase == '|') {
            g |= 0xC0;  // 0x80 & 0x40
        } else {
            g |= 0x80;
        }
    }
    return g;
}

void scidb_binary_writer::put_data(std::string const& data, ENullData nullstatus, EDataType type)
{
    if (nullstatus == eNullable) {
        int8_t nullval = data.empty() ? 0 : -1;
        write(_out, &nullval, sizeof(nullval));
    }
    switch (type) {
    case (eGt8): {
        uint8_t gt = str2gt8(data);
        write(_out, &gt, sizeof(gt));
    } break;
    case (eString): {
        uint32_t sz = data.size()+1;
        write(_out, &sz, sizeof(sz));
        if (sz > 0) write(_out, data.c_str(), sz);
    } break;
    case (eFloat): {
        float d = atof(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eDouble): {
        double d = atof(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eInt8): {
        int8_t d =atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eInt16): {
        int16_t d = atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eInt32): {
        int32_t d = atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eInt64): {
        int64_t d = atol(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eUint8): {
        uint8_t d = atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eUint16): {
        uint16_t d = atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eUint32): {
        uint32_t d = atoi(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    case (eUint64): {
        uint64_t d = atol(data.c_str());
        write(_out, &d, sizeof(d));
    } break;
    }
}

void scidb_binary_writer::put_uint32(uint32_t data)
{
    write(_out, &data, sizeof(data));
}

void scidb_binary_writer::put_int64(int64_t data)
{
    write(_out, &data, sizeof(data));
}
