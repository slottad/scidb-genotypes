/*  $Id: scidb-writers.hpp 394139 2013-03-30 03:26:57Z slottad $
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
 *   Various writer classes for SciDB
 *
 */

#ifndef SCIDB_WRITERS_HPP
#define SCIDB_WRITERS_HPP
#include <string>
#include <map>
#include <stdint.h>

enum ENullData {
    eNullable,
    eNotNullable
};

enum EDataType {
    eString,
    eFloat,
    eDouble,
    eInt8,
    eInt16,
    eInt32,
    eInt64,
    eUint8,
    eUint16,
    eUint32,
    eUint64,
    eGt8
};

class scidb_writer {
public:
    scidb_writer(std::string const& filename);
    virtual ~scidb_writer();

    void set_chrom(std::string const& chrom)  { _chrom = chrom; }
    void set_pos(std::string const& pos) { _pos = pos; }
    virtual void set_var(int64_t var)=0;
    
    virtual void put_prefix()=0;
    virtual void put_separator()=0;
    virtual void put_endrow()=0;

    virtual void put_data(std::string const& data, ENullData nullval, EDataType type)=0;
    virtual void put_uint32(uint32_t data)=0;
    virtual void put_int64(int64_t data)=0;

protected:
    int _out;
    std::string _chrom;
    std::string _pos;
    std::string _prefix;
};

typedef std::map<std::string, int64_t> subjectMap_t;
# define YY_DECL int yylex(scidb_writer& var_writer, scidb_writer& gt_writer, subjectMap_t& subjMap, size_t max_ref_size)
YY_DECL;

class scidb_text_writer: public scidb_writer {
public:
    scidb_text_writer(std::string const& filename, size_t chunksize);
    virtual ~scidb_text_writer();
    virtual void set_var(int64_t var);
    
    virtual void put_prefix();
    virtual void put_separator();
    virtual void put_endrow();

    virtual void put_data(std::string const& data, ENullData nullval, EDataType type);
    virtual void put_uint32(uint32_t data);
    virtual void put_int64(int64_t data);

private:
    size_t _chunksize;
    size_t _rowcount;
    bool _newchunk;
};

class scidb_binary_writer: public scidb_writer {
public:
    scidb_binary_writer(std::string const& filename);
    virtual ~scidb_binary_writer();
    virtual void set_var(int64_t var);
    
    virtual void put_prefix();
    virtual void put_separator();
    virtual void put_endrow();

    virtual void put_data(std::string const& data, ENullData nullval, EDataType type);
    virtual void put_uint32(uint32_t data);
    virtual void put_int64(int64_t data);

private:
    std::string _prefix;
};

#endif // ! SCIDB_WRITERS_HPP
