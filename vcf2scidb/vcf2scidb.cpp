/*  $Id: vcf2scidb.cpp 394139 2013-03-30 03:26:57Z slottad $
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
 *   Utility to convert VCF to a format suitable for loading into SciDB
 *
 */

// Standard includes
#include <iostream>
#include <memory>
#include <fstream>
#include <vector>

// Boost
#include <boost/assign.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>

#include "scidb-writers.hpp"

using namespace std;
using namespace boost;
using namespace boost::program_options;
using namespace boost::algorithm;

extern set<string> chrom_set;
extern int64_t max_pos;
extern int64_t max_var;
extern int64_t max_sampleid;

void read_descriptions(string const& filename, subjectMap_t& subjectMap)
{
    ifstream ifs(filename.c_str(), ifstream::in);
    string line;
    while (getline(ifs,line)) {
        if (line[0] == '#') continue;
        vector<string> cells;
        split(cells, line, is_any_of(","));
        int64_t index = lexical_cast<int64_t>(cells[0]);
        string subject = cells[1];
        subjectMap.insert(make_pair(subject, index));
    }
}

int main( int argc, char** argv)
{
    options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "view help message, then exit")
        ("text,t", "use SciDB text format")
        ("binary,b", "use SciDB binary format")
        ("descriptions,d", value<string>(), "Input a CSV file listing info about the samples")
        ("chunk,c", value<size_t>()->default_value(100000), "loading array chunk size")
        ("var,v", value<string>()->default_value("array_var.scidb"), "variation array output file")
        ("gt,g", value<string>()->default_value("array_gt.scidb"), "genotype array output file")
        ("maxref,m", value<size_t>()->default_value(50), "maximum length ")
    ;
    
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);
    notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    subjectMap_t subjMap;
    if (vm.count("descriptions")) {
        read_descriptions(vm["descriptions"].as<string>(), subjMap);
    }

    size_t chunksize = vm["chunk"].as<size_t>();
    size_t maxref = vm["maxref"].as<size_t>();
    string varfile = vm["var"].as<string>();
    string gtfile = vm["gt"].as<string>();

    unique_ptr<scidb_writer> var_writer;
    unique_ptr<scidb_writer> gt_writer;

    if (vm.count("binary")) {
        var_writer = unique_ptr<scidb_writer>(new scidb_binary_writer(varfile));
        gt_writer = unique_ptr<scidb_writer>(new scidb_binary_writer(gtfile));
    } else {
        var_writer = unique_ptr<scidb_writer>(new scidb_text_writer(varfile, chunksize));
        gt_writer = unique_ptr<scidb_writer>(new scidb_text_writer(gtfile, chunksize));        
    }

    max_var = 0;
    max_pos = 0;
    max_sampleid = 0;
    yylex(*var_writer, *gt_writer, subjMap, maxref);
    cout << chrom_set.size() << " ";
    cout << max_pos << " ";
    cout << max_var << " ";
    cout << max_sampleid << endl;
    return 0;
}
