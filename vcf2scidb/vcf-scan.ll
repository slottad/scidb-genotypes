/*  $Id: vcf-scan.ll 398471 2013-05-06 18:54:57Z slottad $
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
 *   Flex scanner class for vcf files
 *
 */

%{
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp> 
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "scidb-writers.hpp"
using namespace std;
using namespace boost;
size_t colnum = 0;
int64_t cur_var = 0;
string prev_chrom;
string cur_chrom;
string prev_pos;
string cur_pos;
string cur_id;
size_t loc;
int64_t cur_sample;
vector<int64_t> sampleids;
bool skip_row;
set<string> chrom_set;
int64_t max_pos;
int64_t max_var;
int64_t max_sampleid;
%}

%option noyywrap nounput batch

mdline    ^##.*
header    ^#[cC][hH][rR][oO][mM]\t.*
tab       \t
chrom     ^[^#\t\n]+
datum     [^\t\n]+
eol       \n

%%

{mdline}  //cout << "mdline";
{header} {
    vector<string> samples;
    string yystr(yytext);
    split(samples, yystr, is_any_of("\t"));
    samples.erase(samples.begin(), samples.begin()+9);
    sampleids.assign(samples.size(), 0);
    subjectMap_t::iterator smi;
    size_t idx=0;
    for (vector<string>::iterator i=samples.begin(); i != samples.end(); ++i, ++idx) {
        smi = subjMap.find(*i);
        if (smi == subjMap.end()) {
            subjMap[*i] = subjMap.size()+1;
            sampleids[idx] = subjMap.size();
        } else {
            sampleids[idx] = smi->second;
        }
        if (sampleids[idx] > max_sampleid) max_sampleid = sampleids[idx];
        //cout << idx << "<-" << sampleids[idx] << endl;
    }
}
{chrom} {
    colnum = 1;
    cur_chrom = yytext;
    var_writer.set_chrom(cur_chrom);
    gt_writer.set_chrom(cur_chrom);
    chrom_set.insert(cur_chrom);
    skip_row = false;
}
{tab} { ++colnum; }
{datum} {
    string yystr(yytext);
    if (!skip_row) {
        switch (colnum) {
        case 2:  // POS and VAR
        {
            cur_pos = yystr;
            var_writer.set_pos(cur_pos);
            gt_writer.set_pos(cur_pos);
            int64_t pos = atol(cur_pos.c_str());
            if (pos > max_pos) max_pos = pos;
            if ((cur_pos == prev_pos) && (cur_chrom == prev_chrom)) {
                ++cur_var;
            } else {
                prev_chrom = cur_chrom;
                prev_pos = cur_pos;
                cur_var = 1;
            }
            var_writer.set_var(cur_var);
            gt_writer.set_var(cur_var);
            if (cur_var > max_var) max_var = cur_var;
            break;
        }
        case 3: // ID
            if (yystr == ".") yystr.clear();
            cur_id = yystr;
            break;
        case 4: // REF
            if (yystr.size() > max_ref_size) {
                skip_row = true;
            } else {
                var_writer.put_prefix(); 
                var_writer.put_data(cur_id, eNullable, eString);
                var_writer.put_separator();
                var_writer.put_data(yystr, eNotNullable, eString);            
            }
            break;
        case 5: // ALT and Alleles
        {
            var_writer.put_separator();
            var_writer.put_data(yystr, eNotNullable, eString);
            var_writer.put_separator();
            uint32_t alleles = 2 + (uint32_t)count(yystr.begin(),yystr.end(),',');
            var_writer.put_uint32(alleles);
        }
        break;
        case 6: // QUAL
            if (yystr == ".") yystr.clear();
            var_writer.put_separator();
            var_writer.put_data(yystr, eNullable, eFloat);
            break;
        case 7: // FILTER
            if (yystr == ".") yystr.clear();
            var_writer.put_separator();
            var_writer.put_data(yystr, eNullable, eString);
            break;
        case 8: // INFO
            if (yystr == ".") yystr.clear();
            var_writer.put_separator();
            var_writer.put_data(yystr, eNullable, eString);
            break;
        case 9: // FORMAT
            if (yystr == ".") yystr.clear();
            var_writer.put_separator();
            if (yystr.substr(0,2) == "GT") {
                yystr.erase(0,2);
            }
            if (yystr[0] == ':') {
                yystr.erase(0,1);                
            }
            var_writer.put_data(yystr, eNullable, eString);
            var_writer.put_endrow();
            break;
        default: // SAMPLES
        {
            cur_sample=colnum-10;
            string gt;
            string rest;
            size_t div = yystr.find(':');
            if (div == string::npos) {
                gt = yystr;
                rest.clear();
            } else {
                gt = yystr.substr(0,div);
                rest = yystr.substr(div+1, string::npos);
            }
            if ((gt != "./.") && (gt != ".|.")) {
                gt_writer.put_prefix();
                gt_writer.put_int64(sampleids[cur_sample]);
                gt_writer.put_separator();
                gt_writer.put_data(gt, eNotNullable, eGt8);
                gt_writer.put_separator();
                gt_writer.put_data(rest, eNullable, eString);
                gt_writer.put_endrow();
            }
            break;
        }
        }
    }
}
{eol} {
    // mylineno++;
    // cout << endl;
}

