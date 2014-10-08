/*  $Id: vcf2csv.cpp 446804 2014-09-18 13:44:05Z slottad $
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
 *   Utility to convert VCF to two different csv files, suitable for loading into SciDB
 *   Code structure shamelessly stolen from SciDB's splitcsv.c
 *   N.B. This actually uses tabs instead of commas, because SciDB splitcsv cannot otherwise escape commas.
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

// Standard includes
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <map>

// Is 10MB a large enough buffer for a VCF line?
// One hopes, but VCF is pathological
#define RB_SIZE 10485760

using namespace std;

char _line[RB_SIZE];

char* _inputFileName = NULL;
char* _inputSamplesName = NULL;
char* _outputVarName = NULL;
char* _outputGtName = NULL;

FILE* _inputFile = NULL;
FILE* _varFile = NULL;
FILE* _gtFile = NULL;

typedef map<string, string> sampleMap_t;
sampleMap_t _sampleMap;
vector<string> _samples;

string _prevChromPos;
int _curVar = 1;
string _prefix;

void usage()
{
    printf("Utility to split a VCF file into two CSV files.\n"
           "USAGE: vcf2csv <-s SAMPLES> [-i INPUT] file1 file2\n"
           "\t-s SAMPLES\tName of file containing sample descriptions. (REQUIRED)\n"
           "\t-i INPUT\tInput file. (Default = stdin).\n");
}

void haltOnError(const char* errStr)
{
    fprintf(stderr, "ERROR: %s\n", errStr);
    usage();
    exit(EXIT_FAILURE);
}

void parseArgs(int argc, char* argv[]) 
{
    if (argc < 4) {
        haltOnError("Missing some required arguments.\n");
    }

    /* Iterate over the command-line arguments. */
    int i = 0;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            _inputFileName = argv[++i];
        } else if (strcmp(argv[i], "-s") == 0) {
            _inputSamplesName = argv[++i];
        } else 
            break;
    }
    if (_inputSamplesName == NULL) 
        haltOnError("[-s samplefile] is missing, but it is a required argument.");

    if (i >= argc) haltOnError("Missing variation output filename.\n");
    _outputVarName = argv[i];
    if (++i >= argc) haltOnError("Missing genotypes output filename.\n");
    _outputGtName = argv[i];
}

void loadSamples()
{
    FILE* f = fopen(_inputSamplesName, "r");
    if (f == NULL) {
        haltOnError("Failed to open samples file.");
    }
    
    size_t sampleIdx=0;
    size_t sampleCol=0;
    _sampleMap.clear();

    // Read the header
    if (fgets(_line, RB_SIZE, f) == NULL) {
        haltOnError("The sample descriptions file is empty.");
    }
    char* header = strtok(_line, "\t,");
    while ( (header != NULL) && (strcmp(header, "sample") != 0) ) {
        header = strtok(NULL, "\t,");
        ++sampleCol;
    }
    if (header == NULL) {
        haltOnError("The sample descriptions file does not contain a 'sample' column.");
    }

    // Read the samples
    while (fgets(_line, RB_SIZE, f) != NULL) {
        size_t col = 0;
        char* tok = strtok(_line, "\t,");
        while (col < sampleCol) {
            tok = strtok(NULL, "\t,");
            ++col;
        }
        string index(to_string(sampleIdx));
        string sample(tok);
        _sampleMap.insert(make_pair(sample, index));
        //fprintf(stderr, "%s:\t%s\n", sample.c_str(), index.c_str());
        ++sampleIdx;
    }
}


void openFiles()
{
    if (_inputFileName == NULL) {
        _inputFile = stdin;
    } else {
        _inputFile = fopen(_inputFileName, "r");
        if (_inputFile == NULL) {
            haltOnError("Failed to open specified VCF input file.");
        }
    }
    // Resize input buffer
    //setvbuf(_inputFile, NULL, _IOFBF, RB_SIZE*2);

    _varFile = fopen(_outputVarName, "w");
    if (_varFile == NULL) {
        haltOnError("Failed to open variation output file.");
    }

    _gtFile = fopen(_outputGtName, "w");
    if (_gtFile == NULL) {
        haltOnError("Failed to open genotype output file.");
    }
}

void closeFile(FILE* file)
{
    if (file != NULL) {
        fclose(file);
    }
    file = NULL;
}

void closeFiles()
{
    closeFile(_inputFile);
    closeFile(_varFile);
    closeFile(_gtFile);
}

void parseHeader(char* head)
{
    _samples.clear();
    char* pTok = strtok(head, "\t\n");

    size_t col = 0;
    sampleMap_t::iterator smi;
    
    while (pTok != NULL) {
        if (++col > 9) {
            string tok(pTok);
            smi = _sampleMap.find(tok);
            if (smi == _sampleMap.end()) {
                string errStr = "Unknown sample found: " + tok;
                fprintf(stderr, "ERROR: %s\n", errStr.c_str());
                _samples.push_back(string());
            } else {
                _samples.push_back(smi->second);
            }
        }
        pTok = strtok(NULL, "\t\n");
    }
}

size_t allele_count(char* tok) 
{
    size_t count = 2;
    for (size_t i=0; tok[i] != '\0'; i++) {
        if (tok[i] == ',') ++count;
    }
    return count;
}

char* parse_format(char* fmt)
{
    if (fmt[0] == 'G' && fmt[1] == 'T') {
        if (fmt[2] == '\0') return &fmt[2];
        if (fmt[2] == ':') fmt = &fmt[3];
    }
    return fmt;
}

void parseLine(char* line) {
    char* chrom = strtok(line, "\t");
    strtok(NULL, "\t"); // pos

    chrom[strlen(chrom)] = '\t';
    string chromPos(chrom);
    if (chromPos == _prevChromPos) {
        ++_curVar;
    } else {
        _curVar = 1;
    }
    
    _prevChromPos = chromPos;
    _prefix = chromPos + "\t" + to_string(_curVar) + "\t";

    fprintf(_varFile, "%s", _prefix.c_str());

    char* id = strtok(NULL, "\t");
    char* ref = strtok(NULL, "\t");
    char* alt = strtok(NULL, "\t");
    size_t alleles = allele_count(alt);
    
    fprintf(_varFile, "%s\t%s\t%s\t%s\t", id, ref, alt, to_string(alleles).c_str());

    char* qual = strtok(NULL, "\t");
    char* filter = strtok(NULL, "\t");
    char* info = strtok(NULL, "\t\n");
    char* format = strtok(NULL, "\t\n");
    
    if (strcmp(qual,".") == 0) qual[0]='\0';
    if (strcmp(filter,".") == 0) filter[0]='\0';
    if (strcmp(info,".") == 0) info[0]='\0';
    char* new_format = parse_format(format);
    bool has_gt(new_format != format);
   
    fprintf(_varFile, "%s\t%s\t%s\t%s\n", qual, filter, info, new_format);

    size_t idx = 0;
    char* gt = strtok(NULL, "\t");
    while (gt != NULL) {
		if (idx >= _samples.size())
		{
			fprintf(stderr, "gt-index exceeds _samples.size at: chrom_pos=%s, idx=%lu, gt=%s, _samples.size=%lu\n",_prefix.c_str(),idx, gt, _samples.size());
			exit(EXIT_FAILURE);
		}
        if (! _samples[idx].empty() ) {

            fprintf(_gtFile, "%s%s\t", _prefix.c_str(), _samples[idx].c_str());

            if (has_gt) {
                char* pColon = strchr(gt, ':');
                if (pColon) {
                    *pColon = '\t';
                    fprintf(_gtFile, "%s\n", gt);
                } else {
                    fprintf(_gtFile, "%s\t\n", gt);
                }
            } else {
                fprintf(_gtFile, "\t%s\n", gt);
            }
        }
        gt = strtok(NULL, "\t\n");
        ++idx;
    }
}

int main(int argc, char* argv[]) {
    parseArgs(argc, argv);
    loadSamples();
    openFiles();
    _prevChromPos.clear();
    int row = 0;
    while (fgets(_line, RB_SIZE, _inputFile) != NULL) {
        if (strlen(_line) == 1) {
            continue;
        } else if (_line[0] == '#') {
            if (_line[1] == '#') continue;
            parseHeader(_line);
        } else {
            fprintf(stderr, "%d\n", ++row);
            parseLine(_line);
        }
    }
    fprintf(stderr, "\n");
    closeFiles();
    exit(EXIT_SUCCESS);
}
