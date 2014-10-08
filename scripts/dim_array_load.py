#!/opt/python-2.7/bin/python
################################################################################
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Author:  Douglas Slotta
#
################################################################################
import os
import sys
import csv
import argparse
import traceback
import logging as log


def handleException(inst, exitWhenDone, op=None):
    traceback.print_exc()
    if op:
        print >> sys.stderr, "Exception while ", op
    print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
    print >> sys.stderr, "     Exception Value: %r" % inst 
    print >> sys.stderr, ""
    if(exitWhenDone):
        exit(2)

def create_array_def(header, rows, chunksize):
    patterndef = "-t "
    arraydef = "-s '<"
    first = True

    for col in header:
        if first:
            first = False
        else:
            arraydef += ","

        arraydef += col
        if col == 'founder':
            patterndef += "N"
            arraydef += ":bool"
        elif col == 'sex':
            patterndef += "C"
            arraydef += ":char"
        else:
            patterndef += "S"
            arraydef += ":string"

    arraydef += ">[id=0:%s,%s,0]'" % ( rows, chunksize )
    return patterndef + " " + arraydef

def parse_header(filename):
    with open(filename, 'rb') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        header = reader.next()
        if dialect.delimiter == ',':
            delim = "','"
        else:
            delim = "'\\t'"
        return (header, delim)

def get_row_count(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i

def get_chunksize(dimension, rows):
    if dimension == "chrom":
        return 1
    elif dimension == "sample":
        if rows <= 10:
            return rows
        elif rows <= 20:
            return int(rows/2)
        elif rows <= 50:
            return int(rows/3)
        else:
            return int(rows/4)
    else:
        print "The %s dimension is not supported" % dimension
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Load a TSV or CSV file to SciDB')
    parser.add_argument('file', help='filename')
    parser.add_argument('-a', '--array', help='Base array name', required=True)
    parser.add_argument('-d', '--dimension', help='Dimension name [ chrom | sample ]', required=True)
    parser.add_argument('-c', '--host', help='SciDB coordinator host (Default: localhost)', default='localhost')
    parser.add_argument('-p', '--port', help='SciDB host port', type=int, default=1239)
    parser.add_argument('--log_level', help='log output level', type=str, \
        choices=['debug', 'info', 'warning', 'error', 'critical'], default='info')
    args = parser.parse_args()

    log.basicConfig(level='INFO', format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logger = log.getLogger()
    logger.setLevel(args.log_level.upper())

    nrows = get_row_count(args.file)
    chunksize = get_chunksize(args.dimension, nrows)
    (header, delimiter) = parse_header(args.file)
    arrdef = create_array_def(header, nrows, chunksize)
    ver = os.environ['SCIDB_VER']

    loadcmd = "/opt/scidb/%s/bin/loadcsv.py  -x -q -n 1 -d %s -D %s -a %s_%ss %s -i %s -p %s" \
        % (ver, args.host, delimiter, args.array, args.dimension, arrdef, args.file, str(args.port))

    log.info('calling: %s' % loadcmd)
    os.system(loadcmd)

    sys.exit(0) #success

if __name__ == "__main__":
    main()
