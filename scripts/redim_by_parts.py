#!/opt/python-2.7/bin/python
##!/usr/bin/python
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
import argparse
import traceback
import datetime
import scidbutils as sciutil
import scidbapi as scidb


def handleException(inst, exitWhenDone, op=None):
    traceback.print_exc()
    if op:
        print >> sys.stderr, "Exception while ", op
    print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
    print >> sys.stderr, "     Exception Value: %r" % inst 
    print >> sys.stderr, ""
    if(exitWhenDone):
        exit(2)

def get_chunksize(db, array, dim):
    chunksize = sciutil.get_column(db, "project(dimensions(%s),chunk_interval)" % array)[dim]
    return chunksize

def main():
    parser = argparse.ArgumentParser(description='Redimension a large array by parts')
    parser.add_argument('array', help='Base array name')
    parser.add_argument('-c', '--host', help='SciDB coordinator host (Default: localhost)', default='scidb10')
    parser.add_argument('-s', '--size', help='size of chunks to redimension at a time', type=long, default=1000000000)
    parser.add_argument('-b', '--begin', help='row to begin redimension, useful to restart process', type=long, default=0)
    args = parser.parse_args()

    try:
        db = scidb.connect(args.host, 1239)
    except Exception, inst: 
        handleException(inst, True, op="connecting")

    load_array = args.array + "_load"

    start = args.begin
    stop = start + args.size-1
    between = "between(%s,%s,%s)" % (load_array, start, stop)
    items = long(sciutil.get_single_result(db, "count(%s)" % between))
    total = items
    while (items > 0):
        redim = "insert(redimension(%s,%s),%s)" % (between, args.array, args.array)
        sciutil.do_query(db, redim)

        print "range:", start, "-", stop,
        start = stop+1
        stop = stop + args.size
        between = "between(%s,%s,%s)" % (load_array, start, stop)
        items = long(sciutil.get_single_result(db, "count(%s)" % between))
        total += items
        print ",  total:", total, "\r",

    sys.exit(0) #success

if __name__ == "__main__":
    main()
