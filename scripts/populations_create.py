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
import csv
import argparse
import traceback
import shlex
import subprocess as sp
from ScidbQuery import ScidbQuery
import logging as log

log.basicConfig(level='INFO', format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


def handleException(inst, exitWhenDone, op=None):
    traceback.print_exc()
    if op:
        print >> sys.stderr, "Exception while ", op
    print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
    print >> sys.stderr, "     Exception Value: %r" % inst 
    print >> sys.stderr, ""
    if(exitWhenDone):
        exit(2)

def construct_echo(populations):
    if 'global' in populations:
        s = "/bin/echo -e '"
    else:
        s = "/bin/echo -e 'global"
    for pop in populations:
        s += "\\n"
        s += pop
    s += "'"
    return s

def get_populations(baseArray):
    records = query_obj.getRecords('uniq(sort(project(%s_samples,population)))' % baseArray)
    populations = []

    for r in records:
        log.debug('population: %s' % r.population)
        populations.append(r.population)

    if len(populations) < 1:
        raise Exception('no populations found in scidb')

    return populations

def main():
    parser = argparse.ArgumentParser(description='Load a TSV or CSV file to SciDB')
    parser.add_argument('array', help='Base array name')
    parser.add_argument('-c', '--host', help='SciDB coordinator host (Default: localhost)', default='localhost')
    parser.add_argument('-p', '--port', help='SciDB host port', default=1239, type=int)
    parser.add_argument('--log_level', help='log output level', type=str, \
        choices=['debug', 'info', 'warning', 'error', 'critical'], default='info')
    args = parser.parse_args()

    global query_obj; query_obj = ScidbQuery(server=args.host, port=args.port, log_level=args.log_level)

    logger = log.getLogger()
    logger.setLevel(args.log_level.upper())

    populations = get_populations(args.array)
    
    ver = os.environ['SCIDB_VER']

    echo_statement = construct_echo(populations)
    loadpops = "%s | /opt/scidb/%s/bin/loadcsv.py -p %d -x -q -d %s -a %s_populations -s'<population:string>[id=0:%d,%d,0]'" % (echo_statement, ver, args.port, args.host, args.array, len(populations), len(populations)+1)
    os.system(loadpops)
    
    sys.exit(0) #success

if __name__ == "__main__":
    main()
