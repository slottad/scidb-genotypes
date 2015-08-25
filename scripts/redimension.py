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
import traceback
import datetime
from ScidbQuery import ScidbQuery
import gtUtils


def handleException(inst, exitWhenDone, op=None):
    traceback.print_exc()
    if op:
        print >> sys.stderr, "Exception while ", op
    print >> sys.stderr, "     Exception Type: %s" % type(inst)     # the exception instance
    print >> sys.stderr, "     Exception Value: %r" % inst 
    print >> sys.stderr, ""
    if(exitWhenDone):
        exit(2)


def redim_by_parts(array_size):
    log.info('redimensioning by parts - size: %s' % array_size)

    i=0
    j=args.redim_max
    count = 1
    while i < int(array_size):
        redim = """
        insert(
            redimension(
                index_lookup(
                    between({base}_gt_load, {i}, {j}),
                    project({base}_chroms,chrom), 
                    {base}_gt_load.chrom,chromid), 
                {base}_gt), 
            {base}_gt)
        """.format(base = args.array, i = i, j = j)

        query_obj.serverActionOnly(redim)
        i = j + 1
        j += args.redim_max

        if (count % 10) > 0:
            count += 1
            continue
        else:
            # remove earlier array versions to avoid slowing
            latest_version = query_obj.getRecords("aggregate(versions(%s_gt), max(version_id))" % args.array)[0].version_id_max
            log.info('removing gt array versions earlier than version: %s' % latest_version)
            query_obj.serverActionOnly('remove_versions(%s_gt,%s)' % (args.array, latest_version))
            count += 1


def redim_all_at_once():
    log.debug('redimensioning gt array all at once')

    gt_redim = """
        store(
            redimension(
                index_lookup( 
                    {base}_gt_load, 
                    project({base}_chroms, chrom),
                    {base}_gt_load.chrom,chromid),
                {base}_gt), 
            {base}_gt)
    """.format(base = args.array)

    query_obj.serverActionOnly(gt_redim)


def main():
    parser = gtUtils.argparser('host', 'port', 'log_level', description='Load a study to SciDB')
    parser.add_argument('array', help='Base array name')
    parser.add_argument('--nogt', help='skip redimension of the gt array', action='store_true')
    parser.add_argument('--novar', help='skip redimension of the var array', action='store_true')
    parser.add_argument('--redim_max', type=int, help='threshhold for redimensioning by parts', default=1000000000)
    global args; args = parser.parse_args()


    global query_obj; query_obj = ScidbQuery(server=args.host, port=args.port, log_level=args.log_level)
    global log; log = gtUtils.logger(args.log_level)

    # find chunksizes and array limits
    #         min, max, chunk
    # chrom  0, _chroms, 1
    # pos    1, max(_var_load,pos), 200000
    # var    1, max(_var_load,var), max(_var_load,var)
    # sample 0, _samples, _samples

    chrom_high = query_obj.getRecords("project(dimensions(%s_chroms),high)" % args.array)[0].high
    log.debug('chrom_high: %s' % chrom_high)

    pos_high = query_obj.getRecords("aggregate(%s_var_load,max(pos))" % args.array)[0].pos_max
    log.debug('pos_high: %s' % pos_high)

    var_high = query_obj.getRecords("aggregate(%s_var_load,max(var))" % args.array)[0].var_max
    log.debug('var_high: %s' % var_high)

    sample_high = query_obj.getRecords("project(dimensions(%s_samples),high)" % args.array)[0].high
    log.debug('sample_high: %s' % sample_high)

    query = "project(dimensions(%s_samples),chunk_interval)" % args.array
    sample_chunksize = query_obj.getRecords(query)[0].chunk_interval
    log.debug('sample_chunksize: %s' % sample_chunksize)


    # create new arrays
    exists = query_obj.getRecords('show(%s_var)' % args.array)
    if len(exists) > 0:
        log.info('var array exists - will not recreate')
    else:
        var_def = """
        create array {base}_var<id:string null,ref:string,alt:string,alleles:uint32,
        qual:float null,filter:string null,info:string null,format:string null> 
        [chromid=0:{chrom_high},1,0, pos=1:{pos_high},200000,0, var=1:{var_high},{var_high},0]
        """.format(base = args.array, chrom_high = chrom_high, pos_high = pos_high, var_high = var_high)

        query_obj.serverActionOnly(var_def) 

    exists = query_obj.getRecords('show(%s_gt)' % args.array)
    if len(exists) > 0:
        log.info('gt array exists - will not recreate')
    else:
        gt_def = """
        create array {base}_gt <gt:gt16 null,unparsed:string null> 
        [chromid=0:{chrom_high},1,0, 
         pos=1:{pos_high},200000,0, 
         var=1:{var_high},{var_high},0, 
         sampleid=0:{sample_high},{sample_chunksize},0]
        """.format(base = args.array, chrom_high = chrom_high, pos_high = pos_high, 
            var_high = var_high, sample_high = sample_high, sample_chunksize = sample_chunksize)

        query_obj.serverActionOnly(gt_def) 

    # redimension _var_load
    if not args.novar:
        log.info('redimensioning var array')
        var_redim = """
        store(
            redimension(
                index_lookup( 
                    {base}_var_load, 
                    project({base}_chroms, chrom), 
                    {base}_var_load.chrom, 
                    chromid),
                {base}_var),
            {base}_var)
        """.format(base = args.array)

        tstart = datetime.datetime.now()
        query_obj.serverActionOnly(var_redim)
        tstop = datetime.datetime.now()
        tdiff = tstop - tstart
        log.info("finished variation array redimension - time: %s" % tdiff)


    # redimension _gt_load
    if not args.nogt:
        log.info('redimensioning genotype array')

        # large arrays get redimensioned in pieces
        array_size = query_obj.getRecords("aggregate(%s_gt_load,count(pos))" % args.array)[0].pos_count
        tstart = datetime.datetime.now()

        if array_size <= args.redim_max:
            redim_all_at_once()
        else:
            redim_by_parts(array_size)

        tstop = datetime.datetime.now()
        tdiff = tstop - tstart
        log.info("finished genotype array redimension - time: %s" % tdiff)

    sys.exit(0) #success

if __name__ == "__main__":
    main()
