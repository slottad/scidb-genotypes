#!/usr/bin/python
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
sys.path.append(os.getcwd())
sys.path.append('/opt/scidb/13.2' + '/lib')
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

def do_query(db, query):
    result = db.executeQuery(query, "afl")
    db.completeQuery(result.queryID)
    return

def get_single_result(db, query):
    result = db.executeQuery(query, "afl")
    desc = result.array.getArrayDesc()
    attrs = desc.getAttributes()
    attrid = attrs[0].getId()
    attiter = result.array.getConstIterator(attrid)
    chunk = attiter.getChunk()
    chunkiter = chunk.getConstIterator((scidb.swig.ConstChunkIterator.IGNORE_OVERLAPS)|(scidb.swig.ConstChunkIterator.IGNORE_EMPTY_CELLS))
    if not chunkiter.isEmpty(): 
        dataitem = chunkiter.getItem()
        val = scidb.getTypedValue(dataitem, attrs[0].getType())
    else:
        val = "#Error!"
    db.completeQuery(result.queryID)
    return val

def get_column(db, query):
    values = []
    result = db.executeQuery(query, "afl")
    desc = result.array.getArrayDesc()
    attrs = desc.getAttributes()
    attrid = attrs[0].getId()
    attiter = result.array.getConstIterator(attrid)
    while not attiter.end():
        chunk = attiter.getChunk()
        chunkiter = chunk.getConstIterator((scidb.swig.ConstChunkIterator.IGNORE_OVERLAPS)|(scidb.swig.ConstChunkIterator.IGNORE_EMPTY_CELLS))
        while not chunkiter.end():
            if not chunkiter.isEmpty(): 
                dataitem = chunkiter.getItem()
                val = scidb.getTypedValue(dataitem, attrs[0].getType())
                values.append(val)
            else:
                val = "#Error!"
            chunkiter.increment_to_next()
        attiter.increment_to_next()

    db.completeQuery(result.queryID)
    return values

def main():
    parser = argparse.ArgumentParser(description='Compute the allele counts for the given array.')
    parser.add_argument('array', help='The SciDB base array name')
    args = parser.parse_args()

    try: 
        db = scidb.connect("localhost", 1239)
    except Exception, inst: 
        handleException(inst, True, op="connecting")

    max_alleles = get_single_result(db, "aggregate(%s,max(alleles))" % args.array)    
    #print max_alleles

    count = get_single_result(db, "count(%s)" % args.array)    
    #print count

    populations = list(set(get_column(db, "project(%s_samples,population)" % args.array)))
    populations.sort()
    populations.append("global")
    print populations

    idx_base = 0
    idx_max = count * len(populations) * (max_alleles+2) - 1
    loading_array = "create array %s_allele_counts_load<chrom:string,pos:int64,var:int64,count:uint64 null,population:string,allele:int64> [idx=0:*,1000000,0]" % args.array
    do_query(db,loading_array)

    sizes = get_column(db, "project(dimensions(%s),length)" % args.array)
    chunks = get_column(db, "project(dimensions(%s),chunk_interval)" % args.array)

    final_array = "create array %s_allele_counts<count:uint64 null> [chrom(string)=%s,1,0,pos=1:%s,%s,0,var=1:%s,%s,0,population(string)=%s,1,0,allele=-1:%s,%s,0]" % (args.array, sizes[0], sizes[1], chunks[1], sizes[2], chunks[2], len(populations), max_alleles, max_alleles+2)
    #print final_array
    do_query(db,final_array)

    for population in populations:
        if population == "global":
            pop_array = "project(%s_gt,gt)" % args.array
        else:
            pop_array = "cross_join(project(%s_gt,gt),filter(%s_samples,population='%s'),sampleid,row)" % (args.array, args.array, population)

        complete_alleles = []
        for allele in range(-1,max_alleles+1):
            if (allele < 0):
                allele_query = "ploidy(gt)"
            else:
                allele_query = "allele_count(gt,%s)" % allele

            calc_query = "store(apply(unpack(apply(attribute_rename(aggregate(apply(%s,ac,%s),sum(ac),chrom,pos,var),ac_sum,count),population,'%s',allele,%s),row),idx,row+%s),%s_temp1)" % (pop_array, allele_query, population, allele, idx_base, args.array)
            do_query(db,calc_query)

            #bounded_array_query = "create array %s_temp2<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [idx=%s:%s,1000000,0]" % (args.array, idx_base, idx_base+count-1)
            #do_query(db,bounded_array_query)
            staging_array_query = "create array %s_temp2<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [idx=0:*,1000000,0]" % args.array
            do_query(db,staging_array_query)
            
            redim_query = "redimension_store(%s_temp1,%s_temp2)" % (args.array, args.array)
            do_query(db,redim_query)

            insert_query = "insert(%s_temp2,%s_allele_counts_load)" % (args.array, args.array)
            do_query(db,insert_query)

            do_query(db,"remove(%s_temp1)" % args.array)
            do_query(db,"remove(%s_temp2)" % args.array)

            complete_alleles.append(allele)
            print population, "\t", complete_alleles, "\r",
            sys.stdout.flush()
            idx_base += count

        print population, "\t", complete_alleles, " complete"

    do_query(db, "redimension_store(%s_allele_counts_load,%s_allele_counts)" % (args.array, args.array))
    db.disconnect()     #Disconnect from the SciDB server.

    sys.exit(0) #success


if __name__ == "__main__":
    main()


# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(unpack(apply(attribute_rename(aggregate(apply(cross_join(project(interimV3_gt,gt),filter(interimV3_samples,population='CHB'),sampleid,row),ac,ploidy(gt)),sum(ac),chrom,pos,var),ac_sum,count),population,'CHB',allele,-1),row),testG)"
# Query was executed successfully

# real    8m20.829s
# user    0m0.018s
# sys     0m0.008s
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -anq "create array test_stageG<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [row=0:39697779,1000000,0]"
 
# Query was executed successfully
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "redimension_store(testG,test_stageG)"                                                                                                          
# Query was executed successfully

# real    0m9.629s
# user    0m0.017s
# sys     0m0.005s
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(concat(test_stageG,concat(test_stage0,test_stage1)),test_merged_all)"
# ^C
# real    19m42.395s
# user    0m0.017s
# sys     0m0.007s
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(concat(test_stageG,store(concat(test_stage0,test_stage1),tmp)),test_merged_all)"
# Query was executed successfully

# real    3m9.178s
# user    0m0.016s
# sys     0m0.012s
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -q "count(tmp)"
# i,count
# 0,79395560
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -q "count(test_merged_all)"                                                                                                                               
# i,count
# 0,119093340
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -anq "create array test_stageG<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [row=0:39697779,1000000,0]" 
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(unpack(apply(attribute_rename(aggregate(apply(cross_join(project(interimV3_gt,gt),filter(interimV3_samples,population='CHB'),sampleid,row),ac,ploidy(gt)),sum(ac),chrom,pos,var),ac_sum,count),population,'CHB',allele,-1),row),testG)"
# Query was executed successfully

# real    8m20.829s
# user    0m0.018s
# sys     0m0.008s
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -anq "create array test_stageG<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [row=0:39697779,1000000,0]"
 
# Query was executed successfully
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "redimension_store(testG,test_stageG)"                                                                                                          
# Query was executed successfully

# real    0m9.629s
# user    0m0.017s
# sys     0m0.005s
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(concat(test_stageG,concat(test_stage0,test_stage1)),test_merged_all)"
# ^C
# real    19m42.395s
# user    0m0.017s
# sys     0m0.007s
# [gtdev11]/netmnt/gtbridge/vcf_data$ time iquery -anq "store(concat(test_stageG,store(concat(test_stage0,test_stage1),tmp)),test_merged_all)"
# Query was executed successfully

# real    3m9.178s
# user    0m0.016s
# sys     0m0.012s
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -q "count(tmp)"
# i,count
# 0,79395560
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -q "count(test_merged_all)"                                                                                                                               
# i,count
# 0,119093340
# [gtdev11]/netmnt/gtbridge/vcf_data$ iquery -anq "create array test_stageG<chrom:string,pos:int64,var:int64,count:uint64 NULL DEFAULT null,population:string,allele:int64> [row=0:39697779,1000000,0]" 

