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
sys.path.append(os.getcwd())
sys.path.append('/opt/scidb/14.7' + '/lib')
import scidbapi as scidb

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
