#/bin/bash
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

if [ "$2" = "" ]; then
    echo Usage: $0 array filename 
    exit 1
fi

IQUERY=/opt/scidb/$SCIDB_VER/bin/iquery

FILE=$2
FILENAME=$(basename $2)
EXTENSION=${FILENAME##*.}
#ARRAY=${FILENAME%.*}
ARRAY=$1
LOADARRAY=${ARRAY}_samples
GTARRAY=${ARRAY}_gt

SKIPLINES=`grep "^#" $FILE | wc -l`
ROWS=`grep -v "^#" $FILE | wc -l`
POPULATIONS=`grep -v ^# $FILE | cut -d, -f3 | sort -u | wc -l`
CHUNK=`${IQUERY} -a -o csv -q "dimensions(${GTARRAY})" | grep sampleid | cut -f4 -d,`
#echo skip: $SKIPLINES
#echo rows: $ROWS
#echo populations: $POPULATIONS
#echo chunk: $CHUNK

TMPDIR=`mktemp -d /tmp/sdload.XXXXXX` || exit 1
chmod 777 $TMPDIR

echo "Creating arrays"
$IQUERY -anq "create array ${LOADARRAY} <id: int64, sample: string, population: string, founder: bool, sex: char>[row=1:${ROWS}, ${CHUNK},0]"
#$IQUERY -anq "create empty array ${ARRAY} <id: int64, subject: string, founder: bool, sex: char>[population(string)=${POPULATIONS} ,1,0, row=1:${ROWS}, ${ROWS},0]"

echo "Preparing data"
LOADPIPE=$TMPDIR/load_pipe
mkfifo $LOADPIPE
chmod 666 $LOADPIPE
cat $FILE | csv2scidb -s $SKIPLINES -f 1 -c $ROWS -p NSSNC > $LOADPIPE &

echo "Loading data"
$IQUERY -anq "load(${LOADARRAY}, '${LOADPIPE}')"
#cat < $LOADPIPE | less
rm -rf $TMPDIR

#echo "Redimensioning"
#$IQUERY -anq "redimension_store(${LOADARRAY}, ${ARRAY})"
