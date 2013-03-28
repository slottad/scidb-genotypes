#!/bin/bash
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

usage()
{
cat <<EOF
usage: $0 [options] dataset filename(s)

List the variants in a given region. Limited to a single chromosome.

OPTIONS:
   -h/-?   Show this message
   -c      chunksize for load arrays (default: 200,000)
   -m      max reference length to load (default: 50)
   -t      use text mode
   -b      use binary mode (default)
   -x      test mode
   -s      speed test mode
   -d      CVS file containing sample descriptions 
EOF
}

cleanup()
{
    echo "Removing temporary files."
    rm $varloadpipe
    rm $gtloadpipe
    rmdir $tmpdir
}

if [ -z "$SCIDB_VER" ]; then
    echo "SCIDB_VER must be set."
    exit 1
fi

# Parameters
IQUERY=/opt/scidb/$SCIDB_VER/bin/iquery
VCF2SCIDB=vcf2scidb

binary=true
testmode=false
speedtest=false
hasdesc=false
chunksize=200000
refsize=50
while getopts "hc:tbxsd:?" flag
do
    case $flag in
        h)
            echo "Help called"
            usage
            exit 1
            ;;
        c)
            chunksize=$OPTARG
            ;;
        m)
            refsize=$OPTARG
            ;;
        t)
            binary=false
            ;;
        b)
            binary=true
            ;;
        x)
            testmode=true
            ;;           
        s)
            testmode=true
            speedtest=true
            ;;           
        d)
            hasdesc=true
            descriptions=$OPTARG
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

# remaining args
shift $(($OPTIND - 1))
numargs=$#
if [ $numargs -lt 2 ]; then
    echo "Arguments missing"
    usage
    exit 1
else
    array=$1
    files="${@:2}"
fi

var_load_array=${array}_load
gt_load_array=${array}_gt_load
var_array=${array}
gt_array=${array}_gt

if ! $testmode; then
    #echo "Creating load arrays"
    echo -n "Create variation load array: "
    $IQUERY -anq "create array ${var_load_array} < chrom: string, pos: int64, var: int64, id: string null, ref: string, alt: string, alleles: uint32, qual: float null, filter: string null, info: string null, format: string null >[row=1:*,${chunksize},0]"
    echo -n "Create genotype load array: "
    $IQUERY -anq "create array ${gt_load_array} < chrom: string, pos: int64, var: int64, sampleid: int64, gt: gt8, unparsed: string null >[row=1:*,${chunksize},0]"
fi

echo "Preparing data"
if $testmode; then
    if $speedtest; then
        varloadpipe=/dev/null
        gtloadpipe=/dev/null
    else
        varloadpipe=${array}_var.scidb
        gtloadpipe=${array}_gt.scidb
    fi
else
    tmpdir=`mktemp -d /tmp/sdload.XXXXXX` || exit 1
    chmod 777 $tmpdir
    varloadpipe=$tmpdir/varload_pipe
    gtloadpipe=$tmpdir/gtload_pipe
    mkfifo $varloadpipe
    mkfifo $gtloadpipe
    chmod 666 $varloadpipe
    chmod 666 $gtloadpipe
    trap cleanup EXIT
fi

if $binary; then
    options="-b"
else
    options="-t -c ${chunksize}"
fi
if $hasdesc; then
    options="${options} -d ${descriptions}"
fi
options="${options} -m ${refsize} -v ${varloadpipe} -g ${gtloadpipe}"

#echo "Options: ${options}"

case $2 in
    *.gz)
        decompress=zcat
        ;;
    *.bz2)
        decompress=bzcat
        ;;
esac


if ! $testmode; then
    if $binary; then
        echo "Loading binary data"

        $IQUERY -anq "load(${var_load_array}, '${varloadpipe}', -2, '(string, int64, int64, string null, string, string, uint32, float null, string null, string null, string null)')" &
        $IQUERY -anq "load(${gt_load_array}, '${gtloadpipe}', -2, '(string, int64, int64, int64, gt8, string null)')" &
    else
        echo "Loading text data"
        $IQUERY -anq "load(${var_load_array}, '${varloadpipe}')" &
        $IQUERY -anq "load(${gt_load_array}, '${gtloadpipe}')" &
    fi
fi

if [[ -z $decompress ]]; then
    result=(`pv ${files} | ${VCF2SCIDB} ${options}`)
else
    result=(`pv ${files} | ${decompress} | ${VCF2SCIDB} ${options}`)
fi

wait
chroms=${result[0]}
max_pos=${result[1]}
max_var=${result[2]}
max_sampleid=${result[3]}

# echo "Chroms: ${chroms}" 
# echo "max_pos: ${max_pos}" 
# echo "max_var: ${max_var}" 
# echo "max_sampleid: ${max_sampleid}"

if [ $max_sampleid -lt 5 ]; then
    sample_chunksize=$max_sampleid
elif [ $max_sampleid -lt 40 ]; then
    sample_chunksize=$[$max_sampleid/2]
else
    sample_chunksize=$[$max_sampleid/4]
fi

# echo "sample_chunksize: ${sample_chunksize}"

if ! $testmode; then
    echo -n "Create variation array: "
    $IQUERY -anq "create empty array ${var_array} < id: string null, ref: string, alt: string, alleles: uint32, qual: float null, filter: string null, info: string null, format: string null >[ chrom(string)=${chroms},1,0, pos=1:${max_pos},200000,0, var=1:${max_var},${max_var},0 ]"
    echo -n "Create genotype array: "
    $IQUERY -anq "create empty array ${gt_array} < gt: gt8, unparsed: string null >[ chrom(string)=${chroms},1,0, pos=1:${max_pos},200000,0, var=1:${max_var},${max_var},0, sampleid=1:${max_sampleid},${sample_chunksize},0 ]"
    echo "Redimensioning..."
    /usr/bin/time -f "Created variation array, elapsed time: %E\t" $IQUERY -anq "redimension_store(${var_load_array},${var_array})"
    /usr/bin/time -f "Created genotype array, elapsed time: %E\t" $IQUERY -anq "redimension_store(${gt_load_array},${gt_array})"
fi

echo "Finished!"
