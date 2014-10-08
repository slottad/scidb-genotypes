#!/bin/bash

usage()
{
cat <<EOF
usage: $0 [options] dataset filename(s)

Load a study.

OPTIONS:
   -h/-?   Show this message
   -c      chunksize for load arrays (default: 200,000)
   -m      max reference length to load (default: 50)
   -d      SciDB coordinator system (default: localhost)
   -p      SciDB port
   -s      file containing the list of samples (required)
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
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IQUERY=/opt/scidb/$SCIDB_VER/bin/iquery
VCF2CSV=$SCRIPTDIR/vcf2csv
DIMLOAD=$SCRIPTDIR/dim_array_load.py
POPLOAD=$SCRIPTDIR/populations_create.py
LOADCSV=/opt/scidb/$SCIDB_VER/bin/loadcsv.py

chunksize=250000
refsize=50
dbsystem="localhost"
samples=""
port=1239
while getopts "hc:d:s:p:?" flag
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

        d)
            dbsystem=$OPTARG
            ;;
        p)
            port=$OPTARG
            ;;

        s)
            samples=$OPTARG
            ;;
        ?)
            usage
            exit 1
            ;;
    esac
done

if [[ ! $samples ]]; then
    echo "The sample file is a required argument"
    usage
    exit 1
fi

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

echo "calling: $DIMLOAD -a ${array} -d sample -c ${dbsystem} -p ${port} ${samples}"
$DIMLOAD -a ${array} -d sample -c ${dbsystem} -p ${port} ${samples}

echo "calling: $POPLOAD -c ${dbsystem} -p ${port} ${array}"
$POPLOAD -c ${dbsystem} -p ${port} ${array}

var_load_array=${array}_var_load
gt_load_array=${array}_gt_load

var_array=${array}_var
gt_array=${array}_gt

echo "Preparing data"
tmpdir=`mktemp -d /tmp/sdload.XXXXXX` || exit 1
trap cleanup EXIT
chmod 777 $tmpdir
varloadpipe=$tmpdir/varload_pipe
gtloadpipe=$tmpdir/gtload_pipe
mkfifo $varloadpipe
mkfifo $gtloadpipe
chmod 666 $varloadpipe
chmod 666 $gtloadpipe

options="-s ${samples} ${varloadpipe} ${gtloadpipe}"
case $2 in
    *.gz)
        decompress=zcat
        ;;
    *.bz2)
        decompress=bzcat
        ;;
esac


# these loadcsv calls receive data from the fifo written to by vcf2csv
var_load_array_def="<chrom: string, pos: int64, var: int64, id: string null, ref: string, alt: string, alleles: uint32, qual: float null, filter: string null, info: string null, format: string null >[row=0:*,${chunksize},0]"
echo "var_load fifo reader: $LOADCSV -x -q -d ${dbsystem} -p ${port} -D'\t' -a ${var_load_array} -s \"$var_load_array_def\" -i $varloadpipe"
$LOADCSV -x -q -d ${dbsystem} -p ${port} -D'\t' -a ${var_load_array} -s "$var_load_array_def" -i $varloadpipe &


gt_load_array_def="<chrom:string,pos:int64,var:int64,sampleid:int64,gt:gt8 NULL,unparsed:string null>[row=0:*,${chunksize},0]"
echo "gt_load fifo reader: $LOADCSV -x -q -d ${dbsystem} -p ${port} -D'\t' -a ${gt_load_array} -s \"$gt_load_array_def\" -i $gtloadpipe"
$LOADCSV -x -q -d ${dbsystem} -p ${port} -D'\t' -a ${gt_load_array} -s "$gt_load_array_def" -i $gtloadpipe &


# vcf2csv splits the data stream sending output to two named pipes
if [[ -z $decompress ]]; then
    result=(`pv ${files} | ${VCF2CSV} ${options}`)
else
    result=(`pv ${files} | ${decompress} | ${VCF2CSV} ${options}`)
fi

exit $?
