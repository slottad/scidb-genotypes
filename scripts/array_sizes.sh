#FILTER="$1@1"
#echo $FILTER

iquery -olcsv+ -aq "
apply(
project(
 cross_join(
  redimension(
   apply(filter(list('chunk map'), inst=instn), aid, int64(arrid)),
   <total_size: uint64 null>
   [aid= 0:*,1000,0],
   sum(asize) as total_size
  ) as A,
  redimension(
   apply( list('arrays',true), aid, int64(id)),
   <name: string null>
   [aid = 0:*,1000,0]
  ) as B,
  A.aid, B.aid
 ),
 name, total_size
),total_size_human,human_readable(total_size)
)"
