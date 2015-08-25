iquery -olcsv+ -aq "
project(
 cross_join(
  redimension(
   apply(filter(list('chunk map'), inst=instn and attid=0), aid, int64(arrid)),
   <nchunks:uint64 null,
    min_ccnt:uint32 null,
    avg_ccnt:double null,
    max_ccnt:uint32 null,
    total_cnt: uint64 null>
   [aid= 0:*,1000,0],
   count(*) as nchunks,
   min(nelem) as min_ccnt,
   avg(nelem) as avg_ccnt,
   max(nelem) as max_ccnt,
   sum(nelem) as total_cnt
  ) as A,
  redimension(
   apply( list('arrays', true), aid, int64(id)),
   <name: string null>
   [aid = 0:*,1000,0]
  ) as B,
  A.aid, B.aid
 ),
 name, nchunks, min_ccnt, avg_ccnt, max_ccnt, total_cnt
)"
