/*

  Usage:
    sqlite3 gencode.v31_refseq.v109.20190607.sqlite < this.sql

*/

-- Preparation
attach database 'shared/assets/references/grch38/annotations/refseq/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.sqlite' as refseq;
attach database 'shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.sqlite' as gencode;
attach database 'shared/assets/references/grch38/annotations/gencode_refseq/cuffcmp.combined.sqlite' as cuffcmp;

-- create conv table (tx-level and gene-level)
create table conv as 
select cuffcmp.oId, refseq.gene_id as oGId, cuffcmp.nearest_ref, gencode.gene_id as nearest_gid
  from (select distinct oId, nearest_ref from cuffcmp.annotations where oId <> '' and nearest_ref <> '') cuffcmp
    inner join (select gene_id, transcript_id from refseq.annotations) refseq
      on cuffcmp.oId = refseq.transcript_id
      inner join (select distinct gene_id, transcript_id from gencode.annotations where feature = 'exon') gencode
      on cuffcmp.nearest_ref = gencode.transcript_id;

create table annotations as 
-- refseq only
-- 107001
select seqname, source, feature, `start`, `end`, score, strand, frame, ifnull(conv.nearest_gid, gene_id) as gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name, exon_number, '' as oId, ifnull(conv.oGId, '') as oGid, '' as nearest_ref, '' as class_code 
  from
  (select seqname, source, feature, `start`, `end`, score, strand, frame, gene_id, gbkey as gene_type, gene_name, transcript_id, class as transcript_type, transcript_name, exon_number from refseq.annotations
  -- select count(distinct transcript_id) from refseq.annotations
    where transcript_id not in (
    select distinct oId from cuffcmp.annotations
      where class_code = '='
    )) t1
    left join (select oGId, nearest_gid from conv where class_code = '=') conv
    on conv.oGId = t1.gene_id
union all
-- gencode only
-- 173087
select seqname, source, feature, `start`, `end`, score, strand, frame, gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name, exon_number, '' as oId, '' as oGId, '' as nearest_ref, '' as class_code from gencode.annotations
-- select count(distinct transcript_id) from gencode.annotations
  where feature = 'exon' and transcript_id not in (
    select distinct nearest_ref from cuffcmp.annotations
      where class_code = '='
  )
union all
-- intersection
-- 53795 / 53795 -> 1:1
select gencode.seqname, gencode.source, gencode.feature, gencode.`start`, gencode.`end`, gencode.score, gencode.strand, gencode.frame, gencode.gene_id, gencode.gene_type, gencode.gene_name, gencode.transcript_id, gencode.transcript_type, gencode.transcript_name, gencode.exon_number, cuffcmp.oId, '' as oGId, cuffcmp.nearest_ref, cuffcmp.class_code from gencode.annotations gencode
  inner join (select * from cuffcmp.annotations where class_code = '=') cuffcmp
    on gencode.transcript_id = cuffcmp.nearest_ref and gencode.exon_number = cuffcmp.exon_number
    where gencode.feature = 'exon'
;

-- format refseq record matching to gencode record
update annotations set
  score = '.'
  where score is null;

update annotations set
  frame = '0'
  where frame is not null; 

-- merge at gene-level

-- 72557 gene / 333883 tx
select count(distinct transcript_id), count(distinct gene_id) from annotations;

-- gene detail 72557(e:60603, r:11954, i:27326/27236)
select src, sum(cnt) as cnt from
  (select case when prefix = 'ENSG' then 'ensembl' else 'refseq' end as src, cnt from
  (select substr(gene_id, 1, 4) as prefix, count(distinct gene_id) as cnt from annotations
    group by substr(gene_id, 1, 4)) t1
   ) t2
   group by src

-- transcript detail 333883 (e: 226882, r: 107001, i: 53795)
select src, sum(cnt) as cnt from
  (select case when prefix = 'ENST' then 'ensembl' else 'refseq' end as src, cnt from
  (select substr(transcript_id, 1, 4) as prefix, count(distinct transcript_id) as cnt from annotations
    group by substr(transcript_id, 1, 4)) t1
   ) t2
   group by src

-- intersect
select count(distinct oId), count(distinct nearest_ref) from conv where class_code = '=';
select count(distinct oGId), count(distinct nearest_gid) from conv where class_code = '='   

-- 1:n converted genes
select oGId, count(distinct nearest_gid), group_concat(distinct nearest_gid), group_concat(distinct oId) from conv
  where class_code = '='
  group by oGId
  having count(distinct nearest_gid) > 1;

-- n:1 converted genes
select nearest_gid, count(distinct oGId), group_concat(distinct oGId), group_concat(distinct nearest_ref) from conv
  where class_code = '='
  group by nearest_gid
  having count(distinct oGId) > 1;

select count(case when prefix = 'ENSG' then 1 else 0 end) as ensembl_gene, count(case when prefix <> 'ENSG' then 1 else 0 end) as refseq_gene from 
  (select substr(gene_id, 1, 4) as prefix, count(distinct gene_id) from annotations
    group by substr(gene_id, 1, 4)) t1

select case when substr(gene_id, 1, 2) = 'ENSG' then count(gene_id) else count(*) end as count from annotations;

select * from annotations
  where seqname = 'chr1'
  order by seqname, gene_id, transcript_id, cast(exon_number as int)

select * from annotations
  where seqname = 'chr1'
  and transcript_id like 'NM%'
  and gene_id like 'ENSG%'
  order by seqname, gene_id, transcript_id, cast(exon_number as int)
  
select * from 
  (select distinct gene_id, transcript_id from annotations) an
  left join conv
    on an.transcript_id = conv.oId;

-- NOTE: Required SQLite 3.25 above
create table transcripts_ordered as 
select row_number() over(order by seqname, start_min, transcript_id) as `index`, transcript_id, seqname, start_min from
(select seqname, transcript_id, min(start) as start_min from annotations
  group by seqname, transcript_id) t1;

create table annotations_ordered as 
select row_number() over(order by t2.`index`, t1.exon_number) as `index`, t1.* from annotations t1
  left join transcripts_ordered t2
    using(transcript_id)
    order by t2.`index`, t1.exon_number;
