create table annotations_unnormalized as
select
  `index`
   , seqname
   , `source
`
   ,feature
   ,`start`
   ,`
end`
   ,score
   ,strand
   ,frame
   ,a1.gene_id
   ,a2.geneSuperClass
   ,a2.geneClass
   ,a2.geneSubClass
   ,a2.gene_type
   ,a2.gene_name
   ,ifnull
(a3.coding_status, '') as coding_status
   ,ifnull
(a3.cumulative_support, '') as cumulative_support
   ,a2.geneCategory
   ,a2.DHS_type
   ,a1.transcript_id
   ,ifnull
(a3.transcript_type, '') as transcript_type
   ,ifnull
(a3.transcript_name, '') as transcript_name
   ,ifnull
(a3.TIEScore, '') as TIEScore
   ,exon_number
   from annotations a1
     left join
(select distinct gene_id, geneSuperClass, geneClass, geneSubClass, gene_type, gene_name, geneCategory, DHS_type
from annotations
where feature = 'gene')
a2
       using
(gene_id)
       left join
(select distinct transcript_id, coding_status, cumulative_support, transcript_type, transcript_name, '' as TIEScore
from annotations
where feature = 'transcript')
a3
         using
(transcript_id);

alter table annotations rename to annotations_normalized;

create view annotations
as
  select *
  from annotations_unnormalized;

vacuum;

