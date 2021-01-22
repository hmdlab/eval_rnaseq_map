-- Settings
PRAGMA busy_timeout=864000000;
PRAGMA cache_size=48000000;

-- Aggregate
drop table if exists transcript_metrics;
create table if not exists transcript_metrics(qname, recall, precision, f1);

insert into transcript_metrics
  select
    qname,
    tp / (tp + fn) as recall,
    tp / (tp + fp) as precision,
    2 * tp / (2 * tp + fp + fn) as f1
  from
  (
  select
    qname,
    sum(pct_tp) as tp,
    sum(pct_fp) as fp,
    0.0 as tn,
    sum(pct_fn) as fn
  from
    (
    select
      qname,
      tp / total as pct_tp,
      fp / total as pct_fp,
      0.0 as pct_tn,
      fn / total as pct_fn
    from confusion_matrix
    ) t0
  group by qname
  ) t1;

-- XXX: Fix 0 devided NULL value
update transcript_metrics set recall = 0.0 where recall is null;
update transcript_metrics set precision = 0.0 where precision is null;
update transcript_metrics set f1 = 0.0 where f1 is null;

vacuum;
