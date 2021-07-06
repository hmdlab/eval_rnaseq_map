library(tidyverse)

path_sqlite <- here::here("assets/references_tailored/annotations/gencode.v31.annotation_expressed.sqlite")

conn <- RSQLite::SQLite() %>% RSQLite::dbConnect(path_sqlite, synchronous = "off")

annotations <- RSQLite::dbGetQuery(conn, paste0("select * from annotations;"))

RSQLite::dbDisconnect(conn)

lengths_transcript <- annotations %>%
  filter(feature == 'exon') %>%
  mutate(.length = end - start + 1) %>%
  group_by(transcript_id) %>%
  summarize(length = sum(.length))

path_header <- here::here("results/testXX_main_r/header.txt")

header <- read_tsv(path_header) %>%
  set_names(c("cat", "seq", "len"))

header <- header %>%
  mutate(transcript_id = map_chr(seq, ~ gsub("SN:", "", .x))) %>%
  mutate(length = map_chr(len, ~ gsub("LN:", "", .x))) %>%
  select(transcript_id, length)

lengths_transcript %>%
  left_join(header, by = "transcript_id") %>%
  filter(length.x == length.y)