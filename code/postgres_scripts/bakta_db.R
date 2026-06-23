library(tidyverse)
library(DBI)

pg <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),dbname = "glamr_data", host = "localhost", port = "5437", user = "glamr_admin", password = "glamr2023")

bakta_con <- DBI::dbConnect(RSQLite::SQLite(),"~/GLAMR/data/reference/bakta/db/bakta.db")
bakta_db <- RSQLite::dbListTables(bakta_con)

# Make pointers to bakta reference tables
unique_seqs <- tbl(bakta_con, "ups")

unique_seqs %>% collect() %>% dbWriteTable(pg,"unique_uniref_seqs",.,overwrite = TRUE)